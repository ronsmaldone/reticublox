// Reticublox — read 3D structure files from other programs (CrystalMaker, Mercury, Avogadro, VESTA,
// Materials Studio…) into a V2000 molblock that RDKit can read, keeping the real 3D coordinates.
//   .mol / .sdf   bonds and orders from the file (first record)
//   .mol2         bonds and orders from the file (ar = aromatic)
//   .xyz / .pdb   bonds from interatomic distances, bond orders from the geometry
//   .cif          symmetry-expanded, molecules rebuilt across the cell edges, bonds as for .xyz.
//                 Metals are left out, so a MOF CIF falls apart into its organic linkers.
import { v3 } from './nets.js';

export const IMPORT_ACCEPT = '.mol,.sdf,.sd,.mol2,.xyz,.pdb,.ent,.cif';

const COV = { H: 0.31, B: 0.84, C: 0.76, N: 0.71, O: 0.66, F: 0.57, Si: 1.11, P: 1.07, S: 1.05, Cl: 1.02, Se: 1.20, Br: 1.20, I: 1.39 };
const VAL = { H: 1, B: 3, C: 4, N: 3, O: 2, F: 1, Si: 4, P: 3, S: 2, Cl: 1, Se: 2, Br: 1, I: 1 };
const isOrganic = el => el in COV;

/** → { molblock, name, notes[] } — throws with a readable message. */
export function importStructure(text, filename = '') {
  const ext = (filename.match(/\.([a-z0-9]+)$/i)?.[1] || '').toLowerCase();
  const base = filename.replace(/\.[^.]+$/, '') || 'Imported structure';
  if (ext === 'mol' || ext === 'sdf' || ext === 'sd' || (!ext && /V[23]000/.test(text))) {
    const rec = text.split(/^\$\$\$\$/m)[0];
    if (!/M\s+END/.test(rec)) throw new Error('That MOL/SDF file has no “M  END” line.');
    const title = rec.split(/\r?\n/)[0].trim();
    return { molblock: rec.replace(/\r/g, ''), name: title || base, notes: [] };
  }
  let g;
  if (ext === 'mol2') g = readMol2(text);
  else if (ext === 'xyz') g = readXYZ(text);
  else if (ext === 'pdb' || ext === 'ent') g = readPDB(text);
  else if (ext === 'cif') g = readCIF(text);
  else throw new Error(`Unsupported file type “.${ext}”. Use .mol, .sdf, .mol2, .xyz, .pdb or .cif.`);
  const notes = [...(g.notes || [])];
  let { atoms, bonds } = g;
  if (!atoms.length) throw new Error('No atoms found in that file.');

  // metals: SBUs come from the Metal cluster cards; the organic part is what gets imported
  const metals = atoms.filter(a => !isOrganic(a.el));
  if (metals.length) {
    notes.push(`${metals.length} metal atom${metals.length > 1 ? 's' : ''} (${[...new Set(metals.map(a => a.el))].join(', ')}) left out — use a Metal cluster card for the SBU.`);
    const keep = atoms.map(a => isOrganic(a.el));
    const idx = []; let k = 0;
    atoms.forEach((a, i) => { idx[i] = keep[i] ? k++ : -1; });
    atoms = atoms.filter((_, i) => keep[i]);
    if (bonds) bonds = bonds.filter(b => keep[b.a] && keep[b.b]).map(b => ({ ...b, a: idx[b.a], b: idx[b.b] }));
  }
  if (!bonds) bonds = perceiveBonds(atoms);
  // several molecules (solvent, counter-ions, a whole crystal's worth): keep the largest
  ({ atoms, bonds } = largestComponent(atoms, bonds, notes));
  if (bonds.some(b => !b.order)) assignOrders(atoms, bonds);
  return { molblock: writeV2000(atoms, bonds, g.name || base), name: g.name || base, notes };
}

// ---------------------------------------------------------------- readers
const elOf = s => {
  const m = String(s).match(/^[A-Za-z]{1,2}/);
  if (!m) return 'X';
  const e = m[0][0].toUpperCase() + (m[0][1] || '').toLowerCase();
  return (e in COV || METALS.has(e)) ? e : e[0];
};
const METALS = new Set(['Li', 'Na', 'K', 'Mg', 'Ca', 'Al', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Zr', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Hf', 'W', 'Pt', 'Au', 'Hg', 'Pb', 'Ln', 'La', 'Ce', 'Eu', 'Gd', 'Tb', 'Y', 'Sc', 'Ba', 'Sr', 'Cs', 'Rb', 'Bi', 'U', 'Th']);

function readXYZ(text) {
  const L = text.replace(/\r/g, '').split('\n');
  const n = parseInt(L[0]);
  const atoms = [];
  for (const l of L.slice(2, isNaN(n) ? undefined : 2 + n)) {
    const t = l.trim().split(/\s+/);
    if (t.length < 4 || isNaN(+t[1])) continue;
    atoms.push({ el: elOf(t[0]), p: [+t[1], +t[2], +t[3]] });
  }
  return { atoms, bonds: null, name: L[1]?.trim() };
}

function readPDB(text) {
  const atoms = [];
  let name = '';
  for (const l of text.replace(/\r/g, '').split('\n')) {
    if (/^(COMPND|TITLE)/.test(l) && !name) name = l.slice(10).trim();
    if (!/^(ATOM|HETATM)/.test(l)) continue;
    const el = (l.slice(76, 78).trim() || l.slice(12, 16).trim().replace(/\d/g, ''));
    atoms.push({ el: elOf(el), p: [+l.slice(30, 38), +l.slice(38, 46), +l.slice(46, 54)] });
  }
  return { atoms, bonds: null, name };
}

function readMol2(text) {
  const sec = name => (text.split(`@<TRIPOS>${name}`)[1] || '').split('@<TRIPOS>')[0].trim().split(/\r?\n/).filter(Boolean);
  const molL = sec('MOLECULE');
  const ids = new Map(), atoms = [];
  for (const l of sec('ATOM')) {
    const t = l.trim().split(/\s+/);
    ids.set(t[0], atoms.length);
    atoms.push({ el: elOf(t[5]?.split('.')[0] || t[1]), p: [+t[2], +t[3], +t[4]] });
  }
  const bonds = [];
  for (const l of sec('BOND')) {
    const t = l.trim().split(/\s+/);
    const order = { 1: 1, 2: 2, 3: 3, ar: 4, am: 1 }[t[3]] || 1;
    bonds.push({ a: ids.get(t[1]), b: ids.get(t[2]), order });
  }
  return { atoms, bonds, name: molL[0] };
}

// CIF: cell, symmetry operators, atom sites → molecules rebuilt across the cell edges
function readCIF(text) {
  const num = s => parseFloat(String(s).replace(/\(.*\)/, ''));
  const tokens = tokenizeCIF(text);
  const tags = {}, loops = [];
  for (let i = 0; i < tokens.length; i++) {
    const t = tokens[i];
    if (/^loop_$/i.test(t)) {
      const names = [];
      while (tokens[i + 1]?.startsWith('_')) names.push(tokens[++i].toLowerCase());
      const rows = [];
      while (i + 1 < tokens.length && !tokens[i + 1].startsWith('_') && !/^(loop_|data_)/i.test(tokens[i + 1])) {
        const row = [];
        for (let k = 0; k < names.length && i + 1 < tokens.length; k++) row.push(tokens[++i]);
        rows.push(row);
      }
      loops.push({ names, rows });
    } else if (t.startsWith('_')) tags[t.toLowerCase()] = tokens[++i];
  }
  const cell = ['a', 'b', 'c'].map(x => num(tags[`_cell_length_${x}`]));
  const angle = ['alpha', 'beta', 'gamma'].map(x => num(tags[`_cell_angle_${x}`] ?? 90) * Math.PI / 180);
  if (cell.some(isNaN)) throw new Error('This CIF has no unit cell (_cell_length_a/b/c).');
  const M = cellMatrix(cell, angle);
  const loopWith = name => loops.find(l => l.names.includes(name));
  const symL = loopWith('_symmetry_equiv_pos_as_xyz') || loopWith('_space_group_symop_operation_xyz');
  const symCol = symL ? symL.names.findIndex(n => /as_xyz|operation_xyz/.test(n)) : -1;
  const ops = (symL ? symL.rows.map(r => r[symCol]) : ['x,y,z']).map(parseSymop);
  const at = loopWith('_atom_site_fract_x');
  if (!at) throw new Error('This CIF has no fractional atom coordinates (_atom_site_fract_x).');
  const col = n => at.names.indexOf(n);
  const [cx, cy, cz, cl, ct, co] = ['_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z', '_atom_site_label', '_atom_site_type_symbol', '_atom_site_occupancy'].map(col);
  const notes = [];
  // unique atoms in the cell (fractional, wrapped to [0,1))
  const cellAtoms = [];
  let lowOcc = 0;
  for (const r of at.rows) {
    if (co >= 0 && num(r[co]) < 0.5) { lowOcc++; continue; } // minor disorder component
    const el = elOf(ct >= 0 && r[ct] !== '?' && r[ct] !== '.' ? r[ct] : r[cl]);
    const f0 = [num(r[cx]), num(r[cy]), num(r[cz])];
    for (const op of ops) {
      const f = op(f0).map(x => x - Math.floor(x + 1e-9));
      if (!cellAtoms.some(a => a.el === el && cartDist(M, a.f, f) < 0.4)) cellAtoms.push({ el, f });
    }
  }
  if (lowOcc) notes.push(`${lowOcc} disordered site${lowOcc > 1 ? 's' : ''} with occupancy < 0.5 skipped.`);
  // MOF: drop the metals first, so each organic linker is its own finite molecule
  if (cellAtoms.some(a => !isOrganic(a.el))) return readCIFOrganic(cellAtoms, M, notes, tags);
  // bonds across the cell edges (minimum image); unwrap each molecule into one piece
  const n = cellAtoms.length;
  if (n > 3000) throw new Error('This CIF has too many atoms in the cell to rebuild molecules here.');
  const nb = cellAtoms.map(() => []);
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    const [d, shift] = minImage(M, cellAtoms[i].f, cellAtoms[j].f);
    if (bonded(cellAtoms[i].el, cellAtoms[j].el, d)) { nb[i].push([j, shift]); nb[j].push([i, shift.map(s => -s)]); }
  }
  const off = new Array(n).fill(null);
  const comps = [];
  for (let s = 0; s < n; s++) {
    if (off[s]) continue;
    off[s] = [0, 0, 0];
    const comp = [s];
    let infinite = false;
    for (let k = 0; k < comp.length; k++) {
      const i = comp[k];
      for (const [j, sh] of nb[i]) {
        const o = v3.add(off[i], sh);
        if (!off[j]) { off[j] = o; comp.push(j); } else if (v3.len(v3.sub(off[j], o)) > 0.5) infinite = true;
      }
    }
    comps.push({ comp, infinite });
  }
  const organicCount = c => c.comp.filter(i => isOrganic(cellAtoms[i].el)).length;
  const finite = comps.filter(c => !c.infinite && organicCount(c) > 1);
  if (!finite.length) throw new Error(INFINITE);
  finite.sort((a, b) => organicCount(b) - organicCount(a));
  const c = finite[0].comp;
  return {
    atoms: c.map(i => ({ el: cellAtoms[i].el, p: frac2cart(M, v3.add(cellAtoms[i].f, off[i])) })),
    bonds: null, notes,
    name: (tags._chemical_name_common || tags._chemical_name_systematic || '').replace(/^'|'$/g, ''),
  };
}
// MOF CIF: drop the metals first, then each organic linker is its own finite molecule
function readCIFOrganic(cellAtoms, M, notes, tags) {
  const metals = cellAtoms.filter(a => !isOrganic(a.el));
  const org = cellAtoms.filter(a => isOrganic(a.el));
  notes.push(`Metal atoms (${[...new Set(metals.map(a => a.el))].join(', ')}) left out; the largest organic linker was taken.`);
  const n = org.length, nb = org.map(() => []);
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    const [d, shift] = minImage(M, org[i].f, org[j].f);
    if (bonded(org[i].el, org[j].el, d)) { nb[i].push([j, shift]); nb[j].push([i, shift.map(s => -s)]); }
  }
  const off = new Array(n).fill(null);
  let best = null;
  for (let s = 0; s < n; s++) {
    if (off[s]) continue;
    off[s] = [0, 0, 0];
    const comp = [s];
    let infinite = false;
    for (let k = 0; k < comp.length; k++) for (const [j, sh] of nb[comp[k]]) {
      const o = v3.add(off[comp[k]], sh);
      if (!off[j]) { off[j] = o; comp.push(j); } else if (v3.len(v3.sub(off[j], o)) > 0.5) infinite = true;
    }
    if (!infinite && (!best || comp.length > best.length)) best = comp;
  }
  if (!best || best.length < 2) throw new Error(INFINITE);
  return {
    atoms: best.map(i => ({ el: org[i].el, p: frac2cart(M, v3.add(org[i].f, off[i])) })),
    bonds: null, notes,
    name: (tags._chemical_name_common || '').replace(/^'|'$/g, '') || 'Linker from CIF',
  };
}

const INFINITE = 'This CIF is an endless covalent network (e.g. a COF), so there is no single molecule to cut out. In your crystal program, export one building block (with its linking groups) as .mol or .xyz instead.';

function tokenizeCIF(text) {
  const out = [];
  const L = text.replace(/\r/g, '').split('\n');
  for (let i = 0; i < L.length; i++) {
    const line = L[i];
    if (line.startsWith(';')) { // multi-line text field
      let s = line.slice(1);
      while (++i < L.length && !L[i].startsWith(';')) s += '\n' + L[i];
      out.push(s.trim());
      continue;
    }
    const re = /'((?:[^']|'(?=\S))*)'(?=\s|$)|"([^"]*)"(?=\s|$)|(#.*$)|(\S+)/g;
    let m;
    while ((m = re.exec(line))) {
      if (m[3]) break;
      out.push(m[1] ?? m[2] ?? m[4]);
    }
  }
  return out;
}
function parseSymop(s) {
  const parts = String(s).replace(/\s/g, '').toLowerCase().split(',');
  const fns = parts.map(expr => {
    // terms like -x, +y, 1/2, -x+y, 0.5+z
    const terms = expr.match(/[+-]?[^+-]+/g) || [];
    const c = { x: 0, y: 0, z: 0, k: 0 };
    for (const t of terms) {
      const m = t.match(/^([+-]?)(\d*\.?\d*(?:\/\d+)?)\*?([xyz]?)$/);
      if (!m) continue;
      const sg = m[1] === '-' ? -1 : 1;
      let v = m[2] ? (m[2].includes('/') ? m[2].split('/').reduce((a, b) => +a / +b) : +m[2]) : 1;
      if (m[3]) c[m[3]] += sg * v; else c.k += sg * v;
    }
    return f => c.x * f[0] + c.y * f[1] + c.z * f[2] + c.k;
  });
  return f => fns.map(fn => fn(f));
}
function cellMatrix([a, b, c], [al, be, ga]) {
  const cx = c * Math.cos(be), cy = c * (Math.cos(al) - Math.cos(be) * Math.cos(ga)) / Math.sin(ga);
  return [[a, 0, 0], [b * Math.cos(ga), b * Math.sin(ga), 0], [cx, cy, Math.sqrt(Math.max(0, c * c - cx * cx - cy * cy))]];
}
const frac2cart = (M, f) => [0, 1, 2].map(k => f[0] * M[0][k] + f[1] * M[1][k] + f[2] * M[2][k]);
function minImage(M, f1, f2) {
  const d = [0, 1, 2].map(k => f2[k] - f1[k]);
  const sh = d.map(x => -Math.round(x));
  let best = [Infinity, sh];
  // check neighbouring images too (oblique cells)
  for (const dx of [-1, 0, 1]) for (const dy of [-1, 0, 1]) for (const dz of [-1, 0, 1]) {
    const s = [sh[0] + dx, sh[1] + dy, sh[2] + dz];
    const L = v3.len(frac2cart(M, v3.add(d, s)));
    if (L < best[0]) best = [L, s];
  }
  return best;
}
const cartDist = (M, f1, f2) => minImage(M, f1, f2)[0];

// ---------------------------------------------------------------- bonds & orders
const bonded = (e1, e2, d) => isOrganic(e1) && isOrganic(e2) && !(e1 === 'H' && e2 === 'H') &&
  d > 0.4 && d < COV[e1] + COV[e2] + 0.4;
function perceiveBonds(atoms) {
  const bonds = [];
  for (let i = 0; i < atoms.length; i++) for (let j = i + 1; j < atoms.length; j++) {
    const d = v3.len(v3.sub(atoms[i].p, atoms[j].p));
    if (bonded(atoms[i].el, atoms[j].el, d)) bonds.push({ a: i, b: j, order: 0 });
  }
  return bonds;
}
function largestComponent(atoms, bonds, notes) {
  const adj = atoms.map(() => []);
  for (const b of bonds) { adj[b.a].push(b.b); adj[b.b].push(b.a); }
  const comp = new Array(atoms.length).fill(-1);
  const sizes = [];
  for (let s = 0; s < atoms.length; s++) {
    if (comp[s] >= 0) continue;
    const c = sizes.length, q = [s];
    comp[s] = c;
    for (let k = 0; k < q.length; k++) for (const w of adj[q[k]]) if (comp[w] < 0) { comp[w] = c; q.push(w); }
    sizes.push(q.filter(i => atoms[i].el !== 'H').length);
  }
  if (sizes.length === 1) return { atoms, bonds };
  const big = sizes.indexOf(Math.max(...sizes));
  notes.push(`The file holds ${sizes.length} separate molecules or fragments; the largest (${sizes[big]} heavy atoms) was taken.`);
  const idx = []; let k = 0;
  atoms.forEach((a, i) => { idx[i] = comp[i] === big ? k++ : -1; });
  return {
    atoms: atoms.filter((_, i) => comp[i] === big),
    bonds: bonds.filter(b => comp[b.a] === big).map(b => ({ ...b, a: idx[b.a], b: idx[b.b] })),
  };
}

/**
 * Bond orders from the 3D geometry: each atom's spare valence (π electrons it can share) comes from
 * its neighbour count, bond angles and bond lengths, so it works with or without hydrogens in the file.
 * The spare valences are then paired up along bonds, most constrained atoms first (a Kekulé structure).
 */
function assignOrders(atoms, bonds) {
  for (const b of bonds) b.order = 1;
  const nbr = atoms.map(() => []);
  bonds.forEach((b, k) => { nbr[b.a].push([b.b, k]); nbr[b.b].push([b.a, k]); });
  const len = k => v3.len(v3.sub(atoms[bonds[k].a].p, atoms[bonds[k].b].p));
  const angle = (i, j, k) => {
    const u = v3.norm(v3.sub(atoms[j].p, atoms[i].p)), w = v3.norm(v3.sub(atoms[k].p, atoms[i].p));
    return Math.acos(Math.max(-1, Math.min(1, v3.dot(u, w)))) * 180 / Math.PI;
  };
  const spare = atoms.map((a, i) => {
    const N = nbr[i], deg = N.length, el = a.el;
    const minL = deg ? Math.min(...N.map(([, k]) => len(k))) : 9;
    const ang = deg >= 2 ? angle(i, N[0][0], N[1][0]) : 0;
    if (el === 'C') {
      if (deg >= 4) return 0;
      if (deg === 3) {
        const sum = angle(i, N[0][0], N[1][0]) + angle(i, N[0][0], N[2][0]) + angle(i, N[1][0], N[2][0]);
        return sum > 350 ? 1 : 0; // planar → sp2
      }
      if (deg === 2) return ang > 165 ? 2 : (ang < 115 && minL > 1.46) ? 0 : 1;
      if (deg === 1) return minL < 1.24 ? 2 : minL < 1.42 ? 1 : 0;
      return 0;
    }
    if (el === 'N') {
      if (deg >= 3) return 0;
      if (deg === 2) return ang > 165 ? 2 : minL < 1.38 ? 1 : 0;
      if (deg === 1) return minL < 1.20 ? 2 : minL < 1.32 ? 1 : 0;
      return 0;
    }
    if (el === 'O' || el === 'S' || el === 'Se') return deg === 1 && minL < (el === 'O' ? 1.30 : 1.75) ? 1 : 0;
    return 0;
  });
  // pair spare valences along bonds; fewest options first, shortest bond as tie-break
  const opts = i => nbr[i].filter(([j]) => spare[j] > 0);
  for (let guard = 0; guard < atoms.length * 3; guard++) {
    let pick = -1, few = Infinity;
    spare.forEach((s, i) => { if (s > 0) { const o = opts(i).length; if (o && o < few) { few = o; pick = i; } } });
    if (pick < 0) break;
    const [j, k] = opts(pick).sort((x, y) => (opts(x[0]).length - opts(y[0]).length) || (len(x[1]) - len(y[1])))[0];
    bonds[k].order++;
    spare[pick]--; spare[j]--;
  }
  // whatever is left over simply becomes implicit H in RDKit (e.g. the second carboxylate O → OH)
}

function writeV2000(atoms, bonds, title) {
  if (atoms.length > 999) throw new Error('Too many atoms for one building block (max 999).');
  const f = (x, w, d) => x.toFixed(d).padStart(w);
  const i3 = x => String(x).padStart(3);
  const lines = [String(title || '').slice(0, 78), '  Reticublox        3D', '',
    `${i3(atoms.length)}${i3(bonds.length)}  0  0  0  0  0  0  0  0999 V2000`];
  for (const a of atoms) lines.push(`${f(a.p[0], 10, 4)}${f(a.p[1], 10, 4)}${f(a.p[2], 10, 4)} ${a.el.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`);
  for (const b of bonds) lines.push(`${i3(b.a + 1)}${i3(b.b + 1)}${i3(b.order || 1)}  0`);
  lines.push('M  END');
  return lines.join('\n');
}
