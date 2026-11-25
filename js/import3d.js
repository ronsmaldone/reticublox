// Reticublox — read 3D structure files from other programs (CrystalMaker, Mercury, Avogadro, VESTA,
// Materials Studio…) into a V2000 molblock that RDKit can read, keeping the real 3D coordinates.
//   .mol / .sdf   bonds and orders from the file (first record)
//   .mol2         bonds and orders from the file (ar = aromatic)
//   .xyz / .pdb   bonds from interatomic distances, bond orders from the geometry
//   .cif          symmetry-expanded, molecules rebuilt across the cell edges, bonds as for .xyz.
//                 Metals are left out, so a MOF CIF falls apart into its organic linkers.
import { v3, NETS, analyseNet, describeGeometry } from './nets.js';

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
  const { cellAtoms, M, notes, tags } = parseCIFCell(text);
  // MOF: drop the metals first, so each organic linker is its own finite molecule
  if (cellAtoms.some(a => !isOrganic(a.el))) return readCIFOrganic(cellAtoms, M, notes, tags);
  // bonds across the cell edges (minimum image); unwrap each molecule into one piece
  const n = cellAtoms.length;
  if (n > 30000) throw new Error('This CIF has too many atoms in the cell to rebuild molecules here.');
  const nb = cellNeighbours(cellAtoms, M);
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
/** Bonds between the atoms of one cell, each with the cell shift of the partner (minimum image). */
function cellNeighbours(cellAtoms, M) {
  const n = cellAtoms.length, nb = cellAtoms.map(() => []);
  const grid = cellGrid(M, 3.6); // longest bond considered: Cs–O ≈ 3.45 Å
  cellAtoms.forEach((a, i) => grid.add(a.f, i));
  for (let i = 0; i < n; i++) for (const j of grid.near(cellAtoms[i].f)) {
    if (j <= i) continue;
    const [d, shift] = minImage(M, cellAtoms[i].f, cellAtoms[j].f);
    if (bonded(cellAtoms[i].el, cellAtoms[j].el, d)) { nb[i].push([j, shift]); nb[j].push([i, shift.map(s => -s)]); }
  }
  return nb;
}
/**
 * Periodic cell list: fractional positions binned so that every bin is at least `cut` Å across;
 * near(f) lists the atoms in f's bin and the adjacent ones (all of them along short cell edges).
 */
function cellGrid(M, cut) {
  const [a, b, c] = M, V = Math.abs(v3.dot(a, v3.cross(b, c)));
  const width = [V / v3.len(v3.cross(b, c)), V / v3.len(v3.cross(c, a)), V / v3.len(v3.cross(a, b))];
  const nb = width.map(w => Math.max(1, Math.floor(w / cut)));
  const bins = new Map();
  const binOf = f => f.map((x, d) => Math.min(nb[d] - 1, Math.max(0, Math.floor((x - Math.floor(x)) * nb[d]))));
  const key = ([i, j, k]) => i + nb[0] * (j + nb[1] * k);
  const steps = d => (nb[d] >= 3 ? [-1, 0, 1] : [...Array(nb[d]).keys()]);
  return {
    add(f, i) { const k = key(binOf(f)); if (!bins.has(k)) bins.set(k, []); bins.get(k).push(i); },
    near(f) {
      const b0 = binOf(f), keys = new Set(), out = [];
      for (const di of steps(0)) for (const dj of steps(1)) for (const dk of steps(2)) {
        const bb = [di, dj, dk].map((s, d) => (nb[d] >= 3 ? ((b0[d] + s) % nb[d] + nb[d]) % nb[d] : s));
        const k = key(bb);
        if (keys.has(k)) continue;
        keys.add(k);
        const list = bins.get(k);
        if (list) out.push(...list);
      }
      return out;
    },
  };
}
/** CIF text → the unique atoms of one cell (fractional, wrapped), the cell matrix, notes and tags. */
function parseCIFCell(text) {
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
  const grid = cellGrid(M, 1.0); // finds symmetry copies landing on an atom already there
  for (const r of at.rows) {
    if (co >= 0 && num(r[co]) < 0.5) { lowOcc++; continue; } // minor disorder component
    const el = elOf(ct >= 0 && r[ct] !== '?' && r[ct] !== '.' ? r[ct] : r[cl]);
    const f0 = [num(r[cx]), num(r[cy]), num(r[cz])];
    for (const op of ops) {
      const f = op(f0).map(x => x - Math.floor(x + 1e-9));
      if (grid.near(f).some(i => cellAtoms[i].el === el && cartDist(M, cellAtoms[i].f, f) < 0.4)) continue;
      grid.add(f, cellAtoms.length);
      cellAtoms.push({ el, f });
    }
  }
  if (lowOcc) notes.push(`${lowOcc} disordered site${lowOcc > 1 ? 's' : ''} with occupancy < 0.5 skipped.`);
  return { cellAtoms, M, notes, tags };
}
// MOF CIF: drop the metals first, then each organic linker is its own finite molecule
function readCIFOrganic(cellAtoms, M, notes, tags) {
  const metals = cellAtoms.filter(a => !isOrganic(a.el));
  const org = cellAtoms.filter(a => isOrganic(a.el));
  notes.push(`Metal atoms (${[...new Set(metals.map(a => a.el))].join(', ')}) left out; the largest organic linker was taken.`);
  const n = org.length, nb = cellNeighbours(org, M);
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

const INFINITE = 'This CIF is an endless covalent network (a COF), so there is no single molecule to cut out. Switch to “Start from a CIF” at the top of the builder: it splits the framework into its nodes and linkers and finds the net. Or export one building block as .mol or .xyz.';

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
// metal radii (Å, Cordero et al. 2008), for metal–ligand bonds in MOF frameworks
const MCOV = { Li: 1.28, Na: 1.66, K: 2.03, Rb: 2.20, Cs: 2.44, Mg: 1.41, Ca: 1.76, Sr: 1.95, Ba: 2.15, Al: 1.21, Ga: 1.22, In: 1.42,
  Sc: 1.70, Y: 1.90, Ti: 1.60, Zr: 1.75, Hf: 1.75, V: 1.53, Cr: 1.39, Mo: 1.54, W: 1.62, Mn: 1.39, Fe: 1.32, Ru: 1.46, Co: 1.26, Rh: 1.42,
  Ni: 1.24, Pd: 1.39, Pt: 1.36, Cu: 1.32, Ag: 1.45, Au: 1.36, Zn: 1.22, Cd: 1.44, Hg: 1.32, Sn: 1.39, Pb: 1.46, Bi: 1.48,
  La: 2.07, Ce: 2.04, Eu: 1.98, Gd: 1.96, Tb: 1.94, U: 1.96, Th: 2.06 };
const DONOR = new Set(['O', 'N', 'S', 'F', 'Cl', 'Br', 'I', 'P']);
const bonded = (e1, e2, d) => {
  if (d < 0.4) return false;
  const o1 = isOrganic(e1), o2 = isOrganic(e2);
  if (o1 && o2) return !(e1 === 'H' && e2 === 'H') && d < COV[e1] + COV[e2] + 0.4;
  if (!o1 && !o2) return !!(MCOV[e1] && MCOV[e2]) && d < MCOV[e1] + MCOV[e2] + 0.1; // Cu–Cu of a paddlewheel, Zr₆ edges
  const [m, x] = o1 ? [e2, e1] : [e1, e2];
  return !!MCOV[m] && DONOR.has(x) && d < MCOV[m] + COV[x] + 0.35;
};
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

// ================================================================ whole frameworks (COF / MOF CIFs)
/**
 * A COF or MOF crystal structure → its building blocks and net, keeping the computed conformation.
 *   1. rebuild the periodic bond network of one cell (bonds across the cell edges included, metal–ligand bonds too);
 *   2. cut the linkage bonds: COFs at imine C=N, hydrazone, azine, β-ketoenamine C–N; MOFs at each
 *      carboxylate's C–C bond (the carboxylate stays with the metal cluster) and at metal–N bonds;
 *   3. each piece is one copy of a monomer: ≥3 cuts → node, 2 cuts → linker; copies are grouped
 *      into types by formula; the most typical copy of each type is kept, exactly as computed;
 *   4. nodes + the linkers (or direct bonds) between them form the net, matched to the library
 *      by dimension and coordination sequences.
 * Each block comes back as a graph (and, if organic, a molblock) with a * atom where its partner across each
 * cut sat, so its connection points, and the direction of every magnet fitting, are the crystal's own.
 * → { name, cls, dim, net, mode: 'linker'|'direct', types: [{ role, label?, cn, perCell, formula, metal, graph, molblock, spread }], notes }
 */
export function importFramework(text, filename = '') {
  const base = filename.replace(/\.[^.]+$/, '') || 'Framework';
  if (!/\.cif$/i.test(filename) && !/_cell_length_a/.test(text))
    throw new Error('A whole framework needs a CIF (it holds the unit cell). Export the crystal structure as .cif.');
  const { cellAtoms, M, notes } = parseCIFCell(text);
  const isMetal = i => !isOrganic(cellAtoms[i].el);
  const metals = [...new Set(cellAtoms.filter(a => !isOrganic(a.el)).map(a => a.el))];
  const unknown = metals.filter(m => !MCOV[m]);
  if (unknown.length) throw new Error(`Unknown element(s) ${unknown.join(', ')} in this CIF.`);
  const n = cellAtoms.length;
  if (n > 30000) throw new Error(`This cell has ${n} atoms; the framework importer handles up to 30,000.`);
  const el = i => cellAtoms[i].el;
  const nb = cellNeighbours(cellAtoms, M);
  const cartVec = (i, j, sh) => frac2cart(M, v3.add(v3.sub(cellAtoms[j].f, cellAtoms[i].f), sh));
  const heavy = i => nb[i].filter(([j]) => el(j) !== 'H');
  // coordinated solvent (an O on one metal and nothing else heavy, e.g. the axial water of a
  // paddlewheel) is not part of the framework
  const solvent = new Set();
  for (let i = 0; i < n; i++) if (el(i) === 'O') {
    const h = heavy(i);
    if (h.length === 1 && isMetal(h[0][0])) { solvent.add(i); for (const [j] of nb[i]) if (el(j) === 'H') solvent.add(j); }
  }
  if (solvent.size) notes.push('Coordinated solvent (terminal O / water on the metals) left out.');

  // ---- 2. linkage bonds
  const inSmallRing = (i, j, jsh) => { // another path i → j (same image of j) of ≤ 7 bonds?
    const seen = new Set([`${i}|0,0,0`]), q = [[i, [0, 0, 0], 0]];
    for (let k = 0; k < q.length; k++) {
      const [u, o, d] = q[k];
      if (d >= 7) continue;
      for (const [w, sh] of nb[u]) {
        const ow = v3.add(o, sh);
        const atJ = w === j && ow.every((x, m) => x === jsh[m]);
        if (atJ && u === i) continue; // the bond itself
        if (atJ) return true;
        const key = `${w}|${ow}`;
        if (!seen.has(key)) { seen.add(key); q.push([w, ow, d + 1]); }
      }
    }
    return false;
  };
  const cuts = [];
  for (let i = 0; i < n; i++) {
    if (el(i) !== 'C') continue;
    for (const [j, sh] of nb[i]) {
      if (el(j) !== 'N') continue;
      const hc = heavy(i), hn = heavy(j);
      // C: aldehyde-derived (aryl/vinyl C + N, plus H); N: two heavy neighbours (aryl C or N of a hydrazone/azine)
      if (hc.length !== 2 || hn.length !== 2) continue;
      const d = v3.len(cartVec(i, j, sh));
      if (d > 1.38) continue;
      if (inSmallRing(i, j, sh)) continue;
      cuts.push({ a: i, b: j, sh });
    }
  }
  const carboxC = new Set();
  // MOFs: a carboxylate's C–C bond (the carboxylate C is the point of extension and goes with the
  // metal cluster, as in reticular chemistry), and metal–N bonds of pyridyl / azolate linkers
  for (let i = 0; i < n; i++) {
    if (el(i) !== 'C') continue;
    const h = heavy(i), os = h.filter(([j]) => el(j) === 'O'), cs = h.filter(([j]) => el(j) === 'C');
    if (h.length !== 3 || os.length !== 2 || cs.length !== 1) continue;
    // a real carboxylate: both O's bonded only to this C and metals (not an acetal C of a sugar, an ester…)
    if (!os.every(([o]) => heavy(o).every(([x]) => x === i || isMetal(x)))) continue;
    if (!os.some(([o]) => heavy(o).some(([m]) => isMetal(m)))) continue;
    carboxC.add(i);
    const [j, sh] = cs[0];
    if (!inSmallRing(i, j, sh)) cuts.push({ a: i, b: j, sh });
  }
  for (let i = 0; i < n; i++) {
    if (!isMetal(i)) continue;
    for (const [j, sh] of nb[i]) {
      if (el(j) === 'N' && !inSmallRing(i, j, sh)) cuts.push({ a: i, b: j, sh });
      // metal–O of a hydroxyl, ether or alkoxide (e.g. K⁺ on cyclodextrin in CD-MOF), chelating or not;
      // carboxylate O's stay with the metal cluster (their C–C bond is the cut)
      if (el(j) === 'O' && heavy(j).some(([c]) => el(c) === 'C') && !heavy(j).some(([c]) => carboxC.has(c))) cuts.push({ a: i, b: j, sh });
    }
  }
  if (!cuts.length) throw new Error(metals.length
    ? 'No carboxylate or metal–N linkages found to cut this MOF into building blocks. Phosphonate, sulfonate and catecholate MOFs are not supported yet.'
    : 'No imine-type linkages (C=N, C–N of a β-ketoenamine) found to cut this framework into building blocks. Boronate-ester, boroxine and triazine COFs are not supported yet.');
  const cutKey = new Set(cuts.flatMap(c => [`${c.a}>${c.b}|${c.sh}`, `${c.b}>${c.a}|${c.sh.map(s => -s)}`]));
  const isCut = (i, j, sh) => cutKey.has(`${i}>${j}|${sh}`);

  // ---- 3. pieces (each a finite molecule, unwrapped across the cell edges)
  const frag = new Array(n).fill(-1), off = new Array(n).fill(null), frags = [];
  for (let s = 0; s < n; s++) {
    if (frag[s] >= 0 || solvent.has(s)) continue;
    const id = frags.length, atoms = [s];
    frag[s] = id; off[s] = [0, 0, 0];
    for (let k = 0; k < atoms.length; k++) {
      const u = atoms[k];
      for (const [w, sh] of nb[u]) {
        if (isCut(u, w, sh) || solvent.has(w)) continue;
        const o = v3.add(off[u], sh);
        if (frag[w] < 0) { frag[w] = id; off[w] = o; atoms.push(w); }
        else if (v3.len(v3.sub(off[w], o)) > 0.5) throw new Error('After cutting the linkages a piece is still an endless network, so this framework has a linkage Reticublox does not recognise yet.');
      }
    }
    frags.push({ id, atoms, links: [] });
  }
  // links: piece F (copy in cell 0) → piece G in cell `shift`, through cut atom a (in F) – b (in G)
  // several cut bonds between the same two pieces (a K⁺ held by two O's of one cyclodextrin) are one connection
  for (const c of cuts) for (const [x, y, sh] of [[c.a, c.b, c.sh], [c.b, c.a, c.sh.map(s => -s)]]) {
    const F = frags[frag[x]], to = frag[y], shift = v3.sub(v3.add(sh, off[x]), off[y]);
    const pair = { atom: x, partner: y, partnerShift: v3.add(sh, off[x]) };
    const link = F.links.find(l => l.to === to && v3.len(v3.sub(l.shift, shift)) < 0.5);
    if (link) link.pairs.push(pair); else F.links.push({ to, shift, pairs: [pair] });
  }
  const strays = frags.filter(f => f.links.length < 2);
  if (strays.some(f => f.links.length === 1)) notes.push(`${strays.filter(f => f.links.length === 1).length} piece(s) with a single linkage (end caps or defects) left out.`);
  if (strays.some(f => !f.links.length && f.atoms.some(i => el(i) !== 'H'))) notes.push('Unbonded guests (solvent or loose atoms) left out.');
  const real = frags.filter(f => f.links.length >= 2);

  // types: same formula and same number of linkages
  const formula = f => {
    const c = {};
    for (const i of f.atoms) c[el(i)] = (c[el(i)] || 0) + 1;
    return ['C', 'H', ...Object.keys(c).filter(e => e !== 'C' && e !== 'H').sort()].filter(e => c[e]).map(e => e + (c[e] > 1 ? c[e] : '')).join('');
  };
  const typeOf = new Map(), types = [];
  for (const f of real) {
    const key = `${formula(f)}|${f.links.length}`;
    if (!typeOf.has(key)) { typeOf.set(key, types.length); types.push({ key, formula: formula(f), cn: f.links.length, copies: [] }); }
    f.type = typeOf.get(key);
    types[f.type].copies.push(f);
  }

  // ---- 4. the net: nodes are the pieces with ≥3 linkages; 2-c pieces become edges
  const nodes = real.filter(f => f.links.length >= 3);
  if (!nodes.length) throw new Error('No branching building block (3 or more linkages) found: only chains, which Reticublox cannot assemble.');
  const vIdx = new Map(nodes.map((f, k) => [f.id, k]));
  const byId = new Map(real.map(f => [f.id, f]));
  const edges = [];
  let viaLinker = 0, direct = 0;
  for (const f of nodes) for (const l of f.links) {
    let g = byId.get(l.to), shift = l.shift;
    if (!g) continue;
    if (g.links.length === 2) {
      // walk through the linker to the node on its far side
      const back = g.links.findIndex(m => m.to === f.id && v3.len(v3.add(m.shift, shift)) < 0.5);
      const far = g.links[back === 0 ? 1 : 0];
      const h = byId.get(far.to);
      if (!h || h.links.length === 2) throw new Error('Two linkers bonded end to end: this framework’s linkers are not supported yet.');
      shift = v3.add(shift, far.shift);
      g = h; viaLinker++;
    } else direct++;
    edges.push({ i: vIdx.get(f.id), j: vIdx.get(g.id), shift });
  }
  if (viaLinker && direct) throw new Error('This framework mixes node–linker–node and direct node–node linkages, which Reticublox cannot model yet.');
  const mode = viaLinker ? 'linker' : 'direct';
  const adj = nodes.map(() => []);
  for (const e of edges) adj[e.i].push({ j: e.j, off: e.shift });

  // separate nets in the cell (stacked layers, interpenetration): use the one with the first node
  const comp = new Array(nodes.length).fill(-1);
  let nComp = 0;
  for (let s = 0; s < nodes.length; s++) {
    if (comp[s] >= 0) continue;
    const q = [s]; comp[s] = nComp;
    for (let k = 0; k < q.length; k++) for (const { j } of adj[q[k]]) if (comp[j] < 0) { comp[j] = nComp; q.push(j); }
    nComp++;
  }
  const dim = netDimension(adj, comp, 0);
  if (dim === 0) throw new Error('The building blocks found here do not link up into an endless framework (a discrete cage or cluster?), so there is no net to build.');
  if (nComp > 1) notes.push(`The cell holds ${nComp} separate ${dim === 2 ? 'layers' : 'nets (interpenetrated)'}; they are the same framework, so one is modelled.`);

  // match the net: coordination sequences of each node type
  const cs = v => coordSeq(adj, v, 7).join(',');
  const nodeTypes = [...new Set(nodes.map(f => f.type))];
  const csOfType = new Map(nodeTypes.map(t => [t, cs(nodes.findIndex(f => f.type === t && comp[vIdx.get(f.id)] === 0) )]));
  let match = null;
  for (const [key, def] of Object.entries(NETS)) {
    if (def.dim !== dim) continue;
    const an = analyseNet(key);
    const netCS = {};
    an.net.verts.forEach((v, i) => { netCS[v.t] ??= coordSeq(an.nbrs, i, 7).join(','); });
    const map = {};
    for (const t of nodeTypes) {
      const lab = Object.keys(netCS).find(l => netCS[l] === csOfType.get(t) && !Object.values(map).includes(l))
        ?? Object.keys(netCS).find(l => netCS[l] === csOfType.get(t));
      if (!lab) { map.fail = true; break; }
      map[t] = lab;
    }
    if (map.fail || Object.keys(netCS).some(l => !Object.values(map).includes(l))) continue;
    if (mode === 'direct' && !an.bipartite) continue;
    match = { key, map };
    break;
  }
  // not in the library: use the crystal's own net (block centres and links as they sit in the cell)
  const crystal = !match;
  if (crystal) match = { key: 'cif', map: Object.fromEntries(nodeTypes.map((t, k) => [t, String.fromCharCode(65 + k)])) };
  if (!crystal && nodeTypes.length > new Set(Object.values(match.map)).size)
    notes.push(`Two different node molecules sit on equivalent sites of ${match.key}; the first is used for both.`);

  // ---- blocks: the most typical copy of each type, its partners across the cuts as * atoms
  const out = [];
  const usedTypes = mode === 'linker' ? [...nodeTypes, ...new Set(real.filter(f => f.links.length === 2).map(f => f.type))] : nodeTypes;
  const seenLabel = new Set();
  const X = crystal ? crystalGeometry() : null;
  for (const t of usedTypes) {
    const T = types[t];
    const label = T.cn === 2 ? null : match.map[t];
    if (label && seenLabel.has(label)) continue;
    if (label) seenLabel.add(label);
    // a crystal net is built from the first net's nodes, so print one of those
    const pool = crystal && T.cn !== 2 ? T.copies.filter(f => comp[vIdx.get(f.id)] === 0) : T.copies;
    const { copy, spread } = typicalCopy(pool, f => f.atoms.filter(i => el(i) !== 'H').map(i => frac2cart(M, v3.add(cellAtoms[i].f, off[i]))));
    out.push({
      ...(crystal ? { fixed: T.cn === 2 ? X.linkerArms(copy) : X.nodeArms(copy) } : {}),
      role: T.cn === 2 ? 'linker' : 'node', label, cn: T.cn, formula: T.formula,
      perCell: T.copies.length / nComp, spread,
      metal: copy.atoms.some(isMetal),
      ...piece(copy, cellAtoms, M, off, el, `${base} ${T.cn === 2 ? 'linker' : 'node'}`),
    });
  }
  const linkers = out.filter(o => o.role === 'linker');
  if (mode === 'linker' && linkers.length > 1) notes.push(`${linkers.length} different linker molecules found; the first (${linkers[0].formula}) is used.`);
  if (crystal) notes.push(`This net (${nodeTypes.map(t => types[t].cn).join(',')}-connected) is not in the library, so the blocks are placed as in the crystal itself.`);
  return { name: base, cls: metals.length ? 'mof' : 'cof', dim, net: match.key, mode, types: out, notes, ...(crystal ? { crystalNet: X.net(match.map) } : {}) };

  /**
   * The crystal's own net. Each node sits at the centre of its connection points; each link runs from node
   * centre to node centre, and the two magnet faces meet on that line, level with the middle of the cut bond.
   * Arms are given per block as { centre, dirs, t (Å from the centre to the magnet face) } in link order,
   * all in the crystal's frame, so a uniform scale-up of the crystal puts every magnet pair face to face.
   */
  function crystalGeometry() {
    const at = (i, o) => frac2cart(M, v3.add(cellAtoms[i].f, o));
    const mean = P => v3.mul(P.reduce((a, p) => v3.add(a, p), [0, 0, 0]), 1 / P.length);
    const own = l => mean(l.pairs.map(q => at(q.atom, off[q.atom])));
    const cutMid = l => v3.mul(v3.add(own(l), mean(l.pairs.map(q => at(q.partner, q.partnerShift)))), 0.5);
    const centre = f => mean(f.links.map(own));
    const T = sh => frac2cart(M, sh);
    // for a node's link: the node at the far end (through a linker if there is one) and its cell shift
    const far = (f, l) => {
      const g = byId.get(l.to);
      if (!g || g.links.length !== 2) return { h: g, shift: l.shift };
      const back = g.links.findIndex(m => m.to === f.id && v3.len(v3.add(m.shift, l.shift)) < 0.5);
      const fl = g.links[back === 0 ? 1 : 0];
      return { h: byId.get(fl.to), shift: v3.add(l.shift, fl.shift) };
    };
    // distance from node A's centre to where its magnet face sits on the line to B
    const tOf = (A, B, l) => {
      const L = v3.len(v3.sub(B, A)), u = v3.norm(v3.sub(B, A));
      return Math.min(0.95 * L, Math.max(0.05 * L, v3.dot(v3.sub(cutMid(l), A), u)));
    };
    const nodeArms = f => {
      const A = centre(f);
      const arms = f.links.map(l => { const { h, shift } = far(f, l); const B = v3.add(centre(h), T(shift)); return { dir: v3.norm(v3.sub(B, A)), t: tOf(A, B, l) }; });
      return { centre: A, dirs: arms.map(a => a.dir), t: arms.map(a => a.t) };
    };
    const linkerArms = g => {
      // its two nodes, seen from this linker copy
      const ends = g.links.map(l => {
        const f = byId.get(l.to), back = f.links.find(m => m.to === g.id && v3.len(v3.add(m.shift, l.shift)) < 0.5);
        return { P: v3.add(centre(f), T(l.shift)), f, l: back, shift: l.shift };
      });
      const [E0, E1] = ends, mid = v3.mul(v3.add(E0.P, E1.P), 0.5), half = v3.len(v3.sub(E1.P, E0.P)) / 2;
      // the node's own link to this linker, moved into this linker's frame
      const tNode = (E, other) => tOf(v3.sub(E.P, T(E.shift)), v3.sub(other.P, T(E.shift)), E.l);
      return { centre: mid, dirs: [v3.norm(v3.sub(E0.P, mid)), v3.norm(v3.sub(E1.P, mid))], t: [half - tNode(E0, E1), half - tNode(E1, E0)] };
    };
    const net = labelOf => {
      const vs = nodes.filter(f => comp[vIdx.get(f.id)] === 0);
      const vi = new Map(vs.map((f, k) => [f.id, k]));
      const cart = vs.map(centre);
      const nbrs = vs.map((f, k) => f.links.map(l => {
        const { h, shift } = far(f, l);
        return { j: vi.get(h.id), off: shift, vec: v3.sub(v3.add(cart[vi.get(h.id)], T(shift)), cart[k]) };
      }));
      const types = {};
      vs.forEach((f, k) => {
        const lab = labelOf[f.type];
        if (!types[lab]) { const dirs = nbrs[k].map(n => v3.norm(n.vec)); types[lab] = { label: lab, cn: dirs.length, dirs, count: 0, geometry: describeGeometry(dirs) }; }
        types[lab].count++;
      });
      const colour = vs.map(() => -1);
      let bipartite = true;
      for (let s0 = 0; s0 < vs.length; s0++) {
        if (colour[s0] !== -1) continue;
        colour[s0] = 0;
        const q = [s0];
        for (let k = 0; k < q.length; k++) for (const n of nbrs[q[k]]) {
          if (colour[n.j] === -1) { colour[n.j] = 1 - colour[q[k]]; q.push(n.j); } else if (colour[n.j] === colour[q[k]]) bipartite = false;
        }
      }
      return {
        key: 'cif', crystal: true, dim, cell: M, cart, nbrs, types, bipartite, colour,
        net: { verts: vs.map(f => ({ t: labelOf[f.type] })), dim },
        vertsPerCell: vs.length, edgesPerCell: nbrs.reduce((a, x) => a + x.length, 0) / 2,
      };
    };
    return { nodeArms, linkerArms, net };
  }
}

/**
 * One piece as computed: { graph: { atoms [{ el, p }], bonds [{ a, b, order }] } without H, partners
 * across the cuts as * atoms; molblock (organic pieces only, with H) for RDKit's SMILES }.
 */
function piece(f, cellAtoms, M, off, el, title) {
  const atoms = f.atoms.map(i => ({ el: el(i), p: frac2cart(M, v3.add(cellAtoms[i].f, off[i])) }));
  const at = new Map(f.atoms.map((i, k) => [i, k]));
  const bonds = [];
  for (let x = 0; x < atoms.length; x++) for (let y = x + 1; y < atoms.length; y++) {
    if (bonded(atoms[x].el, atoms[y].el, v3.len(v3.sub(atoms[y].p, atoms[x].p)))) bonds.push({ a: x, b: y, order: 0 });
  }
  // partners keep their element while bond orders are worked out (the C=N stays double), then become *
  const partners = [];
  for (const l of f.links) for (const q of l.pairs) {
    const k = atoms.length;
    atoms.push({ el: el(q.partner), p: frac2cart(M, v3.add(cellAtoms[q.partner].f, q.partnerShift)), link: l });
    bonds.push({ a: at.get(q.atom), b: k, order: 0 });
    partners.push(k);
  }
  // bond orders from the organic part only (a carboxylate O bonded to a metal is still C=O / C–O⁻)
  const org = b => isOrganic(atoms[b.a].el) && isOrganic(atoms[b.b].el);
  assignOrders(atoms, bonds.filter(org));
  for (const b of bonds) if (!org(b)) b.order = 1;
  // one * per connection, at the centroid of its partner atoms, bonded to every atom of this piece in it
  for (const l of f.links) {
    const ks = partners.filter(k => atoms[k].link === l);
    const k0 = ks[0];
    atoms[k0].p = v3.mul(ks.reduce((s, k) => v3.add(s, atoms[k].p), [0, 0, 0]), 1 / ks.length);
    for (const k of ks.slice(1)) { atoms[k].gone = true; for (const b of bonds) { if (b.b === k) b.b = k0; if (b.a === k) b.a = k0; } }
    atoms[k0].el = '*';
  }
  // the same own atom bonded to two partners of one connection: keep one bond
  const seenB = new Set();
  for (const b of bonds) { const key = Math.min(b.a, b.b) + '-' + Math.max(b.a, b.b); if (seenB.has(key)) b.dup = true; seenB.add(key); }
  const metal = atoms.some(a => a.el !== '*' && !isOrganic(a.el));
  const keep = atoms.map(a => a.el !== 'H' && !a.gone), idx = [];
  let k = 0;
  atoms.forEach((a, i) => { idx[i] = keep[i] ? k++ : -1; });
  return {
    graph: { atoms: atoms.filter((_, i) => keep[i]).map(a => ({ el: a.el, p: a.p })), bonds: bonds.filter(b => keep[b.a] && keep[b.b] && !b.dup).map(b => ({ a: idx[b.a], b: idx[b.b], order: b.order || 1 })) },
    molblock: metal ? null : molblockOf(atoms, bonds, title),
  };
}

/** V2000 of a piece without the merged-away partner atoms. */
function molblockOf(atoms, bonds, title) {
  const keep = atoms.map(a => !a.gone), idx = [];
  let k = 0;
  atoms.forEach((a, i) => { idx[i] = keep[i] ? k++ : -1; });
  return writeV2000(atoms.filter((_, i) => keep[i]), bonds.filter(b => keep[b.a] && keep[b.b] && !b.dup).map(b => ({ ...b, a: idx[b.a], b: idx[b.b] })), title);
}

/** The copy closest to all the others (shape fingerprint: sorted distances from the centroid). */
function typicalCopy(copies, coords) {
  const fp = copies.map(f => {
    const P = coords(f), c = v3.mul(P.reduce((s, p) => v3.add(s, p), [0, 0, 0]), 1 / P.length);
    return P.map(p => v3.len(v3.sub(p, c))).sort((a, b) => a - b);
  });
  const rms = (a, b) => Math.sqrt(a.reduce((s, x, k) => s + (x - b[k]) ** 2, 0) / a.length);
  let best = 0, bestSum = Infinity, spread = 0;
  fp.forEach((a, i) => {
    let s = 0;
    fp.forEach((b, j) => { if (i !== j) { const r = rms(a, b); s += r; spread = Math.max(spread, r); } });
    if (s < bestSum) { bestSum = s; best = i; }
  });
  return { copy: copies[best], spread };
}

/** Shell sizes of the coordination sequence of vertex v (adjacency with cell shifts). */
function coordSeq(adj, v, shells) {
  const key = (j, o) => `${j}|${o.join(',')}`;
  const seen = new Set([key(v, [0, 0, 0])]);
  let front = [[v, [0, 0, 0]]];
  const out = [];
  for (let s = 0; s < shells; s++) {
    const next = [];
    for (const [u, o] of front) for (const { j, off } of adj[u]) {
      const w = v3.add(o, off), k = key(j, w);
      if (!seen.has(k)) { seen.add(k); next.push([j, w]); }
    }
    out.push(next.length);
    front = next;
  }
  return out;
}

/** Periodic dimension of one connected net: rank of its cycle translation vectors. */
function netDimension(adj, comp, c) {
  const pot = new Array(adj.length).fill(null);
  const s = comp.indexOf(c);
  pot[s] = [0, 0, 0];
  const q = [s], vecs = [];
  for (let k = 0; k < q.length; k++) {
    const u = q[k];
    for (const { j, off } of adj[u]) {
      const p = v3.add(pot[u], off);
      if (!pot[j]) { pot[j] = p; q.push(j); } else { const d = v3.sub(p, pot[j]); if (v3.len(d) > 0.5) vecs.push(d); }
    }
  }
  // rank by Gaussian elimination
  const rows = vecs.map(v => [...v]);
  let rank = 0;
  for (let col = 0; col < 3 && rank < rows.length; col++) {
    const piv = rows.findIndex((r, i) => i >= rank && Math.abs(r[col]) > 1e-6);
    if (piv < 0) continue;
    [rows[rank], rows[piv]] = [rows[piv], rows[rank]];
    for (let i = 0; i < rows.length; i++) if (i !== rank) {
      const f = rows[i][col] / rows[rank][col];
      rows[i] = rows[i].map((x, m) => x - f * rows[rank][m]);
    }
    rank++;
  }
  return rank;
}
