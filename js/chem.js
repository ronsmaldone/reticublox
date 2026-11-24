// Reticublox — chemistry → building block.
// Parses a molecule with RDKit.js, finds the points of extension (connection atoms),
// finds rings, and turns the molecule into a node / linker spec:
//   hub at the core centre, one arm per connection point, a bulge for each ring on an arm.
import { v3, GEOMETRIES } from './nets.js';
import { spacedBulges, atomRadius } from './geometry.js';
import { alignDirs, apply, rotBetween, IDENT } from './rot.js';

let RDKit = null;
export async function initRDKit() {
  if (RDKit) return RDKit;
  RDKit = await window.initRDKitModule({ locateFile: () => 'vendor/RDKit_minimal.wasm' });
  return RDKit;
}

// Linking groups. `atom` = index (within the SMARTS) of the point of extension.
// `ext` = default extension in Å beyond that atom (≈ half of the new bond formed on assembly).
export const GROUPS = {
  dummy:   { label: 'Marked * atoms', smarts: null, ext: 0 },
  cho:     { label: 'Aldehyde –CHO (imine COF)', smarts: '[CX3H1](=O)[#6]', atom: 0, ext: 0.64, drop: true },
  nh2:     { label: 'Amine –NH2 (imine COF)', smarts: '[NX3;H2][#6]', atom: 0, ext: 0.64 },
  cooh:    { label: 'Carboxylic acid (MOF; C = point of extension)', smarts: '[CX3](=O)[OX2H1,OX1-]', atom: 0, ext: 0 },
  boronic: { label: 'Boronic acid –B(OH)2 (boroxine/boronate COF)', smarts: '[BX3]([OX2H1])[OX2H1]', atom: 0, ext: 0.7, drop: true },
  pyridyl: { label: 'Pyridyl N (pillar / M–N)', smarts: '[nX2;r6]', atom: 0, ext: 1.0 },
  nitrile: { label: 'Nitrile –CN (triazine CTF)', smarts: '[CX2]#[NX1]', atom: 0, ext: 0.7, drop: true },
  manual:  { label: 'Click atoms in the drawing', smarts: null, ext: 0.5 },
};

export const PRESETS = {
  'Linkers (2-c)': {
    'Terephthalaldehyde (PDA-CHO)': 'O=Cc1ccc(C=O)cc1',
    'p-Phenylenediamine (PDA)': 'Nc1ccc(N)cc1',
    'Benzidine (BZ)': 'Nc1ccc(-c2ccc(N)cc2)cc1',
    '4,4″-Diamino-p-terphenyl': 'Nc1ccc(-c2ccc(-c3ccc(N)cc3)cc2)cc1',
    'Terephthalic acid (BDC)': 'OC(=O)c1ccc(C(O)=O)cc1',
    'Biphenyl-4,4′-dicarboxylic acid (BPDC)': 'OC(=O)c1ccc(-c2ccc(C(O)=O)cc2)cc1',
    'p-Terphenyl dicarboxylic acid (TPDC)': 'OC(=O)c1ccc(-c2ccc(-c3ccc(C(O)=O)cc3)cc2)cc1',
    '1,4-Benzenediboronic acid (BDBA)': 'OB(O)c1ccc(B(O)O)cc1',
    '4,4′-Bipyridine': 'c1cc(ccn1)-c1ccncc1',
    'Isophthalic acid (bent)': 'OC(=O)c1cccc(C(O)=O)c1',
    '2,6-Naphthalenedicarboxylic acid': 'OC(=O)c1ccc2cc(C(O)=O)ccc2c1',
  },
  'Nodes (3-c)': {
    '1,3,5-Triformylbenzene (TFB)': 'O=Cc1cc(C=O)cc(C=O)c1',
    '1,3,5-Tris(4-aminophenyl)benzene (TAPB)': 'Nc1ccc(cc1)-c1cc(cc(c1)-c1ccc(N)cc1)-c1ccc(N)cc1',
    '1,3,5-Tris(4-formylphenyl)benzene (TFPB)': 'O=Cc1ccc(cc1)-c1cc(cc(c1)-c1ccc(C=O)cc1)-c1ccc(C=O)cc1',
    'Tris(4-aminophenyl)amine (TAPA)': 'Nc1ccc(cc1)N(c1ccc(N)cc1)c1ccc(N)cc1',
    'Trimesic acid (BTC)': 'OC(=O)c1cc(cc(c1)C(O)=O)C(O)=O',
    'BTB (benzene-tribenzoic acid)': 'OC(=O)c1ccc(cc1)-c1cc(cc(c1)-c1ccc(cc1)C(O)=O)-c1ccc(cc1)C(O)=O',
    '1,3,5-Tris(4-aminophenyl)triazine (TAPT)': 'Nc1ccc(cc1)-c1nc(nc(n1)-c1ccc(N)cc1)-c1ccc(N)cc1',
  },
  'Nodes (4-c)': {
    'Tetrakis(4-aminophenyl)methane (TAM, tetrahedral)': 'Nc1ccc(cc1)C(c1ccc(N)cc1)(c1ccc(N)cc1)c1ccc(N)cc1',
    'Tetrakis(4-formylphenyl)methane (tetrahedral)': 'O=Cc1ccc(cc1)C(c1ccc(C=O)cc1)(c1ccc(C=O)cc1)c1ccc(C=O)cc1',
    'Tetrakis(4-aminophenyl)ethylene (ETTA)': 'Nc1ccc(cc1)C(c1ccc(N)cc1)=C(c1ccc(N)cc1)c1ccc(N)cc1',
    '1,3,6,8-Tetrakis(4-formylphenyl)pyrene (TFPPy)': 'O=Cc1ccc(cc1)-c1cc(-c2ccc(C=O)cc2)c2ccc3c(cc(-c4ccc(C=O)cc4)c4ccc1c2c34)-c1ccc(C=O)cc1',
    'Tetrakis(4-aminophenyl)porphyrin (TAPP)': 'Nc1ccc(cc1)-c1c2ccc(n2)c(-c2ccc(N)cc2)c2ccc([nH]2)c(-c2ccc(N)cc2)c2ccc(n2)c(-c2ccc(N)cc2)c2ccc1[nH]2',
    'Biphenyl-3,3′,5,5′-tetracarboxylic acid (nbo)': 'OC(=O)c1cc(cc(c1)C(O)=O)-c1cc(cc(c1)C(O)=O)C(O)=O',
  },
  'Nodes (6-c)': {
    'Hexakis(4-formylphenyl)benzene': 'O=Cc1ccc(cc1)-c1c(-c2ccc(C=O)cc2)c(-c2ccc(C=O)cc2)c(-c2ccc(C=O)cc2)c(-c2ccc(C=O)cc2)c1-c1ccc(C=O)cc1',
  },
};
// MOF-specific organic nodes
PRESETS['Nodes (4-c)']['Tetrakis(4-carboxyphenyl)porphyrin (TCPP)'] =
  'OC(=O)c1ccc(cc1)-c1c2ccc(n2)c(-c2ccc(cc2)C(O)=O)c2ccc([nH]2)c(-c2ccc(cc2)C(O)=O)c2ccc(n2)c(-c2ccc(cc2)C(O)=O)c2ccc1[nH]2';
PRESETS['Nodes (4-c)']['Methanetetrabenzoic acid (MTB, tetrahedral)'] =
  'OC(=O)c1ccc(cc1)C(c1ccc(cc1)C(O)=O)(c1ccc(cc1)C(O)=O)c1ccc(cc1)C(O)=O';

/** Flat preset list tagged with connectivity and framework class (carboxylates / pyridyls → MOF). */
export const PRESET_LIST = Object.entries(PRESETS).flatMap(([group, m]) => {
  const cn = +group.match(/(\d+)-c/)[1];
  return Object.entries(m).map(([name, smiles]) => ({
    name, smiles, cn,
    cls: /C\(O\)=O|C\(=O\)O|OC\(=O\)|ccn|ncc/.test(smiles) ? 'mof' : 'cof',
  }));
});

// Inorganic SBUs (no SMILES): ideal geometry + distance from cluster centre to the
// carboxylate carbon (point of extension), Å. Values are approximate crystallographic distances.
// `el` is the metal printed in molecular-model mode (see sbuCore).
export const SBUS = {
  'Cu₂/Zn₂ paddlewheel (4-c square)': { geom: 'square', r: 2.8, el: 'Cu' },
  'Zn₄O(CO₂)₆ — MOF-5 (6-c octahedral)': { geom: 'octahedral', r: 3.2, el: 'Zn' },
  'M₃O(CO₂)₆ trimer — MIL-88 (6-c trigonal prism)': { geom: 'trigonalPrism', r: 3.6, el: 'Fe' },
  'Zr₆O₄(OH)₄(CO₂)₁₂ — UiO-66 (12-c cuboctahedral)': { geom: 'cuboctahedral', r: 4.5, el: 'Zr' },
  'Zr₆ cluster, 8-c (cubic, e.g. NU-1000 type)': { geom: 'cubic', r: 4.5, el: 'Zr' },
  'Tetrahedral metal / C / Si centre (4-c)': { geom: 'tetrahedral', r: 1.6, el: 'Si' },
};

/**
 * Metal core of each SBU in a fixed frame (Å), with the ideal arm (carboxylate) directions in that
 * same frame. The carboxylate C and O atoms belong to the linker, so only metals (+ μ-O) are printed;
 * each arm's fitting forks from the metals that carboxylate bridges.
 */
function sbuCore(geom, el) {
  const n = v3.norm, s2 = Math.SQRT1_2;
  const az = (deg, z = 0, r = 1) => { const a = deg * Math.PI / 180; return [r * Math.cos(a), r * Math.sin(a), z]; };
  const zr6 = () => [[2.5, 0, 0], [-2.5, 0, 0], [0, 2.5, 0], [0, -2.5, 0], [0, 0, 2.5], [0, 0, -2.5]];
  const edges = pts => pts.flatMap((p, i) => pts.map((q, j) => [i, j]).filter(([a, b]) => b > a && v3.len(v3.sub(pts[a], pts[b])) < 3.6));
  switch (geom) {
    case 'square': // Cu–Cu 2.6 Å, perpendicular to the plane of the four carboxylates
      return { atoms: [[el, [0, 0, 1.3]], [el, [0, 0, -1.3]]], bonds: [[0, 1]], dirs: [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]] };
    case 'octahedral': { // Zn₄ tetrahedron around μ4-O, Zn–O 1.95 Å; carboxylates bridge Zn pairs along ±x, ±y, ±z
      const t = 1.95 / Math.sqrt(3);
      const zn = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].map(p => v3.mul(p, t));
      return { atoms: [['O', [0, 0, 0]], ...zn.map(p => [el, p])], bonds: [1, 2, 3, 4].map(i => [0, i]), dirs: GEOMETRIES.octahedral };
    }
    case 'trigonalPrism': { // M₃ triangle around μ3-O, M–O 1.9 Å; each M pair bridged above and below the plane
      const m = [60, 180, 300].map(a => az(a, 0, 1.9));
      return { atoms: [['O', [0, 0, 0]], ...m.map(p => [el, p])], bonds: [[0, 1], [0, 2], [0, 3]], dirs: [0, 120, 240].flatMap(a => [n(az(a, 1)), n(az(a, -1))]) };
    }
    case 'cuboctahedral': case 'cubic': { // Zr₆ octahedron, Zr–Zr 3.5 Å
      const zr = zr6();
      const dirs = geom === 'cubic' ? [-1, 1].flatMap(a => [-1, 1].flatMap(b => [-1, 1].map(c => n([a, b, c]))))
        : [[1, 1, 0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0], [1, 0, 1], [1, 0, -1], [-1, 0, 1], [-1, 0, -1], [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1]].map(v => v3.mul(v, s2));
      return { atoms: zr.map(p => [el, p]), bonds: edges(zr), dirs };
    }
    default: // single tetrahedral centre: fittings straight from it
      return { atoms: [[el, [0, 0, 0]]], bonds: [], dirs: GEOMETRIES[geom] };
  }
}

/**
 * A 3D molblock from a file → { molblock (heavy atoms, real Å coordinates) | null if flat, smiles }. RDKit checks
 * the chemistry; hydrogens are dropped (they are never printed, and RDKit keeps them as implicit H).
 */
export function prepare3D(molblock) {
  const mol = RDKit.get_mol(molblock, JSON.stringify({ removeHs: true }));
  if (!mol || !mol.is_valid()) {
    mol?.delete();
    throw new Error('RDKit could not make sense of the bonds in that file. Export it as .mol or .sdf (which store the bond orders) and try again.');
  }
  try {
    try { mol.remove_hs_in_place(); } catch { /* older RDKit: removeHs above */ }
    const mb = mol.get_molblock();
    const L = mb.split('\n'), na = parseInt(L[3].slice(0, 3));
    const zs = L.slice(4, 4 + na).map(l => parseFloat(l.slice(20, 30)));
    // a flat file (2D drawing, or a planar molecule exported with z = 0) is handled like a SMILES
    if (zs.every(z => Math.abs(z) < 1e-3)) return { molblock: null, smiles: mol.get_smiles(), flat: true };
    return { molblock: mb, smiles: mol.get_smiles() };
  } finally { mol.delete(); }
}

/**
 * Parse a SMILES string into a 2D molecule graph with coordinates (Å-like units), or — with
 * `molblock` — a 3D structure imported from a file, keeping its real coordinates (mm.is3D).
 */
export function parseMolecule(smiles, molblock = null) {
  const mol = RDKit.get_mol(molblock || smiles.trim());
  if (!mol || !mol.is_valid()) { mol?.delete(); throw new Error(molblock ? 'Could not read the imported structure.' : 'Could not parse that SMILES.'); }
  if (!molblock) try { mol.set_new_coords(true); } catch { try { mol.set_new_coords(); } catch { /* keep coords */ } }
  const mb = mol.get_molblock();
  const lines = mb.split('\n');
  const [na, nb] = [parseInt(lines[3].slice(0, 3)), parseInt(lines[3].slice(3, 6))];
  const atoms = [];
  for (let i = 0; i < na; i++) {
    const l = lines[4 + i];
    atoms.push({ i, x: parseFloat(l.slice(0, 10)), y: parseFloat(l.slice(10, 20)), z: parseFloat(l.slice(20, 30)), el: l.slice(31, 34).trim() });
  }
  const bonds = [];
  for (let k = 0; k < nb; k++) {
    const l = lines[4 + na + k];
    bonds.push({ a: parseInt(l.slice(0, 3)) - 1, b: parseInt(l.slice(3, 6)) - 1, order: parseInt(l.slice(6, 9)) });
  }
  // charges/H counts from the M  CHG lines are not needed for geometry
  const adj = atoms.map(() => []);
  for (const b of bonds) { adj[b.a].push(b.b); adj[b.b].push(b.a); }
  // RDKit 2D depictions use a 1.5 unit bond; scale to ~1.4 Å (aromatic C–C). 3D files are real Å.
  if (!molblock) {
    const f = 1.4 / meanBondLength(atoms, bonds);
    for (const a of atoms) { a.x *= f; a.y *= f; a.z *= f; }
  }
  const mm = { smiles, atoms, bonds, adj, rings: findRings(atoms, adj), rdmol: mol, is3D: !!molblock };
  return mm;
}

function meanBondLength(atoms, bonds) {
  if (!bonds.length) return 1.5;
  return bonds.reduce((s, b) => s + Math.hypot(atoms[b.a].x - atoms[b.b].x, atoms[b.a].y - atoms[b.b].y), 0) / bonds.length;
}

/** Smallest rings (size 3–8) through each bond — an SSSR-like set, adequate for drawing bulges. */
function findRings(atoms, adj) {
  const rings = [];
  const seen = new Set();
  for (let a = 0; a < atoms.length; a++) for (const b of adj[a]) {
    if (b < a) continue;
    // shortest path from b to a not using edge a-b
    const prev = new Map([[b, -1]]);
    const q = [b];
    while (q.length && !prev.has(a)) {
      const u = q.shift();
      for (const w of adj[u]) {
        if ((u === b && w === a) || prev.has(w)) continue;
        prev.set(w, u); q.push(w);
      }
    }
    if (!prev.has(a)) continue;
    const path = [];
    for (let u = a; u !== -1; u = prev.get(u)) path.push(u);
    if (path.length < 3 || path.length > 8) continue;
    const key = [...path].sort((x, y) => x - y).join(',');
    if (seen.has(key)) continue;
    seen.add(key);
    const c = path.reduce((s, i) => [s[0] + atoms[i].x, s[1] + atoms[i].y, s[2] + atoms[i].z], [0, 0, 0]).map(x => x / path.length);
    rings.push({ atoms: path, centre: c, size: path.length });
  }
  return rings;
}

/** Find connection atoms for a group type. Returns [{ atom, ext }]. */
export function findConnections(mm, groupKey) {
  if (groupKey === 'dummy') {
    // a dummy atom stands in for the partner; the point of extension is the atom it is bonded to,
    // extended halfway towards the dummy
    return mm.atoms.filter(a => ['R#', 'R', '*', 'A', 'Q'].includes(a.el)).map(a => ({ atom: mm.adj[a.i][0] ?? a.i, toward: a.i }));
  }
  const g = GROUPS[groupKey];
  if (!g || !g.smarts) return [];
  const q = RDKit.get_qmol(g.smarts);
  let res = [];
  try {
    const raw = JSON.parse(mm.rdmol.get_substruct_matches(q) || '[]');
    // MinimalLib returns "{}" when nothing matches
    res = (Array.isArray(raw) ? raw : []).map(m => ({ atom: m.atoms[g.atom] }));
  } finally { q.delete(); }
  // de-duplicate
  const s = new Set();
  return res.filter(r => (s.has(r.atom) ? false : (s.add(r.atom), true)));
}

export function autoDetectGroup(mm) {
  if (mm.atoms.some(a => ['R#', 'R', '*', 'A', 'Q'].includes(a.el))) return 'dummy';
  for (const k of ['cooh', 'cho', 'nh2', 'boronic', 'nitrile', 'pyridyl']) {
    if (findConnections(mm, k).length >= 2) return k;
  }
  return 'manual';
}

const P2 = a => [a.x, a.y, a.z];

/**
 * Turn a molecule + chosen connection atoms into a block spec (in Å; scale later).
 * opts: { ext (Å), geometryOverride, ringTol (Å) }
 */
export function moleculeToBlock(mm, conns, opts) {
  const n = conns.length;
  if (n < 2) throw new Error('Need at least two connection points to make a block.');
  const pts = conns.map(c => {
    const p = P2(mm.atoms[c.atom]);
    if (c.toward !== undefined) {
      // dummy atom: point of extension is half way along the bond to the dummy
      const d = P2(mm.atoms[c.toward]);
      return v3.add(p, v3.mul(v3.sub(d, p), 0.5));
    }
    return p;
  });
  const centre = findCentre(mm, conns, pts);
  const ext = opts.ext ?? 0;
  const ringTol = opts.ringTol ?? 1.0;

  // rings: centre ring → hub; rings near an arm → bulges
  let hubRing = null;
  const armBulges = pts.map(() => []);
  const unused = [];
  for (const r of mm.rings) {
    const rc = v3.sub(r.centre, centre);
    if (v3.len(rc) < ringTol) { hubRing = r; continue; }
    let placed = false;
    pts.forEach((p, i) => {
      if (placed) return;
      const axis = v3.sub(p, centre);
      const L = v3.len(axis);
      const d = v3.norm(axis);
      const t = v3.dot(rc, d);
      const perp = v3.len(v3.sub(rc, v3.mul(d, t)));
      if (t > 0.5 && t < L && perp < ringTol) { armBulges[i].push({ t, ring: r }); placed = true; }
    });
    if (!placed) unused.push(r);
  }

  let dirs = pts.map(p => v3.norm(v3.sub(p, centre)));
  let lengths = pts.map(p => v3.len(v3.sub(p, centre)) + ext);
  let bulges = armBulges.map(bs => bs.sort((a, b) => a.t - b.t).map(b => b.t));
  dirs = mm.is3D ? orientDirs3D(dirs) : orientDirs(dirs);

  // sp3 centre (C, Si, Ge with four single bonds) → the 2D drawing is a flattened tetrahedron
  // (the depiction centroid can sit an Å or two away from the real centre, so search nearby)
  const isSp3 = a => ['C', 'Si', 'Ge'].includes(a.el) && mm.adj[a.i].length === 4 &&
    mm.bonds.filter(b => b.a === a.i || b.b === a.i).every(b => b.order === 1);
  let ctrAtom = null, best = 2.5;
  for (const a of mm.atoms) {
    const d = v3.len(v3.sub(P2(a), centre));
    if (d < best && isSp3(a)) { best = d; ctrAtom = a; }
  }
  const sp3 = !!ctrAtom;
  let geometry = opts.geometryOverride || '';
  let autoGeometry = false;
  if (!geometry && sp3 && n === 4 && !mm.is3D) { // a 3D file already has the real tetrahedron
    geometry = 'tetrahedral'; autoGeometry = true;
    lengths = pts.map(p => v3.len(v3.sub(p, P2(ctrAtom))) + ext);
  }

  // geometryOverride may also be an explicit list of unit vectors (e.g. the net vertex's arms)
  const ideal = Array.isArray(opts.geometryOverride) ? opts.geometryOverride : GEOMETRIES[geometry];
  if (ideal) {
    if (ideal.length !== n) throw new Error(`That geometry needs ${ideal.length} connection points; this molecule has ${n}.`);
    dirs = ideal;
    // the 2D depiction distorts arm lengths when the real shape is not planar: equalise them
    const mean = lengths.reduce((s, l) => s + l, 0) / n;
    if ((Math.max(...lengths) - Math.min(...lengths)) / mean < 0.4) {
      bulges = bulges.map((bs, i) => bs.map(t => t * mean / lengths[i]));
      lengths = lengths.map(() => mean);
      // same ring count on every arm → use the average ring positions on all arms
      if (bulges.every(b => b.length === bulges[0].length))
        bulges = bulges.map(() => bulges[0].map((_, k) => bulges.reduce((s, b) => s + b[k], 0) / n));
    }
  }
  return { n, centre, dirs, lengths, hubRing: !!hubRing, bulges, unusedRings: unused.length, pts, autoGeometry: autoGeometry ? geometry : '', ctrAtom: ctrAtom ? ctrAtom.i : null, imposed: !!ideal };
}

/** Normal of a set of unit vectors if they are (nearly) coplanar, else null. */
function planeNormal(dirs, tol = 0.05) {
  let n = null;
  for (let a = 0; a < dirs.length && !n; a++) for (let b = a + 1; b < dirs.length; b++) {
    const c = v3.cross(dirs[a], dirs[b]);
    if (v3.len(c) > 0.2) { n = v3.norm(c); break; }
  }
  if (!n) return dirs.length ? (Math.abs(dirs[0][2]) < 0.9 ? v3.norm(v3.cross(dirs[0], [0, 0, 1])) : [1, 0, 0]) : [0, 0, 1];
  return dirs.every(d => Math.abs(v3.dot(d, n)) < tol) ? n : null;
}
/** 3D structures: coplanar arms are turned into the xy plane (arm 1 along +x); others are kept. */
function orientDirs3D(dirs) {
  const n = planeNormal(dirs);
  if (!n) return dirs;
  const R = rotBetween(n, [0, 0, 1]);
  return orientDirs(dirs.map(d => apply(R, d)));
}

/** Planar molecules: put the arms in the xy plane with arm 1 along +x (keeps 2D layout). */
function orientDirs(dirs) {
  const a = Math.atan2(dirs[0][1], dirs[0][0]);
  const c = Math.cos(-a), s = Math.sin(-a);
  return dirs.map(d => [d[0] * c - d[1] * s, d[0] * s + d[1] * c, 0].map(x => (Math.abs(x) < 1e-9 ? 0 : x)));
}

function findCentre(mm, conns, pts) {
  const mid = v3.mul(pts.reduce((s, p) => v3.add(s, p), [0, 0, 0]), 1 / pts.length);
  if (pts.length > 2) {
    // snap to a ring centroid or atom if one sits (almost) at the centroid
    const ring = mm.rings.find(r => v3.len(v3.sub(r.centre, mid)) < 0.8);
    if (ring) return ring.centre;
    return mid;
  }
  // ditopic: straight or bent?
  const [p1, p2] = pts;
  const path = shortestPath(mm, conns[0].atom, conns[1].atom);
  const axis = v3.norm(v3.sub(p2, p1));
  const devP = q => { const r = v3.sub(q, p1); return v3.len(v3.sub(r, v3.mul(axis, v3.dot(r, axis)))); };
  const dev = i => devP(P2(mm.atoms[i]));
  const onPath = new Set(path);
  // straight if every ring the path runs through is centred on the p1–p2 line
  // (para-phenylene path atoms sit ~1.2 Å off-axis, so judge rings by their centroids)
  const pathRings = mm.rings.filter(r => r.atoms.some(a => onPath.has(a)));
  const straight = pathRings.length
    ? pathRings.every(r => devP(r.centre) < 0.8)
    : Math.max(...path.map(dev)) < 1.0;
  if (straight) return mid;
  // bent: the kink is the ring on the path farthest from the p1–p2 line, else the farthest atom
  let best = null, bestD = -1;
  for (const r of mm.rings) {
    if (!r.atoms.some(a => onPath.has(a))) continue;
    const rr = v3.sub(r.centre, p1);
    const d = v3.len(v3.sub(rr, v3.mul(axis, v3.dot(rr, axis))));
    if (d > bestD) { bestD = d; best = r.centre; }
  }
  if (best) return best;
  let far = path[0];
  for (const i of path) if (dev(i) > dev(far)) far = i;
  return P2(mm.atoms[far]);
}

function shortestPath(mm, s, t) {
  const prev = new Map([[s, -1]]);
  const q = [s];
  while (q.length) {
    const u = q.shift();
    if (u === t) break;
    for (const w of mm.adj[u]) if (!prev.has(w)) { prev.set(w, u); q.push(w); }
  }
  const path = [];
  for (let u = t; u !== undefined && u !== -1; u = prev.get(u)) path.push(u);
  return path;
}

/**
 * Convert the Å-scale molecular block into a printable spec (mm).
 * kind: 'node' | 'linker'. Linkers get no hub unless a ring sits at the centre.
 */
export function chemBlockToSpec(cb, { name, kind, scale, P, polarity }) {
  const hubR = kind === 'node' ? P.hubRadius : (cb.hubRing ? P.bulgeRadius : 0);
  return {
    name, kind, polarity, hubRadius: hubR,
    arms: cb.dirs.map((d, i) => ({
      dir: d,
      length: cb.lengths[i] * scale,
      bulges: cb.bulges[i].map(t => ({ t: t * scale, r: P.bulgeRadius })),
    })),
  };
}

/**
 * Molecular (stick) model: the real atoms and bonds, with a magnet fitting at each connection point.
 * The fittings stand for the bonds made on assembly (imine C=N, B–O, M–O₂C coordination…), so
 * leaving atoms (aldehyde O, boronic OH, nitrile N, dummy atoms, H) are removed; carboxylate O's stay.
 * Arms keep the block convention (dir + length centre → magnet face) so fit checks, the assembly
 * and the bill of materials work unchanged; the atoms are posed so each arm tip lies on its direction.
 * targets: unit vectors the arms must follow (the net's, or an ideal geometry) or null (as drawn).
 */
export function moleculeModelSpec(mm, conns, cb, { group, ext = 0, scale, P, targets = null, name, kind, polarity }) {
  const n = conns.length;
  const c0 = cb.ctrAtom != null ? P2(mm.atoms[cb.ctrAtom]) : cb.centre;
  let pos = mm.atoms.map(a => v3.sub(P2(a), c0));
  let pts = cb.pts.map(p => v3.sub(p, c0));
  const u = pts.map(p => v3.norm(p));
  const tset = targets || cb.dirs;
  // A conjugated (2D-drawn) molecule stays rigid and flat; only a drawn sp³ centre (e.g. TAM) is
  // opened up into 3D. A flat molecule meeting a non-flat vertex (BTC in tbo) sits in the plane
  // normal to the arms' resultant, and its magnet fittings bend to meet the net instead.
  // A 3D structure from a file is already its real shape: it is only turned as a whole.
  const rigid = cb.ctrAtom == null || mm.is3D;
  let aim = tset, nt = null;
  const sum = tset.reduce((a, t) => v3.add(a, t), [0, 0, 0]);
  if (rigid && planeNormal(u) && tset.length >= 3 && v3.len(sum) > 0.05 && tset.some(t => Math.abs(v3.dot(t, v3.norm(sum))) > 0.02)) {
    nt = v3.norm(sum);
    aim = tset.map(t => v3.norm(v3.sub(t, v3.mul(nt, v3.dot(t, nt)))));
  }
  const al = alignDirs(u, aim, { up: [0, 0, 1] }) || { R: IDENT, perm: u.map((_, i) => i) };
  const tgt = u.map((_, i) => tset[al.perm[i]]);
  pos = pos.map(p => apply(al.R, p));
  pts = pts.map(p => apply(al.R, p));
  if (nt) {
    // slide the flat molecule along the cone axis so its connection atoms lie on the arm lines
    // (as in HKUST-1, where the flat BTC sits off the ideal vertex point)
    const hs = conns.map((c, i) => {
      const d = tgt[i], dn = v3.dot(d, nt), dp = v3.len(v3.sub(d, v3.mul(nt, dn)));
      const r = v3.len(v3.sub(pos[c.atom], v3.mul(nt, v3.dot(pos[c.atom], nt))));
      return dp > 1e-6 ? r * dn / dp : 0;
    });
    const h = hs.reduce((a, b) => a + b, 0) / hs.length;
    pos = pos.map(p => v3.add(p, v3.mul(nt, h)));
    pts = pts.map(p => v3.add(p, v3.mul(nt, h)));
  }

  // core (stays put) and branches (each rotated rigidly so its point of extension lies on its arm)
  let core = [];
  if (cb.ctrAtom != null) core = [cb.ctrAtom];
  else {
    const ring = mm.rings.find(r => v3.len(v3.sub(r.centre, cb.centre)) < 1.0);
    if (ring) core = ring.atoms;
    else {
      const at = mm.atoms.find(a => v3.len(v3.sub(P2(a), cb.centre)) < 0.35);
      if (at) core = [at.i];
    }
  }
  if (core.length && !rigid) {
    const inCore = new Set(core), seen = new Set(core);
    for (const start of mm.atoms.map(a => a.i)) {
      if (seen.has(start)) continue;
      const comp = [start]; seen.add(start);
      for (let k = 0; k < comp.length; k++) for (const w of mm.adj[comp[k]]) if (!seen.has(w)) { seen.add(w); comp.push(w); }
      const set = new Set(comp);
      const ci = conns.map((c, i) => (set.has(c.atom) ? i : -1)).filter(i => i >= 0);
      if (ci.length !== 1) continue;
      const i = ci[0];
      const pivotIdx = comp.flatMap(a => mm.adj[a]).find(w => inCore.has(w));
      if (pivotIdx === undefined) continue;
      const pv = pos[pivotIdx];
      // the point on the arm line at the same distance from the pivot (rigid rotation keeps that distance)
      const r = v3.len(v3.sub(pts[i], pv)), tp = v3.dot(tgt[i], pv);
      const disc = tp * tp - v3.dot(pv, pv) + r * r;
      const want = v3.mul(tgt[i], disc > 0 ? tp + Math.sqrt(disc) : v3.len(pts[i]));
      const Q = rotBetween(v3.norm(v3.sub(pts[i], pv)), v3.norm(v3.sub(want, pv)));
      for (const a of comp) pos[a] = v3.add(pv, apply(Q, v3.sub(pos[a], pv)));
      pts[i] = v3.add(pv, apply(Q, v3.sub(pts[i], pv)));
    }
  }

  // atoms to leave out
  const drop = new Set(mm.atoms.filter(a => a.el === 'H' || ['R#', 'R', '*', 'A', 'Q'].includes(a.el)).map(a => a.i));
  for (const c of conns) if (c.toward !== undefined) drop.add(c.toward);
  if (GROUPS[group]?.drop) for (const c of conns)
    for (const w of mm.adj[c.atom]) if (mm.adj[w].length === 1 && mm.atoms[w].el !== 'C') drop.add(w);

  const clean = x => (Math.abs(x) < 1e-6 ? 0 : x);
  const idx = new Map();
  const atoms = [];
  for (const a of mm.atoms) {
    if (drop.has(a.i)) continue;
    idx.set(a.i, atoms.length);
    atoms.push({ p: v3.mul(pos[a.i], scale).map(clean), el: a.el, src: a.i });
  }
  const bonds = mm.bonds.filter(b => idx.has(b.a) && idx.has(b.b)).map(b => [idx.get(b.a), idx.get(b.b), b.order]);

  // arm (fitting) lengths: chemistry, but long enough for the flare to clear the atom
  let lengths = conns.map((c, i) => {
    const A = atoms[idx.get(c.atom)].p, d = tgt[i];
    const f = P.flareLength + atomRadius(mm.atoms[c.atom].el, P) * 0.7 + 1;
    const dA = v3.dot(d, A), disc = f * f - v3.dot(A, A) + dA * dA;
    const Lmin = dA + (disc > 0 ? Math.sqrt(disc) : f);
    return Math.max((v3.len(pts[i]) + ext) * scale, Lmin);
  });
  if (targets || Math.max(...lengths) - Math.min(...lengths) < 0.1 * Math.max(...lengths)) {
    const L = Math.max(...lengths); lengths = lengths.map(() => L);
  }
  return {
    name, kind, polarity, hubRadius: 0, style: 'model',
    arms: tgt.map((d, i) => ({ dir: d.map(clean), length: lengths[i], bulges: [] })),
    model: {
      atoms, bonds, fittings: conns.map(c => idx.get(c.atom)),
      // how far each fitting's rod bends away from its net edge (degrees)
      bend: Math.round(Math.max(...tgt.map((d, i) => {
        const A = atoms[idx.get(conns[i].atom)].p, F = v3.sub(v3.mul(d, lengths[i]), v3.mul(d, P.flareLength));
        return Math.acos(Math.max(-1, Math.min(1, v3.dot(v3.norm(v3.sub(F, A)), d)))) * 180 / Math.PI;
      }))),
    },
  };
}

/** SBU (inorganic cluster) → spec in mm. */
export const sbuCN = key => GEOMETRIES[SBUS[key].geom].length;

export function sbuToSpec(key, { scale, P, polarity, ext = 0, dirs = null }) {
  const s = SBUS[key];
  return {
    name: key.split(' ')[0], kind: 'node', polarity, hubRadius: P.hubRadius,
    arms: (dirs || GEOMETRIES[s.geom]).map(d => ({ dir: d, length: (s.r + ext) * scale, bulges: [] })),
    minLengthNote: true,
  };
}

/**
 * SBU as a molecular model: the metal core posed so its arms lie along `dirs` (the net's, or the
 * ideal cluster shape), with a magnet fitting per carboxylate that forks from the metals it bridges.
 */
const clean = x => (Math.abs(x) < 1e-6 ? 0 : x);
export function sbuModelSpec(key, { scale, P, polarity, ext = 0, dirs = null }) {
  const s = SBUS[key];
  const core = sbuCore(s.geom, s.el);
  const tgt = dirs || GEOMETRIES[s.geom];
  const al = alignDirs(core.dirs, tgt, { improper: false });
  const R = al && al.err < 0.2 ? al.R : alignDirs(core.dirs, GEOMETRIES[s.geom], { improper: false }).R;
  const arms = al && al.err < 0.2 ? tgt : GEOMETRIES[s.geom];
  const atoms = core.atoms.map(([el, p]) => ({ p: apply(R, p).map(x => clean(x * scale)), el }));
  // each arm's fitting starts at the metal(s) nearest its direction (ties → bridging, e.g. both Cu)
  const fittings = arms.map(d => {
    const metals = atoms.map((a, i) => i).filter(i => atoms[i].el === s.el);
    const score = i => v3.len(atoms[i].p) < 1e-6 ? 0 : v3.dot(v3.norm(atoms[i].p), d);
    const best = Math.max(...metals.map(score));
    return metals.filter(i => score(i) > best - 0.05);
  });
  const length = arms.map((d, k) => {
    // tip at least `f` from every fitted atom so the flare clears it (as for molecules)
    const f = P.flareLength + atomRadius(s.el, P) * 0.7 + 1;
    const need = Math.max(...fittings[k].map(i => {
      const A = atoms[i].p, dA = v3.dot(A, d);
      return dA + Math.sqrt(Math.max(0, f * f - v3.dot(A, A) + dA * dA));
    }));
    return Math.max((s.r + ext) * scale, need);
  });
  const L = Math.max(...length); // all carboxylates are alike: one length
  return {
    name: key.split(' ')[0] + ' (model)', kind: 'node', polarity, hubRadius: P.hubRadius, style: 'model',
    arms: arms.map(d => ({ dir: d.map(clean), length: L, bulges: [] })),
    model: { atoms, bonds: core.bonds, fittings, bend: 0, sbu: true },
  };
}

export { spacedBulges };
