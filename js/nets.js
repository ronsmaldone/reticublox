// Reticublox — net (topology) library.
// Each net is a periodic embedding: cell vectors, vertices in fractional coordinates
// with a "shape class" label, and optional rules for which vertex classes may bond.
// Edges are found geometrically (shortest allowed distance, within a tolerance) and
// every net is checked against its expected coordination numbers at load time.

const S3 = Math.sqrt(3);
const HEX = [[1, 0, 0], [-0.5, S3 / 2, 0], [0, 0, 1]];

/** All positions generated from `p` by the Fm-3m point operations + F-centring (fractional, deduplicated). */
function fmOrbit(p) {
  const perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
  const cent = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]];
  const seen = new Map();
  for (const pm of perms) for (let s = 0; s < 8; s++) for (const c of cent) {
    const q = pm.map((k, i) => {
      const v = (((s >> i) & 1 ? -1 : 1) * p[k] + c[i]) % 1;
      const w = v < 0 ? v + 1 : v;
      return Math.abs(w - 1) < 1e-9 ? 0 : w;
    });
    seen.set(q.map(x => x.toFixed(6)).join(','), q);
  }
  return [...seen.values()];
}

/** Apply a list of symmetry operations (fractional) + centring translations, deduplicated mod 1. */
function orbitOps(p, ops, cent = [[0, 0, 0]]) {
  const seen = new Map();
  for (const op of ops) {
    const q = op(p);
    for (const c of cent) {
      const w = q.map((v, i) => { let r = (v + c[i]) % 1; if (r < 0) r += 1; return Math.abs(r - 1) < 1e-9 ? 0 : r; });
      seen.set(w.map(x => x.toFixed(6)).join(','), w);
    }
  }
  return [...seen.values()];
}
// P-43m: all coordinate permutations with an even number of sign changes
const P_43M = [];
for (const pm of [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]])
  for (const s of [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
    P_43M.push(p => pm.map((k, i) => s[i] * p[k]));
// I-43d (No. 220), International Tables general position; body centring added separately
const h = 0.5, q4 = 0.25, t4 = 0.75;
const I_43D = [
  ([x, y, z]) => [x, y, z], ([x, y, z]) => [-x + h, -y, z + h], ([x, y, z]) => [-x, y + h, -z + h], ([x, y, z]) => [x + h, -y + h, -z],
  ([x, y, z]) => [z, x, y], ([x, y, z]) => [z + h, -x + h, -y], ([x, y, z]) => [-z + h, -x, y + h], ([x, y, z]) => [-z, x + h, -y + h],
  ([x, y, z]) => [y, z, x], ([x, y, z]) => [-y, z + h, -x + h], ([x, y, z]) => [y + h, -z + h, -x], ([x, y, z]) => [-y + h, -z, x + h],
  ([x, y, z]) => [y + q4, x + q4, z + q4], ([x, y, z]) => [-y + q4, -x + t4, z + t4], ([x, y, z]) => [y + t4, -x + q4, -z + t4], ([x, y, z]) => [-y + t4, x + t4, -z + q4],
  ([x, y, z]) => [x + q4, z + q4, y + q4], ([x, y, z]) => [-x + t4, z + t4, -y + q4], ([x, y, z]) => [-x + q4, -z + t4, y + t4], ([x, y, z]) => [x + t4, -z + q4, -y + t4],
  ([x, y, z]) => [z + q4, y + q4, x + q4], ([x, y, z]) => [z + t4, -y + q4, -x + t4], ([x, y, z]) => [-z + t4, y + t4, -x + q4], ([x, y, z]) => [-z + q4, -y + t4, x + t4],
];
const I_CENT = [[0, 0, 0], [0.5, 0.5, 0.5]];

export const NETS = {
  sql: {
    name: 'sql — square lattice', dim: 2,
    examples: 'COF-366, 2D porphyrin COFs, Cu-paddlewheel layers',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    verts: [{ t: 'A', f: [0, 0, 0] }],
    cn: { A: 4 },
  },
  hcb: {
    name: 'hcb — honeycomb', dim: 2,
    examples: 'COF-1, COF-5, TAPB–PDA COF, most [3+2] and [3+3] 2D COFs',
    cell: HEX,
    verts: [{ t: 'A', f: [1 / 3, 2 / 3, 0] }, { t: 'A', f: [2 / 3, 1 / 3, 0] }],
    cn: { A: 3 },
  },
  kgm: {
    name: 'kgm — kagome', dim: 2,
    examples: 'Dual-pore [4+2] COFs (e.g. ETTA-based)',
    cell: HEX,
    verts: [{ t: 'A', f: [0.5, 0, 0] }, { t: 'A', f: [0, 0.5, 0] }, { t: 'A', f: [0.5, 0.5, 0] }],
    cn: { A: 4 },
  },
  hxl: {
    name: 'hxl — hexagonal (6-c)', dim: 2,
    examples: '6-connected planar nodes, e.g. hexaphenylbenzene / HATN cores',
    cell: HEX,
    verts: [{ t: 'A', f: [0, 0, 0] }],
    cn: { A: 6 },
  },
  kgd: {
    name: 'kgd — kagome dual (3,6-c)', dim: 2,
    examples: '[6+3] 2D COFs and MOFs',
    cell: HEX,
    verts: [{ t: 'H', f: [0, 0, 0] }, { t: 'T', f: [1 / 3, 2 / 3, 0] }, { t: 'T', f: [2 / 3, 1 / 3, 0] }],
    pairs: [['H', 'T']],
    cn: { H: 6, T: 3 },
  },
  dia: {
    name: 'dia — diamond', dim: 3,
    examples: 'COF-300, COF-320, ZIF-like tetrahedral frameworks',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    verts: [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0],
            [0.25, 0.25, 0.25], [0.25, 0.75, 0.75], [0.75, 0.25, 0.75], [0.75, 0.75, 0.25]]
      .map(f => ({ t: 'A', f })),
    cn: { A: 4 },
  },
  pcu: {
    name: 'pcu — primitive cubic', dim: 3,
    examples: 'MOF-5 / IRMOF series (Zn4O + BDC), pillared-paddlewheel MOFs',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    verts: [{ t: 'A', f: [0, 0, 0] }],
    cn: { A: 6 },
  },
  bcu: {
    name: 'bcu — body-centred cubic (8-c)', dim: 3,
    examples: '8-connected cluster MOFs',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    verts: [{ t: 'A', f: [0, 0, 0] }, { t: 'A', f: [0.5, 0.5, 0.5] }],
    cn: { A: 8 },
  },
  fcu: {
    name: 'fcu — face-centred cubic (12-c)', dim: 3,
    examples: 'UiO-66 / UiO-67 (Zr6 cluster)',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    verts: [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]].map(f => ({ t: 'A', f })),
    cn: { A: 12 },
  },
  nbo: {
    name: 'nbo — NbO (square planar)', dim: 3,
    examples: 'MOF-505, NOTT-100 (Cu paddlewheel + tetracarboxylate)',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    verts: [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0], [0, 0, 0.5], [0, 0.5, 0], [0.5, 0, 0]]
      .map(f => ({ t: 'A', f })),
    cn: { A: 4 },
  },
  tbo: {
    name: 'tbo — twisted boracite (3,4-c)', dim: 3,
    examples: 'HKUST-1 / Cu3(BTC)2, MOF-399. With straight edges a square paddlewheel (S) forces the BTC vertex (T) to be slightly pyramidal (117°); in the crystal the bend sits at the carboxylate instead. Use direct node–node edges: BTC node ↔ paddlewheel node',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    // x = (2 − √2)/4 makes the 4-c vertex exactly square planar with all edges equal
    verts: [
      ...fmOrbit(Array(3).fill((2 - Math.SQRT2) / 4)).map(f => ({ t: 'T', f })),
      ...fmOrbit([0, 0.25, 0.25]).map(f => ({ t: 'S', f })),
    ],
    pairs: [['T', 'S']],
    cn: { T: 3, S: 4 },
  },
  ctn: {
    name: 'ctn — C3N4 (3,4-c)', dim: 3,
    examples: 'COF-102, COF-103 (tetrahedral node + triangular boroxine / HHTP). Near-ideal: flat triangles (2° off) and tetrahedra (2° off)',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    // I-43d: 3-c on 16c (x,x,x), x tuned so both vertices are within ~2.3° of ideal; 4-c on 12a
    verts: [
      ...orbitOps([0.21458, 0.21458, 0.21458], I_43D, I_CENT).map(f => ({ t: 'T', f })),
      ...orbitOps([0.375, 0, 0.25], I_43D, I_CENT).map(f => ({ t: 'Q', f })),
    ],
    pairs: [['T', 'Q']],
    cn: { T: 3, Q: 4 },
  },
  bor: {
    name: 'bor — boracite (3,4-c)', dim: 3,
    examples: 'COF-105, COF-108 (tetrahedral node + triangle). The tetrahedral node (Q) is kept ideal, the same block as in dia; straight edges then make the 3-c node (T) pyramidal (arms 19.5° out of plane)',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    // P-43m: 3-c on 4e (x,x,x) with x = 1/4 (ideal tetrahedra), 4-c on 3c (0,½,½)
    verts: [
      ...orbitOps([0.25, 0.25, 0.25], P_43M).map(f => ({ t: 'T', f })),
      ...orbitOps([0, 0.5, 0.5], P_43M).map(f => ({ t: 'Q', f })),
    ],
    pairs: [['T', 'Q']],
    cn: { T: 3, Q: 4 },
  },
  pts: {
    name: 'pts — PtS (4,4-c), square kept ideal', dim: 3,
    examples: 'MOF-11 (Cu paddlewheel + tetrahedral ATC), 3D porphyrin COFs (TAPP/TFPP + tetrahedral). c/a = 2: the square node (S) is exact, the tetrahedral node (Q) is flattened (90°, 90°, 120°)',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 2]],
    // P4₂/mmc: square 4-c on 2c, tetrahedral 4-c on 2e
    verts: [
      { t: 'S', f: [0, 0.5, 0] }, { t: 'S', f: [0.5, 0, 0.5] },
      { t: 'Q', f: [0, 0, 0.25] }, { t: 'Q', f: [0, 0, 0.75] },
    ],
    pairs: [['S', 'Q']],
    cn: { S: 4, Q: 4 },
  },
  pts_t: {
    name: 'pts — PtS (4,4-c), tetrahedron kept ideal', dim: 3,
    examples: 'Same net as pts, with c/a = √2: the tetrahedral node (Q) is exact (the same block as in dia, bor and ctn) and the square node (S) becomes a 70.5° × 109.5° rectangle',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, Math.SQRT2]],
    verts: [
      { t: 'S', f: [0, 0.5, 0] }, { t: 'S', f: [0.5, 0, 0.5] },
      { t: 'Q', f: [0, 0, 0.25] }, { t: 'Q', f: [0, 0, 0.75] },
    ],
    pairs: [['S', 'Q']],
    cn: { S: 4, Q: 4 },
  },
  stp: {
    name: 'stp — (4,6-c) trigonal prism + square', dim: 3,
    examples: 'PCN-600 (Fe₃O(CO₂)₆ trimer + TCPP porphyrin). With c = a/√3 both nodes are exact: an ideal square and a trigonal prism',
    cell: [[1, 0, 0], [-0.5, S3 / 2, 0], [0, 0, 1 / S3]],
    // P6/mmm: 6-c prism on 2c, square 4-c on 3g
    verts: [
      { t: 'P', f: [1 / 3, 2 / 3, 0] }, { t: 'P', f: [2 / 3, 1 / 3, 0] },
      { t: 'S', f: [0.5, 0, 0.5] }, { t: 'S', f: [0, 0.5, 0.5] }, { t: 'S', f: [0.5, 0.5, 0.5] },
    ],
    pairs: [['P', 'S']],
    cn: { P: 6, S: 4 },
  },
  srs: {
    name: 'srs — SrSi2 (3-c, chiral)', dim: 3,
    examples: 'Chiral 3-c MOFs/COFs; gyroid-like',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    verts: (() => {
      const base = [[1 / 8, 1 / 8, 1 / 8], [3 / 8, 7 / 8, 5 / 8], [7 / 8, 5 / 8, 3 / 8], [5 / 8, 3 / 8, 7 / 8]];
      const all = [];
      for (const b of base) {
        all.push({ t: 'A', f: b });
        all.push({ t: 'A', f: b.map(x => (x + 0.5) % 1) });
      }
      return all;
    })(),
    cn: { A: 3 },
  },
  acs: {
    name: 'acs — trigonal prism (6-c)', dim: 3,
    examples: 'MIL-88, MOF-235 (M3O trimer + linear dicarboxylate)',
    cell: [[1, 0, 0], [-0.5, S3 / 2, 0], [0, 0, 1.15]],
    verts: [{ t: 'A', f: [1 / 3, 2 / 3, 0.25] }, { t: 'A', f: [2 / 3, 1 / 3, 0.75] }],
    pairs: [['A', 'A']],
    onlyBetweenSites: true, // bonds only between the two layers (prisms), not within a layer
    cn: { A: 6 },
  },
  sod: {
    name: 'sod — sodalite', dim: 3,
    examples: 'ZIF-8, ZIF-90 (Zn + imidazolate); truncated-octahedron cages',
    cell: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    // Im-3m, vertex at (¼, 0, ½): the 6 positions of one orbit + body centring = 12 per cubic cell
    verts: [[0.25, 0, 0.5], [0.75, 0, 0.5], [0.5, 0.25, 0], [0.5, 0.75, 0], [0, 0.5, 0.25], [0, 0.5, 0.75]]
      .flatMap(f => [f, f.map(x => (x + 0.5) % 1)]).map(f => ({ t: 'A', f })),
    cn: { A: 4 },
  },
};

// Per-net metadata for the builder: space group (plane group for 2D layers), a short plain
// name, and real MOF / COF examples ('' = no common examples of that class).
const META = {
  sql:   { sg: 'p4mm',     title: 'square lattice',          mof: 'Cu-paddlewheel square-grid layers, Zn₂(bdc)₂ sheets', cof: 'COF-366 (porphyrin), 2D square imine COFs' },
  hcb:   { sg: 'p6mm',     title: 'honeycomb',               mof: 'Cu₃(HHTP)₂ conductive 2D MOFs', cof: 'COF-1, COF-5, TAPB–PDA, most [3+2] and [3+3] COFs' },
  kgm:   { sg: 'p6mm',     title: 'kagome',                  mof: 'Kagome paddlewheel layers', cof: 'Dual-pore [4+2] COFs (ETTA-based)' },
  hxl:   { sg: 'p6mm',     title: 'hexagonal 6-c',           mof: '', cof: '6-c planar cores (HATN, hexaphenylbenzene)' },
  kgd:   { sg: 'p6mm',     title: 'kagome dual (3,6-c)',     mof: '2D [6+3] MOFs', cof: '[6+3] COFs' },
  dia:   { sg: 'Fd-3m',    title: 'diamond',                 mof: 'Diamondoid Cu(I)/Zn MOFs, ZIF-like nets', cof: 'COF-300, COF-320 (tetrahedral + linear)' },
  pcu:   { sg: 'Pm-3m',    title: 'primitive cubic',         mof: 'MOF-5 / IRMOF series, pillared paddlewheels (DMOF)', cof: '' },
  bcu:   { sg: 'Im-3m',    title: 'body-centred cubic 8-c',  mof: '8-c Zr₆ cluster MOFs', cof: '' },
  fcu:   { sg: 'Fm-3m',    title: 'face-centred cubic 12-c', mof: 'UiO-66, UiO-67 (Zr₆ cluster)', cof: '' },
  nbo:   { sg: 'Im-3m',    title: 'NbO (square planar)',     mof: 'MOF-505, NOTT-100', cof: '' },
  tbo:   { sg: 'Fm-3m',    title: 'twisted boracite (3,4-c)', mof: 'HKUST-1 / Cu₃(BTC)₂, MOF-399', cof: '' },
  ctn:   { sg: 'I-43d',    title: 'C₃N₄ (3,4-c)',            mof: '', cof: 'COF-102, COF-103' },
  bor:   { sg: 'P-43m',    title: 'boracite (3,4-c)',        mof: '', cof: 'COF-105, COF-108' },
  pts:   { sg: 'P4₂/mmc',  title: 'PtS (4,4-c) · square ideal', mof: 'MOF-11', cof: '3D porphyrin COFs' },
  pts_t: { sg: 'P4₂/mmc',  title: 'PtS (4,4-c) · tetrahedron ideal', mof: 'MOF-11', cof: '3D porphyrin COFs' },
  stp:   { sg: 'P6/mmm',   title: 'trigonal prism + square (4,6-c)', mof: 'PCN-600 (Fe₃O trimer + TCPP)', cof: '' },
  srs:   { sg: 'I4₁32',    title: 'SrSi₂ (3-c, chiral)',     mof: 'Chiral 3-c MOFs', cof: 'srs COFs (rare)' },
  acs:   { sg: 'P6₃/mmc',  title: 'trigonal prism 6-c',      mof: 'MIL-88, MOF-235', cof: '' },
  sod:   { sg: 'Im-3m',    title: 'sodalite',                mof: 'ZIF-8, ZIF-90 (Zn imidazolates)', cof: '' },
};
for (const [k, m] of Object.entries(META)) Object.assign(NETS[k], m, { note: NETS[k].examples });

// ---------- small vector helpers ----------
export const v3 = {
  add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
  sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
  mul: (a, s) => [a[0] * s, a[1] * s, a[2] * s],
  dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
  cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
  len: a => Math.hypot(a[0], a[1], a[2]),
  norm: a => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; },
};

function frac2cart(cell, f) {
  return [0, 1, 2].map(k => f[0] * cell[0][k] + f[1] * cell[1][k] + f[2] * cell[2][k]);
}

/**
 * Resolve a net: find edges, validate coordination, normalise so edge length = 1,
 * work out vertex shape classes, bipartite colouring and canonical arm directions.
 */
export function analyseNet(key) {
  const net = NETS[key];
  if (!net) throw new Error(`Unknown net ${key}`);
  const cell = net.cell;
  const nV = net.verts.length;
  const cart = net.verts.map(v => frac2cart(cell, v.f));
  const range = net.dim === 2 ? [[-1, 0, 1], [-1, 0, 1], [0]] : [[-1, 0, 1], [-1, 0, 1], [-1, 0, 1]];
  const allowed = (ta, tb) => !net.pairs || net.pairs.some(([p, q]) => (p === ta && q === tb) || (p === tb && q === ta));

  // candidate pairs
  const cands = [];
  for (let i = 0; i < nV; i++) for (let j = 0; j < nV; j++) {
    if (!allowed(net.verts[i].t, net.verts[j].t)) continue;
    if (net.onlyBetweenSites && net.verts[i].f[2] === net.verts[j].f[2]) continue;
    for (const a of range[0]) for (const b of range[1]) for (const c of range[2]) {
      if (i === j && a === 0 && b === 0 && c === 0) continue;
      const pj = v3.add(cart[j], frac2cart(cell, [a, b, c]));
      const d = v3.len(v3.sub(pj, cart[i]));
      cands.push({ i, j, off: [a, b, c], d });
    }
  }
  const dmin = Math.min(...cands.map(c => c.d));
  const edgesDir = cands.filter(c => c.d < dmin * 1.03);

  // neighbour vectors per vertex (normalised so edge length = 1)
  const nbrs = Array.from({ length: nV }, () => []);
  for (const e of edgesDir) {
    const pj = v3.add(cart[e.j], frac2cart(cell, e.off));
    nbrs[e.i].push({ j: e.j, off: e.off, vec: v3.mul(v3.sub(pj, cart[e.i]), 1 / dmin) });
  }

  // validate coordination
  const problems = [];
  net.verts.forEach((v, i) => {
    if (net.cn[v.t] !== nbrs[i].length) problems.push(`vertex ${i} (${v.t}) has ${nbrs[i].length} neighbours, expected ${net.cn[v.t]}`);
  });

  // undirected unique edge list (i<j or same with positive offset)
  const edges = [];
  for (const e of edgesDir) {
    const key2 = `${e.j},${e.i},${e.off.map(x => -x).join(',')}`;
    if (edges.some(x => x.key === key2)) continue;
    edges.push({ ...e, key: `${e.i},${e.j},${e.off.join(',')}` });
  }

  // bipartite 2-colouring (periodic: colour must also be consistent through offsets — we just
  // check the quotient graph and require no odd cycles there; good enough for these nets)
  const colour = new Array(nV).fill(-1);
  let bipartite = true;
  for (let s = 0; s < nV; s++) {
    if (colour[s] !== -1) continue;
    colour[s] = 0;
    const q = [s];
    while (q.length) {
      const u = q.shift();
      for (const n of nbrs[u]) {
        if (colour[n.j] === -1) { colour[n.j] = 1 - colour[u]; q.push(n.j); }
        else if (colour[n.j] === colour[u]) bipartite = false;
      }
    }
  }

  // shape classes: one canonical arm set per vertex type label
  const types = {};
  net.verts.forEach((v, i) => {
    if (!types[v.t]) types[v.t] = { label: v.t, cn: nbrs[i].length, dirs: canonicalDirs(nbrs[i].map(n => v3.norm(n.vec))), count: 0 };
    types[v.t].count++;
  });
  for (const t of Object.values(types)) t.geometry = describeGeometry(t.dirs);

  const scaledCell = cell.map(r => v3.mul(r, 1 / dmin));
  return {
    key, net, dim: net.dim, cell: scaledCell, cart: cart.map(c => v3.mul(c, 1 / dmin)),
    nbrs, edges, types, bipartite, colour, problems,
    edgesPerCell: edges.length, vertsPerCell: nV,
  };
}

/** Rotate a set of unit vectors into a print-friendly canonical orientation. */
export function canonicalDirs(dirs) {
  let sum = dirs.reduce((s, d) => v3.add(s, d), [0, 0, 0]);
  // within ~3° of planar (e.g. ctn's triangles): flatten — magnets absorb that little strain
  if (dirs.length >= 3 && v3.len(sum) > 1e-6) {
    const ax = v3.norm(sum);
    if (dirs.every(d => Math.abs(v3.dot(d, ax)) < 0.05)) {
      dirs = dirs.map(d => v3.norm(v3.sub(d, v3.mul(ax, v3.dot(d, ax)))));
      sum = dirs.reduce((s, d) => v3.add(s, d), [0, 0, 0]);
    }
  }
  // planar? find normal via cross products
  let normal = null;
  for (let a = 0; a < dirs.length && !normal; a++) for (let b = a + 1; b < dirs.length; b++) {
    const c = v3.cross(dirs[a], dirs[b]);
    if (v3.len(c) > 0.2) { normal = v3.norm(c); break; }
  }
  let up;
  const planar = !normal || dirs.every(d => Math.abs(v3.dot(d, normal)) < 1e-3);
  if (planar) up = normal || perpendicular(dirs[0]);
  else {
    // choose an arm as "up" that leaves the most arms horizontal (best for printing)
    let best = -1;
    for (const d of dirs) {
      const score = dirs.filter(e => Math.abs(v3.dot(e, d)) < 1e-3).length;
      if (score > best) { best = score; up = d; }
    }
    // pyramidal arm sets (e.g. the 117° 3-c vertex of tbo): stand it on its arm tips, hub up
    if (best === 0 && v3.len(sum) > 0.1) up = v3.mul(sum, -1);
  }
  // build rotation taking `up` -> z, first arm's horizontal component -> x
  const z = v3.norm(up);
  let ref = dirs.find(d => Math.abs(v3.dot(d, z)) < 0.99) || [1, 0, 0];
  const x = v3.norm(v3.sub(ref, v3.mul(z, v3.dot(ref, z))));
  const y = v3.cross(z, x);
  return dirs.map(d => clean([v3.dot(d, x), v3.dot(d, y), v3.dot(d, z)]));
}

function perpendicular(a) {
  const t = Math.abs(a[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  return v3.norm(v3.cross(a, t));
}
const clean = a => a.map(x => (Math.abs(x) < 1e-9 ? 0 : x));

export function describeGeometry(dirs) {
  const n = dirs.length;
  const planar = dirs.every(d => Math.abs(d[2]) < 1e-3);
  const angles = [];
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++)
    angles.push(Math.round(Math.acos(Math.max(-1, Math.min(1, v3.dot(dirs[i], dirs[j])))) * 180 / Math.PI));
  const uniq = [...new Set(angles)].sort((a, b) => a - b);
  const names = {
    2: planar && uniq.includes(180) ? 'linear' : 'bent',
    // angle sum ≥ 354°: a slightly buckled node from a crystal structure is still trigonal planar
    3: planar || angles.reduce((x, y) => x + y, 0) >= 354 ? 'trigonal planar' : 'trigonal pyramidal',
    4: planar ? (uniq[0] === 90 ? 'square planar' : 'planar 4-c (rectangle/kagome)') : (uniq[0] === 109 ? 'tetrahedral' : 'distorted tetrahedral'),
    6: planar ? 'hexagonal planar' : (uniq.includes(180) && uniq[0] === 90 ? 'octahedral' : 'trigonal prismatic'),
    8: 'cubic',
    12: 'cuboctahedral',
  };
  return { name: names[n] || `${n}-c`, planar, angles: uniq };
}

/** Ideal geometries used by the custom-block and chemistry geometry override. */
export const GEOMETRIES = {
  linear: [[1, 0, 0], [-1, 0, 0]],
  bent120: [[1, 0, 0], [-0.5, S3 / 2, 0]],
  trigonal: [0, 120, 240].map(a => polar(a)),
  square: [0, 90, 180, 270].map(a => polar(a)),
  tetrahedral: canonicalDirs([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].map(v3.norm)),
  hexagonal: [0, 60, 120, 180, 240, 300].map(a => polar(a)),
  octahedral: [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
  trigonalPrism: canonicalDirs((() => {
    const h = 1; const out = []; // arms 45° above/below the plane — the exact prism of stp (≈ acs)
    for (const a of [0, 120, 240]) for (const z of [h, -h]) { const p = polar(a); out.push(v3.norm([p[0], p[1], z])); }
    return out;
  })()),
  cubic: canonicalDirs([-1, 1].flatMap(a => [-1, 1].flatMap(b => [-1, 1].map(c => v3.norm([a, b, c]))))),
  cuboctahedral: canonicalDirs([[1, 1, 0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0], [1, 0, 1], [1, 0, -1], [-1, 0, 1], [-1, 0, -1], [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1]].map(v3.norm)),
};
function polar(deg) { const r = deg * Math.PI / 180; return clean([Math.cos(r), Math.sin(r), 0]); }

/**
 * Build a finite patch of the net for the assembly preview.
 * Returns vertex instances (position + neighbour directions) and edges inside the patch.
 */
export function buildPatch(an, reps) {
  const R = an.dim === 2 ? [reps, reps, 1] : [reps, reps, reps];
  const id = (a, b, c, i) => `${a},${b},${c},${i}`;
  const mk = (a, b, c, i) => ({ pos: v3.add(an.cart[i], frac2cart(an.cell, [a, b, c])), i, type: an.net.verts[i].t, dirs: an.nbrs[i].map(n => v3.norm(n.vec)) });
  const nbrKeys = k => { const [a, b, c, i] = k.split(',').map(Number); return an.nbrs[i].map(n => id(a + n.off[0], b + n.off[1], c + n.off[2], n.j)); };
  // core: the requested block of unit cells
  const verts = new Map();
  for (let a = 0; a < R[0]; a++) for (let b = 0; b < R[1]; b++) for (let c = 0; c < R[2]; c++)
    an.cart.forEach((_, i) => verts.set(id(a, b, c, i), mk(a, b, c, i)));
  // close rings that the cell boundary cuts: add outside vertices (one cell deep) bonded to ≥2 patch vertices
  const lo = an.dim === 2 ? [-1, -1, 0] : [-1, -1, -1], hi = an.dim === 2 ? [R[0], R[1], 0] : R;
  // (single pass against the core only, so dense nets don't snowball)
  const core = new Set(verts.keys());
  for (let a = lo[0]; a <= hi[0]; a++) for (let b = lo[1]; b <= hi[1]; b++) for (let c = lo[2]; c <= hi[2]; c++)
    an.cart.forEach((_, i) => {
      const k = id(a, b, c, i);
      if (!core.has(k) && nbrKeys(k).filter(x => core.has(x)).length >= 2) verts.set(k, mk(a, b, c, i));
    });
  // prune nodes that would hang by a single bond, then keep the largest connected piece
  for (let changed = true; changed;) {
    changed = false;
    for (const k of [...verts.keys()]) if (nbrKeys(k).filter(x => verts.has(x)).length < 2 && verts.size > 2) { verts.delete(k); changed = true; }
  }
  {
    const seen = new Set(); let bestComp = [];
    for (const s of verts.keys()) {
      if (seen.has(s)) continue;
      const comp = [s]; seen.add(s);
      for (let qi = 0; qi < comp.length; qi++) for (const x of nbrKeys(comp[qi])) if (verts.has(x) && !seen.has(x)) { seen.add(x); comp.push(x); }
      if (comp.length > bestComp.length) bestComp = comp;
    }
    const keep = new Set(bestComp);
    for (const k of [...verts.keys()]) if (!keep.has(k)) verts.delete(k);
  }
  const keys = [...verts.keys()];
  const index = new Map(keys.map((k, n) => [k, n]));
  const all = keys.map(k => verts.get(k));
  const edges = [];
  const adj = all.map(() => []);
  for (const [k, v] of verts) {
    const [a, b, c] = k.split(',').map(Number);
    for (const n of an.nbrs[v.i]) {
      const k2 = id(a + n.off[0], b + n.off[1], c + n.off[2], n.j);
      if (!verts.has(k2)) continue;
      const ia = index.get(k), ib = index.get(k2);
      adj[ia].push(ib);
      if (k < k2) edges.push({ a: ia, b: ib });
    }
  }
  // centre the patch
  const ctr = v3.mul(all.reduce((s, v) => v3.add(s, v.pos), [0, 0, 0]), 1 / all.length);
  for (const v of all) v.pos = v3.sub(v.pos, ctr);
  // 2-colour the finite patch (for magnet polarity when nodes connect directly)
  let bipartite = true;
  all.forEach(v => (v.colour = -1));
  for (let s = 0; s < all.length; s++) {
    if (all[s].colour !== -1) continue;
    all[s].colour = 0;
    const q = [s];
    while (q.length) {
      const u = q.shift();
      for (const w of adj[u]) {
        if (all[w].colour === -1) { all[w].colour = 1 - all[u].colour; q.push(w); }
        else if (all[w].colour === all[u].colour) bipartite = false;
      }
    }
  }
  return { verts: all, edges: edges.map(e => [all[e.a].pos, all[e.b].pos, e.a, e.b]), bipartite };
}
