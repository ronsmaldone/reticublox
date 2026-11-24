// Reticublox — small 3×3 rotation helpers (row-major [[r00,r01,r02],[r10..],[r20..]]).
import { v3 } from './nets.js';

export const apply = (R, p) => [
  R[0][0] * p[0] + R[0][1] * p[1] + R[0][2] * p[2],
  R[1][0] * p[0] + R[1][1] * p[1] + R[1][2] * p[2],
  R[2][0] * p[0] + R[2][1] * p[1] + R[2][2] * p[2],
];
export const mulM = (A, B) => A.map((row, i) => [0, 1, 2].map(j => row[0] * B[0][j] + row[1] * B[1][j] + row[2] * B[2][j]));
export const IDENT = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const ang = (a, b) => Math.acos(Math.max(-1, Math.min(1, v3.dot(a, b))));

/** Smallest rotation taking unit vector u onto unit vector w. */
export function rotBetween(u, w) {
  const c = v3.dot(u, w);
  if (c > 1 - 1e-12) return IDENT;
  let k = v3.cross(u, w);
  if (v3.len(k) < 1e-9) { // antiparallel: 180° about any perpendicular
    k = v3.norm(Math.abs(u[0]) < 0.9 ? v3.cross(u, [1, 0, 0]) : v3.cross(u, [0, 1, 0]));
    return axisAngle(k, Math.PI);
  }
  return axisAngle(v3.norm(k), Math.acos(Math.max(-1, Math.min(1, c))));
}
export function axisAngle(k, t) {
  const [x, y, z] = k, c = Math.cos(t), s = Math.sin(t), C = 1 - c;
  return [
    [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
    [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
    [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
  ];
}
/** Orthonormal frame (as matrix columns) from two non-parallel vectors. */
function frame(a, b, flip = 1) {
  const x = v3.norm(a);
  const y = v3.mul(v3.norm(v3.sub(b, v3.mul(x, v3.dot(b, x)))), flip);
  const z = v3.cross(x, y);
  return [[x[0], y[0], z[0]], [x[1], y[1], z[1]], [x[2], y[2], z[2]]];
}
const T = M => [0, 1, 2].map(i => [0, 1, 2].map(j => M[j][i]));

/**
 * Best rotation R (proper unless `improper` allowed) and assignment perm so that R·A[i] ≈ B[perm[i]].
 * `up` (optional) fixes the roll when every A is collinear: R maps [0,0,1] towards `up`.
 * Returns { R, perm, err (max angle, rad) }.
 */
export function alignDirs(A, B, { improper = true, up = null } = {}) {
  const n = A.length;
  const i1 = A.findIndex(a => Math.abs(v3.dot(a, A[0])) < 0.97);
  if (i1 < 0 || B.every(b => Math.abs(v3.dot(b, B[0])) >= 0.97)) {
    // collinear set (linear linker): A[0] → B[0], roll from `up`
    let R = rotBetween(A[0], B[0]);
    const z0 = apply(R, [0, 0, 1]);
    let want = up || [0, 0, 1];
    want = v3.sub(want, v3.mul(B[0], v3.dot(want, B[0])));
    if (v3.len(want) < 1e-6) want = Math.abs(B[0][0]) < 0.9 ? v3.cross(B[0], [1, 0, 0]) : v3.cross(B[0], [0, 1, 0]);
    want = v3.norm(want);
    const zp = v3.norm(v3.sub(z0, v3.mul(B[0], v3.dot(z0, B[0]))));
    const roll = Math.atan2(v3.dot(v3.cross(zp, want), B[0]), v3.dot(zp, want));
    R = mulM(axisAngle(B[0], roll), R);
    return { R, perm: assign(A.map(a => apply(R, a)), B), err: 0 };
  }
  const a01 = ang(A[0], A[i1]);
  let best = null;
  for (const flip of improper ? [1, -1] : [1]) {
    const FA = frame(A[0], A[i1]);
    for (let j0 = 0; j0 < n; j0++) for (let j1 = 0; j1 < n; j1++) {
      if (j0 === j1 || Math.abs(ang(B[j0], B[j1]) - a01) > 0.35) continue;
      const R = mulM(frame(B[j0], B[j1], flip), T(FA));
      const RA = A.map(a => apply(R, a));
      const perm = assign(RA, B);
      if (!perm) continue;
      const err = Math.max(...RA.map((a, i) => ang(a, B[perm[i]])));
      if (!best || err < best.err - 1e-9) best = { R, perm, err };
    }
    if (best && best.err < 0.03) break; // a proper rotation fits — don't mirror
  }
  return best;
}
/** Greedy one-to-one nearest assignment. */
function assign(RA, B) {
  const used = new Set(), perm = [];
  const pairs = [];
  RA.forEach((a, i) => B.forEach((b, j) => pairs.push([v3.dot(a, b), i, j])));
  pairs.sort((x, y) => y[0] - x[0]);
  for (const [, i, j] of pairs) if (perm[i] === undefined && !used.has(j)) { perm[i] = j; used.add(j); }
  return perm.length === RA.length && !perm.includes(undefined) ? perm : null;
}
