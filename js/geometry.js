// Reticublox — solid modelling of node / linker blocks with Manifold (robust mesh booleans).
//
// A block "spec" is the single data model shared by every way of making blocks:
// {
//   name, kind: 'node' | 'linker',
//   hubRadius: number | 0,            // sphere at the centre (mm), 0 = none
//   arms: [{ dir:[x,y,z] (unit), length: mm from centre to magnet face,
//            bulges: [{ t: mm from centre, r: mm }] }],
//   polarity: 'N' | 'S'               // which magnet pole faces outwards
// }
import Module from '../vendor/manifold.js';
import { v3 } from './nets.js';

let wasm = null;
export async function initManifold() {
  if (wasm) return wasm;
  wasm = await Module({ locateFile: () => new URL('../vendor/manifold.wasm', import.meta.url).href });
  wasm.setup();
  return wasm;
}

export const DEFAULT_PRINT = {
  rodRadius: 4.5,      // mm
  flareRadius: 7.0,    // mm (end "trumpet" that holds the magnet)
  flareLength: 9.0,    // mm
  hubRadius: 8.0,      // default node hub sphere (mm)
  bulgeRadius: 6.5,    // ring / bulge sphere (mm)
  magnetDiameter: 6.0, // mm
  magnetThickness: 3.0,// mm
  clearance: 0.15,     // mm added to pocket radius and depth
  flatBase: 0.7,       // mm shaved off the underside of planar blocks (0 = off)
  atomRadius: 5.5,     // mm, atom spheres in molecular models (bonds use the rod radius)
  segments: 64,
};

/** Rotation (column-major 4x4 for Manifold) taking +z onto `d`, then translating by `t`. */
function alignZ(d, t = [0, 0, 0]) {
  const z = v3.norm(d);
  const helper = Math.abs(z[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
  const x = v3.norm(v3.cross(helper, z));
  const y = v3.cross(z, x);
  return [x[0], x[1], x[2], 0, y[0], y[1], y[2], 0, z[0], z[1], z[2], 0, t[0], t[1], t[2], 1];
}

/** Check a spec against print parameters; returns a list of human-readable warnings. */
export function checkSpec(spec, P) {
  const w = [];
  const pocketR = P.magnetDiameter / 2 + P.clearance;
  const flareR = effectiveFlareRadius(P);
  if (flareR > P.flareRadius + 1e-6)
    w.push(`Flare radius raised to ${flareR.toFixed(1)} mm so the magnet pocket keeps a 1 mm wall.`);
  if (P.flareLength < P.magnetThickness + P.clearance + 1)
    w.push('Flare is shorter than the magnet pocket — increase flare length.');
  const hub = spec.hubRadius || P.rodRadius;
  if (spec.model) {
    spec.arms.forEach((a, i) => {
      for (const A of fitAtoms(spec.model, i)) {
        const f = v3.len(v3.sub(armTip(a, a.length - P.flareLength), A.p));
        if (f < atomRadius(A.el, P) * 0.5) { w.push(`Magnet fitting ${i + 1} runs into its atom — press “Rebuild model & parts” after changing print settings.`); break; }
      }
    });
    const { atoms, bonds } = spec.model;
    const minBond = Math.min(...bonds.map(([i, j]) => v3.len(v3.sub(atoms[i].p, atoms[j].p))));
    if (bonds.length && minBond < 2 * (P.atomRadius ?? 5.5) + 1)
      w.push(`Atoms touch (bonds ${minBond.toFixed(1)} mm long): rings print nearly closed. Molecular models look best at 9–12 mm per Å, or use a smaller atom radius.`);
  } else spec.arms.forEach((a, i) => {
    const straight = a.length - P.flareLength;
    if (straight < hub * 0.6) w.push(`Arm ${i + 1} is short (${a.length.toFixed(1)} mm): the flare runs into the hub. Increase the scale or arm length.`);
    for (const b of a.bulges) {
      if (b.t + b.r > a.length - P.flareLength * 0.4) w.push(`A bulge on arm ${i + 1} overlaps the magnet flare.`);
    }
  });
  if (isPlanar(spec) && P.flatBase > 0 && !spec.model) {
    const floor = P.rodRadius - P.flatBase;
    if (floor < pocketR + 0.6) w.push(`Flat-base cut (${P.flatBase} mm) would thin or open the magnet pocket floor. Use ≤ ${Math.max(0, P.rodRadius - pocketR - 0.6).toFixed(1)} mm, or a thicker rod.`);
  }
  if (spec.model && !isPlanar(spec) && flatCutZ(spec, P) !== null) w.push('Flat molecule: prints lying flat without supports; the magnet fittings tilt slightly to meet the net.');
  else if (!isPlanar(spec)) w.push(spec.model ? '3D molecular model: print with tree supports.' : '3D block: print with supports (tree supports work well). One arm points straight up.');
  return w;
}

export function effectiveFlareRadius(P) {
  return Math.max(P.flareRadius, P.magnetDiameter / 2 + P.clearance + 1.0, P.rodRadius);
}

/** Point at distance t along an arm. Arms start at the centre, or (crystal pieces) at their atom `at`. */
export const armTip = (a, t = a.length) => v3.add(a.at || [0, 0, 0], v3.mul(a.dir, t));

export const isPlanar = spec => spec.arms.every(a => Math.abs(a.dir[2]) < 1e-3 && Math.abs(armTip(a)[2]) < 1e-3) &&
  (!spec.model || spec.model.atoms.every(a => Math.abs(a.p[2]) < 1e-3));

/**
 * Height of the flat-base cut, or null if the piece can't lie flat. Blocks: planar ones only.
 * Molecular models: whenever every atom lies in one horizontal plane (a flat molecule) and the cut
 * leaves every magnet pocket closed, even if the fittings tilt out of that plane (BTC in tbo).
 */
export function flatCutZ(spec, P) {
  const r = flatPrint(spec, P);
  return r ? r.cut : null;
}
/** { cut, flip }: flip = print upside down (rotate 180° about x) so tilted fittings point up. */
export function flatPrint(spec, P) {
  if (P.flatBase <= 0) return null;
  if (!spec.model) return isPlanar(spec) ? { cut: -(P.rodRadius - P.flatBase), flip: false } : null;
  const zs0 = spec.model.atoms.map(a => a.p[2]);
  if (Math.max(...zs0) - Math.min(...zs0) > 1e-3) return null;
  const flip = spec.arms.reduce((t, a) => t + armTip(a)[2], 0) < zs0[0] * spec.arms.length - 1e-6;
  const sg = flip ? -1 : 1;
  const z0 = sg * zs0[0];
  const cut = z0 - (P.rodRadius - P.flatBase);
  const pocketR = P.magnetDiameter / 2 + P.clearance, pocketH = P.magnetThickness + P.clearance;
  for (const a of spec.arms) {
    const rad = pocketR * Math.sqrt(Math.max(0, 1 - a.dir[2] ** 2));
    if (Math.min(sg * armTip(a)[2], sg * armTip(a, a.length - pocketH)[2]) - rad < cut + 0.6) return null;
  }
  return { cut, flip };
}

/** Atoms a fitting starts from: one connection atom, or several metals bridged by one carboxylate. */
export const fitAtoms = (model, k) => [].concat(model.fittings[k]).map(i => model.atoms[i]);

/** Atom sphere radius in molecular models: second-row atoms at the set size, heavier ones a bit larger. */
export const atomRadius = (el, P) => (P.atomRadius ?? 5.5) * (['C', 'N', 'O', 'B', 'F', 'H'].includes(el) ? 1 : 1.15);

/** Build the printable solid. Returns a Manifold (caller owns it). */
export function buildSolid(spec, P) {
  const { Manifold } = wasm;
  const seg = P.segments;
  const flareR = effectiveFlareRadius(P);
  const pocketR = P.magnetDiameter / 2 + P.clearance;
  const pocketH = P.magnetThickness + P.clearance;
  const parts = [];
  if (spec.model) return finish(buildModelParts(spec, P, parts, seg, flareR), spec, P, pocketR, pocketH, seg);

  parts.push(Manifold.sphere(Math.max(spec.hubRadius || 0, P.rodRadius), seg));
  for (const arm of spec.arms) {
    const straight = Math.max(0.5, arm.length - P.flareLength);
    // rod runs 1 mm into the flare so the union is one watertight shell (no face-to-face contact)
    parts.push(Manifold.cylinder(straight + 1, P.rodRadius, P.rodRadius, seg).transform(alignZ(arm.dir)));
    const flareStart = v3.mul(arm.dir, straight);
    parts.push(Manifold.cylinder(arm.length - straight, P.rodRadius, flareR, seg).transform(alignZ(arm.dir, flareStart)));
    for (const b of arm.bulges) parts.push(Manifold.sphere(b.r, seg).translate(v3.mul(arm.dir, b.t)));
  }
  return finish(parts, spec, P, pocketR, pocketH, seg);
}

/** Molecular model: atom spheres, bond rods, and a rod + flare from each connection atom to its arm tip. */
function buildModelParts(spec, P, parts, seg, flareR) {
  const { Manifold } = wasm;
  const { atoms, bonds } = spec.model;
  const aseg = Math.max(24, Math.round(seg * 0.75));
  for (const a of atoms) parts.push(Manifold.sphere(Math.max(atomRadius(a.el, P), P.rodRadius), aseg).translate(a.p));
  for (const [i, j] of bonds) {
    const d = v3.sub(atoms[j].p, atoms[i].p);
    parts.push(Manifold.cylinder(v3.len(d), P.rodRadius, P.rodRadius, aseg).transform(alignZ(d, atoms[i].p)));
  }
  // fitting: rod from the atom to the flare start, then a flare along the net edge (so magnet faces
  // meet squarely); a sphere rounds the elbow when the rod has to bend to reach the edge
  spec.arms.forEach((arm, k) => {
    const F = armTip(arm, arm.length - P.flareLength);
    for (const { p: A } of fitAtoms(spec.model, k)) {
      const r = v3.sub(F, A);
      parts.push(Manifold.cylinder(v3.len(r) + 0.5, P.rodRadius, P.rodRadius, seg).transform(alignZ(r, A)));
    }
    parts.push(Manifold.sphere(P.rodRadius, seg).translate(F));
    parts.push(Manifold.cylinder(P.flareLength, P.rodRadius, flareR, seg).transform(alignZ(arm.dir, F)));
  });
  return parts;
}

/** Union the parts, cut the magnet pockets, flatten the base of planar pieces, sit on z = 0. */
function finish(parts, spec, P, pocketR, pocketH, seg) {
  const { Manifold } = wasm;
  let solid = Manifold.union(parts);
  parts.forEach(p => p.delete());

  // magnet pockets, along each arm (net edge) direction
  const pockets = spec.arms.map(arm =>
    Manifold.cylinder(pocketH + 1, pocketR, pocketR, seg).transform(alignZ(arm.dir, armTip(arm, arm.length - pocketH))));
  const pocketU = Manifold.union(pockets);
  pockets.forEach(p => p.delete());
  let s2 = solid.subtract(pocketU);
  solid.delete(); pocketU.delete(); solid = s2;

  // flat base for pieces that lie flat
  const fp = flatPrint(spec, P);
  if (fp) {
    if (fp.flip) { s2 = solid.rotate([180, 0, 0]); solid.delete(); solid = s2; }
    s2 = solid.trimByPlane([0, 0, 1], fp.cut);
    solid.delete(); solid = s2;
  }
  // sit on the build plate
  const bb = solid.boundingBox();
  s2 = solid.translate([0, 0, -bb.min[2]]);
  solid.delete();
  return s2;
}

/** Manifold -> flat triangle arrays usable by three.js and the STL writer. */
export function solidToMesh(solid) {
  const m = solid.getMesh();
  const np = m.numProp;
  const pos = new Float32Array(m.triVerts.length * 3);
  for (let i = 0; i < m.triVerts.length; i++) {
    const v = m.triVerts[i];
    pos[i * 3] = m.vertProperties[v * np];
    pos[i * 3 + 1] = m.vertProperties[v * np + 1];
    pos[i * 3 + 2] = m.vertProperties[v * np + 2];
  }
  return { positions: pos, volume: solid.volume(), triangles: m.triVerts.length / 3, bbox: solid.boundingBox() };
}

/** Binary STL from a flat triangle position array. */
export function meshToSTL(positions, name = 'reticublox') {
  const nTri = positions.length / 9;
  const buf = new ArrayBuffer(84 + nTri * 50);
  const dv = new DataView(buf);
  const header = `Reticublox ${name}`.slice(0, 79);
  for (let i = 0; i < header.length; i++) dv.setUint8(i, header.charCodeAt(i));
  dv.setUint32(80, nTri, true);
  let o = 84;
  for (let t = 0; t < nTri; t++) {
    const p = i => [positions[t * 9 + i * 3], positions[t * 9 + i * 3 + 1], positions[t * 9 + i * 3 + 2]];
    const a = p(0), b = p(1), c = p(2);
    const n = v3.norm(v3.cross(v3.sub(b, a), v3.sub(c, a)));
    for (const x of [n, a, b, c]) { dv.setFloat32(o, x[0], true); dv.setFloat32(o + 4, x[1], true); dv.setFloat32(o + 8, x[2], true); o += 12; }
    dv.setUint16(o, 0, true); o += 2;
  }
  return buf;
}

// ---------- spec factories ----------

/** Evenly spaced bulges along a segment [t0, t1]. */
export function spacedBulges(count, t0, t1, r) {
  const out = [];
  for (let k = 1; k <= count; k++) out.push({ t: t0 + (t1 - t0) * k / (count + 1), r });
  return out;
}

/** Generic node from a set of unit directions. */
export function makeNodeSpec({ name, dirs, armLength, hubRadius, bulgesPerArm = 0, bulgeRadius, polarity = 'N', flareLength }) {
  return {
    name, kind: 'node', hubRadius, polarity,
    arms: dirs.map(d => ({ dir: d, length: armLength, bulges: spacedBulges(bulgesPerArm, hubRadius, armLength - flareLength, bulgeRadius) })),
  };
}

/** Straight ditopic linker of total magnet-face-to-face length L with n ring bulges. */
export function makeLinkerSpec({ name, length, rings = 1, bulgeRadius, polarity = 'S', spacing }) {
  const half = length / 2;
  // ring centres spaced like a para-oligophenylene: `spacing` mm apart, centred on the middle
  const sp = spacing ?? Math.min(bulgeRadius * 2.2, length / (rings + 1));
  const ts = [];
  for (let k = 0; k < rings; k++) ts.push((k - (rings - 1) / 2) * sp);
  const centre = ts.find(t => Math.abs(t) < 1e-6) !== undefined;
  const arm = sign => ({
    dir: [sign, 0, 0], length: half,
    bulges: ts.filter(t => sign * t > 1e-6).map(t => ({ t: Math.abs(t), r: bulgeRadius })),
  });
  return { name, kind: 'linker', hubRadius: centre ? bulgeRadius : 0, polarity, arms: [arm(1), arm(-1)] };
}
