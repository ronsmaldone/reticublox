// Reticublox — three.js viewer: exact printable mesh of one block, or an assembled patch of a net.
import * as THREE from 'three';
import { OrbitControls } from '../vendor/OrbitControls.js';
import { v3 } from './nets.js';
import { alignDirs, apply } from './rot.js';

export class Viewer {
  constructor(host) {
    this.host = host;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    host.appendChild(this.renderer.domElement);
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(40, 1, 0.5, 20000);
    this.camera.up.set(0, 0, 1);
    this.camera.position.set(120, -160, 140);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x445566, 1.6));
    const d1 = new THREE.DirectionalLight(0xffffff, 1.6); d1.position.set(1, -2, 3); this.scene.add(d1);
    const d2 = new THREE.DirectionalLight(0xffffff, 0.6); d2.position.set(-2, 1, 1); this.scene.add(d2);
    this.content = new THREE.Group();
    this.scene.add(this.content);
    this.grid = null;
    this.setTheme();
    new ResizeObserver(() => this.resize()).observe(host);
    this.resize();
    const loop = () => { this.controls.update(); this.renderer.render(this.scene, this.camera); requestAnimationFrame(loop); };
    loop();
  }

  setTheme() {
    const dark = matchMedia('(prefers-color-scheme: dark)').matches;
    this.scene.background = new THREE.Color(dark ? 0x1b1f24 : 0xf3f4f6);
    this.gridColours = dark ? [0x3a4048, 0x2a2f36] : [0xc9ced6, 0xe1e4e8];
  }

  resize() {
    const w = this.host.clientWidth, h = this.host.clientHeight;
    if (!w || !h) return;
    this.renderer.setSize(w, h);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
  }

  clear() {
    this.content.traverse(o => { o.geometry?.dispose(); });
    this.content.clear();
    if (this.grid) { this.scene.remove(this.grid); this.grid = null; }
  }

  addGrid(size, z = 0) {
    const n = Math.max(10, Math.round(size / 10));
    this.grid = new THREE.GridHelper(n * 10, n, ...this.gridColours);
    this.grid.rotation.x = Math.PI / 2;
    this.grid.position.z = z;
    this.scene.add(this.grid);
  }

  frame(box, keepView = false) {
    const c = box.getCenter(new THREE.Vector3());
    const r = box.getSize(new THREE.Vector3()).length() / 2 || 50;
    this.controls.target.copy(c);
    if (!keepView) {
      const dir = new THREE.Vector3(0.55, -0.75, 0.65).normalize();
      this.camera.position.copy(c).addScaledVector(dir, r * 3.3);
    }
    this.camera.near = r / 100; this.camera.far = r * 50;
    this.camera.updateProjectionMatrix();
  }

  /** Show the exact printable mesh (flat triangle array, mm). */
  showMesh(positions, colour = 0x3d4b5c) {
    this.clear();
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    g.computeVertexNormals();
    const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ color: colour, roughness: 0.55, metalness: 0.05, flatShading: false }));
    this.content.add(m);
    g.computeBoundingBox();
    const box = g.boundingBox;
    this.addGrid(Math.max(box.max.x - box.min.x, box.max.y - box.min.y) * 1.6, 0);
    this.frame(box);
  }

  /**
   * Assembled patch.
   * patch: { verts:[{pos, dirs, type, colour}], edges:[[p,q]] } in net units (edge = 1)
   * look: { edge: mm, nodeLook: {type: {arm, hub, bulges:[t], colour}}, linker: {len, bulges:[t]} | null, P }
   */
  showAssembly(patch, look) {
    this.clear();
    const P = look.P;
    const seg = 24;
    const mats = {};
    const mat = c => (mats[c] ??= new THREE.MeshStandardMaterial({ color: c, roughness: 0.5 }));
    const unitCyl = new THREE.CylinderGeometry(1, 1, 1, seg, 1, true); // along y
    const sphere = new THREE.SphereGeometry(1, seg, seg / 2);
    const cone = new THREE.CylinderGeometry(1, P.rodRadius / P.flareR, 1, seg, 1, false); // top r=1 (flare), bottom r=rod/flare
    const Y = new THREE.Vector3(0, 1, 0);
    const place = (geom, m, from, dir, len, rScale) => {
      const mesh = new THREE.Mesh(geom, m);
      const d = new THREE.Vector3(...dir);
      mesh.quaternion.setFromUnitVectors(Y, d);
      mesh.scale.set(rScale, len, rScale);
      mesh.position.set(...v3.add(from, v3.mul(dir, len / 2)));
      this.content.add(mesh);
    };
    const ball = (m, at, r) => { const s = new THREE.Mesh(sphere, m); s.position.set(...at); s.scale.setScalar(r); this.content.add(s); };
    const armMesh = (m, origin, dir, L, bulges, flip = false) => {
      const straight = Math.max(0.5, L - P.flareLength);
      place(unitCyl, m, origin, dir, straight, P.rodRadius);
      // flare
      const fStart = v3.add(origin, v3.mul(dir, straight));
      place(cone, m, fStart, dir, L - straight, P.flareR);
      for (const t of bulges) ball(m, v3.add(origin, v3.mul(dir, t)), P.bulgeRadius);
    };
    // molecular model posed at `origin` by rotation R: atoms, bonds, and fittings out to the arm tips
    const drawModel = (m, model, origin, R) => {
      const at = p => v3.add(origin, apply(R, p));
      const A = model.atoms.map(a => at(a.p));
      model.atoms.forEach((a, k) => ball(m, A[k], Math.max(a.r, P.rodRadius)));
      for (const [i, j] of model.bonds) {
        const d = v3.sub(A[j], A[i]);
        place(unitCyl, m, A[i], v3.norm(d), v3.len(d), P.rodRadius);
      }
      for (const arm of model.arms) {
        const d = apply(R, arm.dir), tip = at(arm.tip);
        const F = v3.sub(tip, v3.mul(d, P.flareLength));
        for (const from of arm.from) { // several when one carboxylate bridges two metals
          const f = at(from), r = v3.sub(F, f);
          place(unitCyl, m, f, v3.norm(r), v3.len(r), P.rodRadius);
        }
        ball(m, F, P.rodRadius);
        place(cone, m, F, d, P.flareLength, P.flareR);
      }
    };
    const pose = (model, dirs) => {
      const al = alignDirs(model.arms.map(a => a.dir), dirs, { up: [0, 0, 1] });
      return al && al.err < 0.2 ? al.R : null;
    };
    const E = look.edge;
    for (const v of patch.verts) {
      const nl = look.nodeLook[v.type];
      if (!nl) continue;
      const pos = v3.mul(v.pos, E);
      const m = mat(look.direct && v.colour === 1 ? nl.colourB : nl.colour);
      const R = nl.model && pose(nl.model, v.dirs);
      if (R) { drawModel(m, nl.model, pos, R); continue; }
      ball(m, pos, Math.max(nl.hub, P.rodRadius));
      for (const d of v.dirs) armMesh(m, pos, d, nl.arm, nl.bulges);
    }
    if (look.linker) {
      const m = mat(look.linker.colour);
      for (const [p, q] of patch.edges) {
        const a = v3.mul(p, E), b = v3.mul(q, E);
        const mid = v3.mul(v3.add(a, b), 0.5);
        const d = v3.norm(v3.sub(b, a));
        const half = look.linker.len / 2;
        const R = look.linker.model && pose(look.linker.model, [d, v3.mul(d, -1)]);
        if (R) { drawModel(m, look.linker.model, mid, R); continue; }
        if (look.linker.hub) ball(m, mid, look.linker.hub);
        else ball(m, mid, P.rodRadius);
        armMesh(m, mid, d, half, look.linker.bulges);
        armMesh(m, mid, v3.mul(d, -1), half, look.linker.bulges);
      }
    }
    const box = new THREE.Box3().setFromObject(this.content);
    if (look.dim === 2) this.addGrid(Math.max(box.max.x - box.min.x, box.max.y - box.min.y) * 1.3, box.min.z);
    this.frame(box);
  }

  snapshot() { return this.renderer.domElement.toDataURL('image/png'); }
}
