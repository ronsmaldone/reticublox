// Reticublox — application: a four-step builder
//   1. MOF or COF   2. 2D or 3D   3. topology (by net or space group)   4. building blocks
// Step 4 has one card per node type (+ the linker). Each card is Generic (shape from the
// topology), Molecule (preset / SMILES / drawn) or SBU (metal cluster, MOFs only).
import { NETS, analyseNet, buildPatch, describeGeometry, GEOMETRIES, v3 } from './nets.js';
import {
  initManifold, DEFAULT_PRINT, buildSolid, solidToMesh, meshToSTL, checkSpec, armTip,
  makeNodeSpec, makeLinkerSpec, effectiveFlareRadius, atomRadius,
} from './geometry.js';
import {
  initRDKit, GROUPS, PRESET_LIST, SBUS, sbuCN, parseMolecule, findConnections, autoDetectGroup,
  moleculeToBlock, chemBlockToSpec, sbuToSpec, sbuModelSpec, moleculeModelSpec, prepare3D, moleculeFromGraph,
} from './chem.js';
import { importStructure, importFramework, IMPORT_ACCEPT } from './import3d.js';
import { Viewer } from './viewer.js';

const $ = s => document.querySelector(s);
const $$ = s => [...document.querySelectorAll(s)];
const NODE_COLOURS = [0x5a6a80, 0x7a4fb5, 0x2f7d5b, 0xb5612f];
const NODE_COLOURS_B = [0x8b97a8, 0xb79be0, 0x86c7a8, 0xe0a47a];
const LINKER_COLOUR = 0x4f86d6;
const GEOM_NOTE_NETS = ['tbo', 'ctn', 'bor', 'pts', 'pts_t', 'stp'];
const GEOM_LABELS = {
  linear: 'Linear', bent120: 'Bent 120°', trigonal: 'Trigonal planar', square: 'Square planar', tetrahedral: 'Tetrahedral',
  hexagonal: 'Hexagonal planar', octahedral: 'Octahedral', trigonalPrism: 'Trigonal prism', cubic: 'Cubic (8-c)', cuboctahedral: 'Cuboctahedral (12-c)',
};

// ---------------- state ----------------
const store = {
  get(k, d) { try { const v = localStorage.getItem('reticublox:' + k); return v ? JSON.parse(v) : d; } catch { return d; } },
  set(k, v) { try { localStorage.setItem('reticublox:' + k, JSON.stringify(v)); } catch { /* private mode */ } },
};
const freshWiz = () => ({ cls: null, dim: null, net: null, browse: 'net', step: 1, mode: 'linker', scale: 6, slots: {}, open: null });
const state = {
  P: { ...DEFAULT_PRINT, ...store.get('print', {}) },
  parts: store.get('parts', []),
  selected: null,
  preview: null,
  view: 'block',
  asm: store.get('asm', null),
  wiz: { ...freshWiz(), ...store.get('wiz', {}) },
  built: false,
  builder: store.get('builder', 'scratch'), // 'scratch' (steps 1–4) | 'cif' (from a crystal structure)
  cif: store.get('cif', null),               // { name, fw, scale, sel }
  cifError: '',
};
const W = state.wiz;
let nextId = Math.max(0, ...state.parts.map(p => p.id)) + 1;
const meshCache = new Map();
const netCache = new Map();
const molCache = new Map();
let viewer, rdkitReady = false;

const save = () => { store.set('parts', state.parts); store.set('print', state.P); store.set('asm', state.asm); store.set('wiz', W); store.set('builder', state.builder); store.set('cif', state.cif); };
const status = t => ($('#status').textContent = t);
const fmt = (x, d = 1) => Number(x).toFixed(d);
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const net = k => { if (!netCache.has(k)) netCache.set(k, analyseNet(k)); return netCache.get(k); };

// ---------------- printable-ness ----------------
function ensurePrintable(spec) {
  if (spec.model) return spec; // molecular models size their fittings themselves
  const P = state.P;
  const notes = [];
  const minLen = Math.max(spec.hubRadius || 0, P.rodRadius) + P.flareLength + 2;
  for (const a of spec.arms) {
    if (a.length < minLen) { notes.push(`Arm lengthened from ${fmt(a.length)} to ${fmt(minLen)} mm (printable minimum).`); a.length = minLen; }
    a.bulges = a.bulges.filter(b => b.t + b.r * 0.5 < a.length - P.flareLength * 0.5);
  }
  spec.notes = [...new Set([...(spec.notes || []), ...notes])];
  return spec;
}
function restyle(spec) {
  const P = state.P;
  if (spec.kind === 'node' && spec.hubRadius) spec.hubRadius = P.hubRadius;
  if (spec.kind === 'linker' && spec.hubRadius) spec.hubRadius = P.bulgeRadius;
  for (const a of spec.arms) for (const b of a.bulges) b.r = P.bulgeRadius;
}
function meshFor(spec, key) {
  const k = key + JSON.stringify(state.P);
  if (key && meshCache.has(k)) return meshCache.get(k);
  const solid = buildSolid(spec, state.P);
  const m = solidToMesh(solid);
  solid.delete();
  if (key) meshCache.set(k, m);
  return m;
}
const allAngles = dirs => {
  const out = [];
  for (let i = 0; i < dirs.length; i++) for (let j = i + 1; j < dirs.length; j++)
    out.push(Math.acos(Math.max(-1, Math.min(1, v3.dot(dirs[i], dirs[j])))) * 180 / Math.PI);
  return out.sort((a, b) => a - b).map(a => Math.round(a));
};

// =====================================================================
//  STEP 4 — slots (one per vertex type, plus the edge / linker)
// =====================================================================
const slotKeys = () => {
  if (!W.net) return [];
  const keys = Object.keys(net(W.net).types).map(l => 'v:' + l);
  if (W.mode === 'linker') keys.push('edge');
  return keys;
};
function slotInfo(key) {
  if (key === 'edge') return { cn: 2, dirs: [[1, 0, 0], [-1, 0, 0]], geometry: { name: 'linear', planar: true }, count: net(W.net).edgesPerCell, kind: 'linker' };
  const t = net(W.net).types[key.slice(2)];
  return { ...t, kind: 'node' };
}
function slotTitle(key) {
  if (key === 'edge') return 'Linker';
  const n = Object.keys(net(W.net).types).length;
  return n > 1 ? `Node ${key.slice(2)}` : 'Node';
}
function slot(key) {
  W.slots[key] ??= key === 'edge'
    ? { source: 'generic', generic: { length: 54, rings: 1 } }
    : { source: 'generic', generic: { arm: null, bulges: 0 }, sbuGeom: 'net' };
  return W.slots[key];
}
const defaultArm = () => (W.mode === 'direct' ? 45 : 28);

function getMol(key) {
  const s = slot(key);
  if (!s.mol?.smiles || !rdkitReady) return null;
  const src = s.mol.molblock || s.mol.smiles;
  const c = molCache.get(key);
  if (c && c.src === src) return c.mm;
  const mm = parseMolecule(s.mol.smiles, s.mol.molblock || null);
  molCache.set(key, { src, mm });
  return mm;
}

/** Build the (unsaved) spec for a slot and judge how well it fits the net. */
function slotSpec(key) {
  const s = slot(key), t = slotInfo(key), P = state.P;
  const label = key === 'edge' ? '' : key.slice(2);
  const multi = Object.keys(net(W.net).types).length > 1;
  let spec = null, cb = null;
  try {
    if (s.source === 'generic') {
      spec = key === 'edge'
        ? makeLinkerSpec({ name: `${W.net} linker (${s.generic.rings} ring${s.generic.rings === 1 ? '' : 's'})`, length: +s.generic.length, rings: +s.generic.rings, bulgeRadius: P.bulgeRadius })
        : makeNodeSpec({ name: `${W.net} node${multi ? ' ' + label : ''} (${t.cn}-c)`, dirs: t.dirs, armLength: +(s.generic.arm ?? defaultArm()), hubRadius: P.hubRadius, bulgesPerArm: +s.generic.bulges, bulgeRadius: P.bulgeRadius, flareLength: P.flareLength });
    } else if (s.source === 'sbu') {
      if (!s.sbu) return { error: 'Choose a metal cluster.' };
      const fits = sbuCN(s.sbu) === t.cn;
      const dirs = s.sbuGeom === 'net' && fits ? t.dirs : null;
      spec = s.sbuStyle === 'model' ? sbuModelSpec(s.sbu, { scale: W.scale, P, dirs }) : sbuToSpec(s.sbu, { scale: W.scale, P, dirs });
    } else {
      if (!s.mol?.smiles) return { error: 'Choose, type or draw a molecule.' };
      if (!rdkitReady) return { error: 'Chemistry toolkit still loading…' };
      const mm = getMol(key);
      if (s.mol.conns.length < 2) return { error: 'Mark at least two connection points (click atoms in the drawing).' };
      const g = s.mol.geom;
      const override = g === 'net' ? (s.mol.conns.length === t.cn ? t.dirs : '') : g;
      cb = moleculeToBlock(mm, s.mol.conns, { ext: +s.mol.ext, geometryOverride: override });
      spec = s.mol.style === 'model'
        ? moleculeModelSpec(mm, s.mol.conns, cb, { group: s.mol.group, ext: +s.mol.ext, scale: W.scale, P, targets: cb.imposed ? cb.dirs : null, name: shortName(s.mol.name) + ' (model)', kind: t.kind })
        : chemBlockToSpec(cb, { name: shortName(s.mol.name), kind: t.kind, scale: W.scale, P });
    }
  } catch (e) { return { error: e.message }; }
  spec.kind = t.kind;
  ensurePrintable(spec);
  return { spec, fit: judgeFit(spec, t, key === 'edge'), cb };
}
/** How well a block's arms match a net site (t: { cn, dirs }). */
function judgeFit(spec, t, isEdge, hint = ' — use “Match the net”') {
  if (spec.arms.length !== t.cn) return { ok: false, msg: `has ${spec.arms.length} connections, needs ${t.cn}` };
  const radial = spec.arms.map(a => v3.norm(armTip(a)));
  const A = allAngles(radial), B = allAngles(t.dirs);
  // pairwise angles miss "flat vs pyramidal", so also compare how far arms tilt out of plane
  const dev = Math.max(...A.map((x, i) => Math.abs(x - B[i])), Math.round(Math.abs(tilt(radial) - tilt(t.dirs))));
  if (spec.model?.bend > 20) return { ok: false, msg: `molecule shape is far from the net — magnet fittings bend by ${spec.model.bend}°` };
  if (dev > 8) return { ok: false, msg: isEdge ? 'bent — the net needs a straight linker' : `angles off by up to ${dev}°${hint}` };
  if (dev > 2) return { ok: true, msg: `fits (angles within ${dev}°)` };
  if (spec.model?.bend > 3) return { ok: true, msg: `fits — molecule kept flat, fittings bend ${spec.model.bend}° to meet the net` };
  return { ok: true, msg: 'fits the net' };
}
/** Mean out-of-plane tilt (degrees) of arms relative to the plane normal to their resultant. */
function tilt(dirs) {
  const sum = dirs.reduce((s, d) => v3.add(s, d), [0, 0, 0]);
  const meanTilt = n => dirs.reduce((s, d) => s + Math.asin(Math.min(1, Math.abs(v3.dot(d, n)))), 0) / dirs.length * 180 / Math.PI;
  if (v3.len(sum) >= 0.3) return meanTilt(v3.norm(sum)); // clearly pyramidal / 3D
  // (nearly) flat: tilt out of the best plane through the arms (a slightly buckled node from a crystal)
  let best = 0;
  for (let i = 0; i < dirs.length; i++) for (let j = i + 1; j < dirs.length; j++) {
    const c = v3.cross(dirs[i], dirs[j]);
    if (v3.len(c) < 0.2) continue;
    const t = meanTilt(v3.norm(c));
    if (best === 0 || t < best) best = t;
  }
  return best;
}
const shortName = s => { const m = String(s).match(/\(([^)]+)\)\s*$/); return (m ? m[1] : String(s)).slice(0, 40); };

// =====================================================================
//  Wizard rendering
// =====================================================================
const CLS = { mof: 'MOF', cof: 'COF' };
const reachable = n => n === 1 || (n === 2 && W.cls) || (n === 3 && W.cls && W.dim) || (n === 4 && W.cls && W.dim && W.net);

function stepShell(n, title, summary, body) {
  const open = W.step === n, can = reachable(n);
  return `<section class="step ${open ? 'open' : ''} ${can ? '' : 'locked'}" data-step="${n}">
    <header data-act="goto" data-step="${n}"><span class="num">${n}</span><span class="stitle">${title}</span>
      <span class="ssum">${open ? '' : (summary || '')}</span>${!open && can ? '<span class="chg">Change</span>' : ''}</header>
    ${open ? `<div class="sbody">${body}</div>` : ''}</section>`;
}

function step1() {
  const card = (v, t, d) => `<button class="card ${W.cls === v ? 'on' : ''}" data-act="cls" data-v="${v}"><b>${t}</b><span>${d}</span></button>`;
  return stepShell(1, 'Framework type', W.cls ? CLS[W.cls] : '', `<div class="cards">
    ${card('mof', 'MOF', 'Metal–organic framework. Metal clusters (SBUs) joined by organic linkers through coordination bonds.')}
    ${card('cof', 'COF', 'Covalent organic framework. Organic nodes and linkers joined by covalent bonds (imine, boronate ester, boroxine…).')}
  </div>`);
}
function step2() {
  const card = (v, t, d) => `<button class="card ${W.dim === v ? 'on' : ''}" data-act="dim" data-v="${v}"><b>${t}</b><span>${d}</span></button>`;
  const ex = W.cls === 'mof' ? ['Cu₃(HHTP)₂', 'MOF-5, HKUST-1, UiO-66'] : ['COF-1, COF-5, TAPB–PDA', 'COF-300, COF-102, COF-105'];
  return stepShell(2, 'Dimensionality', W.dim ? `${W.dim}D` : '', `<div class="cards">
    ${card(2, '2D', `Layers and sheets, e.g. ${ex[0]}. Blocks print flat without supports.`)}
    ${card(3, '3D', `Extended frameworks, e.g. ${ex[1]}. Non-planar nodes print with supports.`)}
  </div>`);
}
function cnLabel(k) { return [...new Set(Object.values(net(k).types).map(t => t.cn))].sort((a, b) => a - b).join(',') + '-c'; }
function netRow(k) {
  const n = NETS[k], ex = n[W.cls];
  return `<button class="netrow ${W.net === k ? 'on' : ''} ${ex ? '' : 'dim'}" data-act="net" data-v="${k}">
    <span class="nk">${k.replace('_t', '')}</span><span class="nt">${n.title}</span><span class="nsg">${n.sg} · ${cnLabel(k)}</span>
    <small>${ex ? esc(ex) : `No common ${CLS[W.cls]} examples — still buildable`}</small></button>`;
}
function step3() {
  const keys = Object.keys(NETS).filter(k => NETS[k].dim === W.dim);
  let list;
  if (W.browse === 'sg') {
    const groups = {};
    for (const k of keys) (groups[NETS[k].sg] ??= []).push(k);
    list = Object.entries(groups).sort().map(([sg, ks]) => `<div class="sghead">${W.dim === 2 ? 'Plane group' : 'Space group'} ${sg}</div>${ks.map(netRow).join('')}`).join('');
  } else {
    const sorted = [...keys].sort((a, b) => (!!NETS[b][W.cls]) - (!!NETS[a][W.cls]));
    list = sorted.map(netRow).join('');
  }
  return stepShell(3, 'Topology', W.net ? `${W.net.replace('_t', '')} · ${NETS[W.net].title}` : '', `
    <div class="seg" id="browse"><button data-act="browse" data-v="net" class="${W.browse === 'net' ? 'on' : ''}">By net</button><button data-act="browse" data-v="sg" class="${W.browse === 'sg' ? 'on' : ''}">By space group</button></div>
    <div class="netlist">${list}</div>`);
}

function step4() {
  if (!W.net) return stepShell(4, 'Building blocks', '', '');
  const an = net(W.net);
  const bip = buildPatch(an, 2).bipartite;
  const note = GEOM_NOTE_NETS.includes(W.net) ? `<p class="hint note">${esc(NETS[W.net].note)}</p>` : '';
  const anyChem = slotKeys().some(k => ['molecule', 'sbu'].includes(slot(k).source));
  const body = `${note}
    <div class="lbl">How do nodes connect?</div>
    <div class="seg"><button data-act="mode" data-v="linker" class="${W.mode === 'linker' ? 'on' : ''}">Through a linker</button><button data-act="mode" data-v="direct" class="${W.mode === 'direct' ? 'on' : ''}" ${bip ? '' : 'disabled title="This net is not bipartite, so direct magnet bonds would repel somewhere"'}>Directly node–node</button></div>
    <p class="hint">${W.mode === 'linker'
      ? (W.cls === 'mof' ? 'SBU ⟷ linker ⟷ SBU, e.g. MOF-5 (Zn₄O + BDC), UiO-66.' : 'Node ⟷ linker ⟷ node, e.g. [3+2] imine COFs such as TAPB–PDA.')
      : (W.cls === 'mof' ? 'Node ⟷ node, e.g. HKUST-1 (paddlewheel ⟷ BTC), PCN-600.' : 'Node ⟷ node, e.g. [3+3] COFs, COF-102 (tetrahedral ⟷ boroxine).')}</p>
    ${anyChem ? `<div class="scale-row">Model scale <input type="number" data-act="scale" value="${W.scale}" min="1" step="0.5"> mm per Å</div>` : ''}
    ${slotKeys().map(slotCard).join('')}
    <button class="primary" data-act="build">${state.built ? 'Rebuild model & parts' : 'Build model & parts'}</button>`;
  return stepShell(4, 'Building blocks', '', body);
}

function slotCard(key) {
  const s = slot(key), t = slotInfo(key), open = W.open === key;
  const r = slotSpec(key);
  const fitHtml = r.error ? `<span class="fit bad">${esc(r.error)}</span>` : `<span class="fit ${r.fit.ok ? 'ok' : 'bad'}">${r.fit.ok ? '✓' : '⚠'} ${esc(r.fit.msg)}</span>`;
  const sources = key === 'edge' ? [['generic', 'Generic'], ['molecule', 'Molecule']]
    : [['generic', 'Generic'], ['molecule', 'Molecule'], ...(W.cls === 'mof' ? [['sbu', 'Metal cluster (SBU)']] : [])];
  let summary = '';
  if (r.spec) {
    const L = r.spec.arms.map(a => a.length);
    summary = key === 'edge' ? `${esc(r.spec.name)} · ${fmt(L[0] + L[1])} mm` : `${esc(r.spec.name)} · arms ${fmt(Math.min(...L))}${Math.max(...L) - Math.min(...L) > 0.5 ? '–' + fmt(Math.max(...L)) : ''} mm`;
    if (r.spec.model) summary += ` · ${r.spec.model.atoms.length} atoms`;
  }
  const sub = key === 'edge' ? `straight, 2-c · ${t.count} per cell` : `${t.cn}-c ${t.geometry.name} · ${t.count} per cell`;
  return `<div class="slot ${open ? 'open' : ''}" data-slot="${key}">
    <div class="slot-head" data-act="openslot"><div><b>${slotTitle(key)}</b><small>${sub}</small></div><span class="caret">${open ? '▾' : '▸'}</span></div>
    <div class="seg small">${sources.map(([v, l]) => `<button data-act="src" data-v="${v}" class="${s.source === v ? 'on' : ''}">${l}</button>`).join('')}</div>
    <div class="slot-sum">${summary}</div>${fitHtml}
    ${open ? `<div class="editor">${editor(key, r)}</div>` : ''}
  </div>`;
}

/** The chosen connection points as removable chips (a sure way to undo a mis-click). */
function connChips(key, m) {
  let mm = null;
  try { mm = getMol(key); } catch { /* drawing shows the error */ }
  const conns = m.conns || [];
  if (!conns.length) return '<div class="conns"><span class="hint">No connection points yet: click atoms in the drawing.</span></div>';
  const lab = c => `${mm?.atoms[c.atom]?.el ?? ''}${c.atom + 1}`;
  return `<div class="conns"><span class="lbl">Connection points (${conns.length})</span>
    ${conns.map(c => `<span class="chip">${esc(lab(c))}<button data-act="rmconn" data-v="${c.atom}" title="Remove this connection point">×</button></span>`).join('')}
    <button class="ghost small" data-act="clearconns">Clear all</button></div>`;
}

function numField(label, path, value, attrs = '') {
  return `<label>${label}<input type="number" data-field="${path}" value="${value}" ${attrs}></label>`;
}

function editor(key, r) {
  const s = slot(key), t = slotInfo(key);
  if (s.source === 'generic') {
    return key === 'edge'
      ? `<div class="grid2">${numField('Length, magnet face to face <span class="unit">mm</span>', 'generic.length', s.generic.length, 'min="30" step="2"')}
         ${numField('Rings (bulges)', 'generic.rings', s.generic.rings, 'min="0" max="6"')}</div>
         <p class="hint">A stand-in for any straight ditopic linker: one bulge per aromatic ring.</p>`
      : `<div class="grid2">${numField('Arm length, centre → magnet <span class="unit">mm</span>', 'generic.arm', s.generic.arm ?? defaultArm(), 'min="15" step="1"')}
         ${numField('Bulges per arm', 'generic.bulges', s.generic.bulges, 'min="0" max="3"')}</div>
         <p class="hint">Arms point exactly along the ${W.net} edges, so the model closes without strain.</p>`;
  }
  if (s.source === 'sbu') {
    const opts = Object.keys(SBUS).map(k => `<option value="${esc(k)}" ${s.sbu === k ? 'selected' : ''} ${sbuCN(k) === t.cn ? '' : 'disabled'}>${esc(k)}${sbuCN(k) === t.cn ? '' : ` (${sbuCN(k)}-c)`}</option>`).join('');
    return `<label>Metal cluster<select data-field="sbu"><option value="">— choose a ${t.cn}-c cluster —</option>${opts}</select></label>
      <label>Arm directions<select data-field="sbuGeom"><option value="net" ${s.sbuGeom === 'net' ? 'selected' : ''}>Match the net (recommended)</option><option value="ideal" ${s.sbuGeom === 'ideal' ? 'selected' : ''}>Ideal cluster shape</option></select></label>
      ${s.sbu ? `<div class="lbl">Print as</div>
      <div class="seg small"><button data-act="style" data-v="block" class="${s.sbuStyle === 'model' ? '' : 'on'}">Reticular block</button><button data-act="style" data-v="model" class="${s.sbuStyle === 'model' ? 'on' : ''}">Molecular model</button></div>` : ''}
      <p class="hint">${s.sbuStyle === 'model'
        ? 'The metal core (e.g. the Cu–Cu pair of a paddlewheel, Zn₄O, Zr₆) as atoms and bonds. Each magnet fitting stands for one carboxylate and forks from the metals it bridges; the carboxylate C and O atoms come with the linker. Best at 9–12 mm per Å.'
        : 'Arm length = cluster centre → carboxylate C, × scale. Short arms are lengthened to the printable minimum.'}</p>
      ${s.sbuStyle === 'model' && r?.spec?.model ? `<div class="info"><table>
        <tr><th>Printed</th><td>${r.spec.model.atoms.length} atoms, ${r.spec.model.bonds.length} bonds</td></tr>
        <tr><th>Fittings</th><td>${r.spec.arms.length}, each from ${[...new Set(r.spec.model.fittings.map(f => f.length))].join('/')} metal atom${r.spec.model.fittings[0].length > 1 ? 's' : ''}</td></tr>
        <tr><th>Centre → magnet</th><td>${fmt(r.spec.arms[0].length, 0)} mm</td></tr></table></div>` : ''}`;
  }
  // molecule
  const m = s.mol || {};
  const suggested = PRESET_LIST.filter(p => p.cn === t.cn && p.cls === W.cls);
  const other = PRESET_LIST.filter(p => p.cn === t.cn && p.cls !== W.cls);
  const opt = p => `<option value="${esc(p.smiles)}" ${m.smiles === p.smiles ? 'selected' : ''}>${esc(p.name)}</option>`;
  const groupOpts = Object.entries(GROUPS).map(([k, g]) => `<option value="${k}" ${m.group === k ? 'selected' : ''}>${g.label}</option>`).join('');
  const geomOpts = [['net', 'Match the net (recommended)'], ['', 'From the 2D structure'], ...Object.entries(GEOM_LABELS)]
    .map(([v, l]) => `<option value="${v}" ${(m.geom ?? 'net') === v ? 'selected' : ''}>${l}</option>`).join('');
  const cb = r?.cb;
  const isModel = m.style === 'model';
  const styleSeg = `<div class="lbl">Print as</div>
    <div class="seg small"><button data-act="style" data-v="block" class="${isModel ? '' : 'on'}">Reticular block</button><button data-act="style" data-v="model" class="${isModel ? 'on' : ''}">Molecular model</button></div>
    <p class="hint">${isModel
      ? `The atoms and bonds as drawn (H left out), with a magnet fitting at each connection point standing in for the ${W.cls === 'mof' ? 'metal–carboxylate coordination bond' : 'new covalent (dynamic) bond'}. ${GROUPS[m.group]?.drop ? 'Atoms lost on bond formation are greyed out and not printed.' : ''} Best at 9–12 mm per Å.`
      : 'Hub, arms and ring bulges: the shape of the molecule reduced to its topology.'}</p>`;
  const info = cb && isModel && r.spec?.model ? `<table>
      <tr><th>Connections</th><td>${cb.n} magnet fittings${cb.autoGeometry ? ` · sp³ centre → ${cb.autoGeometry}` : ''}</td></tr>
      <tr><th>Printed</th><td>${r.spec.model.atoms.length} atoms, ${r.spec.model.bonds.length} bonds</td></tr>
      <tr><th>Fitting bend</th><td>${r.spec.model.bend > 1 ? `${r.spec.model.bend}° (the bond, not the molecule, bends to meet the net)` : 'none (straight)'}</td></tr>
      <tr><th>Centre → magnet</th><td>${[...new Set(r.spec.arms.map(a => fmt(a.length, 0)))].join(', ')} mm</td></tr></table>`
    : cb ? `<table>
      <tr><th>Connections</th><td>${cb.n}${cb.autoGeometry ? ` · sp³ centre → ${cb.autoGeometry}` : ''}</td></tr>
      <tr><th>Arms</th><td>${cb.lengths.map(l => fmt(l, 1)).join(', ')} Å → ${r.spec.arms.map(a => fmt(a.length, 0)).join(', ')} mm</td></tr>
      <tr><th>Ring bulges</th><td>${cb.bulges.reduce((a, b) => a + b.length, 0) + (cb.hubRing ? 1 : 0)}${cb.hubRing ? ' (one at the centre)' : ''}${cb.unusedRings ? ` · ${cb.unusedRings} off-axis ring(s) not shown` : ''}</td></tr></table>` : '';
  return `<label>Example ${t.cn}-connected molecules<select data-act="preset"><option value="">— choose —</option>
      ${suggested.length ? `<optgroup label="Typical for ${CLS[W.cls]}s">${suggested.map(opt).join('')}</optgroup>` : ''}
      ${other.length ? `<optgroup label="Other ${t.cn}-c molecules">${other.map(opt).join('')}</optgroup>` : ''}</select></label>
    <label>SMILES<div class="row"><input type="text" data-act="smiles" value="${esc(m.smiles || '')}" placeholder="e.g. O=Cc1ccc(C=O)cc1" spellcheck="false"><button data-act="loadsmiles">Load</button></div></label>
    <div class="row"><button class="ghost" data-act="draw">✎ Draw a structure…</button><button class="ghost" data-act="import" title="${IMPORT_ACCEPT}">⤓ Import 3D structure…</button></div>
    ${m.molblock ? `<p class="hint">3D structure imported from a file: its real shape (twisted rings, bond angles) is kept. The drawing below is a projection; click atoms to change connection points. Typing a SMILES or drawing replaces it with a flat 2D structure.</p>` : ''}
    <div class="molview" id="molView">${m.smiles ? '' : '<p class="hint">Pick, type or draw a molecule. Click atoms to toggle connection points.</p>'}</div>
    <div class="legend">${isModel ? '<span class="dot conn"></span>connection point → magnet fitting' : '<span class="dot conn"></span>connection point <span class="dot ring"></span>ring → bulge <span class="dot hub"></span>ring → hub'} · click an atom to add, click again to remove</div>
    ${m.smiles ? connChips(key, m) : ''}
    ${m.smiles ? `${styleSeg}<div class="grid2">
      <label>Linking group<select data-field="mol.group" data-act="group">${groupOpts}</select></label>
      ${numField('Extension past atom <span class="unit">Å</span>', 'mol.ext', m.ext ?? 0, 'step="0.05" min="0"')}</div>
      <label>Arm directions<select data-field="mol.geom">${geomOpts}</select></label>
      <div class="info">${info}</div>` : ''}`;
}

function renderWizard() {
  if (W.step > 1 && !reachable(W.step)) W.step = [4, 3, 2, 1].find(reachable);
  if (state.builder === 'cif') { $('#wizard').innerHTML = builderSwitch() + renderCif(); return; }
  $('#wizard').innerHTML = builderSwitch() + `<div class="wtitle">Build a framework model</div>` + step1() + step2() + step3() + step4();
  if (W.open && slot(W.open).source === 'molecule') drawMolecule(W.open);
}

function previewSlot(key) {
  const r = slotSpec(key);
  state.preview = r.spec || null;
  state.view = 'block';
  refreshView();
}

// ---------------- wizard events ----------------
function selectNet(k) {
  if (W.net !== k) {
    W.net = k; W.slots = {}; W.open = null; molCache.clear();
    const an = net(k);
    W.mode = Object.keys(an.types).length > 1 && buildPatch(an, 2).bipartite ? 'direct' : 'linker';
    state.built = false;
  }
  W.step = 4;
}

$('#wizard').addEventListener('click', e => {
  const el = e.target.closest('[data-act]');
  if (!el || el.disabled) return;
  const act = el.dataset.act, v = el.dataset.v;
  const key = el.closest('[data-slot]')?.dataset.slot;
  switch (act) {
    case 'goto': { const n = +el.dataset.step; if (reachable(n)) W.step = n; break; }
    case 'cls': W.cls = v; W.step = W.dim ? (W.net ? 4 : 3) : 2; break;
    case 'dim': if (W.dim !== +v) { W.dim = +v; W.net = null; W.slots = {}; } W.step = 3; break;
    case 'browse': W.browse = v; break;
    case 'net': selectNet(v); break;
    case 'mode': W.mode = v; W.open = null; state.built = false; break;
    case 'openslot': W.open = W.open === key ? null : key; if (W.open) previewSlot(key); break;
    case 'src': slot(key).source = v; W.open = key; previewSlot(key); break;
    case 'loadsmiles': loadSlotMolecule(key, el.closest('.row').querySelector('input').value, 'Custom molecule'); break;
    case 'draw': openDraw(key); return;
    case 'rmconn': toggleConn(key, +v); return;
    case 'clearconns': { const s = slot(key); if (s.mol) { s.mol.conns = []; s.mol.group = 'manual'; } previewSlot(key); break; }
    case 'import': importFor = key; fileIn.value = ''; fileIn.click(); return;
    case 'builder': state.builder = v; state.preview = null; if (v === 'cif' && state.cif?.sel != null) previewCif(state.cif.sel); else refreshView(); break;
    case 'cifopen': fwIn.value = ''; fwIn.click(); return;
    case 'cifblock': state.cif.sel = +v; previewCif(+v); break;
    case 'cifbuild': buildCif(); break;
    case 'style': {
      const s = slot(key);
      if (s.source === 'sbu') s.sbuStyle = v;
      else if (!s.mol) return;
      else s.mol.style = v;
      if (v === 'model' && W.scale < 8) W.scale = 10; // rings stay open; the scale row shows the change
      previewSlot(key);
      break;
    }
    case 'build': buildModel(); break;
    default: return;
  }
  save(); renderWizard();
});

$('#wizard').addEventListener('change', e => {
  const el = e.target;
  const key = el.closest('[data-slot]')?.dataset.slot;
  if (el.dataset.act === 'cifscale') { state.cif.scale = +el.value || 10; if (state.cif.sel != null) previewCif(state.cif.sel); }
  else if (el.dataset.act === 'scale') { W.scale = +el.value || 6; }
  else if (el.dataset.act === 'preset') { if (el.value) loadSlotMolecule(key, el.value, el.selectedOptions[0].textContent); }
  else if (el.dataset.act === 'smiles') { loadSlotMolecule(key, el.value, 'Custom molecule'); }
  else if (el.dataset.field && key) {
    const s = slot(key);
    const [a, b] = el.dataset.field.split('.');
    const val = el.type === 'number' ? +el.value : el.value;
    if (b) s[a][b] = val; else s[a] = val;
    if (el.dataset.act === 'group' && s.mol && getMol(key)) {
      s.mol.conns = s.mol.group === 'manual' ? s.mol.conns : findConnections(getMol(key), s.mol.group);
      s.mol.ext = GROUPS[s.mol.group].ext;
    }
  } else return;
  if (key) previewSlot(key);
  save(); renderWizard();
});

function loadSlotMolecule(key, smiles, name, molblock = null, notes = []) {
  smiles = (smiles || '').trim();
  if (!smiles || !key) return;
  if (!rdkitReady) { setWarnings(['The chemistry toolkit is still loading — try again in a moment.']); return; }
  const s = slot(key);
  try {
    molCache.delete(key);
    s.mol = { smiles, name, geom: s.mol?.geom ?? (key === 'edge' ? '' : 'net'), style: s.mol?.style ?? 'block' };
    if (molblock) s.mol.molblock = molblock; // imported 3D structure: real coordinates kept
    const mm = getMol(key);
    s.mol.group = autoDetectGroup(mm);
    s.mol.conns = findConnections(mm, s.mol.group);
    s.mol.ext = GROUPS[s.mol.group].ext;
    s.source = 'molecule';
    W.open = key;
    previewSlot(key);
    if (notes.length) setWarnings(notes);
  } catch (err) {
    s.mol = null;
    setWarnings([err.message]);
  }
  save(); renderWizard();
}

// ---------------- molecule drawing (SVG, clickable atoms) ----------------
function drawMolecule(key) {
  const host = $('#molView');
  const mm = host && getMol(key);
  if (!mm) return;
  const s = slot(key);
  const r = slotSpec(key);
  const cb = r.cb;
  const W2 = 320, pad = 22;
  // 3D structures are drawn projected onto their best-fit plane; 2D ones as laid out
  const pr = projector(mm);
  const A2 = mm.atoms.map(a => pr([a.x, a.y, a.z]));
  const xs = A2.map(q => q[0]), ys = A2.map(q => q[1]);
  const [x0, x1, y0, y1] = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
  const k = Math.min((W2 - 2 * pad) / Math.max(x1 - x0, 1), 26);
  const H = Math.max(120, (y1 - y0) * k + 2 * pad);
  const ox = (W2 - (x1 - x0) * k) / 2;
  const X = x => ox + (x - x0) * k, Y = y => H - pad - (y - y0) * k;
  const conn = new Set(s.mol.conns.map(c => c.atom));
  const isModel = !!r.spec?.model;
  const kept = isModel ? new Set(r.spec.model.atoms.map(a => a.src)) : null;
  const faded = i => kept && !kept.has(i);
  let svg = `<svg viewBox="0 0 ${W2} ${H}" xmlns="http://www.w3.org/2000/svg" font-family="system-ui" font-size="11">`;
  const C2 = cb && pr(cb.centre);
  if (cb) for (const p3 of cb.pts) { const p = pr(p3); svg += `<line x1="${X(C2[0])}" y1="${Y(C2[1])}" x2="${X(p[0])}" y2="${Y(p[1])}" style="stroke:var(--conn)" stroke-width="6" stroke-opacity=".13" stroke-linecap="round"/>`; }
  for (const ring of isModel ? [] : mm.rings) {
    const st = !cb ? 'off' : v3.len(v3.sub(ring.centre, cb.centre)) < 1.0 ? 'hub' : isBulge(ring, cb) ? 'ring' : 'off';
    const col = st === 'hub' ? 'var(--hub)' : st === 'ring' ? 'var(--ring)' : 'var(--muted)';
    const rc = pr(ring.centre);
    svg += `<circle cx="${X(rc[0])}" cy="${Y(rc[1])}" r="${k * 0.62}" fill="none" style="stroke:${col}" stroke-width="${st === 'off' ? 1 : 2}" ${st === 'off' ? 'stroke-dasharray="3 3" opacity=".5"' : 'opacity=".8"'}/>`;
  }
  for (const b of mm.bonds) {
    const a = A2[b.a], c = A2[b.b];
    const [xa, ya, xc, yc] = [X(a[0]), Y(a[1]), X(c[0]), Y(c[1])];
    const L = Math.hypot(xc - xa, yc - ya) || 1, nx = -(yc - ya) / L * 3, ny = (xc - xa) / L * 3;
    const op = faded(b.a) || faded(b.b) ? ' opacity=".3"' : '';
    const line = o => `<line x1="${xa + nx * o}" y1="${ya + ny * o}" x2="${xc + nx * o}" y2="${yc + ny * o}" style="stroke:var(--ink)" stroke-width="1.4"${op}/>`;
    svg += b.order === 2 ? line(-0.5) + line(0.5) : b.order === 3 ? line(-1) + line(0) + line(1) : line(0);
  }
  // connection points drawn last so they sit on top of crowded neighbours
  for (const a of [...mm.atoms].sort((p, q) => conn.has(p.i) - conn.has(q.i))) {
    const [ax, ay] = A2[a.i];
    const lab = a.el === 'C' ? '' : (['R#', 'R', 'A', 'Q', '*'].includes(a.el) ? '*' : a.el);
    svg += `<g class="atom${conn.has(a.i) ? ' on' : ''}" data-i="${a.i}" data-x="${X(ax).toFixed(1)}" data-y="${Y(ay).toFixed(1)}"${faded(a.i) ? ' opacity=".35"' : ''}><title>${conn.has(a.i) ? 'Connection point: click to remove' : 'Click to make a connection point'}</title><circle class="hit" cx="${X(ax)}" cy="${Y(ay)}" r="${Math.min(9, Math.max(5, k * 0.45)).toFixed(1)}" fill="transparent"/>`;
    if (lab) svg += `<circle cx="${X(ax)}" cy="${Y(ay)}" r="7" style="fill:var(--bg)"/><text x="${X(ax)}" y="${Y(ay) + 4}" text-anchor="middle" style="fill:var(--ink)">${lab}</text>`;
    if (conn.has(a.i)) svg += `<circle cx="${X(ax)}" cy="${Y(ay)}" r="5.5" style="fill:var(--conn)" opacity=".85"/>`;
    svg += '</g>';
  }
  if (cb) svg += `<circle cx="${X(C2[0])}" cy="${Y(C2[1])}" r="3" style="fill:var(--hub)"/>`;
  host.innerHTML = svg + '</svg>';
}
/** 2D view of a molecule: identity for 2D drawings; best-fit plane (PCA) for 3D structures. */
function projector(mm) {
  if (!mm.is3D) return p => [p[0], p[1]];
  const P3 = mm.atoms.map(a => [a.x, a.y, a.z]);
  const c = v3.mul(P3.reduce((s, p) => v3.add(s, p), [0, 0, 0]), 1 / P3.length);
  const C = [0, 1, 2].map(i => [0, 1, 2].map(j => P3.reduce((s, p) => s + (p[i] - c[i]) * (p[j] - c[j]), 0)));
  const mul = v => [0, 1, 2].map(i => C[i][0] * v[0] + C[i][1] * v[1] + C[i][2] * v[2]);
  const top = (start, not) => {
    let v = start;
    for (let k = 0; k < 60; k++) {
      v = mul(v);
      if (not) v = v3.sub(v, v3.mul(not, v3.dot(v, not)));
      v = v3.len(v) > 1e-12 ? v3.norm(v) : start;
    }
    return v;
  };
  const e1 = top([1, 0.3, 0.2]);
  const e2 = top(v3.norm(v3.cross(e1, Math.abs(e1[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0])), e1);
  return p => { const d = v3.sub(p, c); return [v3.dot(d, e1), v3.dot(d, e2)]; };
}
function isBulge(r, cb) {
  return cb.pts.some(p => {
    const axis = v3.sub(p, cb.centre), L = v3.len(axis), d = v3.norm(axis);
    const rc = v3.sub(r.centre, cb.centre), t = v3.dot(rc, d);
    return t > 0.5 && t < L && v3.len(v3.sub(rc, v3.mul(d, t))) < 1.0;
  });
}
/** Toggle atom i as a connection point: first click adds it, second click removes it. */
function toggleConn(key, i) {
  const s = slot(key);
  if (!s.mol) return;
  s.mol.group = 'manual';
  const at = s.mol.conns.findIndex(c => c.atom === i);
  if (at >= 0) s.mol.conns.splice(at, 1); else s.mol.conns.push({ atom: i });
  previewSlot(key);
  save(); renderWizard();
}
$('#wizard').addEventListener('click', e => {
  const svg = e.target.closest('#molView svg');
  if (!svg || !W.open) return;
  // pick the atom nearest the click, not whichever hit-circle happens to be on top; a marked
  // connection point wins a near-tie so clicking its red dot always removes it
  const pt = new DOMPoint(e.clientX, e.clientY).matrixTransform(svg.getScreenCTM().inverse());
  let best = null, bestD = 14;
  for (const g of svg.querySelectorAll('.atom')) {
    const d = Math.hypot(+g.dataset.x - pt.x, +g.dataset.y - pt.y) - (g.classList.contains('on') ? 3 : 0);
    if (d < bestD) { bestD = d; best = g; }
  }
  if (best) toggleConn(W.open, +best.dataset.i);
});

// ---------------- JSME drawing dialog ----------------
let jsme = null, drawFor = null;
function loadJSME() {
  return new Promise((resolve, reject) => {
    if (jsme) return resolve(jsme);
    window.jsmeOnLoad = () => { jsme = new window.JSApplet.JSME('jsmeHost', '100%', '420px', { options: 'oldlook,star,nocanonize,noquery' }); resolve(jsme); };
    const sc = document.createElement('script');
    sc.src = 'vendor/jsme/jsme.nocache.js';
    sc.onerror = () => reject(new Error('Could not load the JSME editor.'));
    document.head.appendChild(sc);
  });
}
// ---------------- 3D structure import (.mol/.sdf/.mol2/.xyz/.pdb/.cif) ----------------
let importFor = null;
const fileIn = Object.assign(document.createElement('input'), { type: 'file', accept: IMPORT_ACCEPT, hidden: true });
document.body.append(fileIn);
fileIn.addEventListener('change', async () => {
  const f = fileIn.files[0];
  if (!f || !importFor) return;
  if (!rdkitReady) { setWarnings(['The chemistry toolkit is still loading — try again in a moment.']); return; }
  try {
    const imp = importStructure(await f.text(), f.name);
    const { molblock, smiles, flat } = prepare3D(imp.molblock);
    const name = imp.name && imp.name !== 'Imported structure' ? imp.name : f.name.replace(/\.[^.]+$/, '');
    if (flat) imp.notes.push('The file is flat (all z = 0), so it was read as a 2D structure.');
    loadSlotMolecule(importFor, smiles, flat ? name : `${name} · 3D`, molblock, imp.notes);
  } catch (e) { setWarnings([`Could not import ${f.name}: ${e.message}`]); }
});

// ---------------- whole-framework import (COF CIF → net + node and linker blocks) ----------------
const fwIn = Object.assign(document.createElement('input'), { type: 'file', accept: '.cif', hidden: true });
document.body.append(fwIn);
fwIn.addEventListener('change', async () => {
  const f = fwIn.files[0];
  if (!f) return;
  status(`Reading ${f.name}…`);
  try {
    const fw = importFramework(await f.text(), f.name);
    state.cif = { name: fw.name, fw, scale: state.cif?.scale ?? 10, sel: 0 };
    state.cifError = '';
    cifCache.clear(); cifK = { key: '', k: 0 };
    registerCifNet();
    previewCif(0);
  } catch (e) { state.cifError = `Could not read ${f.name}: ${e.message}`; }
  status(rdkitReady ? 'Ready' : 'Loading chemistry toolkit…');
  setWarnings([]);
  save(); renderWizard();
});

function builderSwitch() {
  const b = (v, l) => `<button data-act="builder" data-v="${v}" class="${state.builder === v ? 'on' : ''}">${l}</button>`;
  return `<div class="seg builder">${b('scratch', 'Build from scratch')}${b('cif', 'Start from a CIF')}</div>`;
}

// Each block of an imported framework: its molecule exactly as computed, connection points at the
// cut linkage bonds (the * atoms), fittings along the crystal's own bond directions.
const cifCache = new Map();
function cifBlock(i, k = null) {
  const C = state.cif, t = C.fw.types[i];
  if (!t.graph) return { error: 'Saved from an older version: choose the CIF again.' };
  try {
    let m = cifCache.get(i);
    if (!m) {
      const mm = moleculeFromGraph(t.graph);
      m = { mm, conns: findConnections(mm, 'dummy') };
      cifCache.set(i, m);
    }
    // SMILES for organic pieces (RDKit); metal clusters are shown by formula
    if (m.smiles === undefined && rdkitReady) { try { m.smiles = t.molblock ? prepare3D(t.molblock).smiles : ''; } catch { m.smiles = ''; } }
    const kind = t.role === 'linker' ? 'linker' : 'node';
    const multi = C.fw.types.filter(x => x.role === 'node').length > 1;
    const name = `${C.name} ${t.metal ? 'SBU' : kind}${multi && kind === 'node' ? ' ' + t.label : ''}`;
    if (t.fixed) {
      // a net not in the library: the block sits exactly as in the crystal, its arms reaching the magnet
      // faces on the lines to its neighbours, all scaled by one factor (k mm per Å) so the model closes
      const kk = k ?? cifArmScale();
      const cb = moleculeToBlock(m.mm, m.conns, { ext: 0, geometryOverride: '' });
      const spec = moleculeModelSpec(m.mm, m.conns, cb, { group: 'dummy', ext: 0, scale: C.scale, P: state.P, kind, maxBend: 35, name,
        fixed: { centre: t.fixed.centre, dirs: t.fixed.dirs, lengths: t.fixed.t.map(x => x * kk) } });
      spec.kind = kind;
      const b = spec.model.bend;
      return { spec, smiles: m.smiles, fit: { ok: true, msg: `placed as in the crystal${b > 3 ? `; fittings bend up to ${b}°` : ''}` } };
    }
    // the molecule stays exactly as in the crystal (rigid); its magnet fittings point along the edges of
    // the ideal net with equal arms, so every copy fits every site and the model closes
    const site = kind === 'linker' ? { cn: 2, dirs: [[1, 0, 0], [-1, 0, 0]] } : net(C.fw.net).types[t.label];
    const cb = moleculeToBlock(m.mm, m.conns, { ext: 0, geometryOverride: site.dirs });
    const spec = moleculeModelSpec(m.mm, m.conns, cb, { group: 'dummy', ext: 0, scale: C.scale, P: state.P, targets: cb.dirs, kind, maxBend: 20, name });
    spec.kind = kind;
    const b = spec.model.bend;
    return { spec, smiles: m.smiles, fit: { ok: true, msg: b > 3 ? `molecule as in the crystal; fittings bend up to ${b}° to meet the net` : 'molecule as in the crystal; fittings straight' } };
  } catch (e) { return { error: e.message }; }
}
/**
 * Crystal nets: the one scale (mm per Å) for every arm, at least the model scale and large enough that
 * every magnet fitting clears its atoms (a short metal–O link may need it raised).
 */
let cifK = { key: '', k: 0 };
function cifArmScale() {
  const C = state.cif, key = JSON.stringify([C.name, C.scale, state.P]);
  if (cifK.key === key) return cifK.k;
  let k = C.scale;
  C.fw.types.forEach((t, i) => {
    if (!t.fixed) return;
    const r = cifBlock(i, C.scale);
    if (r.spec) t.fixed.t.forEach((x, a) => { k = Math.max(k, r.spec.model.minLen[a] / x); });
  });
  cifK = { key, k: Math.ceil(k * 10) / 10 };
  return cifK.k;
}
/** Make the crystal's own net (when the library has none that fits) available to the assembly. */
function registerCifNet() {
  if (state.cif?.fw?.crystalNet) netCache.set('cif', state.cif.fw.crystalNet);
}
const netTitle = k => (k === 'cif' ? 'net taken from the crystal' : `${k.replace('_t', '')} · ${NETS[k].title}`);
function previewCif(i) {
  const r = cifBlock(i);
  state.preview = r.spec || null;
  state.view = 'block';
  refreshView();
}
const cifKey = t => (t.role === 'linker' ? 'edge' : 'v:' + t.label);

function renderCif() {
  const C = state.cif;
  const section = (n, title, body) => `<section class="step open"><header><span class="num">${n}</span><span class="stitle">${title}</span></header><div class="sbody">${body}</div></section>`;
  const pick = `<div class="row"><button class="${C ? 'ghost' : 'primary'}" data-act="cifopen">⤓ ${C ? 'Choose another CIF…' : 'Choose a CIF…'}</button>${C ? `<span class="cifname">${esc(C.name)}.cif</span>` : ''}</div>`;
  let html = `<div class="wtitle">Model a framework from its crystal structure</div>` + section(1, 'Crystal structure', `
    <p class="hint">A COF or MOF unit cell (.cif) from Materials Studio, Mercury, VESTA, CrystalMaker or similar. The framework is cut into its building blocks at the bonds that link them, and the magnet connections go exactly where those bonds were:</p>
    <ul class="hint"><li><b>COFs:</b> imine, hydrazone, azine and β-ketoenamine bonds.</li><li><b>MOFs:</b> each carboxylate's C–C bond (the carboxylate stays with the metal cluster, the SBU) and metal–N bonds of pyridyl and azolate linkers.</li></ul>
    <p class="hint">The net is worked out from how the blocks connect. Every block prints as a molecular model in the exact conformation from the file.</p>
    ${pick}${state.cifError ? `<p class="fit bad">${esc(state.cifError)}</p>` : ''}`);
  if (!C) return html;
  const fw = C.fw;
  const blocks = fw.types.map((t, i) => [t, cifBlock(i), i]);
  const card = ([t, r, i]) => {
    const L = r.spec?.arms.map(reach) || [];
    const size = !r.spec ? '' : t.role === 'linker' ? `${fmt(L[0] + L[1], 0)} mm magnet to magnet` : `arms ${fmt(Math.min(...L), 0)}${Math.max(...L) - Math.min(...L) > 0.5 ? '–' + fmt(Math.max(...L), 0) : ''} mm`;
    return `<div class="slot cifblock ${C.sel === i ? 'open' : ''}" data-act="cifblock" data-v="${i}">
      <div class="slot-head"><div><b>${t.role === 'linker' ? 'Linker' : (t.metal ? 'Metal cluster (SBU)' : 'Node') + (fw.types.filter(x => x.role === 'node').length > 1 ? ` ${t.label}` : '')}</b><small>${t.formula} · ${t.cn}-c · ${t.perCell} per cell</small></div></div>
      <div class="slot-sum">${r.spec ? `${r.spec.model.atoms.length} atoms printed · ${size}` : ''}</div>
      ${r.error ? `<span class="fit bad">${esc(r.error)}</span>` : `<span class="fit ${r.fit.ok ? 'ok' : 'bad'}">${r.fit.ok ? '✓' : '⚠'} ${esc(r.fit.msg)}</span>`}
      ${r.smiles ? `<div class="smiles" title="SMILES (* = the partner across each cut linkage)">${esc(r.smiles)}</div>` : ''}
      ${t.spread > 0.05 ? `<p class="hint">The ${t.perCell} copies in the cell differ by up to ${fmt(t.spread, 2)} Å; the most typical one is used.</p>` : ''}
    </div>`;
  };
  html += section(2, 'What was found', `
    <div class="info"><table>
      <tr><th>Net</th><td>${esc(netTitle(fw.net))} (${fw.dim}D)</td></tr>
      <tr><th>Type</th><td>${fw.cls === 'mof' ? 'MOF' : 'COF'} · ${fw.mode === 'linker' ? 'node ⟷ linker ⟷ node' : 'node ⟷ node, directly'}</td></tr></table></div>
    ${blocks.map(card).join('')}
    ${fw.notes.length ? `<div class="hint">${fw.notes.map(esc).join('<br>')}</div>` : ''}
    <p class="hint">Click a block to see it. Magnet fittings start at the atoms of the cut bonds (forking when one link has several) and run to ${fw.crystalNet ? 'the line between neighbouring blocks, as they sit in the crystal' : 'the edges of the ideal net'}.</p>`);
  // size: node centre → node centre at this scale
  const node = blocks.find(([t]) => t.role === 'node')?.[1].spec, link = blocks.find(([t]) => t.role === 'linker')?.[1].spec;
  const E = node ? 2 * meanArm(node) + (link ? reach(link.arms[0]) + reach(link.arms[1]) : 0) : 0;
  html += section(3, 'Size and build', `
    <div class="scale-row">Model scale <input type="number" data-act="cifscale" value="${C.scale}" min="1" step="0.5"> mm per Å</div>
    ${fw.crystalNet ? `<p class="hint">Molecules at ${C.scale} mm per Å; the framework itself (block positions and arms) is the crystal scaled by ${fmt(cifArmScale(), 1)} mm per Å${cifArmScale() > C.scale ? ', a little more so every magnet fitting clears its atoms' : ''}.</p>`
      : E ? `<p class="hint">Node centre to node centre: ${fmt(E, 0)} mm. Molecular models print best at 9–12 mm per Å; lower it for smaller pieces.</p>` : ''}
    <button class="primary" data-act="cifbuild">Build model & parts</button>`);
  return html;
}
function buildCif() {
  const C = state.cif;
  if (!C) return;
  const results = C.fw.types.map((t, i) => [cifKey(t), cifBlock(i)]);
  const errors = results.filter(([, r]) => r.error).map(([k, r]) => `${k === 'edge' ? 'Linker' : 'Node'}: ${r.error}`);
  if (errors.length) { setWarnings(errors); return; }
  assemble(C.fw.net, C.fw.mode, results.map(([k, r]) => [k, r.spec]), C.fw.crystalNet ? cifArmScale() : null);
}

async function openDraw(key) {
  drawFor = key;
  $('#drawDlg').showModal();
  try {
    const j = await loadJSME();
    const cur = slot(key).mol?.smiles;
    if (cur) j.readGenericMolecularInput(cur); else j.reset();
  } catch (e) { $('#jsmeHost').textContent = e.message; }
}
$('#drawCancel').addEventListener('click', () => $('#drawDlg').close());
$('#drawUse').addEventListener('click', () => {
  const smi = jsme?.smiles();
  $('#drawDlg').close();
  if (smi && drawFor) loadSlotMolecule(drawFor, smi, 'Drawn molecule');
});

// =====================================================================
//  Build the model: turn every slot into a part, set up the assembly
// =====================================================================
function buildModel() {
  const keys = slotKeys();
  const results = keys.map(k => [k, slotSpec(k)]);
  const errors = results.filter(([, r]) => r.error).map(([k, r]) => `${slotTitle(k)}: ${r.error}`);
  if (errors.length) { setWarnings(errors); return; }
  assemble(W.net, W.mode, results.map(([k, r]) => [k, r.spec]));
  state.built = true;
  save(); renderWizard();
}
/** Turn block specs (keyed 'v:<label>' / 'edge') into parts and set up the assembled model. */
function assemble(netKey, mode, specs, crystalScale = null) {
  state.parts = state.parts.filter(p => p.source !== 'wizard');
  const an = net(netKey);
  const patch = buildPatch(an, 2);
  const nodes = {};
  let linker = null;
  for (const [k, s] of specs) {
    const spec = structuredClone(s);
    if (k === 'edge') { spec.polarity = 'S'; linker = addPart(spec).id; continue; }
    const label = k.slice(2);
    if (mode === 'linker') spec.polarity = 'N';
    else {
      const cols = new Set(patch.verts.filter(v => v.type === label).map(v => v.colour));
      spec.polarity = cols.size > 1 ? 'N+S' : cols.has(0) ? 'N' : 'S';
    }
    nodes[label] = addPart(spec).id;
  }
  state.asm = { net: netKey, nodes, linker, reps: an.dim === 2 ? 3 : (an.vertsPerCell > 20 ? 1 : 2), ...(crystalScale ? { k: crystalScale } : {}) };
  $('#asmReps').value = state.asm.reps;
  state.preview = null;
  state.view = 'assembly';
  save(); renderParts(); refreshView();
}
function addPart(spec) {
  const part = { id: nextId++, spec, source: 'wizard' };
  state.parts.push(part);
  return part;
}

// =====================================================================
//  Viewer, parts list, assembly, export
// =====================================================================
function setWarnings(list) { $('#warnings').innerHTML = list.map(w => `<div>${esc(w)}</div>`).join(''); }

function showBlock(spec, key) {
  if (!spec) { viewer.clear(); $('#viewInfo').textContent = ''; setWarnings([]); return; }
  try {
    const m = meshFor(spec, key);
    viewer.showMesh(m.positions, spec.kind === 'linker' ? LINKER_COLOUR : NODE_COLOURS[0]);
    const b = m.bbox; const sz = [0, 1, 2].map(i => b.max[i] - b.min[i]);
    $('#viewInfo').textContent = `${spec.name} — ${fmt(sz[0], 0)} × ${fmt(sz[1], 0)} × ${fmt(sz[2], 0)} mm · ${spec.arms.length} magnets · ≈${fmt(m.volume / 1000 * 1.24, 0)} g PLA`;
    setWarnings([...(spec.notes || []), ...checkSpec(spec, state.P)]);
  } catch (e) { console.error(e); setWarnings([`Could not build this block: ${e.message}`]); }
}

function refreshView() {
  $$('#viewMode button').forEach(b => b.classList.toggle('on', b.dataset.v === state.view));
  $('#repsWrap').style.display = state.view === 'assembly' ? '' : 'none';
  if (state.view === 'assembly') return showAssembly();
  if (state.preview) return showBlock(state.preview, null);
  const p = state.parts.find(p => p.id === state.selected);
  showBlock(p?.spec, p ? 'part' + p.id : null);
}

function describeArms(spec) {
  const L = spec.arms.map(reach);
  if (spec.kind === 'linker') return `${fmt(L[0] + L[1])} mm long`;
  return Math.max(...L) - Math.min(...L) < 0.5 ? `arms ${fmt(L[0])} mm` : `arms ${fmt(Math.min(...L))}–${fmt(Math.max(...L))} mm`;
}
const polTag = pol => pol === 'N+S' ? '<span class="tag NS" title="Print some with N-out and some with S-out magnets — see the bill of materials">N &amp; S out</span>' : `<span class="tag ${pol}" title="Magnet pole facing outwards">${pol} out</span>`;

function renderParts() {
  const ul = $('#partsList');
  if (!state.parts.length) { ul.innerHTML = '<p class="hint">Parts appear here after step 4 → “Build model & parts”.</p>'; return; }
  ul.innerHTML = state.parts.map(p => {
    const s = p.spec;
    const g = (s.kind === 'linker' ? 'linker' : `${s.arms.length}-c ${describeGeometry(s.arms.map(a => a.dir)).name}`) + (s.model ? ' · molecular model' : '');
    return `<li data-id="${p.id}" class="${p.id === state.selected && state.view === 'block' && !state.preview ? 'sel' : ''}">
      <span class="pname" title="${esc(s.name)}">${esc(s.name)}</span>
      <span class="pmeta">${g} · ${describeArms(s)} · ${polTag(s.polarity)}</span>
      <span class="pact"><button data-act="stl" title="Download STL">STL</button><button data-act="del" title="Remove">✕</button></span>
    </li>`;
  }).join('');
}
$('#partsList').addEventListener('click', e => {
  const li = e.target.closest('li'); if (!li) return;
  const id = +li.dataset.id, part = state.parts.find(p => p.id === id), act = e.target.dataset.act;
  if (act === 'del') {
    state.parts = state.parts.filter(p => p.id !== id);
    if (state.selected === id) state.selected = null;
    save(); renderParts(); refreshView(); return;
  }
  if (act === 'stl') return downloadSTL(part);
  state.selected = id; state.preview = null; state.view = 'block';
  renderParts(); refreshView();
});

const safeName = s => s.replace(/[^\w.-]+/g, '_').replace(/_+/g, '_').slice(0, 60);
const fileName = p => `${safeName(p.spec.name)}_${p.spec.polarity.replace('+', '')}out.stl`;
const stlBuffer = part => meshToSTL(meshFor(part.spec, 'part' + part.id).positions, part.spec.name);
function downloadBlob(blob, name) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}
const downloadSTL = part => downloadBlob(new Blob([stlBuffer(part)], { type: 'model/stl' }), fileName(part));
$('#dlAll').addEventListener('click', async () => {
  if (!state.parts.length) return;
  status('Building STL files…');
  const zip = new JSZip();
  for (const p of state.parts) zip.file(fileName(p), stlBuffer(p));
  zip.file('ASSEMBLY_GUIDE.txt', assemblyGuide());
  downloadBlob(await zip.generateAsync({ type: 'blob' }), 'reticublox_parts.zip');
  status('Ready');
});
$('#clearParts').addEventListener('click', () => {
  if (!state.parts.length || !confirm('Remove all parts?')) return;
  state.parts = []; state.selected = null; state.asm = null; state.built = false;
  $('#bom').innerHTML = '';
  save(); renderParts(); renderWizard(); refreshView();
});

/** What the assembly preview needs to draw a molecular model: atoms, bonds, and each arm's atom + tip. */
const modelLook = spec => spec.model ? {
  atoms: spec.model.atoms.map(a => ({ p: a.p, r: atomRadius(a.el, state.P) })),
  bonds: spec.model.bonds,
  arms: spec.arms.map((a, k) => ({ dir: a.dir, from: [].concat(spec.model.fittings[k]).map(i => spec.model.atoms[i].p), tip: armTip(a) })),
} : null;
/** Centre → magnet face distance of an arm (arms of crystal pieces start at their atom, not the centre). */
const reach = a => v3.len(armTip(a));
const meanArm = spec => spec.arms.reduce((s, a) => s + reach(a), 0) / spec.arms.length;
function showAssembly() {
  const a = state.asm;
  if (!a) { viewer.clear(); setWarnings(['Finish step 4 and press “Build model & parts” to see the assembled framework.']); $('#viewInfo').textContent = ''; return; }
  const an = net(a.net);
  const part = id => state.parts.find(p => p.id === id);
  const nodeLook = {};
  let i = 0;
  for (const t of Object.values(an.types)) {
    const p = part(a.nodes?.[t.label]);
    if (!p) continue;
    nodeLook[t.label] = { arm: meanArm(p.spec), hub: p.spec.hubRadius, bulges: p.spec.arms[0].bulges.map(b => b.t), colour: NODE_COLOURS[i % 4], colourB: NODE_COLOURS_B[i % 4], partId: p.id, model: modelLook(p.spec) };
    i++;
  }
  const lp = part(a.linker);
  const patch = buildPatch(an, a.reps || 2);
  if (!Object.keys(nodeLook).length || !patch.edges.length) { viewer.clear(); setWarnings(['Nothing to show — rebuild the model.']); return; }
  const [, , ia, ib] = patch.edges[0];
  const linkLen = lp ? reach(lp.spec.arms[0]) + reach(lp.spec.arms[1]) : 0;
  // a crystal net is in Å, scaled as a whole; a library net has unit edges, scaled by the arm lengths
  const E = an.crystal ? (a.k || 10) : (nodeLook[patch.verts[ia].type]?.arm ?? 0) + (nodeLook[patch.verts[ib].type]?.arm ?? 0) + linkLen;
  viewer.showAssembly(patch, {
    edge: E, nodeLook, dim: an.dim, direct: !lp,
    P: { ...state.P, flareR: effectiveFlareRadius(state.P) },
    linker: lp ? { len: linkLen, hub: lp.spec.hubRadius, bulges: lp.spec.arms[0].bulges.map(b => b.t), colour: LINKER_COLOUR, model: modelLook(lp.spec) } : null,
  });
  setWarnings([]);
  const counts = {};
  let magnets = 0;
  for (const v of patch.verts) {
    const nl = nodeLook[v.type]; if (!nl) continue;
    const pol = !lp && v.colour === 1 ? 'S' : 'N';
    const key = `${part(nl.partId).spec.name}|${pol}`;
    counts[key] = (counts[key] || 0) + 1;
    magnets += an.types[v.type].cn;
  }
  if (lp) magnets += patch.edges.length * 2;
  const reps = a.reps || 2;
  $('#bom').innerHTML = `<b>Bill of materials</b> <span class="hint">(${an.dim === 3 ? `${reps}×${reps}×${reps}` : `${reps}×${reps}`} cells)</span>
    <table>${Object.entries(counts).map(([k, n]) => { const [nm, pol] = k.split('|'); return `<tr><td>${esc(nm)}</td><td>${polTag(pol)}</td><td>× ${n}</td></tr>`; }).join('')}
    ${lp ? `<tr><td>${esc(lp.spec.name)}</td><td>${polTag('S')}</td><td>× ${patch.edges.length}</td></tr>` : ''}
    <tr><td>Magnets</td><td></td><td>× ${magnets}</td></tr></table>
    <div class="hint">${an.crystal ? `Placed as in the crystal, scaled ${fmt(E, 1)} mm per Å.` : `Edge (centre–centre) ${fmt(E, 0)} mm.`} ${lp ? 'Node magnets N-out, linker magnets S-out: nodes only stick to linkers.' : 'Alternate N-out / S-out nodes (dark / light in the preview).'}</div>`;
  $('#viewInfo').textContent = `${a.net === 'cif' ? 'crystal net' : a.net.replace('_t', '')} — ${patch.verts.length} nodes, ${patch.edges.length} ${lp ? 'linkers' : 'direct bonds'}`;
}
$('#asmReps').addEventListener('change', e => { if (state.asm) { state.asm.reps = +e.target.value; save(); refreshView(); } });

function assemblyGuide() {
  const P = state.P;
  return [
    'RETICUBLOX — assembly guide', '',
    state.asm ? `Framework: ${state.builder === 'cif' && state.cif ? `from ${state.cif.name}.cif, ` : CLS[W.cls] ? CLS[W.cls] + ' ' : ''}${net(state.asm.net).dim}D, ${state.asm.net === 'cif' ? `net taken from the crystal, scaled ${state.asm.k} mm per Å` : `net ${state.asm.net.replace('_t', '')} (${NETS[state.asm.net].title}, ${NETS[state.asm.net].sg})`}` : '',
    `Magnets: ${P.magnetDiameter} x ${P.magnetThickness} mm discs (pocket ${(P.magnetDiameter + 2 * P.clearance).toFixed(2)} x ${(P.magnetThickness + P.clearance).toFixed(2)} mm).`, '',
    'POLARITY (important):',
    '  "N out": insert every magnet with its NORTH pole facing outwards. "S out": SOUTH pole outwards.',
    '  "N & S out": print this block twice as often as one colour class and fill half with N-out, half with S-out',
    '  magnets (counts in the bill of materials). Tip: stack all magnets in one column and mark the exposed face.',
    '  Press-fit, or add a drop of CA glue.', '',
    'PRINTING: planar blocks and flat molecular models have a flat underside (no supports). 3D pieces print with tree supports.',
    '  Molecular models: the magnet fittings stand for the bonds formed on assembly (imine, boronate, metal-carboxylate...).', '',
    'PARTS:', ...state.parts.map(p => `  - ${p.spec.name} [${p.spec.kind}, ${p.spec.arms.length} arms, ${p.spec.polarity}-out] ${describeArms(p.spec)}`),
    ...($('#bom').innerText ? ['', 'BILL OF MATERIALS:', $('#bom').innerText] : []),
  ].join('\n');
}

// ---------------- print settings ----------------
const PRINT_FIELDS = [
  ['rodRadius', 'Rod radius', 'mm'], ['flareRadius', 'Flare radius', 'mm'], ['flareLength', 'Flare length', 'mm'],
  ['hubRadius', 'Node hub radius', 'mm'], ['bulgeRadius', 'Ring bulge radius', 'mm'],
  ['magnetDiameter', 'Magnet diameter', 'mm'], ['magnetThickness', 'Magnet thickness', 'mm'],
  ['clearance', 'Pocket clearance', 'mm'], ['flatBase', 'Flat-base cut', 'mm'], ['atomRadius', 'Atom radius (molecular models)', 'mm'], ['segments', 'Smoothness', 'segments'],
];
$('#printForm').innerHTML = PRINT_FIELDS.map(([k, l, u]) =>
  `<label>${l} <span class="unit">${u}</span><input type="number" step="${k === 'segments' ? 8 : 0.05}" data-k="${k}" value="${state.P[k]}"></label>`).join('');
$('#applyPrint').addEventListener('click', () => {
  $$('#printForm input').forEach(i => (state.P[i.dataset.k] = +i.value));
  state.parts.forEach(p => { restyle(p.spec); p.spec.notes = []; ensurePrintable(p.spec); });
  meshCache.clear();
  save(); renderParts(); renderWizard(); refreshView();
});

$('#viewMode').addEventListener('click', e => {
  const b = e.target.closest('button'); if (!b) return;
  state.view = b.dataset.v;
  if (state.view === 'block' && !state.preview && !state.selected && state.parts[0]) state.selected = state.parts[0].id;
  renderParts(); refreshView();
});

// ---------------- boot ----------------
(async function boot() {
  viewer = new Viewer($('#viewer'));
  registerCifNet();
  if (W.net && !NETS[W.net]) Object.assign(W, freshWiz());
  renderWizard();
  renderParts();
  if (state.asm) $('#asmReps').value = state.asm.reps;
  await initManifold();
  status('Loading chemistry toolkit…');
  state.view = state.asm ? 'assembly' : 'block';
  refreshView();
  try {
    await initRDKit();
    rdkitReady = true;
    status('Ready');
    renderWizard();
  } catch (e) {
    console.error(e);
    status('Chemistry toolkit failed to load — generic blocks still work');
  }
  window.__reticublox = { state, W, analyseNet, buildPatch, slotSpec, GEOMETRIES };
})();
