// Reticublox desktop app (Electron). Serves the web app from inside the app bundle on a private
// app:// scheme (WebAssembly and ES modules need a real origin; file:// won't do), in its own window.
const { app, BrowserWindow, protocol, shell, Menu } = require('electron');
const path = require('path');
const fs = require('fs/promises');

const ROOT = __dirname; // index.html, css/, js/, vendor/ sit next to this file inside the bundle
const TYPES = {
  '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css',
  '.wasm': 'application/wasm', '.json': 'application/json', '.svg': 'image/svg+xml', '.png': 'image/png',
  '.gif': 'image/gif', '.jpg': 'image/jpeg', '.md': 'text/markdown', '.txt': 'text/plain',
};

protocol.registerSchemesAsPrivileged([{
  scheme: 'app',
  privileges: { standard: true, secure: true, supportFetchAPI: true, corsEnabled: true, stream: true },
}]);

function createWindow() {
  const win = new BrowserWindow({
    width: 1440, height: 900, minWidth: 900, minHeight: 600,
    title: 'Reticublox', backgroundColor: '#14171c', show: false,
    webPreferences: { contextIsolation: true, sandbox: true },
  });
  if (process.env.RETICUBLOX_SELFTEST) selfTest(win);
  else win.once('ready-to-show', () => win.show());
  // links to the web open in the normal browser, never inside the app
  win.webContents.setWindowOpenHandler(({ url }) => {
    if (/^https?:/.test(url)) shell.openExternal(url);
    return { action: 'deny' };
  });
  win.webContents.on('will-navigate', (e, url) => {
    if (!url.startsWith('app://')) { e.preventDefault(); if (/^https?:/.test(url)) shell.openExternal(url); }
  });
  win.loadURL('app://reticublox/index.html');
}

// RETICUBLOX_SELFTEST=1 …/Reticublox: load the app hidden, exercise the WebAssembly libraries and
// the structure editor, print the result and quit (used to check a fresh build).
function selfTest(win) {
  win.webContents.on('console-message', (_e, level, msg) => { if (level >= 3) console.log('page error:', msg); });
  win.webContents.once('did-finish-load', async () => {
    try {
      const r = await win.webContents.executeJavaScript(`(async () => {
        const until = async (f, ms = 30000) => { const t = Date.now(); while (!(await f())) { if (Date.now() - t > ms) return false; await new Promise(r => setTimeout(r, 200)); } return true; };
        const status = await until(() => /Ready/.test(document.querySelector('#status')?.textContent || ''));
        const chem = await import('./js/chem.js'), geo = await import('./js/geometry.js'), imp = await import('./js/import3d.js');
        const mm = chem.parseMolecule('OC(=O)c1cc(C(O)=O)cc(C(O)=O)c1');
        const g = chem.autoDetectGroup(mm);
        const P = { flareLength: 10, rodRadius: 3, flareRadius: 5, hubRadius: 6, atomRadius: 5.5, magnetDiameter: 6, magnetThickness: 3, clearance: 0.15, flatBase: 0.8, segments: 32, bulgeRadius: 6 };
        const spec = chem.sbuModelSpec('Cu₂/Zn₂ paddlewheel (4-c square)', { scale: 10, P });
        const s = geo.buildSolid(spec, P); const solid = s.decompose().length + ' shell(s), genus ' + s.genus(); s.delete();
        const xyz = '3\\nwater-free test\\nC 0 0 0\\nO 1.21 0 0\\nC -0.7 1.2 0.3\\n';
        const i3 = imp.importStructure(xyz, 't.xyz').molblock.includes('M  END');
        const sc = document.createElement('script'); sc.src = 'vendor/jsme/jsme.nocache.js'; document.body.append(sc);
        const jsme = await until(() => !!window.JSApplet, 20000);
        const builders = [...document.querySelectorAll('.seg.builder button')].map(b => b.textContent).join(' | ');
        return { status, group: g, solid, import3d: i3, cifBuilder: typeof imp.importFramework === 'function' && /CIF/.test(builders), jsme, zip: !!window.JSZip };
      })()`);
      console.log('SELFTEST', JSON.stringify(r));
    } catch (e) { console.log('SELFTEST failed', e.message); }
    app.quit();
  });
}

app.whenReady().then(() => {
  protocol.handle('app', async req => {
    const rel = decodeURIComponent(new URL(req.url).pathname);
    const file = path.normalize(path.join(ROOT, rel));
    if (!file.startsWith(ROOT + path.sep)) return new Response('Not found', { status: 404 });
    try {
      const body = await fs.readFile(file);
      return new Response(body, { headers: { 'content-type': TYPES[path.extname(file).toLowerCase()] || 'application/octet-stream' } });
    } catch {
      return new Response('Not found', { status: 404 });
    }
  });
  Menu.setApplicationMenu(Menu.buildFromTemplate([
    { role: 'appMenu' }, { role: 'fileMenu' }, { role: 'editMenu' },
    { label: 'View', submenu: [{ role: 'reload' }, { role: 'toggleDevTools' }, { type: 'separator' }, { role: 'resetZoom' }, { role: 'zoomIn' }, { role: 'zoomOut' }, { type: 'separator' }, { role: 'togglefullscreen' }] },
    { role: 'windowMenu' },
    { role: 'help', submenu: [{ label: 'Reticublox README', click: () => shell.openPath(path.join(ROOT, 'README.md')) }] },
  ]));
  createWindow();
  app.on('activate', () => { if (!BrowserWindow.getAllWindows().length) createWindow(); });
});
app.on('window-all-closed', () => { if (process.platform !== 'darwin') app.quit(); });
