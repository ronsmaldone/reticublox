# Reticublox Builder

Design 3D-printable **node** and **linker** blocks that snap together with magnets into physical models of
covalent organic frameworks (COFs) and metal–organic frameworks (MOFs).

## Use it online
**https://ronsmaldone.github.io/reticublox/** runs in any modern browser (Chrome, Edge, Firefox or Safari) on Mac, Windows or Linux, with nothing to install. Your parts and settings are saved in that browser.

## Download
**Mac app:** download `Reticublox-mac-arm64.zip` (Apple Silicon, M1 and later) or `Reticublox-mac-x64.zip` (Intel Macs) from the [latest release](https://github.com/ronsmaldone/reticublox/releases/latest), unzip it and drag Reticublox to Applications. The first time, right-click the app and choose **Open**. Any other computer can run it from the source code; see below.

## Run it

```bash
cd reticublox
python3 serve.py
```

This opens http://localhost:8765 in your browser. Everything runs locally. The libraries are in `vendor/`, so you don't need an internet connection. Opening `index.html` directly (`file://`) won't work, because WebAssembly modules need to be served over HTTP.

## Desktop app (Mac)
`desktop/build-mac.sh` builds a standalone Mac app. It doesn't need Python, a server or an internet connection: drag it to Applications and open it like any other app. There are two builds, each with a zip for sharing:

| Mac | App | Zip |
|---|---|---|
| Apple Silicon (M1 and later) | `dist/mac-arm64/Reticublox.app` | `dist/Reticublox-mac-arm64.zip` |
| Intel | `dist/mac-x64/Reticublox.app` | `dist/Reticublox-mac-x64.zip` |

To see which one a Mac needs, open  → About This Mac: **Chip** means Apple Silicon and **Processor** means Intel. The Intel build also runs on Apple Silicon through Rosetta, but more slowly.

- **First launch:** the app is signed ad hoc, not with a paid Apple Developer ID, so macOS says it's from an unidentified developer. Right-click the app and choose **Open**, then **Open** again. You only need to do this once per Mac. If macOS says the app "is damaged", run `xattr -cr /Applications/Reticublox.app` once; this happens to apps downloaded through some browsers or mail clients.
- **Saved work:** the app keeps its own parts and settings, separate from the browser version.
- **Rebuilding after changes to the web app:** run `desktop/build-mac.sh`, or `desktop/build-mac.sh x64` for Intel Macs.
  - The first build downloads Electron's prebuilt Mac app (about 105 MB) from GitHub into `desktop/.cache/`. Node.js isn't needed.
  - To check a fresh build without opening a window, run `RETICUBLOX_SELFTEST=1 dist/mac-arm64/Reticublox.app/Contents/MacOS/Reticublox`.
- **Files:** `desktop/main.cjs` is the app window. It serves the files from inside the app on a private `app://` address, because WebAssembly needs a real origin. `desktop/icon.svg` and `desktop/icon.icns` are the app icon.

## Building a model in four steps

The left panel walks you through four steps. Each finished step collapses to a one-line summary; click it to change your choice.

1. **Framework type: MOF or COF.** This sets which example structures are suggested. The metal-cluster (SBU) option appears only for MOFs.
2. **Dimensionality: 2D or 3D.**
3. **Topology.** Browse **by net** or **by space group** (plane group for 2D layers). Each net shows its space group, connectivity and real MOF or COF examples. Nets with examples for your framework type are listed first.
   - 2D nets: sql, hcb, kgm, hxl, kgd
   - 3D nets: dia, pcu, bcu, fcu, nbo, tbo, ctn, bor, pts (two variants), stp, srs, acs
4. **Building blocks.**
   - First choose whether nodes connect **through a linker** or **directly node–node**. Direct bonding is only offered for nets where it works, and the app sets a sensible default.
   - There is one card per node type, plus one for the linker. Each card can be:
     - **Generic:** the shape comes from the topology, with arms exactly along the net's edges. You set the arm length (or linker length) and the number of bulges.
     - **Molecule:** pick one of the suggested building blocks, type a SMILES, draw one in the built-in editor, or **import a 3D structure** from another program (see below). Reactive groups (CHO, NH₂, COOH, B(OH)₂, nitrile, pyridyl, `*` atoms) are detected automatically. Click an atom in the drawing to make it a connection point, and click it again to remove it. The chosen points are also listed under the drawing, each with an × to remove it, plus **Clear all**.
     - **Metal cluster (SBU), MOFs only:** paddlewheel, Zn₄O, M₃O trimer, Zr₆ (12-c or 8-c), or a tetrahedral centre. Only clusters with the right connectivity can be chosen.
   - A Molecule card can be printed two ways (**Print as**):
     - **Reticular block:** hub, arms and ring bulges, which is the molecule reduced to its topology.
     - **Molecular model:** the actual atoms and bonds of the drawn structure, with a magnet fitting at each connection point. See below.
   - Every card says whether its block fits the net. For molecules and SBUs, **Arm directions → Match the net** (the default) snaps the arms onto the net's exact angles. The molecule still sets the arm lengths and ring bulges.
   - Press **Build model & parts**. This makes the parts list, the assembled 3D preview, and a bill of materials with magnet polarities.

## How molecules become blocks
- **Arm length** is the distance from the core centre to the point of extension, plus an *extension*. The extension is half of the bond made on assembly; for example, 0.64 Å for an imine. Scale is set in mm per Å (default 6).
- **Bulges**: a ring whose centre lies on an arm becomes a sphere at that position. A ring at the centre becomes the hub.
- **3D cores**: 2D drawings flatten 3D molecules. When a quaternary sp³ C or Si sits at the centre (as in TAM), the app switches to tetrahedral geometry and equalises the arms. For any other shape, use the **Arm directions** setting.
- If an arm is too short to hold a magnet, it's lengthened to the printable minimum, and the app tells you.

## Molecular models
Set a Molecule card's **Print as** to **Molecular model** to print the molecule itself as a stick model. You can still use Generic cards (and SBU cards) for the other parts, so molecules and topology blocks can be mixed in one model.
- **Atoms** are spheres (*Atom radius* under Print settings, 5.5 mm by default; heavier atoms are 15% larger). **Bonds** are rods with the same radius as the block arms, so a flat molecule's underside is cut flat and prints without supports. Bond orders aren't shown and hydrogens are left out.
- **Magnet fittings** replace the reactive groups. Each fitting is a rod and flared magnet socket, and it stands for the bond made on assembly:
  - imine C=N for aldehydes and amines;
  - B–O for boronic acids;
  - triazine C–N for nitriles;
  - M–N for pyridyl groups;
  - metal–carboxylate coordination for carboxylic acids.
- Atoms that leave when the bond forms (the aldehyde O, the boronic acid OH groups, the nitrile N, `*` atoms) are greyed out in the drawing and not printed. Carboxylate O atoms stay, with the fitting between them where the metal binds.
- **Arm directions → Match the net**: the magnet faces sit exactly where a generic block's would, so molecules, generic blocks and SBUs all fit together.
  - **Flat (conjugated) molecules stay flat.** When the net's vertex is slightly pyramidal (BTC in tbo, 10°; TAPB in bor, 19.5°), the flat molecule sits a little off the vertex point, as it does in the real crystal. Only the magnet sockets tilt to meet the net's edges. The model still prints lying flat without supports: it's turned over so the sockets tilt upwards, and the pockets stay closed.
  - **A drawn sp³ centre (TAM, tetraphenylsilane) is opened up into a real tetrahedron.** Each arm is turned about the central atom. These models print with supports.
  - If a flat molecule can't reach the net (for example a square porphyrin in the rectangular *pts — tetrahedron kept ideal* variant), the card shows ⚠ and how far the fittings would have to bend.
- Every fitting is long enough for its flare to clear its atom. Because of that, molecular arms are usually a little longer than the chemistry alone would give.
- **Scale:** at 6 mm per Å, the atoms touch and rings close up. When you switch a card to Molecular model, the scale goes up to 10 mm per Å if it was lower. 9–12 mm per Å works well, and the card warns if atoms touch.
- The assembled-model preview shows each molecule posed at its site.
- **Metal clusters (SBU cards) can be molecular models too** (Print as → Molecular model). The metal core is printed as atoms and bonds: the Cu–Cu pair of a paddlewheel, Zn₄ around μ₄-O, M₃ around μ₃-O, the Zr₆ octahedron, or a single tetrahedral centre. Each magnet fitting stands for one carboxylate and forks from the metals it bridges (both Cu atoms in a paddlewheel, a Zn pair in Zn₄O). The carboxylate C and O atoms come with the linker model, so a BTC model plus a paddlewheel model assemble into the full Cu₂(O₂C)₄ unit. Cluster models are 3D and print with tree supports.

## Importing 3D structures
SMILES and drawings are flat, so molecules that are really twisted come out flat: in hexaphenylbenzene, for example, the six outer rings stand at 65–90° like a propeller. To keep the real shape, use **⤓ Import 3D structure…** on a Molecule card and choose a file exported from CrystalMaker, Mercury, Avogadro, VESTA, Materials Studio or similar.

| File | What is read |
|---|---|
| `.mol`, `.sdf` | Atoms, bonds and bond orders (the first molecule in the file). This is the most reliable format. |
| `.mol2` | Atoms, bonds and bond orders. |
| `.xyz`, `.pdb` | Atoms only. Bonds are found from the distances between atoms, and bond orders from the bond lengths and angles. Hydrogens are optional. |
| `.cif` | The atoms are expanded by the space-group symmetry and molecules are rebuilt across the cell edges. Bonds are then found as for `.xyz`. Disordered sites with occupancy below 0.5 are skipped. |

- **The real 3D shape is kept,** including twisted rings and true bond angles.
  - **Match the net** turns the whole molecule onto its site without bending it. If its connection points don't quite match the net, the magnet fittings bend instead, and the card shows by how much.
  - Molecular models of non-flat molecules print with tree supports.
- **Connection points:** the reactive groups are found as usual. You can still click atoms in the drawing, which shows the molecule projected onto its best-fit plane.
- **Hydrogens** are left out, as for SMILES.
- **Metal atoms** are left out; use a Metal cluster card for the SBU. Importing a **MOF CIF** therefore gives you its organic linker.
  - A CIF of a COF, or any other endless covalent network, can't be split into one molecule here. Use **Start from a CIF** (below) instead.
- **Several molecules in one file** (for example solvent): the largest is used.
- **Flat files** (all z = 0) are read as a 2D structure, like a SMILES.
- Typing a SMILES, picking a preset or drawing replaces the imported structure.

## Starting from a CIF (whole COF or MOF)
The builder has two modes, switched at the top of the panel:

- **Build from scratch** (the four steps above) makes geometric, generic blocks from a chosen net, or from molecules you pick or draw.
- **Start from a CIF** takes a COF or MOF crystal structure, for example a geometry-optimised cell from Materials Studio or a MOF-5 CIF from CrystalMaker. It places the magnet connections exactly where the blocks were bonded in the crystal, and works out everything else itself:

1. **It rebuilds the bond network of the cell,** including the bonds that cross the cell edges and the metal–ligand bonds.
2. **It cuts the linkage bonds.**
   - **COFs:** imine C=N (also hydrazone and azine) and the C–N bond of β-ketoenamines.
   - **MOFs:** each carboxylate's C–C bond, so the carboxylate stays with the metal cluster as the point of extension (MOF-5 gives a Zn₄O(CO₂)₆ SBU and a C₆H₄ linker). Also metal–N bonds of pyridyl and azolate linkers, and metal–O bonds to hydroxyl and ether oxygens (CD-MOF-1 gives K⁺ and γ-cyclodextrin blocks). Coordinated solvent (water on the metals) is left out.
   - Several bonds between the same two pieces (a K⁺ held by two oxygens of one cyclodextrin) make one connection, with a forked fitting.
3. **Each piece is one copy of a building block.** A piece with three or more linkages is a node (an SBU if it holds metals); one with two is a linker. Copies with the same formula are grouped, and the most typical copy of each is used.
4. **It finds the net** from how the nodes connect (its dimension and coordination sequences), and picks it from the library, for example hcb. If no library net matches (CD-MOF-1's 4,8-connected net, say), it uses the crystal's own net instead: see below.

It then shows what it found (the net, each node and linker with its formula, SMILES and fit), lets you set the model scale, and builds the parts and the assembled model. The two builders keep their own settings, so switching between them loses nothing.

Every block is a molecular model that **keeps the conformation from the file**: ring twists, side chains (alkyl and alkoxy chains are printed as computed) and any buckling of the layer. The magnet fittings start at the atoms of the cut bonds (for example N on the amine side and C on the aldehyde side of an imine, the carboxylate C and the ring C in a MOF, or Zn and N in ZIF-8) and run exactly along the edges of the ideal net, with equal arms. That is what makes the pieces fit together: a crystal cell's copies of the same block are rarely identical (in an optimised COF they are often near mirror images of each other), so one printed piece can only fill every site if its magnets follow the ideal net rather than one copy's own bond directions. The molecule itself is turned as a whole to fit its site as closely as possible; where an atom sits off its edge, its fitting rod bends to meet the edge (at most 20°, by lengthening the arm if needed).

- **Nets not in the library:** the blocks are placed exactly as they sit in the crystal, and each fitting runs from its atoms to the line joining the two block centres, level with the middle of the cut bond. Block positions and arms are the crystal scaled up by one factor (at least the model scale, more if a fitting needs room to clear its atoms; the panel shows it), so every magnet pair meets face to face. This needs a symmetric crystal structure (an experimental CIF); fittings bend by up to 35°.
- **Large cells:** CIFs that list the whole cell, or high-symmetry space groups, are fine (CD-MOF-1 has about 2,100 atoms per cell).
- **Several copies** of the same node or linker in the cell are never exactly alike after optimisation. The note after the import says by how much they differ.
- **Stacked layers** or interpenetrated nets in one cell are the same framework, so one is modelled.
- **Scale:** 10 mm per Å to start with. A COF with large pores makes big pieces at that scale; the panel shows the node-to-node distance so you can pick a smaller one.
- **Not supported yet:** boronate-ester, boroxine and triazine COFs; phosphonate, sulfonate and catecholate MOFs; and nets that aren't in the library.
- Use the **CIF**, not a MOL2 or PDB export of a supercell. Those files lose the bonds that cross the cell edges, so the framework falls apart.

## Magnets and polarity
- Each arm ends in a flared tip with a blind pocket for a disc magnet. The default is 6 × 3 mm with 0.15 mm clearance; you can change this under **Print settings**.
- **Node** magnets go in **N-out** and **linker** magnets go in **S-out**. Nodes then attract linkers and repel each other.
- In direct node–node mode, the nodes alternate between N-out and S-out. The bill of materials lists how many of each to make.
- Each STL filename includes its polarity, and `ASSEMBLY_GUIDE.txt` is included in the ZIP.

## HKUST-1 (tbo)
1. Choose **MOF → 3D → tbo**. **Directly node–node** is selected automatically.
2. Set node T to **Molecule → Trimesic acid (BTC)** and node S to **Metal cluster → Cu₂/Zn₂ paddlewheel**, or leave both as Generic. Set both to **Print as → Molecular model** for an atom-and-bond model. Then press **Build model & parts**. BTC goes in N-out and the paddlewheel S-out.
3. One unit cell takes 32 BTC nodes, 24 paddlewheels and 192 magnets.

In tbo, the arms can't all be straight and equal while the paddlewheel is square *and* BTC is flat. The generated BTC block is therefore slightly pyramidal (117° between arms, each arm about 10° out of plane). It prints standing on its arm tips. In the real crystal, the bend is at the carboxylate instead.

## 3D COFs with tetrahedral nodes (dia, ctn, bor)
Use **Node – node (direct)**. The triangular node is N-out and the tetrahedral node is S-out.

| Net | Examples | Node geometry |
|---|---|---|
| dia | COF-300 (with a linker) | Ideal tetrahedron |
| ctn | COF-102, COF-103 | Nearly ideal: triangles flattened to exactly planar, tetrahedra within 2° |
| bor | COF-105, COF-108 | Tetrahedral node kept ideal (the same block as in dia), so the 3-c node is pyramidal (arms 19.5° out of plane) |

| pts | MOF-11, 3D porphyrin COFs | Comes in two variants; see below |
| stp | PCN-600 (Fe₃O trimer + TCPP) | Both exact: square, and a trigonal prism with arms 45° above and below the plane. This is the same prism as in acs, so one M₃O-trimer block works in both MIL-88 (acs) and PCN-600 (stp). |

**pts:** with straight equal edges, the square node and the tetrahedral node can't both be ideal.
- **pts — square kept ideal** (c/a = 2): the tetrahedral node becomes 90°/90°/120°. Use this with porphyrin or paddlewheel squares.
- **pts — tetrahedron kept ideal** (c/a = √2): the square becomes a 70.5° × 109.5° rectangle. Use this to reuse your dia/bor/ctn tetrahedral block.

In bor's cubic cell, straight sticks can't make the triangle flat *and* the tetrahedron ideal. If a molecule card's **Arm directions** is set to *From the 2D structure*, the card warns that a flat triangle doesn't fit. Keep it on *Match the net*.

The Assembly preview closes rings that cross the cell boundary and removes nodes that would hang by a single bond. The patch it shows, and its bill of materials, is therefore one connected model you can build.

## Printing
- Planar blocks get a thin flat cut on the underside and print without supports.
- 3D blocks (tetrahedral, octahedral and others) print with one arm pointing straight up. Use tree supports.

## Code map
- `js/nets.js`: net library, edge detection, coordination checks, patch generation and bipartite colouring
- `js/geometry.js`: solids using Manifold CSG (block rods, flares and bulges; molecular-model atoms, bonds and fittings; magnet pockets; flat base) and the STL writer
- `js/chem.js`: RDKit.js parsing, ring and connection-group detection, molecule-to-block and molecule-to-model conversion, presets and SBUs
- `js/import3d.js`: 3D file import (MOL/SDF, MOL2, XYZ, PDB, CIF with symmetry expansion), bond and bond-order perception; whole-framework import (cut the linkages, group the monomers, match the net)
- `desktop/`: the Mac app (Electron shell, icon, build script)
- `js/rot.js`: rotation helpers (matching a set of arm directions to a site in the net)
- `js/viewer.js`: three.js viewer for single blocks and the assembled patch
- `js/main.js`: the four-step builder, parts list, assembly preview and export

To add a net, add an entry to `NETS` with its cell, vertex fractional coordinates and expected coordination numbers. At load time the app checks every vertex against those coordination numbers.

## Third-party libraries (in `vendor/`)
- three.js (MIT)
- manifold-3d (Apache-2.0)
- RDKit.js (BSD-3)
- JSME (BSD-3)
- JSZip (MIT)
