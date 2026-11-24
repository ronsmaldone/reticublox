# Third-party software

Reticublox bundles the following libraries in `vendor/`, unmodified. Each is distributed under its own licence, summarised below; the full texts are at the links given.

| Library | Version | Licence | Copyright |
|---|---|---|---|
| [three.js](https://github.com/mrdoob/three.js) (`three.module.js`, `OrbitControls.js`) | r170 | MIT | © 2010–2024 three.js authors |
| [manifold-3d](https://github.com/elalish/manifold) (`manifold.js`, `manifold.wasm`) | 3.x | [Apache-2.0](https://www.apache.org/licenses/LICENSE-2.0) | © The Manifold Authors |
| [RDKit.js](https://github.com/rdkit/rdkit-js) (`RDKit_minimal.js`, `RDKit_minimal.wasm`) | MinimalLib | BSD-3-Clause | © Rational Discovery LLC, Greg Landrum, Julie Penzotti and others |
| [JSME Molecule Editor](https://jsme-editor.github.io/) (`jsme/`) | 2024 | BSD-3-Clause | © Peter Ertl and Bruno Bienfait |
| [JSZip](https://github.com/Stuk/jszip) (`jszip.min.js`) | 3.10.1 | MIT (dual-licensed MIT or GPLv3; used under MIT) | © 2009–2016 Stuart Knightley; includes pako (MIT) |

The Mac app (`desktop/`, built by `desktop/build-mac.sh`) also contains [Electron](https://github.com/electron/electron) (MIT, © Electron contributors, © GitHub Inc.), which includes Chromium and other components under their own licences. Their notices ship inside the app at `Reticublox.app/Contents/Resources/LICENSE` and `LICENSES.chromium.html`.

## MIT licence (three.js, JSZip, Electron)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## BSD 3-Clause licence (RDKit, JSME)

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Apache License 2.0 (manifold-3d)

manifold-3d is licensed under the Apache License, Version 2.0. You may obtain a copy of the licence at https://www.apache.org/licenses/LICENSE-2.0. It is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
