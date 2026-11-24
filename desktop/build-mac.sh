#!/bin/bash
# Build Reticublox.app — a standalone Mac app — from the web app in the folder above.
# No Node.js needed: it takes Electron's prebuilt Mac app, puts Reticublox inside, renames and
# re-signs it (ad hoc), and zips it for sharing.
#
#   desktop/build-mac.sh            Apple Silicon (arm64) build
#   desktop/build-mac.sh x64        Intel Mac build
#
# Output: dist/Reticublox.app and dist/Reticublox-mac-<arch>.zip
set -euo pipefail

ELECTRON_VERSION="${ELECTRON_VERSION:-33.2.1}"
ARCH="${1:-arm64}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
CACHE="$HERE/.cache"
DIST="$ROOT/dist"
ZIP="electron-v${ELECTRON_VERSION}-darwin-${ARCH}.zip"
URL="https://github.com/electron/electron/releases/download/v${ELECTRON_VERSION}/${ZIP}"
VERSION="$(sed -n 's/.*"version": *"\([^"]*\)".*/\1/p' "$HERE/package.json")"

mkdir -p "$CACHE" "$DIST"
if [ ! -f "$CACHE/$ZIP" ]; then
  echo "Downloading Electron ${ELECTRON_VERSION} (${ARCH}) from GitHub…"
  curl -fL --progress-bar -o "$CACHE/$ZIP.part" "$URL"
  mv "$CACHE/$ZIP.part" "$CACHE/$ZIP"
fi

APP="$DIST/Reticublox.app"
rm -rf "$APP" "$DIST/Electron.app"
ditto -x -k "$CACHE/$ZIP" "$DIST"
mv "$DIST/Electron.app" "$APP"
# Electron's and Chromium's licence notices travel inside the app
mv "$DIST/LICENSE" "$DIST/LICENSES.chromium.html" "$APP/Contents/Resources/"
rm -f "$DIST/version"

# the web app goes in Contents/Resources/app (Electron loads it instead of its default app)
RES="$APP/Contents/Resources"
rm -f "$RES/default_app.asar"
mkdir -p "$RES/app"
cp -R "$ROOT/index.html" "$ROOT/css" "$ROOT/js" "$ROOT/vendor" "$ROOT/README.md" "$ROOT/LICENSE" "$ROOT/THIRD_PARTY_NOTICES.md" "$RES/app/"
cp "$HERE/main.cjs" "$HERE/package.json" "$RES/app/"
cp "$HERE/icon.icns" "$RES/electron.icns"

# name, identity and executable
PL="$APP/Contents/Info.plist"
mv "$APP/Contents/MacOS/Electron" "$APP/Contents/MacOS/Reticublox"
plutil -replace CFBundleExecutable -string Reticublox "$PL"
plutil -replace CFBundleName -string Reticublox "$PL"
plutil -replace CFBundleDisplayName -string Reticublox "$PL"
plutil -replace CFBundleIdentifier -string org.reticublox.app "$PL"
plutil -replace CFBundleShortVersionString -string "$VERSION" "$PL"
plutil -replace CFBundleVersion -string "$VERSION" "$PL"
plutil -replace NSHumanReadableCopyright -string "Reticublox — magnetic COF/MOF model blocks" "$PL"
plutil -replace LSApplicationCategoryType -string public.app-category.education "$PL"

# changing the bundle breaks Electron's signature: re-sign ad hoc (required on Apple Silicon)
xattr -cr "$APP"
codesign --force --deep --sign - "$APP"
codesign --verify --deep --strict "$APP"

OUT="$DIST/Reticublox-mac-${ARCH}.zip"
rm -f "$OUT"
ditto -c -k --sequesterRsrc --keepParent "$APP" "$OUT"
echo "Built $APP"
echo "Zipped for sharing: $OUT ($(du -h "$OUT" | cut -f1))"
