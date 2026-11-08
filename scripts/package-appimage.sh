#!/usr/bin/env bash
# Builds a Linux AppImage for the xtalforge CLI.
#
# Usage: scripts/package-appimage.sh [OUTPUT_DIR]
# Requires target/release/xtalforge and curl; downloads appimagetool on first run.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VERSION=$(sed -n 's/^version = "\([^"]*\)"/\1/p' "$ROOT/Cargo.toml" | head -n 1)
ARCH=${ARCH:-x86_64}
OUT=${1:-"$ROOT/artifacts"}
BIN="$ROOT/target/release/xtalforge"
TOOL="$ROOT/target/appimagetool-$ARCH.AppImage"

[ -x "$BIN" ] || { echo "package-appimage: build $BIN first" >&2; exit 2; }
[ -f "$ROOT/assets/xtalforge.desktop" ] || { echo "package-appimage: missing desktop file" >&2; exit 2; }
[ -f "$ROOT/assets/xtalforge.png" ] || { echo "package-appimage: missing icon" >&2; exit 2; }

mkdir -p "$OUT"
appdir=$(mktemp -d)
chmod 755 "$appdir"
trap 'rm -rf "$appdir"' EXIT

install -Dm755 "$BIN" "$appdir/usr/bin/xtalforge"
install -Dm644 "$ROOT/assets/xtalforge.desktop" \
  "$appdir/usr/share/applications/xtalforge.desktop"
install -Dm644 "$ROOT/assets/xtalforge.png" \
  "$appdir/usr/share/icons/hicolor/256x256/apps/xtalforge.png"
install -Dm644 "$ROOT/assets/io.github.jigonzalez930209.xtalforge.appdata.xml" \
  "$appdir/usr/share/metainfo/io.github.jigonzalez930209.xtalforge.appdata.xml"
ln -s usr/share/applications/xtalforge.desktop "$appdir/xtalforge.desktop"
ln -s usr/share/icons/hicolor/256x256/apps/xtalforge.png "$appdir/xtalforge.png"

cat > "$appdir/AppRun" <<'APPRUN'
#!/bin/sh
HERE=$(dirname "$(readlink -f "$0")")
exec "$HERE/usr/bin/xtalforge" "$@"
APPRUN
chmod +x "$appdir/AppRun"

if [ ! -x "$TOOL" ]; then
  curl --fail --location --retry 3 --output "$TOOL" \
    "https://github.com/AppImage/appimagetool/releases/download/continuous/appimagetool-$ARCH.AppImage"
  chmod +x "$TOOL"
fi

APPIMAGE_EXTRACT_AND_RUN=1 ARCH="$ARCH" "$TOOL" \
  "$appdir" "$OUT/xtalforge-$VERSION-$ARCH.AppImage"
echo "wrote $OUT/xtalforge-$VERSION-$ARCH.AppImage"
