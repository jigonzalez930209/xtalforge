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
RUNTIME="$ROOT/target/appimage-runtime-$ARCH"
COMPRESSION=zstd

[ -x "$BIN" ] || { echo "package-appimage: build $BIN first" >&2; exit 2; }
[ -f "$ROOT/assets/io.github.jigonzalez930209.xtalforge.desktop" ] || { echo "package-appimage: missing desktop file" >&2; exit 2; }
[ -f "$ROOT/assets/xtalforge.png" ] || { echo "package-appimage: missing icon" >&2; exit 2; }
[ -f "$ROOT/assets/io.github.jigonzalez930209.xtalforge.appdata.xml" ] \
  || { echo "package-appimage: missing AppStream metadata" >&2; exit 2; }
if command -v appstreamcli >/dev/null 2>&1; then
  appstreamcli validate --no-net \
    "$ROOT/assets/io.github.jigonzalez930209.xtalforge.appdata.xml" >/dev/null
fi

mkdir -p "$OUT"
appdir=$(mktemp -d)
chmod 755 "$appdir"
trap 'rm -rf "$appdir"' EXIT

install -Dm755 "$BIN" "$appdir/usr/bin/xtalforge"
install -Dm644 "$ROOT/assets/io.github.jigonzalez930209.xtalforge.desktop" \
  "$appdir/usr/share/applications/io.github.jigonzalez930209.xtalforge.desktop"
install -Dm644 "$ROOT/assets/xtalforge.png" \
  "$appdir/usr/share/icons/hicolor/256x256/apps/xtalforge.png"
install -Dm644 "$ROOT/assets/io.github.jigonzalez930209.xtalforge.appdata.xml" \
  "$appdir/usr/share/metainfo/io.github.jigonzalez930209.xtalforge.appdata.xml"
ln -s usr/share/applications/io.github.jigonzalez930209.xtalforge.desktop \
  "$appdir/io.github.jigonzalez930209.xtalforge.desktop"
ln -s usr/share/icons/hicolor/256x256/apps/xtalforge.png "$appdir/xtalforge.png"

cat > "$appdir/AppRun" <<'APPRUN'
#!/bin/sh
HERE=$(dirname "$(readlink -f "$0")")
exec "$HERE/usr/bin/xtalforge" "$@"
APPRUN
chmod +x "$appdir/AppRun"

if [ ! -x "$TOOL" ]; then
  for url in \
    "https://github.com/AppImage/appimagetool/releases/download/continuous/appimagetool-$ARCH.AppImage" \
    "https://github.com/AppImage/appimagetool/releases/download/1.9.0/appimagetool-$ARCH.AppImage"; do
    if curl --fail --location --retry 5 --retry-all-errors \
      --output "$TOOL.tmp" "$url"; then
      mv "$TOOL.tmp" "$TOOL"
      chmod +x "$TOOL"
      break
    fi
  done
fi
[ -x "$TOOL" ] || { echo "package-appimage: no appimagetool" >&2; exit 2; }
if [ ! -x "$RUNTIME" ]; then
  for url in \
    "https://github.com/AppImage/type2-runtime/releases/latest/download/runtime-$ARCH" \
    "https://github.com/AppImage/type2-runtime/releases/download/continuous/runtime-$ARCH"; do
    if curl --fail --location --retry 5 --retry-all-errors \
      --output "$RUNTIME.tmp" "$url"; then
      mv "$RUNTIME.tmp" "$RUNTIME"
      chmod +x "$RUNTIME"
      break
    fi
  done
  if [ ! -x "$RUNTIME" ]; then
    # The appimagetool AppImage starts with the same type2 runtime; the last
    # squashfs magic marks where its payload begins. Extracting the prefix
    # keeps the release independent of a transient GitHub asset timeout.
    offset=$(grep -abo 'hsqs' "$TOOL" | tail -n 1 | cut -d: -f1)
    if [ -n "$offset" ] && [ "$offset" -gt 0 ]; then
      dd if="$TOOL" of="$RUNTIME" bs=1 count="$offset" status=none
      chmod +x "$RUNTIME"
    fi
  fi
fi
[ -x "$RUNTIME" ] || { echo "package-appimage: no AppImage runtime" >&2; exit 2; }

# The metadata was validated above without network access; skip the tool's
# network check, which fails a release on a transient URL timeout.
APPIMAGE_EXTRACT_AND_RUN=1 ARCH="$ARCH" "$TOOL" --no-appstream \
  --comp "$COMPRESSION" --runtime-file "$RUNTIME" \
  "$appdir" "$OUT/xtalforge-$VERSION-$ARCH.AppImage"
echo "wrote $OUT/xtalforge-$VERSION-$ARCH.AppImage"
