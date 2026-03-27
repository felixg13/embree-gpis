#!/usr/bin/env bash
set -e

# ---------------------------------------------------------------------------
# Usage: ./render.sh [scene] [extra gpis_hair flags...]
#
# Scenes:
#   instant    — one curl, 400x300, 1spp  (instant feedback)
#   test       — one curl, 600x400, 8spp  (quick quality check)
#   dual       — curl + straight hair side by side (default)
#   curls      — four curls in a row with spacing
#   hair       — straight hair asset alone, close up
#   compare    — same curl twice, for BSDF comparison
#   profile    — two hair meshes side by side; right one rotated 90° for side view
#   gpis       — one curl, raymarching mode, 400x300, 4spp
# ---------------------------------------------------------------------------

SCENE="${1:-dual}"
shift || true   # remaining args forwarded to gpis_hair

cmake --build build 2>&1 | grep -v "^\[" || true

OUT="renders/${SCENE}.ppm"
PNG="renders/${SCENE}.png"

case "$SCENE" in
  instant)
    ./build/gpis_hair \
      assets/curl.m3hair \
      --spp 1 --width 400 --height 300 \
      -o "$OUT" "$@"
    ;;

  test)
    ./build/gpis_hair \
      assets/curl.m3hair \
      --spp 8 --width 400 --height 400 \
      -o "$OUT" "$@"
    ;;

  test_gpis)
    ./build/gpis_hair \
      assets/curl.m3hair \
      --spp 8 --width 400 --height 400 \
      --mode raymarching \
      -o "$OUT" "$@"
    ;;

  dual)
    ./build/gpis_hair \
      assets/curl.m3hair assets/hair.m3hair \
      --spacing 13 \
      --spp 16 --width 1200 --height 600 \
      -o "$OUT" "$@"
    ;;

  curls)
    ./build/gpis_hair \
      assets/curl.m3hair assets/curl.m3hair \
      assets/curl.m3hair assets/curl.m3hair \
      --spacing 2.5 \
      --spp 16 --width 1200 --height 600 \
      -o "$OUT" "$@"
    ;;

  hair)
    ./build/gpis_hair \
      assets/hair.m3hair \
      --spacing 0 \
      --spp 16 --width 800 --height 800 \
      -o "$OUT" "$@"
    ;;

  compare)
    ./build/gpis_hair \
      assets/curl.m3hair assets/curl.m3hair \
      --spacing 3.0 \
      --spp 64 --width 1000 --height 600 \
      -o "$OUT" "$@"
    ;;

  profile)
    ./build/gpis_hair \
      assets/hair.m3hair \
      --rotate 0:90 \
      --spp 16 --width 400 --height 400 \
      -o "$OUT" "$@"
    ;;

  profile_gpis)
    ./build/gpis_hair \
      assets/hair.m3hair \
      --rotate 0:90 \
      --spp 16 --width 400 --height 400 \
      --mode raymarching \
      -o "$OUT" "$@"
    ;;

  gpis)
    ./build/gpis_hair \
      assets/curl.m3hair \
      --spp 4 --width 800 --height 600 \
      --mode raymarching \
      -o "$OUT" "$@"
    ;;

  *)
    echo "Unknown scene: $SCENE"
    echo "Available: instant | test | dual | curls | hair | compare | profile | gpis"
    exit 1
    ;;
esac

magick "$OUT" "$PNG"
xdg-open "$PNG" 2>/dev/null || open "$PNG" 2>/dev/null || true
