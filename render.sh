#!/usr/bin/env bash
set -e

SCENE="${1:-dual}"
shift || true

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
      --spp 16 --width 1000 --height 1000 \
      -o "$OUT" "$@"
    ;;

  profile_gpis)
    ./build/gpis_hair \
      assets/hair.m3hair \
      --rotate 0:90 \
      --spp 16 --width 1000 --height 1000 \
      --mode raymarching \
      -o "$OUT" "$@"
    ;;

  curl_gpis)
    ./build/gpis_hair \
      assets/curl.m3hair \
      --spp 16 --width 1000 --height 1000 \
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
    echo "Available: instant | test | test_gpis | dual | curls | hair | compare | profile | profile_gpis | curl_gpis | gpis"
    exit 1
    ;;
esac

magick "$OUT" "$PNG"
xdg-open "$PNG" 2>/dev/null || open "$PNG" 2>/dev/null || true
