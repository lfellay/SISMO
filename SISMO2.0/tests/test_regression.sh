#!/bin/sh
set -eu

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
  printf 'Usage: %s SISMO MODEL REFERENCE [CONFIG]\n' "$0" >&2
  exit 2
fi

SISMO=$1
MODEL=$2
REFERENCE=$3
CONFIG=${4:-}

[ -f "$MODEL" ] || {
  printf 'missing regression model: %s\n' "$MODEL" >&2
  exit 2
}
[ -f "$REFERENCE" ] || {
  printf 'missing regression reference: %s\n' "$REFERENCE" >&2
  exit 2
}

TMP_BASE=${TMPDIR:-/tmp}
TEST_DIR=$(mktemp -d "$TMP_BASE/sismo-regression.XXXXXX")
trap 'rm -rf "$TEST_DIR"' EXIT HUP INT TERM
OUT="$TEST_DIR/result"

if [ -n "$CONFIG" ]; then
  "$SISMO" "$MODEL" "$OUT" "$CONFIG"
else
  "$SISMO" "$MODEL" "$OUT"
fi

cmp "$OUT.sismo" "$REFERENCE"
printf 'SISMO numerical regression passed: output is byte-identical to %s\n' "$REFERENCE"
