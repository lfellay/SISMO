#!/bin/sh
set -eu

SISMO=${1:-./sismo}
case "$SISMO" in
  /*) ;;
  *) SISMO=$(CDPATH= cd -- "$(dirname -- "$SISMO")" && pwd)/$(basename -- "$SISMO") ;;
esac

TMP_BASE=${TMPDIR:-/tmp}
TEST_DIR=$(mktemp -d "$TMP_BASE/sismo-cli-test.XXXXXX")
trap 'rm -rf "$TEST_DIR"' EXIT HUP INT TERM

fail() {
  printf 'test_cli.sh: %s\n' "$1" >&2
  exit 1
}

expect_status() {
  expected=$1
  shift
  set +e
  "$@" >"$TEST_DIR/stdout" 2>"$TEST_DIR/stderr"
  status=$?
  set -e
  if [ "$status" -ne "$expected" ]; then
    printf 'expected status %s, got %s for:' "$expected" "$status" >&2
    printf ' %s' "$@" >&2
    printf '\nstdout:\n' >&2
    sed -n '1,20p' "$TEST_DIR/stdout" >&2
    printf 'stderr:\n' >&2
    sed -n '1,20p' "$TEST_DIR/stderr" >&2
    exit 1
  fi
}

expect_status 0 "$SISMO" --help
[ ! -s "$TEST_DIR/stderr" ] || fail '--help wrote to stderr'

expect_status 0 "$SISMO" missing.mod --help
expect_status 0 "$SISMO" --help ignored-extra-argument

expect_status 1 "$SISMO" missing.mod output --scan-points
grep -q 'runtime options have moved' "$TEST_DIR/stdout" ||
  fail 'obsolete option did not produce the migration message'

cat >"$TEST_DIR/unknown.conf" <<'EOF'
unknown_setting = 1
EOF
expect_status 2 "$SISMO" missing.mod output "$TEST_DIR/unknown.conf"
grep -q 'unknown configuration key' "$TEST_DIR/stderr" ||
  fail 'unknown configuration key was not diagnosed'

cat >"$TEST_DIR/duplicate.conf" <<'EOF'
scan_points = 2000
SCAN_POINTS = 2001
EOF
expect_status 2 "$SISMO" missing.mod output "$TEST_DIR/duplicate.conf"
grep -q 'duplicate key' "$TEST_DIR/stderr" ||
  fail 'duplicate configuration key was not diagnosed'

cat >"$TEST_DIR/nan.conf" <<'EOF'
tolerance = NaN
EOF
expect_status 2 "$SISMO" missing.mod output "$TEST_DIR/nan.conf"
grep -q 'valid real number' "$TEST_DIR/stderr" ||
  fail 'non-finite configuration value was not rejected'

SISMO_CONFIG="$TEST_DIR/does-not-exist.conf" \
  expect_status 2 "$SISMO" missing.mod output
grep -q 'cannot open file' "$TEST_DIR/stderr" ||
  fail 'SISMO_CONFIG did not take precedence'

printf 'SISMO CLI tests passed\n'
