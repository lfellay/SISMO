#!/bin/sh
set -eu

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
MESA2SISMO=${1:-"$ROOT_DIR/mesa2SISMO/mesa2SISMO"}
INTSISMO=${2:-"$ROOT_DIR/intSISMO/intSISMO"}

case "$MESA2SISMO" in
  /*) ;;
  *) MESA2SISMO=$(CDPATH= cd -- "$(dirname -- "$MESA2SISMO")" && pwd)/$(basename -- "$MESA2SISMO") ;;
esac
case "$INTSISMO" in
  /*) ;;
  *) INTSISMO=$(CDPATH= cd -- "$(dirname -- "$INTSISMO")" && pwd)/$(basename -- "$INTSISMO") ;;
esac

TMP_BASE=${TMPDIR:-/tmp}
TEST_DIR=$(mktemp -d "$TMP_BASE/sismo-model-test.XXXXXX")
trap 'rm -rf "$TEST_DIR"' EXIT HUP INT TERM

fail() {
  printf 'test_model_pipeline.sh: %s\n' "$1" >&2
  exit 1
}

expect_failure() {
  label=$1
  shift
  set +e
  "$@" >"$TEST_DIR/failure.out" 2>"$TEST_DIR/failure.err"
  status=$?
  set -e
  if [ "$status" -eq 0 ]; then
    fail "$label unexpectedly succeeded"
  fi
}

expect_cksum() {
  file=$1
  expected=$2
  actual=$(cksum "$file" | awk '{print $1 " " $2}')
  [ "$actual" = "$expected" ] ||
    fail "legacy-format checksum changed for $(basename -- "$file"): $actual"
}

write_minimal_profile() {
  output=$1
  awk 'BEGIN {
    print "1 2 3 4 5 6 7 8"
    print "model_number num_zones star_mass photosphere_r photosphere_L Teff star_age initial_z"
    print "1 20 1.0 1.0 1.0 6000 1.0e8 0.02"
    print ""
    print "1 2 3 4 5 6 7 8"
    print "zone mass_grams radius_cm rho pressure gamma1 brunt_N2 temperature"
    for (j = 1; j <= 20; ++j) {
      mass = (21-j)*1.0e32
      radius = (21-j)*1.0e10
      printf "%d %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n", \
        j, mass, radius, j+0.0, j*1.0e15, 1.6666666666666667, j*1.0e-6, j*1.0e6
    }
  }' >"$output"
}

write_full_profile() {
  output=$1
  variant=$2
  atmosphere=$3
  awk -v variant="$variant" -v atmosphere="$atmosphere" 'BEGIN {
    print "1 2 3 4 5 6 7 8"
    print "model_number num_zones star_mass photosphere_r photosphere_L Teff star_age initial_z"
    print "1 20 1.0 1.0 1.0 6000 1.0e8 0.02"
    print ""
    print "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24"
    print "zone mass_grams radius_cm temperature rho pressure grada gamma1 chiRho chiT cp cv luminosity opacity gradT gradr x z Cprho CpT Qrho QT tau brunt_N2"
    scale = (variant == "changed") ? 2.0 : 1.0
    for (j = 1; j <= 20; ++j) {
      mass = (21-j)*1.0e32
      radius = (21-j)*1.0e10
      tau = 1.0
      if (atmosphere == "yes") {
        if (j == 1) tau = 0.01
        else if (j == 2) tau = 0.10
        else if (j == 3) tau = 0.50
      }
      printf "%d %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n", \
        j, mass, radius, j*1.0e6, j+0.0, j*1.0e15, 0.4, 1.6666666666666667, \
        1.1, 0.9, scale*1.0e8, scale*6.0e7, 1.0, scale*0.3, 0.25, 0.30, \
        0.70, 0.02, scale*0.10, scale*0.20, scale*0.30, scale*0.40, tau, j*1.0e-6
    }
  }' >"$output"
}

write_minimal_profile "$TEST_DIR/minimal.data"
write_full_profile "$TEST_DIR/full.data" normal no
write_full_profile "$TEST_DIR/full-changed.data" changed no
write_full_profile "$TEST_DIR/full-atmosphere.data" normal yes

"$MESA2SISMO" --help >"$TEST_DIR/mesa-help"
grep -q -- '--nonad' "$TEST_DIR/mesa-help" ||
  fail 'mesa2SISMO help does not advertise --nonad'

"$MESA2SISMO" "$TEST_DIR/minimal.data" "$TEST_DIR/minimal-default.madmod"
"$MESA2SISMO" "$TEST_DIR/minimal.data" "$TEST_DIR/minimal-explicit.madmod" --adiabatic
cmp -s "$TEST_DIR/minimal-default.madmod" "$TEST_DIR/minimal-explicit.madmod" ||
  fail 'default and explicit adiabatic models differ'

magic=$(dd if="$TEST_DIR/minimal-default.madmod" bs=8 count=1 2>/dev/null)
[ "$magic" = SISMOAD2 ] || fail 'default model does not use the compact SISMOAD2 format'

expect_failure 'minimal non-adiabatic conversion' \
  "$MESA2SISMO" "$TEST_DIR/minimal.data" "$TEST_DIR/minimal-nonad.madmod" --nonad
grep -qi 'missing.*column' "$TEST_DIR/failure.out" ||
  fail 'minimal non-adiabatic failure did not identify a missing full-model field'

expect_failure 'conflicting mesa2SISMO mode flags' \
  "$MESA2SISMO" "$TEST_DIR/minimal.data" "$TEST_DIR/conflict.madmod" --adiabatic --nonad

"$MESA2SISMO" "$TEST_DIR/full.data" "$TEST_DIR/full-default.madmod"
"$MESA2SISMO" "$TEST_DIR/full-changed.data" "$TEST_DIR/full-changed-default.madmod"
cmp -s "$TEST_DIR/full-default.madmod" "$TEST_DIR/full-changed-default.madmod" ||
  fail 'non-adiabatic source fields changed the default compact model'

"$MESA2SISMO" "$TEST_DIR/full.data" "$TEST_DIR/full-nonad.madmod" --nonad
"$MESA2SISMO" "$TEST_DIR/full-changed.data" "$TEST_DIR/full-changed-nonad.madmod" --nonad
expect_cksum "$TEST_DIR/full-nonad.madmod" '1445830678 8472'
if cmp -s "$TEST_DIR/full-nonad.madmod" "$TEST_DIR/full-changed-nonad.madmod"; then
  fail 'non-adiabatic mode did not preserve changed full-model fields'
fi

legacy_magic=$(dd if="$TEST_DIR/full-nonad.madmod" bs=8 count=1 2>/dev/null)
[ "$legacy_magic" != SISMOAD2 ] || fail 'non-adiabatic mode unexpectedly wrote the compact format'

"$INTSISMO" --help >"$TEST_DIR/int-help"
grep -q -- '--nonad' "$TEST_DIR/int-help" ||
  fail 'intSISMO help does not advertise --nonad'

"$INTSISMO" "$TEST_DIR/minimal-default.madmod" 64 1 radial
"$INTSISMO" "$TEST_DIR/minimal-explicit.madmod" 64 1 radial --adiabatic
cmp -s "$TEST_DIR/minimal-default.osc.mod" "$TEST_DIR/minimal-explicit.osc.mod" ||
  fail 'default and explicit adiabatic remeshing differ'

cp "$TEST_DIR/minimal-default.osc.mod" "$TEST_DIR/minimal-before-failure.osc.mod"
cp "$TEST_DIR/minimal-default.grid.d" "$TEST_DIR/minimal-before-failure.grid.d"
expect_failure 'compact model in non-adiabatic intSISMO mode' \
  "$INTSISMO" "$TEST_DIR/minimal-default.madmod" 64 --nonad
cmp -s "$TEST_DIR/minimal-before-failure.osc.mod" "$TEST_DIR/minimal-default.osc.mod" ||
  fail 'failed non-adiabatic read replaced an existing OSC model'
cmp -s "$TEST_DIR/minimal-before-failure.grid.d" "$TEST_DIR/minimal-default.grid.d" ||
  fail 'failed non-adiabatic read replaced existing grid metadata'

expect_failure 'conflicting intSISMO mode flags' \
  "$INTSISMO" "$TEST_DIR/minimal-default.madmod" 64 --adiabatic --nonad

"$INTSISMO" "$TEST_DIR/full-nonad.madmod" 64 1 radial --nonad
cmp -s "$TEST_DIR/minimal-default.osc.mod" "$TEST_DIR/full-nonad.osc.mod" ||
  fail 'compact and full models produced different adiabatic OSC structures'
grep -Eq 'radial[[:space:]]+nonadiabatic' "$TEST_DIR/full-nonad.grid.d" ||
  fail 'grid metadata did not record non-adiabatic input mode'

"$MESA2SISMO" "$TEST_DIR/full-atmosphere.data" "$TEST_DIR/full-atmosphere.madmod" --nonad
expect_cksum "$TEST_DIR/full-atmosphere.madmod" '1711662951 8964'
"$INTSISMO" --nonad "$TEST_DIR/full-atmosphere.madmod" 64 1 radial
[ -s "$TEST_DIR/full-atmosphere.osc.mod" ] || fail 'atmosphere model did not produce OSC output'

printf '%s\n' 'Adiabatic-first model pipeline tests passed'
