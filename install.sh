#!/bin/sh
set -eu

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
BIN_DIR="$ROOT_DIR/bin"
MAKE_CMD=${MAKE:-make}
FC_CMD=${FC:-gfortran}
INSTALL_CMD=${INSTALL:-install}

usage() {
  cat <<EOF
Usage: ./install.sh [--bin-dir DIR]

Build SISMO 2.0, mesa2SISMO, and intSISMO, then install them and the
SISMO configuration in DIR.
The default destination is:
  $ROOT_DIR/bin

Environment variables:
  FC       GNU Fortran compiler (default: gfortran)
  MAKE     make command (default: make)
  INSTALL  install command (default: install)
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --bin-dir)
      if [ "$#" -lt 2 ] || [ -z "$2" ]; then
        printf '%s\n' 'install.sh: --bin-dir requires a directory' >&2
        exit 2
      fi
      BIN_DIR=$2
      shift 2
      ;;
    --bin-dir=*)
      BIN_DIR=${1#*=}
      if [ -z "$BIN_DIR" ]; then
        printf '%s\n' 'install.sh: --bin-dir requires a directory' >&2
        exit 2
      fi
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'install.sh: unknown option: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

for command_name in "$MAKE_CMD" "$FC_CMD" "$INSTALL_CMD" mktemp; do
  if ! command -v "$command_name" >/dev/null 2>&1; then
    printf 'install.sh: required command not found: %s\n' "$command_name" >&2
    exit 1
  fi
done

build_program() {
  label=$1
  source_dir=$2
  shift 2

  printf '\nBuilding %s\n' "$label"
  "$MAKE_CMD" -j1 -C "$source_dir" clean
  "$MAKE_CMD" -j1 -C "$source_dir" FC="$FC_CMD" "$@"
}

build_program 'SISMO 2.0' "$ROOT_DIR/SISMO2.0"
build_program 'mesa2SISMO' "$ROOT_DIR/mesa2SISMO"
build_program 'intSISMO' "$ROOT_DIR/intSISMO" debug=no

TMP_BASE=${TMPDIR:-/tmp}
STAGE_DIR=$(mktemp -d "$TMP_BASE/sismo-install.XXXXXX")
COMMIT_DIR=
BACKUP_DIR=
LOCK_DIR=
LOCK_HELD=0
COMMIT_ACTIVE=0
INSTALLED_NAMES=
BACKED_UP_NAMES=

cleanup_and_exit() {
  status=$1
  trap - EXIT HUP INT TERM
  set +e

  if [ "$COMMIT_ACTIVE" -eq 1 ]; then
    for name in $INSTALLED_NAMES; do
      rm -f "$BIN_DIR/$name"
    done
    for name in $BACKED_UP_NAMES; do
      mv "$BACKUP_DIR/$name" "$BIN_DIR/$name"
    done
  fi

  if [ -n "$COMMIT_DIR" ]; then
    rm -rf "$COMMIT_DIR"
  fi
  rm -rf "$STAGE_DIR"
  if [ "$LOCK_HELD" -eq 1 ]; then
    rmdir "$LOCK_DIR"
  fi
  exit "$status"
}

trap 'cleanup_and_exit $?' EXIT
trap 'exit 1' HUP INT TERM

"$INSTALL_CMD" -m 755 "$ROOT_DIR/SISMO2.0/sismo" "$STAGE_DIR/sismo"
"$INSTALL_CMD" -m 755 "$ROOT_DIR/mesa2SISMO/mesa2SISMO" "$STAGE_DIR/mesa2SISMO"
"$INSTALL_CMD" -m 755 "$ROOT_DIR/intSISMO/intSISMO" "$STAGE_DIR/intSISMO"
"$INSTALL_CMD" -m 644 "$ROOT_DIR/SISMO2.0/sismo.conf" "$STAGE_DIR/sismo.conf"

"$STAGE_DIR/sismo" --help >/dev/null
"$STAGE_DIR/mesa2SISMO" --help >/dev/null
"$STAGE_DIR/intSISMO" --help >/dev/null

mkdir -p "$BIN_DIR"
LOCK_DIR="$BIN_DIR/.sismo-install.lock"
if ! mkdir "$LOCK_DIR" 2>/dev/null; then
  printf 'install.sh: another installation is active in %s\n' "$BIN_DIR" >&2
  exit 1
fi
LOCK_HELD=1

COMMIT_DIR=$(mktemp -d "$BIN_DIR/.sismo-commit.XXXXXX")
BACKUP_DIR="$COMMIT_DIR/backup"
mkdir "$BACKUP_DIR"

"$INSTALL_CMD" -m 755 "$STAGE_DIR/sismo" "$COMMIT_DIR/sismo"
"$INSTALL_CMD" -m 755 "$STAGE_DIR/mesa2SISMO" "$COMMIT_DIR/mesa2SISMO"
"$INSTALL_CMD" -m 755 "$STAGE_DIR/intSISMO" "$COMMIT_DIR/intSISMO"
"$INSTALL_CMD" -m 644 "$STAGE_DIR/sismo.conf" "$COMMIT_DIR/sismo.conf.default"
COMMIT_NAMES='sismo mesa2SISMO intSISMO sismo.conf.default'

if [ -e "$BIN_DIR/sismo.conf" ] || [ -L "$BIN_DIR/sismo.conf" ]; then
  CONFIG_STATUS=preserved
else
  "$INSTALL_CMD" -m 644 "$STAGE_DIR/sismo.conf" "$COMMIT_DIR/sismo.conf"
  COMMIT_NAMES="$COMMIT_NAMES sismo.conf"
  CONFIG_STATUS=created
fi

# Every new artifact is complete and on the destination filesystem before
# replacement begins. If any rename fails or the process is interrupted, the
# EXIT trap removes new files and restores the previous complete tool set.
COMMIT_ACTIVE=1
for name in $COMMIT_NAMES; do
  if [ -d "$BIN_DIR/$name" ] && [ ! -L "$BIN_DIR/$name" ]; then
    printf 'install.sh: destination is a directory: %s\n' "$BIN_DIR/$name" >&2
    exit 1
  fi
  if [ -e "$BIN_DIR/$name" ] || [ -L "$BIN_DIR/$name" ]; then
    BACKED_UP_NAMES="$BACKED_UP_NAMES $name"
    mv "$BIN_DIR/$name" "$BACKUP_DIR/$name"
  fi
  INSTALLED_NAMES="$INSTALLED_NAMES $name"
  mv "$COMMIT_DIR/$name" "$BIN_DIR/$name"
done
COMMIT_ACTIVE=0

printf '\nInstalled SISMO 2.0 tools in %s\n' "$BIN_DIR"
printf '  %s\n' "$BIN_DIR/sismo" "$BIN_DIR/mesa2SISMO" "$BIN_DIR/intSISMO"
printf 'Configuration (%s):\n  %s\n' "$CONFIG_STATUS" "$BIN_DIR/sismo.conf"
printf 'Default configuration template (updated):\n  %s\n' "$BIN_DIR/sismo.conf.default"
