#!/bin/bash
# Build a variant of ascot5_main from a clean tree and store it as
# variants/ascot5_main.<tag>, with the compiler log next to it.
#
# Usage: tools/build_variant.sh <tag> [extra make variables...]
#   e.g. tools/build_variant.sh vlen64 GPU_VLEN=64
#        tools/build_variant.sh reg128 MAXREGCOUNT=128 MINFO=1
#        tools/build_variant.sh cpu CPU=1      (CPU build of the same source)
set -e

TAG=${1:?usage: $0 <tag> [make vars]}
shift

CPU=0
ARGS=()
for a in "$@"; do
    if [ "$a" == "CPU=1" ]; then CPU=1; else ARGS+=("$a"); fi
done

module purge
if [ $CPU == 1 ]; then
    # Same environment as install_develop_mp_hpc.sh
    source /opt/cray/pe/cpe/24.07/restore_lmod_system_defaults.sh
    module load PrgEnv-gnu
    module load cray-hdf5-parallel
    module load craype-x86-milan
    # endcond.c passes the wrong RNG pointer type to flr_losses_eval (dead
    # code, FLR losses are never enabled); newer GCC makes that an error.
    MAKEARGS=(MPI=1 CC=cc
              "FLAGS=-O3 -march=native -fopenmp -fPIC -std=c11 -D_XOPEN_SOURCE=700 -Wno-error=incompatible-pointer-types")
else
    source /opt/cray/pe/cpe/25.09/restore_lmod_system_defaults.sh
    module load PrgEnv-nvidia
    module load cudatoolkit
    module load cray-hdf5-parallel
    module load craype-x86-milan
    module load craype-accel-nvidia80
    MAKEARGS=(GPU=1 ACC=1 MPI=1 CC=cc)
fi

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

# Top-level clean removes build/, so keep finished variants aside.
VARDIR="$ROOT/variants"
mkdir -p "$VARDIR"
make clean > /dev/null
make ascot5_main "${MAKEARGS[@]}" "${ARGS[@]}" \
    > "$VARDIR/build.$TAG.log" 2>&1 \
    || { echo "Build failed, see $VARDIR/build.$TAG.log"; exit 1; }
cp build/ascot5_main "$VARDIR/ascot5_main.$TAG"
echo "Built $VARDIR/ascot5_main.$TAG ($(git rev-parse --short HEAD) $*)"
