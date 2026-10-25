#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
export JULIA_DEPOT_PATH="$ROOT/.julia_depot"
exec julia --startup-file=no --project="$ROOT" "$ROOT/scripts/benchmark.jl" "$@"
