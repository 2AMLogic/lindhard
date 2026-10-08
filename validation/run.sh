#!/usr/bin/env bash
# Regenerate the results tables in docs/validation.md (one command).
#
#   validation/run.sh                 level 1 (full statistics, release build),
#                                     plus the committed level-2 summaries and
#                                     level-3 results
#   validation/run.sh --experiments   also rerun the level-3 comparisons
#   validation/run.sh --oracles       also rerun the level-2 oracles that are
#                                     configured (RUSTBCA_BIN, OPENTRIM_BIN)
#   validation/run.sh --backscatter   also rerun the level-3 electron
#                                     backscatter comparison (about an hour
#                                     on two cores)
#
# Exits non-zero if a level-1 check fails; the docs are then left untouched.
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

run_experiments=0
run_oracles=0
run_backscatter=0
for arg in "$@"; do
  case "$arg" in
    --experiments) run_experiments=1 ;;
    --oracles) run_oracles=1 ;;
    --backscatter) run_backscatter=1 ;;
    -h | --help)
      sed -n '2,15p' "$0"
      exit 0
      ;;
    *)
      echo "unknown argument: $arg" >&2
      exit 2
      ;;
  esac
done

# Every level-3 dataset must carry its provenance (validation/data/README.md).
python3 validation/experiments/run.py --check
python3 validation/experiments/backscatter.py --check
python3 validation/experiments/se_yield.py --datasets
python3 validation/experiments/test_run.py
python3 validation/experiments/test_se_yield.py

level1="$(mktemp "${TMPDIR:-/tmp}/lindhard-level1.XXXXXX")"
trap 'rm -f "$level1"' EXIT

LINDHARD_VALIDATION=full LINDHARD_VALIDATION_OUT="$level1" \
  cargo test --release -p lindhard --test validation

if [ "$run_oracles" = 1 ]; then
  python3 validation/oracles/run.py
fi
if [ "$run_experiments" = 1 ]; then
  python3 validation/experiments/run.py
fi
if [ "$run_backscatter" = 1 ]; then
  python3 validation/experiments/backscatter.py
fi

python3 validation/update_docs.py --level1 "$level1"
