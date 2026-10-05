#!/usr/bin/env bash
# Regenerate the results tables in docs/validation.md (one command).
#
#   validation/run.sh                 level 1 (full statistics, release build),
#                                     plus the committed level-2 summaries and
#                                     level-3 results
#   validation/run.sh --experiments   also rerun the level-3 comparisons
#   validation/run.sh --oracles       also rerun the level-2 oracles that are
#                                     configured (RUSTBCA_BIN, OPENTRIM_BIN)
#
# Exits non-zero if a level-1 check fails; the docs are then left untouched.
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

run_experiments=0
run_oracles=0
for arg in "$@"; do
  case "$arg" in
    --experiments) run_experiments=1 ;;
    --oracles) run_oracles=1 ;;
    -h | --help)
      sed -n '2,12p' "$0"
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

python3 validation/update_docs.py --level1 "$level1"
