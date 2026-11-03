#!/usr/bin/env bash
# =============================================================================
# tools/update_expected.sh - refresh the published reference results from a run.
#
#   tools/update_expected.sh <run dir> <expected dir> [--write]
#
# The expected_*/ folders are what a reader compares their own run against, so
# they have to be the output of the run the article reports - no older file
# left behind, and nothing added without being noticed.
#
# A wholesale copy does neither safely. A run writes files the reference set
# does not keep (the multi-gigabyte .rds, the log, the MedDRA mapping audit),
# and the reference set holds files a given run may not have rewritten. So
# this walks what is already tracked, takes each file from the run, and says
# plainly what it found: updated, unchanged, missing from the run, or new in
# the run and not yet tracked.
#
# Nothing is copied without --write. Read the report first.
# =============================================================================
set -euo pipefail

RUN=${1:?usage: update_expected.sh <run dir> <expected dir> [--write]}
EXP=${2:?usage: update_expected.sh <run dir> <expected dir> [--write]}
WRITE=${3:-}

[ -d "$RUN" ] || { echo "no such run directory: $RUN" >&2; exit 1; }
[ -d "$EXP" ] || { echo "no such expected directory: $EXP" >&2; exit 1; }

echo "run      $RUN"
echo "expected $EXP"
[ "$WRITE" = "--write" ] && echo "mode     WRITING" || echo "mode     dry run"
echo

# The tracked list comes from git, which a reader who downloaded the zip does
# not have. Without a repository there is nothing to compare "tracked" against,
# so fall back to every file already in the expected folder - the same set, for
# a clean checkout.
list_tracked() {
  if git -C "$(dirname "$EXP")" rev-parse --git-dir >/dev/null 2>&1; then
    git ls-files "$EXP"
  else
    find "$EXP" -type f | sort
  fi
}

updated=0; same=0; absent=0; added=0

# every file git already tracks under the expected folder
while IFS= read -r tracked; do
  rel=${tracked#"$EXP"/}
  src="$RUN/$rel"
  if [ ! -f "$src" ]; then
    printf '  MISSING FROM THE RUN  %s\n' "$rel"
    absent=$((absent + 1))
    continue
  fi
  if cmp -s "$src" "$tracked"; then
    same=$((same + 1))
    continue
  fi
  printf '  update                %s\n' "$rel"
  updated=$((updated + 1))
  [ "$WRITE" = "--write" ] && cp -p "$src" "$tracked"
done < <(list_tracked)

# anything the run produced that the reference set does not yet carry
while IFS= read -r -d '' f; do
  rel=${f#"$RUN"/}
  case "$rel" in
    *.rds|pipeline_log.txt|llt_mappings_applied.csv|master_*) continue ;;
  esac
  if [ ! -e "$EXP/$rel" ]; then
    printf '  NEW IN THE RUN        %s\n' "$rel"
    added=$((added + 1))
    if [ "$WRITE" = "--write" ]; then
      mkdir -p "$(dirname "$EXP/$rel")"
      cp -p "$f" "$EXP/$rel"
    fi
  fi
done < <(find "$RUN" -type f -print0 | sort -z)

echo
echo "  $updated updated, $same already identical, $absent missing from the run, $added new"
if [ "$WRITE" != "--write" ]; then
  echo
  echo "  Dry run. Check the MISSING lines above - a tracked file the run did not"
  echo "  produce would otherwise be left behind at its old values. Then re-run"
  echo "  with --write."
fi
