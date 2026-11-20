#!/usr/bin/env bash
# Path-limited commit in prFBA with index.lock retry (another agent commits other paths concurrently).
#   usage: pcommit.sh <msgfile> <path> [path ...]      (new paths are `git add`ed by name first)
set -u
R=/home/freiburger/Documents/prFBA
msg=$1; shift
cd "$R"
for i in $(seq 1 20); do
  if [ -e .git/index.lock ]; then sleep 5; continue; fi
  new=()
  for p in "$@"; do git ls-files --error-unmatch -- "$p" >/dev/null 2>&1 || new+=("$p"); done
  if [ ${#new[@]} -gt 0 ]; then git add -- "${new[@]}" 2>/tmp/pcommit.err || { grep -q index.lock /tmp/pcommit.err && { sleep 5; continue; }; cat /tmp/pcommit.err; exit 1; }; fi
  if git commit -q -F "$msg" -- "$@" 2>/tmp/pcommit.err; then git log --oneline -1; exit 0; fi
  grep -q index.lock /tmp/pcommit.err && { sleep 5; continue; }
  cat /tmp/pcommit.err; exit 1
done
echo "pcommit: gave up after retries" >&2; exit 1
