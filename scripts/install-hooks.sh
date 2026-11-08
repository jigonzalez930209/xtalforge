#!/usr/bin/env bash
# Installs the optional local Git hooks (Task 0.2.2).
#
# Usage: scripts/install-hooks.sh
#
# The pre-commit hook runs the rule check, `cargo fmt --check`, the
# one-directory commit check and clippy on the crates touched by the commit.
# The commit-msg hook enforces the feat/ref subject and the task list in the
# body (AGENTS.md section 2.1). Neither hook replaces CI.
set -eu

root=$(git rev-parse --show-toplevel)
cd "$root"
hooks=$(git rev-parse --git-path hooks)

for hook in pre-commit commit-msg; do
  target="$hooks/$hook"
  if [ -e "$target" ] && ! grep -q 'xtalforge-managed hook' "$target"; then
    echo "Refusing to overwrite $target." >&2
    echo "Remove it or merge the xtalforge hook manually." >&2
    exit 1
  fi
done

mkdir -p "$hooks"
for hook in pre-commit commit-msg; do
  cp "scripts/hooks/$hook" "$hooks/$hook"
  chmod +x "$hooks/$hook"
  echo "installed $hooks/$hook"
done

echo "Hooks are optional; CI enforces the same rules on every push."
