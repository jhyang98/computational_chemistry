#!/usr/bin/env bash
# Synchronize only this project directory to computational_chemistry/soft_east.
set -euo pipefail

SYNC_CLONE=${SYNC_CLONE:-/tmp/computational_chemistry_soft_east_upload}
SYNC_REMOTE=${SYNC_REMOTE:-https://github.com/jhyang98/computational_chemistry.git}
SYNC_BRANCH=${SYNC_BRANCH:-master}
COMMIT_MESSAGE=${1:-Update soft_east results}

project_dir=$(pwd)
if [ ! -f "$project_dir/AGENTS.md" ] || [ ! -f "$project_dir/STATUS.md" ]; then
    echo "error: run this script from the soft_east project root" >&2
    exit 1
fi

if [ ! -d "$SYNC_CLONE/.git" ]; then
    git clone "$SYNC_REMOTE" "$SYNC_CLONE"
fi

if [ "$(git -C "$SYNC_CLONE" remote get-url origin)" != "$SYNC_REMOTE" ]; then
    echo "error: unexpected origin in $SYNC_CLONE" >&2
    exit 1
fi

unrelated=$(git -C "$SYNC_CLONE" status --porcelain | awk '$2 !~ /^soft_east\// { print }')
if [ -n "$unrelated" ]; then
    echo "error: the sync clone has changes outside soft_east/" >&2
    printf '%s\n' "$unrelated" >&2
    exit 1
fi

git -C "$SYNC_CLONE" pull --rebase origin "$SYNC_BRANCH"
mkdir -p "$SYNC_CLONE/soft_east"
rsync -a --delete --exclude='.git/' "$project_dir/" "$SYNC_CLONE/soft_east/"
git -C "$SYNC_CLONE" add -- soft_east

if git -C "$SYNC_CLONE" diff --cached --quiet; then
    echo "GitHub sync: no changes"
    exit 0
fi

git -C "$SYNC_CLONE" commit -m "$COMMIT_MESSAGE"
git -C "$SYNC_CLONE" push origin "$SYNC_BRANCH"
echo "GitHub sync complete: $COMMIT_MESSAGE"
