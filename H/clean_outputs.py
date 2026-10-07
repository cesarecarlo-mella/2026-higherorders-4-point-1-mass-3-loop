#!/usr/bin/env python3
"""Erase the contents of every folder under this directory, keeping the folder tree itself.

Deletes all files recursively inside the folder this script lives in (i.e. "C"),
but leaves every subfolder in place (empty). The script itself is never deleted,
and neither is solution.nb.

Usage:
    python3 clean_outputs.py           # asks for confirmation, then deletes
    python3 clean_outputs.py --dry-run # only lists what would be deleted
    python3 clean_outputs.py --yes     # skip confirmation prompt
"""

import argparse
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="list files without deleting them")
    parser.add_argument("--yes", "-y", action="store_true", help="do not ask for confirmation")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    self_path = Path(__file__).resolve()
    keep = {self_path, root / "solution.nb"}

    files = sorted(p for p in root.rglob("*") if p.is_file() and p not in keep)

    if not files:
        print("Nothing to delete.")
        return

    print(f"Found {len(files)} file(s) under {root}:")
    for f in files:
        print(f"  {f.relative_to(root)}")

    if args.dry_run:
        print("\nDry run: no files were deleted.")
        return

    if not args.yes:
        answer = input(f"\nDelete all {len(files)} file(s) listed above? [y/N] ").strip().lower()
        if answer != "y":
            print("Aborted.")
            return

    deleted = 0
    for f in files:
        try:
            f.unlink()
            deleted += 1
        except OSError as e:
            print(f"  Failed to delete {f}: {e}", file=sys.stderr)

    print(f"Deleted {deleted} file(s). Folder structure preserved.")


if __name__ == "__main__":
    main()
