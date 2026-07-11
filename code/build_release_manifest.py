from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "metadata" / "FILE_MANIFEST_SHA256.csv"
EXCLUDED_PARTS = {".git", "out", "__pycache__", ".pytest_cache"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def included(path: Path) -> bool:
    relative = path.relative_to(ROOT)
    if path == OUTPUT:
        return False
    return not any(part in EXCLUDED_PARTS for part in relative.parts)


def main() -> None:
    files = sorted(path for path in ROOT.rglob("*") if path.is_file() and included(path))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["path", "bytes", "sha256"])
        for path in files:
            writer.writerow([path.relative_to(ROOT).as_posix(), path.stat().st_size, sha256(path)])
    print(f"Wrote {OUTPUT} with {len(files)} entries")


if __name__ == "__main__":
    main()
