from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = ""
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                records.append((name, "".join(chunks)))
            name = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if name:
        records.append((name, "".join(chunks)))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-gap-fraction", type=float, default=0.5)
    args = parser.parse_args()
    records = read_fasta(args.input)
    if not records:
        raise ValueError("empty alignment")
    length = len(records[0][1])
    if any(len(sequence) != length for _, sequence in records):
        raise ValueError("unaligned sequences")
    keep = [
        index
        for index in range(length)
        if sum(sequence[index] in "-." for _, sequence in records) / len(records)
        <= args.max_gap_fraction
    ]
    with args.output.open("w", encoding="ascii") as handle:
        for name, sequence in records:
            trimmed = "".join(sequence[index] for index in keep)
            handle.write(f">{name}\n{trimmed}\n")
    print(f"sequences={len(records)} original_sites={length} retained_sites={len(keep)}")


if __name__ == "__main__":
    main()
