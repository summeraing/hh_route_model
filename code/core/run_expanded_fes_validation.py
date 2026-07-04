from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import requests


UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

GROUPS = {
    "eukaryota": "taxonomy_id:2759",
    "alphaproteobacteria": "taxonomy_id:28211",
    "archaea": "taxonomy_id:2157",
    "nonalpha_bacteria": "taxonomy_id:2 AND NOT taxonomy_id:28211",
}

AA_ALLOWED = set("ACDEFGHIKLMNPQRSTVWYBXZJ")

CURATED_TERMS = {
    "IBA57": '(gene_exact:IBA57 OR gene:IBA57 OR protein_name:"IBA57" OR protein_name:"CAF17")',
    "ISCA1": '(gene:ISCA1 OR protein_name:"ISCA1" OR protein_name:"iron-sulfur cluster assembly protein IscA")',
    "ISCA2": '(gene:ISCA2 OR protein_name:"ISCA2" OR protein_name:"iron-sulfur cluster assembly protein IscA")',
    "ISCU": '(gene:ISCU OR gene:iscU OR gene:sufU OR gene:nifU OR protein_name:"iron-sulfur cluster assembly scaffold" OR protein_name:"iron-sulfur cluster assembly protein IscU")',
    "NFS1": '(gene:NFS1 OR gene:iscS OR gene:nifS OR gene:sufS OR protein_name:"cysteine desulfurase")',
    "FDX1_FDX2": '(gene:FDX1 OR gene:FDX2 OR protein_name:"ferredoxin")',
    "GLRX5": '(gene:GLRX5 OR gene:grx5 OR protein_name:"glutaredoxin")',
    "ABCB7_ATM1": '(gene:ABCB7 OR gene:ATM1 OR protein_name:"ABC transporter" OR protein_name:"ATP-binding cassette")',
    "NFU1": '(gene:NFU1 OR gene:nfuA OR protein_name:"NFU1" OR protein_name:"iron-sulfur cluster scaffold NFU")',
    "BOLA3": '(gene:BOLA3 OR gene:bolA OR protein_name:"BolA")',
    "NUBPL": '(gene:NUBPL OR gene:IND1 OR protein_name:"NADH dehydrogenase assembly factor")',
    "LYRM4_ISD11": '(gene:LYRM4 OR gene:ISD11 OR protein_name:"ISD11" OR protein_name:"LYR motif-containing protein 4")',
    "CIAO1": '(gene:CIAO1 OR protein_name:"CIAO1" OR protein_name:"cytosolic iron-sulfur assembly component 1")',
    "CIAO2B": '(gene:CIAO2B OR gene:FAM96B OR protein_name:"CIAO2B" OR protein_name:"FAM96B")',
    "CIAO3": '(gene:CIAO3 OR gene:NARFL OR protein_name:"CIAO3" OR protein_name:"NARFL")',
    "MMS19": '(gene:MMS19 OR protein_name:"MMS19")',
    "NUBP1_NBP35": '(gene:NUBP1 OR gene:NBP35 OR protein_name:"NBP35" OR protein_name:"nucleotide-binding protein 1")',
    "NUBP2_CFD1": '(gene:NUBP2 OR gene:CFD1 OR gene:NAR1 OR protein_name:"CFD1" OR protein_name:"nucleotide-binding protein 2")',
}


@dataclass
class Record:
    header: str
    seq: str


def tool(name: str) -> str | None:
    return shutil.which(name)


def safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.|-]+", "_", text)[:180]


def clean_protein_sequence(seq: str) -> str:
    seq = seq.upper().replace("*", "")
    seq = seq.replace("U", "X").replace("O", "X")
    return "".join(ch if ch in AA_ALLOWED else "X" for ch in seq)


def fasta_records(text: str) -> list[Record]:
    records: list[Record] = []
    header = None
    chunks: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header and chunks:
                records.append(Record(header, "".join(chunks)))
            header = line[1:]
            chunks = []
        else:
            chunks.append(line)
    if header and chunks:
        records.append(Record(header, "".join(chunks)))
    return records


def marker_query(query_terms: str) -> str:
    terms = [t.strip() for t in query_terms.split("|") if t.strip()]
    expanded = []
    for term in terms:
        if " " in term:
            expanded.append(f'"{term}"')
        else:
            expanded.extend([f"gene:{term}", f'protein_name:"{term}"', term])
    return "(" + " OR ".join(expanded) + ")"


def count_fasta(path: Path) -> int:
    if not path.exists():
        return 0
    return path.read_text(encoding="utf-8").count(">")


def fetch_group(
    marker_id: str,
    query_terms: str,
    group: str,
    max_records: int,
    out_fasta: Path,
    timeout: int,
    retries: int,
) -> int:
    base_query = CURATED_TERMS.get(marker_id, marker_query(query_terms))
    query = f"{base_query} AND ({GROUPS[group]})"
    params = {
        "format": "fasta",
        "query": query,
        "size": max(220, max_records * 6),
    }
    last_exc: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(UNIPROT_SEARCH, params=params, timeout=timeout)
            r.raise_for_status()
            break
        except Exception as exc:  # noqa: BLE001
            last_exc = exc
            if attempt >= retries:
                raise
            time.sleep(min(20, 2 * attempt))
    else:
        raise RuntimeError(f"fetch failed without exception: {last_exc}")
    records = fasta_records(r.text)

    selected: list[Record] = []
    seen_seq: set[str] = set()
    for rec in records:
        seq = clean_protein_sequence(rec.seq)
        if len(seq) < 50 or len(seq) > 1800:
            continue
        if seq in seen_seq:
            continue
        seen_seq.add(seq)
        selected.append(Record(rec.header, seq))
        if len(selected) >= max_records:
            break

    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w", encoding="utf-8") as fh:
        for i, rec in enumerate(selected, 1):
            fh.write(f">{marker_id}|{group}|{i}|{safe(rec.header)}\n")
            for j in range(0, len(rec.seq), 80):
                fh.write(rec.seq[j : j + 80] + "\n")
    return len(selected)


def combine_fastas(paths: list[Path], out: Path, per_file_cap: int = 0) -> int:
    n = 0
    with out.open("w", encoding="utf-8") as dest:
        for path in paths:
            if not path.exists():
                continue
            text = path.read_text(encoding="utf-8")
            records = fasta_records(text)
            if per_file_cap and per_file_cap > 0:
                records = records[:per_file_cap]
            for rec in records:
                seq = clean_protein_sequence(rec.seq)
                if not seq:
                    continue
                n += 1
                dest.write(f">{rec.header}\n")
                for j in range(0, len(seq), 80):
                    dest.write(seq[j : j + 80] + "\n")
    return n


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("RUN", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def run_tree(marker_dir: Path, marker_id: str, threads: str, bootstrap: int, fixed: bool) -> Path | None:
    mafft = tool("mafft")
    iqtree = tool("iqtree2") or tool("iqtree")
    if not mafft or not iqtree:
        print(f"SKIP {marker_id}: mafft or iqtree not found")
        return None

    combined = marker_dir / f"{marker_id}.combined.faa"
    aligned = marker_dir / f"{marker_id}.aligned.faa"
    prefix = marker_dir / ("fixed_LGFR6" if fixed else "model_finder")
    tree = prefix.with_suffix(".treefile")
    if tree.exists() and tree.stat().st_size > 0:
        print(f"REUSE_TREE {marker_id} {tree}")
        return tree
    if not aligned.exists() or aligned.stat().st_size == 0:
        with aligned.open("w", encoding="utf-8") as out:
            subprocess.run(
                [mafft, "--auto", "--thread", str(threads), str(combined)],
                stdout=out,
                check=True,
            )

    model = "LG+F+R6" if fixed else "MFP"
    cmd = [
        iqtree,
        "-s",
        str(aligned),
        "-m",
        model,
        "-alrt",
        str(bootstrap),
        "-B",
        str(bootstrap),
        "-T",
        str(threads),
        "--prefix",
        str(prefix),
        "--redo",
    ]
    run(cmd)
    return tree if tree.exists() else None


def nearest_summary(tree_path: Path, marker_id: str, run_label: str) -> dict[str, object]:
    try:
        from Bio import Phylo
    except Exception as exc:
        return {
            "marker_id": marker_id,
            "run_label": run_label,
            "status": f"biopython_missing:{exc}",
        }

    tree = Phylo.read(str(tree_path), "newick")
    terms = tree.get_terminals()
    euk = [t for t in terms if "|eukaryota|" in t.name]
    non = [t for t in terms if "|eukaryota|" not in t.name]
    counts = {g: 0 for g in GROUPS if g != "eukaryota"}
    if not euk or not non:
        return {"marker_id": marker_id, "run_label": run_label, "status": "not_estimable"}

    for terminal in euk:
        best = min(non, key=lambda other: tree.distance(terminal, other))
        group = best.name.split("|")[1]
        counts[group] = counts.get(group, 0) + 1

    total = len(euk)
    row: dict[str, object] = {
        "marker_id": marker_id,
        "run_label": run_label,
        "status": "ok",
        "n_eukaryotic_terminals": total,
    }
    for group, value in counts.items():
        row[f"nearest_{group}_n"] = value
        row[f"nearest_{group}_fraction"] = value / total if total else 0
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", default="AUTO")
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--download-only", action="store_true")
    ap.add_argument("--reuse-fastas", action="store_true", help="Reuse existing FASTA files in the output directory when present.")
    ap.add_argument("--skip-fetch", action="store_true", help="Do not contact UniProt; only use existing FASTA files.")
    ap.add_argument("--fetch-timeout", type=int, default=120)
    ap.add_argument("--fetch-retries", type=int, default=3)
    ap.add_argument(
        "--tree-mode",
        choices=["both", "model_finder", "fixed"],
        default="both",
        help="Which tree analyses to run after FASTA preparation.",
    )
    ap.add_argument("--per-group-cap", type=int, default=0, help="Optional cap on records from each group FASTA for tree screening.")
    ap.add_argument("--limit", type=int, default=0, help="Optional smoke-test limit on marker rows.")
    args = ap.parse_args()

    manifest = pd.read_csv(args.manifest)
    if args.limit and args.limit > 0:
        manifest = manifest.head(args.limit).copy()
    out = Path(args.out)
    raw = out / "raw_fastas"
    marker_root = out / "markers"
    out.mkdir(parents=True, exist_ok=True)

    retrieval_rows = []
    summary_rows = []

    for _, marker in manifest.iterrows():
        marker_id = marker["marker_id"]
        max_per_group = int(marker.get("max_per_group", 60))
        marker_dir = marker_root / marker_id
        marker_dir.mkdir(parents=True, exist_ok=True)

        fastas = []
        group_counts: dict[str, int] = {}
        for group in GROUPS:
            fasta = raw / marker_id / f"{group}.faa"
            status = "ok"
            existing_n = count_fasta(fasta)
            if args.skip_fetch:
                n = existing_n
                status = "reused" if existing_n > 0 else ("existing_empty" if fasta.exists() else "missing_existing_fasta")
            elif args.reuse_fastas and existing_n > 0:
                n = existing_n
                status = "reused"
            else:
                try:
                    n = fetch_group(
                        marker_id,
                        marker["query_terms"],
                        group,
                        max_per_group,
                        fasta,
                        timeout=args.fetch_timeout,
                        retries=args.fetch_retries,
                    )
                except Exception as exc:  # noqa: BLE001
                    n = 0
                    status = f"fetch_failed:{type(exc).__name__}"
                    fasta.parent.mkdir(parents=True, exist_ok=True)
                    fasta.write_text("", encoding="utf-8")
            retrieval_rows.append(
                {
                    "marker_id": marker_id,
                    "group": group,
                    "selected_records": n,
                    "max_per_group": max_per_group,
                    "status": status,
                }
            )
            group_counts[group] = n
            fastas.append(fasta)
            time.sleep(0.4)

        n_total = combine_fastas(fastas, marker_dir / f"{marker_id}.combined.faa", per_file_cap=args.per_group_cap)
        if args.download_only or n_total < 12:
            summary_rows.append(
                {
                    "marker_id": marker_id,
                    "run_label": "not_run",
                    "status": "download_only_or_too_few_records",
                    "n_total_records": n_total,
                }
            )
            continue
        n_euk = min(group_counts.get("eukaryota", 0), args.per_group_cap) if args.per_group_cap else group_counts.get("eukaryota", 0)
        n_non_euk = 0
        for group in GROUPS:
            if group == "eukaryota":
                continue
            group_n = group_counts.get(group, 0)
            n_non_euk += min(group_n, args.per_group_cap) if args.per_group_cap else group_n
        if n_euk == 0 or n_non_euk == 0:
            summary_rows.append(
                {
                    "marker_id": marker_id,
                    "run_label": "not_run",
                    "status": "not_estimable_no_eukaryote_or_no_reference",
                    "n_total_records": n_total,
                    "compartment_class": marker["compartment_class"],
                    "route_expectation": marker["route_expectation"],
                }
            )
            continue

        tree_modes: list[tuple[bool, str]] = []
        if args.tree_mode in {"both", "model_finder"}:
            tree_modes.append((False, "model_finder"))
        if args.tree_mode in {"both", "fixed"}:
            tree_modes.append((True, "fixed_LGFR6"))

        for fixed, label in tree_modes:
            tree = run_tree(marker_dir, marker_id, str(args.threads), args.bootstrap, fixed=fixed)
            if tree:
                row = nearest_summary(tree, marker_id, label)
                row["n_total_records"] = n_total
                row["compartment_class"] = marker["compartment_class"]
                row["route_expectation"] = marker["route_expectation"]
                summary_rows.append(row)

    pd.DataFrame(retrieval_rows).to_csv(out / "expanded_fes_retrieval_summary.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(out / "expanded_fes_nearest_neighbor_summary.csv", index=False)
    manifest.to_csv(out / "expanded_fes_frozen_marker_manifest.csv", index=False)

    with (out / "run_environment.txt").open("w", encoding="utf-8") as fh:
        fh.write(f"mafft={tool('mafft')}\n")
        fh.write(f"iqtree2={tool('iqtree2')}\n")
        fh.write(f"iqtree={tool('iqtree')}\n")
        fh.write(f"bootstrap={args.bootstrap}\n")
        fh.write(f"threads={args.threads}\n")


if __name__ == "__main__":
    main()
