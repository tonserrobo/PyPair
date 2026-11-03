#!/usr/bin/env python3
"""Parse a BWA-MEM SAM file and a GNU time log into one benchmark row."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

TIME_RE = re.compile(r"^\s*(\d+):(\d+(?:\.\d+)?)$")
MAX_RSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
ELAPSED_RE = re.compile(r"Elapsed \(wall clock\) time.*?:\s*([0-9:.]+)")


def parse_elapsed(text: str) -> float:
    match = TIME_RE.match(text.strip())
    if not match:
        return 0.0
    minutes, seconds = match.groups()
    return int(minutes) * 60 + float(seconds)


def parse_time_log(path: Path) -> tuple[float, float]:
    wall = 0.0
    max_rss_mb = 0.0
    if not path.exists():
        return wall, max_rss_mb
    for line in path.read_text(errors="replace").splitlines():
        match = ELAPSED_RE.search(line)
        if match:
            wall = parse_elapsed(match.group(1))
        match = MAX_RSS_RE.search(line)
        if match:
            max_rss_mb = int(match.group(1)) / 1024
    return wall, max_rss_mb


def read_truth(path: Path) -> dict[str, tuple[str, int]]:
    with path.open() as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {row["read_id"]: (row["chrom"], int(row["pos0"])) for row in rows}


def parse_sam(sam: Path, truth: dict[str, tuple[str, int]]) -> dict:
    mapped = exact = within10 = 0
    seen: set[str] = set()
    with sam.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            read_id = fields[0]
            flag = int(fields[1])
            if flag & 0x4 or flag & 0x900:
                continue
            if read_id in seen:
                continue
            seen.add(read_id)
            if read_id not in truth:
                continue
            mapped += 1
            chrom, true_pos = truth[read_id]
            if fields[2] == chrom:
                pos = int(fields[3]) - 1
                if pos == true_pos:
                    exact += 1
                if abs(pos - true_pos) <= 10:
                    within10 += 1

    n_reads = len(truth)
    return {
        "n_reads": n_reads,
        "mapped": mapped,
        "exact": exact,
        "within_10bp": within10,
        "mapped_rate": mapped / n_reads if n_reads else 0.0,
        "exact_rate": exact / n_reads if n_reads else 0.0,
        "within_10bp_rate": within10 / n_reads if n_reads else 0.0,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Parse BWA-MEM SAM + time log")
    parser.add_argument("--sam", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--time-log", type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    truth = read_truth(args.truth)
    metrics = parse_sam(args.sam, truth)
    wall, max_rss_mb = parse_time_log(args.time_log) if args.time_log else (0.0, 0.0)

    row = {
        "label": args.label,
        **metrics,
        "wall_seconds": wall,
        "max_rss_mb": max_rss_mb,
        "ms_per_read": (wall * 1000 / metrics["n_reads"]) if metrics["n_reads"] else 0.0,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_header = not args.out.exists()
    with args.out.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        if write_header:
            writer.writeheader()
        writer.writerow(row)

    print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
