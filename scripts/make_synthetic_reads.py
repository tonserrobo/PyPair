#!/usr/bin/env python3
"""Generate deterministic single-end reads and truth positions for BWA tests."""

from __future__ import annotations

import argparse
import random
from pathlib import Path

BASES = "ACGT"


def read_fasta(path: Path) -> tuple[str, str]:
    name = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
            else:
                chunks.append(line)
    if name is None:
        raise ValueError(f"no FASTA header in {path}")
    return name, "".join(chunks).upper()


def generate(
    reference: Path,
    out_prefix: Path,
    n_reads: int,
    read_length: int,
    seed: int,
    error_rate: float,
) -> None:
    name, seq = read_fasta(reference)
    if len(seq) <= read_length:
        raise ValueError("reference is shorter than the read length")

    rng = random.Random(seed)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fastq_path = out_prefix.with_name(out_prefix.name + "_reads.fastq")
    truth_path = out_prefix.with_name(out_prefix.name + "_truth.tsv")

    with fastq_path.open("w") as fastq, truth_path.open("w") as truth:
        truth.write("read_id\tchrom\tpos0\n")
        for index in range(n_reads):
            while True:
                pos = rng.randrange(0, len(seq) - read_length + 1)
                window = seq[pos:pos + read_length]
                if "N" not in window:
                    break

            read = list(window)
            if error_rate:
                for offset in range(len(read)):
                    if rng.random() < error_rate:
                        read[offset] = rng.choice(
                            [base for base in BASES if base != read[offset]]
                        )
            read_text = "".join(read)
            read_id = f"read_{index:06d}"
            fastq.write(f"@{read_id}\n{read_text}\n+\n{'I' * len(read_text)}\n")
            truth.write(f"{read_id}\t{name}\t{pos}\n")

    print(f"{reference.name}: {n_reads} reads -> {fastq_path}, {truth_path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate synthetic reads for BWA-MEM")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, required=True)
    parser.add_argument("--n-reads", type=int, default=10_000)
    parser.add_argument("--read-length", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--error-rate", type=float, default=0.0)
    args = parser.parse_args(argv)

    generate(
        args.reference,
        args.out_prefix,
        args.n_reads,
        args.read_length,
        args.seed,
        args.error_rate,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
