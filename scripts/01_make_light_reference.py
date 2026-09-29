#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys


# ============================================================
# PROJECT PATHS
# ============================================================

PROJECT = Path(__file__).resolve().parents[1]

HG38 = PROJECT.parent / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
PROVIRUS = PROJECT.parent / "pTG vector LTR to LTR.txt"

OUTPUT = PROJECT / "reference" / "selected_chr_plus_provirus.fa"


# ============================================================
# CHROMOSOMES TO INCLUDE
# IMPORTANT:
# Your hg38 FASTA uses "1", "2", "3"... NOT "chr1", "chr2"
# ============================================================

CHROMOSOMES = [
    "12",
    "2",
    "6",
    "9",
    "11",
    "8",
    "5",
    "7",
    "22",
    "10",
]


# ============================================================
# CHECK INPUT FILES
# ============================================================

if not HG38.exists():
    print("ERROR: hg38 FASTA not found:")
    print(HG38)
    sys.exit(1)

if not PROVIRUS.exists():
    print("ERROR: provirus FASTA not found:")
    print(PROVIRUS)
    sys.exit(1)


# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

OUTPUT.parent.mkdir(parents=True, exist_ok=True)


# ============================================================
# PRINT INFORMATION
# ============================================================

print("Creating lightweight reference...")
print()

print("Human chromosomes included:")

for chrom in CHROMOSOMES:
    print("   ", chrom)

print("   provirus")
print()


# ============================================================
# EXTRACT SELECTED HUMAN CHROMOSOMES
# ============================================================

cmd = [
    "samtools",
    "faidx",
    str(HG38),
    *CHROMOSOMES
]

print("Extracting chromosomes...")

try:

    with open(OUTPUT, "w") as out:

        subprocess.run(
            cmd,
            stdout=out,
            stderr=None,
            check=True
        )

except subprocess.CalledProcessError:

    print()
    print("ERROR: samtools could not extract the chromosomes.")
    print()
    print("Please check that the chromosome names in the FASTA are:")
    print(CHROMOSOMES)
    sys.exit(1)


# ============================================================
# READ PROVIRUS SEQUENCE
# ============================================================

print("Reading provirus sequence...")

sequence = []

with open(PROVIRUS) as f:

    for line in f:

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            continue

        sequence.append(line)


sequence = "".join(sequence).upper()


# ============================================================
# CHECK PROVIRUS LENGTH
# ============================================================

print("Provirus length:", len(sequence))

if len(sequence) != 3683:

    print()
    print("WARNING:")
    print("Expected provirus length = 3683 bp")
    print("Actual provirus length   =", len(sequence))
    print()


# ============================================================
# APPEND PROVIRUS TO REFERENCE
# ============================================================

print("Adding provirus to reference...")

with open(OUTPUT, "a") as out:

    out.write(">provirus\n")

    for i in range(0, len(sequence), 80):

        out.write(sequence[i:i + 80] + "\n")


# ============================================================
# FINISHED
# ============================================================

print()
print("Reference successfully created:")
print(OUTPUT)

print()
print("Reference contains:")
print("  - 10 selected human chromosomes")
print("  - 1 provirus sequence")
print()
print("Done.")
