#!/usr/bin/env python3

from pathlib import Path
import re
import sys
import subprocess

SAM = Path(sys.argv[1])
LTR_FILE = Path("reference/canonical_LTR.fa")

# Human chromosomes included in this project
HUMAN_CHR = {
    "12", "2", "6", "9", "11",
    "8", "5", "7", "22", "10"
}

# Minimum exact soft-clip length.
# We start with 20 bp to avoid meaningless 1–4 bp matches.
MIN_CLIP = 20


def read_fasta(path):
    return "".join(
        line.strip()
        for line in path.read_text().splitlines()
        if not line.startswith(">")
    )


def reverse_complement(seq):
    table = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )
    return seq.translate(table)[::-1]


def exact_ltr_match(seq, ltr):
    """
    Return:
        (True, orientation, LTR_start, LTR_end)
    if seq is an exact substring of the canonical LTR.
    Otherwise return:
        (False, "", None, None)
    """

    if len(seq) < MIN_CLIP:
        return False, "", None, None

    # Forward orientation
    p = ltr.find(seq)

    if p >= 0:
        return True, "forward", p + 1, p + len(seq)

    # Reverse-complement orientation
    rc = reverse_complement(seq)
    p = ltr.find(rc)

    if p >= 0:
        return True, "reverse_complement", p + 1, p + len(seq)

    return False, "", None, None


def reference_length(cigar):
    """
    Count bases consumed on the reference.
    M, D, N, =, X consume reference.
    I, S, H, P do not.
    """
    total = 0

    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)

        if op in "MDN=X":
            total += n

    return total


ltr = read_fasta(LTR_FILE)

print("Canonical LTR length:", len(ltr))
print("Minimum exact soft-clip length:", MIN_CLIP)
print()

sam_text = subprocess.check_output(
    ["samtools", "view", str(SAM)],
    text=True
)

results = []

for line in sam_text.splitlines():

    f = line.split("\t")

    if len(f) < 11:
        continue

    read_id = f[0]
    flag = int(f[1])
    chrom = f[2]
    pos = int(f[3])
    mapq = int(f[4])
    cigar = f[5]
    seq = f[9]

    # Only primary alignments
    if flag & 256:
        continue

    if flag & 2048:
        continue

    # Only selected human chromosomes
    if chrom not in HUMAN_CHR:
        continue

    # Need a soft clip
    left_match = re.match(r"^(\d+)S", cigar)
    right_match = re.search(r"(\d+)S$", cigar)

    clips = []

    if left_match:
        n = int(left_match.group(1))
        clips.append(("left", n, seq[:n]))

    if right_match:
        n = int(right_match.group(1))
        clips.append(("right", n, seq[-n:]))

    if not clips:
        continue

    # Reference end of alignment
    ref_len = reference_length(cigar)
    ref_end = pos + ref_len - 1

    for side, clip_len, clip_seq in clips:

        matched, orientation, ltr_start, ltr_end = exact_ltr_match(
            clip_seq,
            ltr
        )

        if not matched:
            continue

        # Breakpoint:
        # left soft clip  -> genomic alignment starts at POS
        # right soft clip -> genomic alignment ends at ref_end
        if side == "left":
            breakpoint = pos
        else:
            breakpoint = ref_end

        results.append({
            "read_id": read_id,
            "chr": chrom,
            "mapq": mapq,
            "cigar": cigar,
            "side": side,
            "softclip_length": clip_len,
            "breakpoint": breakpoint,
            "ltr_orientation": orientation,
            "ltr_start": ltr_start,
            "ltr_end": ltr_end
        })


# Remove duplicate read/side records
unique = {}
for r in results:
    key = (r["read_id"], r["side"])
    unique[key] = r

results = list(unique.values())

print("Exact LTR soft-clip matches:", len(results))
print()

for r in sorted(results, key=lambda x: (x["chr"], x["breakpoint"], x["read_id"])):

    print(
        r["read_id"],
        "\tCHR:", r["chr"],
        "\tMAPQ:", r["mapq"],
        "\tSIDE:", r["side"],
        "\tCLIP:", r["softclip_length"],
        "\tBREAKPOINT:", r["breakpoint"],
        "\tLTR:", r["ltr_start"], "-", r["ltr_end"],
        "\tORIENTATION:", r["ltr_orientation"]
    )
