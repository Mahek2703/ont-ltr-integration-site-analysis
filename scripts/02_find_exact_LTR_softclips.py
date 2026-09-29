from pathlib import Path
import re

sam_file = Path("results/test_chunk.sam")
ltr_file = Path("reference/canonical_LTR.fa")

# Read canonical LTR
ltr = "".join(
    line.strip()
    for line in ltr_file.read_text().splitlines()
    if not line.startswith(">")
)

def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

def extract_softclip(seq, cigar):
    m = re.match(r"^(\d+)S", cigar)
    if m:
        n = int(m.group(1))
        return seq[:n], "left"

    m = re.search(r"(\d+)S$", cigar)
    if m:
        n = int(m.group(1))
        return seq[-n:], "right"

    return None, None

matches = []

with sam_file.open() as f:
    for line in f:
        if line.startswith("@"):
            continue

        fields = line.rstrip().split("\t")

        read_id = fields[0]
        flag = int(fields[1])
        chrom = fields[2]
        pos = int(fields[3])
        mapq = int(fields[4])
        cigar = fields[5]
        seq = fields[9]

        # Skip unmapped reads
        if flag & 4:
            continue

        softclip, side = extract_softclip(seq, cigar)

        if softclip is None:
            continue

        if softclip in ltr:
            orientation = "forward"
        elif revcomp(softclip) in ltr:
            orientation = "reverse_complement"
        else:
            continue

        matches.append(
            (read_id, chrom, pos, mapq, cigar, side,
             len(softclip), orientation, softclip)
        )

print("Canonical LTR length:", len(ltr))
print("Exact LTR soft-clipped matches:", len(matches))
print()

for row in matches:
    print(
        "read_id:", row[0],
        "| chr:", row[1],
        "| position:", row[2],
        "| MAPQ:", row[3],
        "| CIGAR:", row[4],
        "| side:", row[5],
        "| softclip_length:", row[6],
        "| orientation:", row[7]
    )
