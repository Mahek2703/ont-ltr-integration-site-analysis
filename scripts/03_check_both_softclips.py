from pathlib import Path
import re

sam_file = Path("results/test_chunk.sam")
ltr_file = Path("reference/canonical_LTR.fa")

ltr = "".join(
    line.strip()
    for line in ltr_file.read_text().splitlines()
    if not line.startswith(">")
)

def revcomp(seq):
    return seq.translate(
        str.maketrans("ACGTNacgtn", "TGCANtgcan")
    )[::-1]

def check_clip(seq):
    if seq in ltr:
        return "forward"
    if revcomp(seq) in ltr:
        return "reverse_complement"
    return None

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

        if flag & 4:
            continue

        # Left soft clip
        m_left = re.match(r"^(\d+)S", cigar)
        if m_left:
            n = int(m_left.group(1))
            clip = seq[:n]
            orientation = check_clip(clip)

            if orientation:
                matches.append(
                    (read_id, chrom, pos, mapq, cigar,
                     "left", n, orientation)
                )

        # Right soft clip
        m_right = re.search(r"(\d+)S$", cigar)
        if m_right:
            n = int(m_right.group(1))
            clip = seq[-n:]
            orientation = check_clip(clip)

            if orientation:
                matches.append(
                    (read_id, chrom, pos, mapq, cigar,
                     "right", n, orientation)
                )

print("Canonical LTR length:", len(ltr))
print("Exact soft-clipped matches:", len(matches))
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
