from pathlib import Path
import csv
import re

INPUT = Path("results/all_LTR_matches.txt")
OUTPUT = Path("results/final_integration_sites.csv")

pattern = re.compile(
    r'^([0-9a-f-]+)\s+'
    r'CHR:\s*(\S+)\s+'
    r'MAPQ:\s*(\d+)\s+'
    r'SIDE:\s*(\S+)\s+'
    r'CLIP:\s*(\d+)\s+'
    r'BREAKPOINT:\s*(\d+)\s+'
    r'LTR:\s*(\d+)\s*-\s*(\d+)\s+'
    r'ORIENTATION:\s*(\S+)'
)

records = []

for line in INPUT.read_text().splitlines():
    match = pattern.match(line)
    if not match:
        continue

    (
        read_id,
        chrom,
        mapq,
        side,
        clip,
        breakpoint,
        ltr_start,
        ltr_end,
        orientation,
    ) = match.groups()

    records.append({
        "read_id": read_id,
        "chr": f"chr{chrom}",
        "mapq": int(mapq),
        "position": int(breakpoint),
        "cluster": f"chr{chrom}:{breakpoint}",
        "softclipped_matched": "exact",
        "side": side,
        "clip_length": int(clip),
        "ltr_start": int(ltr_start),
        "ltr_end": int(ltr_end),
        "orientation": orientation,
    })

columns = [
    "read_id",
    "chr",
    "mapq",
    "position",
    "cluster",
    "softclipped_matched",
    "side",
    "clip_length",
    "ltr_start",
    "ltr_end",
    "orientation",
]

with OUTPUT.open("w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns)
    writer.writeheader()
    writer.writerows(records)

print(f"Records written: {len(records)}")
print(f"Output: {OUTPUT}")
