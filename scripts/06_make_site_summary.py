import csv
from collections import defaultdict
from statistics import mean

INPUT = "results/final_integration_sites.csv"
OUTPUT = "results/integration_site_summary.csv"

groups = defaultdict(list)

with open(INPUT, newline="") as f:
    for row in csv.DictReader(f):
        key = (row["chr"], int(row["position"]))
        groups[key].append(row)

with open(OUTPUT, "w", newline="") as f:
    writer = csv.writer(f)

    writer.writerow([
        "chr",
        "position",
        "supporting_reads",
        "min_mapq",
        "max_mapq",
        "mean_clip_length",
        "min_clip_length",
        "max_clip_length",
        "qc_flag"
    ])

    for (chrom, position), rows in sorted(groups.items()):
        mapqs = [int(r["mapq"]) for r in rows]
        clips = [int(r["clip_length"]) for r in rows]

        if len(rows) > 1 and min(mapqs) >= 60:
            qc_flag = "multi_read_MAPQ60"
        elif len(rows) == 1 and min(mapqs) >= 60:
            qc_flag = "single_read_MAPQ60"
        elif len(rows) == 1 and max(mapqs) == 0:
            qc_flag = "single_read_MAPQ0"
        else:
            qc_flag = "review"

        writer.writerow([
            chrom,
            position,
            len(rows),
            min(mapqs),
            max(mapqs),
            f"{mean(clips):.1f}",
            min(clips),
            max(clips),
            qc_flag
        ])

print(f"Site summaries written: {len(groups)}")
print(f"Output: {OUTPUT}")
