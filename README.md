# ONT LTR-Based Integration Site Analysis

## Overview

A reproducible wet-lab and computational workflow for identifying candidate genomic integration sites of an LTR-containing construct using Oxford Nanopore (ONT) sequencing.

## Workflow

### Wet Lab

Genomic DNA
→ LTR-targeted linear PCR
→ Streptavidin bead capture
→ Second-strand synthesis
→ Linker ligation
→ Final PCR
→ ONT sequencing

### Dry Lab

ONT reads
→ minimap2 alignment
→ LTR soft-clip detection
→ Exact LTR matching
→ Breakpoint identification
→ Candidate integration sites

## Methods

- Reference: GRCh38 + provirus sequence
- Chromosomes analysed: chr2, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr22
- Alignment: minimap2
- Analysis: Python 3
- LTR criterion: ≥20 bp soft clip with an exact match to the canonical LTR
- Primary alignments were analysed and both LTR orientations were considered.

## Experimental Primers

Biotinylated LTR primer:

Biotin-GAGCTCTCTGGCTAACTAGG

Final PCR LTR primer:

GAGCTCTCTGGCTAACTAGG

Linker cassette primer:

GACCCGGGAGATCTGAATTC

## Results

A total of 20 exact LTR-supporting reads were identified across 6 candidate genomic positions.

| Candidate site | Reads | MAPQ |
|---|---:|---:|
| chr12:4,083,263 | 1 | 60 |
| chr12:115,354,578 | 1 | 60 |
| chr12:120,865,490 | 15 | 60 |
| chr2:66,315,975 | 1 | 60 |
| chr7:21,552,438 | 1 | 60 |
| chr9:42,120,176 | 1 | 0 |

chr12:120,865,490 showed the strongest repeated read-level support, accounting for 15/20 (75%) of the exact LTR-supported reads.

## Repository

scripts/    → analysis scripts
reference/  → reference sequences
results/    → final result tables

Main outputs:

- results/final_integration_sites.csv
- results/integration_site_summary.csv
- results/all_LTR_matches.txt

Large sequencing and alignment files are excluded from the repository.

## Conclusion

This workflow identifies candidate LTR–genome junctions from ONT reads.

Under the applied computational criteria, chr12:120,865,490 is the most strongly supported candidate in this dataset.

Candidate sites require independent experimental validation for confirmation.
