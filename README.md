[English](README.md) | [Русский](README.ru.md)

# GRDH - guide RNA design helper

**A desktop application for searching, selecting, and analyzing guide RNAs (gRNAs) for the CRISPR/Cas9 system.**

Runs locally, with no server connection or execution queues.

---

## Table of Contents

- [Features](#features)
- [Comparison with Other Tools](#comparison-with-other-tools)
- [Installation](#installation)
- [User Guide](#user-guide)
- [Roadmap](#roadmap)
- [Authors](#authors)
- [License](#license)

---

## Features

### Input Sequence Sources

- `.fasta` file with optional region boundaries
- Custom sequence entered manually
- NCBI nucleotide locus ID with optional region boundaries
- CRISPOR output file (`.xls`) - analyzed on equal footing with local search

### gRNA Scoring & Evaluation

- **Doench 16 (Azimuth)** - cleavage efficiency; the model was extracted from the original Azimuth and reproduces its results bit-exactly
- **OTS (Our Total Score)** - a custom score accounting for:
  - GC content
  - oligonucleotide repeats
  - nucleotides at key spacer positions (3, 14, 16, 18, 20)
  - presence of GCC at positions 16–20
  - spacer accessibility in the secondary structure (positions 18–20 and 51–53)
  - 8 or more consecutive paired nucleotides and more than 12 paired nucleotides in total
  - N-position of the PAM
- **Restriction site search** near the cut site, with counts within 100, 250, and 500 bp radii

### gRNA Secondary Structure

- Dot-bracket structure and MFE calculated by the [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) library
- 2D structure rendering to `.png` - custom implementation (not the standard ViennaRNA visualizer)
- Accessibility calculation for key spacer positions
- Structure image saved for each gRNA

### Interactive Genomic Map

- Coordinate-based zoom - from full-sequence overview down to individual nucleotides
- Arrows indicating gRNA direction (sense/antisense), scaled with zoom
- Red Cas9 cut-site line; PAM displayed alongside the spacer
- Automatic lane packing for overlapping gRNAs
- Hover tooltips with gRNA parameters

### gRNA Table

- Bidirectional sync: selecting a gRNA on the map highlights the corresponding table row and vice versa
- Delete gRNAs from the table with automatic map update
- Cell highlighting by quality criteria: Doench < 50, T in the N-position of PAM, oligo-T, GCC at 16–20, 8+ consecutive paired nucleotides, 12+ total paired nucleotides
- Sort by any column

### Paired gRNA Selection Mode

- Select an anchor gRNA; the rest are automatically filtered by cut-site distance (50–2000 bp)
- Suitable for paired-nickase and deletion experiments
- Reset button to restore the full list

### Export

- Export remaining table gRNAs to `.xlsx` with structure images
- All structure images from the current session are saved in the `output` folder

---

## Comparison with Other Tools

Compared against **CRISPOR**, **CRISPR-P 2.0**, **CRISPRon**, **CRISPR DB**, **sgDesigner**, and **RNAfold** (ViennaRNA).

### Unique GRDH Features

| Feature | Description |
|---|---|
| Paired gRNA selection | Filter by cut-site distance (50–2000 bp) - not available in any of the listed tools |
| Bidirectional map ↔ table sync | Selecting a gRNA on the map highlights the table row and vice versa; deletion from the table updates the map |
| Structure + map + scoring in one tool | Secondary structure, interactive map, and custom scoring combined in a single application |
| CRISPOR import | Import CRISPOR `.xls` output, then add secondary structure, OTS, and map visualization |
| Local GUI application | CRISPOR and CRISPR-P 2.0 are web-only; sgDesigner and CRISPRon are command-line only |

### Strengths

| Feature | GRDH | CRISPOR | CRISPR-P 2.0 | CRISPRon | CRISPR DB | sgDesigner | RNAfold |
|---|---|---|---|---|---|---|---|
| Doench 16 (Azimuth) | + | + | - (Doench 14) | - | - | - | - |
| gRNA secondary structure | + custom rendering | - | + ViennaRNA visualizer (dot-bracked structure only) | ⚠️ folding energy, no rendering | - | - | + (structure only) |
| Restriction sites | + 100/250/500 bp | + | + | - | - | - | - |
| Interactive map | + zoom to nucleotides | - table only | + | - | - | - | - |
| Off-target | - (planned) | + MIT + CFD | + | + CRISPRoff | - | - | - |

### Limitations

- **Off-target analysis is not implemented** (planned). CRISPOR uses MIT and CFD scores; CRISPR-P 2.0 uses its own scoring. For off-target evaluation in this version, we recommend using CRISPOR, whose output GRDH can import
- **Only SpCas9 is supported.** CRISPR-P 2.0 and CRISPOR support other Cas systems (planned)
- **No cloning primer design** - available in CRISPOR
- **No microhomology score** - available in CRISPR-P 2.0
- **OTS is a heuristic score**, not trained on experimental data. CRISPRon and CRISPR DB likely predict cleavage efficiency more accurately
- **No genome context** - GRDH works with arbitrary sequences. CRISPOR and CRISPR-P 2.0 search for gRNAs in the context of a selected genome

---

## User Guide

1. Select a gRNA search mode: `.fasta` / custom sequence / NCBI / CRISPOR.
2. **a.** For `.fasta`, specify the file path. Optionally set region boundaries.\
   **b.** For a custom sequence, paste it into the input field (ATGC only, max 5 kbp).\
   **c.** For NCBI, enter the nucleotide locus ID. Optionally set region boundaries.\
   **d.** For CRISPOR, specify the path to the CRISPOR `.xls` results file.
3. For local search, the target sequence length must be between 50 bp and 5 kbp. If an error occurs, adjust the boundary values.
4. Click the "Analyze" button. On successful start, a progress bar, a gRNA counter, and the current spacer sequence will appear.
5. After analysis completes, a map and a results table will appear.
   - Zooming the map with the mouse wheel reveals individual nucleotides. The red line marks the Cas9 cut site; the PAM is displayed next to the spacer. Hovering over a spacer shows its parameters.
   - Clicking a gRNA on the map highlights the corresponding table row. Selecting a table row highlights the gRNA on the map.
6. The table contains key gRNA information. It is initially sorted by Doench 16; values below 50 are highlighted in red.
   - Sequence and PAM cells are highlighted when: T is in the N-position of the PAM (NGG), 8+ consecutive paired nucleotides or 12+ total paired nucleotides in the spacer, GCC at positions 16–20, or four consecutive T's.
   - You can delete selected gRNAs or all red-flagged ones - the map updates automatically.
7. For paired gRNA selection, click the "Pair selection" button, then click a gRNA on the map. Only gRNAs with a cut site within 50–2000 bp of the selected one remain. The "Reset filter" button restores the full list.
8. Export to `.xlsx` includes only the gRNAs remaining in the table. All structure images from the current session are stored in the `output` folder.

---

## Roadmap
- [ ] Building an application for Windows 11 and Archlinux
- [ ] English translation and locale switching
- [ ] Read target sequences from `.gb` files and add `.gb` annotations to the genomic map
- [ ] Off-target search against a selected genome (`.fasta`, `.gb`, NCBI ID)
- [ ] Include Doench 14 in scoring
- [ ] Expand the list of supported Cas systems and their scoring
- [ ] Optimize and validate OTS scores

---

## Authors

| Name | Role |
|---|---|
| Alexander D. Lukin | Development and maintenance |
| Tatyana A. Frankevich | Testing, biological expertise |
| Natalya V. Permyakova | Scientific supervision, concept, validation of selection criteria, biological expertise, implementation requirements |

**Laboratory of Plant Bioengineering, Institute of Cytology and Genetics, SB RAS (Novosibirsk, Russia)**

The authors thank Sophia A. Markova for providing the icons.

---

## License

GNU Affero General Public License v3.0 - see the [LICENSE](LICENSE) file.
