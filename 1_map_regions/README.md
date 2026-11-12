# Sequence Homology Mapping for FUNCODE

This step maps regulatory elements (e.g., DHS peaks) between species using pairwise sequence alignment. FUNCODE performs **bidirectional mapping** to create a comprehensive set of sequence-homologous regulatory element pairs (REPs).

## Overview

FUNCODE's mapping approach is highly similar to UCSC's liftOver tool - both use the same underlying netted LASTZ alignments. The key distinction is that FUNCODE maps summit positions (single nucleotide) and then extends to fixed-width regions, ensuring mapped regions are centered on orthologous regulatory element summits.

## Bidirectional Mapping Strategy

FUNCODE maps DHSs in **both directions** and combines them to create the final core set of REPs:

```
Forward Mapping (hg38 → mm10):
  Human DHS summits → Aligned mouse positions → 200bp windows

Reverse Mapping (mm10 → hg38):
  Mouse DHS midpoints → Aligned human positions → 200bp windows

Final Core Set = Forward Set + Non-duplicate Reverse Mappings
```

### Duplicate Removal
A pair from the reverse mapping is considered a duplicate if both its human and mouse elements overlap by >50bp with any pair in the forward set. After removing duplicates, the remaining reverse-mapped pairs are combined with the forward set.

## Input Data

### 1. Regulatory Element Peaks
- **Human**: DHS peaks from ENCODE (e.g., ENCFF503GCK)
- **Mouse**: DHS peaks from ENCODE (e.g., ENCFF910SRW)

### 2. Pairwise Sequence Alignments
Download from UCSC Genome Browser:

| Direction | Alignment File | Source |
|-----------|---------------|--------|
| hg38 → mm10 | `hg38.mm10.net.axt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm10/ |
| mm10 → hg38 | `mm10.hg38.net.axt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsHg38/ |

### 3. Chromosome Size Files
- Human hg38: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/
- Mouse mm10: https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/

## Pipeline Scripts

### Forward Mapping (Human → Mouse)

```bash
# 1. Prepare human DHS summit coordinates
Rscript R/make_region_file.R

# 2. Split axt.net file by chromosome
# (preprocessing step)

# 3. Map summits to mouse genome
python py/base_run_hg38tomm10.py

# 4. Calculate sequence identity
python py/two_region_hg38tomm10.py
```

### Reverse Mapping (Mouse → Human)

```bash
# 1. Map mouse DHS midpoints to human genome
python py/base_run_mm10tohg38.py

# 2. Calculate sequence identity
python py/two_region_mm10tohg38.py
```

### Zebrafish Extension
For human-zebrafish mapping:
```bash
python py/base_run_hg38todanRer10.py
```

### Pig Extension
For human-pig mapping (run from the project root; the single genome-wide `hg38.susScr11.net.axt` from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsSusScr11/ is searched on the hg38 side for forward mapping and on the susScr11 side for reverse mapping):
```bash
python 1_map_regions/py/base_run_hg38tosusScr11.py
python 1_map_regions/py/base_run_susScr11tohg38.py
Rscript 1_map_regions/R/make_region_file_hg38_susScr11.R
python 1_map_regions/py/two_region_hg38tosusScr11.py
```

## Output
Each mapping produces:
- Paired coordinates (200bp windows in each species)
- Alignment block IDs
- Strand orientation
- Gap status
- Percent sequence homology (PSH)

## Comparison with liftOver

| Feature | FUNCODE | liftOver |
|---------|---------|----------|
| Underlying alignment | UCSC netted LASTZ | UCSC netted LASTZ |
| Mapping unit | Summit/midpoint → 200bp window | Arbitrary interval |
| Bidirectional | Yes (computed both ways) | One direction per run |
| Sequence identity | Computed and stored | Not provided |
| Alignment metadata | Block IDs, strand, gaps | Minimal |

## References

- UCSC Genome Browser: https://genome.ucsc.edu/
- LASTZ documentation: https://lastz.github.io/lastz/
- ENCODE Portal: https://www.encodeproject.org/

## File Structure

```
1_map_regions/
├── py/
│   ├── base_module.py           # Core mapping functions
│   ├── region_module.py         # Region handling utilities
│   ├── base_run_hg38tomm10.py   # Human → Mouse mapping
│   ├── base_run_mm10tohg38.py   # Mouse → Human mapping
│   ├── base_run_hg38todanRer10.py  # Human → Zebrafish mapping
│   ├── base_run_hg38tosusScr11.py  # Human → Pig mapping
│   ├── base_run_susScr11tohg38.py  # Pig → Human mapping
│   ├── two_region_hg38tomm10.py    # Sequence identity (H→M)
│   ├── two_region_mm10tohg38.py    # Sequence identity (M→H)
│   └── two_region_hg38tosusScr11.py  # Sequence identity (H→P)
├── R/
│   ├── make_region_file.R       # Prepare input region files
│   └── make_region_file_hg38_susScr11.R  # Human-pig core set
└── README.md
```
