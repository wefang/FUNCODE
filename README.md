## Functional Conservation of DNA Elements (FUNCODE)

Scripts for computing FUNCODE scores from ENCODE DNase-seq, ATAC-seq and Histone ChIP-seq data.

### UCSC Browser track hubs:
[hg38](https://genome.ucsc.edu/cgi-bin/hgTracks?hubUrl=https://raw.githubusercontent.com/wefang/funcode/main/track_hubs/hub/hub.txt&genome=hg38&position=chr14:53920000-53990000)
[mm10](https://genome.ucsc.edu/cgi-bin/hgTracks?hubUrl=https://raw.githubusercontent.com/wefang/funcode/main/track_hubs/hub/hub.txt&genome=mm10&position=chr14:46352000-46422000)

![image](https://github.com/user-attachments/assets/6ccbdfd9-aa9e-4140-970c-fffc46db1371)

![image](https://github.com/user-attachments/assets/46fcaad6-dfc8-42f4-a8f7-6479f83e06b1)


### Output types on ENCODE data portal
These scripts produce data corresponding to two output types on the ENCODE data portal:

(1) "functional conservation quantification"

(2) "functional conservation mapping".

- "functional conservation mapping" is relevant if the user would like to map regulatory elements across species.

The following datasets contain "functional conservation mapping":
| Accession   | Description                                                                                                                                                                                                               |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| ENCSR174NMF | FUNCODE conservation scores of context-dependent variable acitivity (CO-V) for chromatin accessibility at human(hg38)-mouse(mm10) DHS pairs derived from nearby gene homology                                             |
| ENCSR079JCW | FUNCODE conservation scores of context-dependent variable acitivity (CO-V) for H3K27ac histone modification at human(hg38)-mouse(mm10) DHS pairs derived from nearby gene homology                                        |
| ENCSR799PIS | FUNCODE conservation scores of context-dependent variable acitivity (CO-V) for H3K4me1 histone modification at human(hg38)-mouse(mm10) DHS pairs derived from nearby gene homology                                        |
| ENCSR241JXZ | FUNCODE conservation scores of context-dependent variable acitivity (CO-V) for H3K4me3 histone modification at human(hg38)-mouse(mm10) DHS pairs derived from nearby gene homology                                        |
| ENCSR705SNU | FUNCODE conservation scores of baseline (CO-B) and context-dependent variable activity (CO-V) for chromatin accessibility, H3K27ac, H3K4me1, H3K4me3 at human(hg38)-mouse(mm10) DHS pairs derived from sequence alignment |


- "functional conservation quantification" is relevant if the user would like to visualize the functional conservation signal in either human or mouse genome. Both '.bed' and '.bigWig' files are provided.

The following datasets contain "functional conservation quantification":
| Accession   | Description                                                                                                                                                                |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| ENCSR666RQQ | FUNCODE conservation scores of baseline activity (CO-B) for chromatin accessibility at human DHSs (hg38) paired to mouse using sequence alignment                          |
| ENCSR383DSZ | FUNCODE conservation scores of baseline activity (CO-B) for chromatin accessibility at mouse DHSs (mm10) paired to human using sequence alignment                          |
| ENCSR794MIL | FUNCODE conservation scores of baseline activity (CO-B) for H3K27ac histone modification at human DHSs (hg38) paired to mouse using sequence alignment                     |
| ENCSR468BFL | FUNCODE conservation scores of baseline activity (CO-B) for H3K27ac histone modification at mouse DHSs (mm10) paired to human using sequence alignment                     |
| ENCSR166FHJ | FUNCODE conservation scores of baseline activity (CO-B) for H3K4me1 histone modification at human DHSs (hg38) paired to mouse using sequence alignment                     |
| ENCSR312OMI | FUNCODE conservation scores of baseline activity (CO-B) for H3K4me1 histone modification at mouse DHSs (mm10) paired to human using sequence alignment                     |
| ENCSR623QIK | FUNCODE conservation scores of baseline activity (CO-B) for H3K4me3 histone modification at human DHSs (hg38) paired to mouse using sequence alignment                     |
| ENCSR840YRP | FUNCODE conservation scores of baseline activity (CO-B) for H3K4me3 histone modification at mouse DHSs (mm10) paired to human using sequence alignment                     |
| ENCSR289SXF | FUNCODE conservation scores of context-dependent variable activity (CO-V) for chromatin accessibility at human DHSs (hg38) paired to mouse using sequence alignment        |
| ENCSR989MYS | FUNCODE conservation scores of context-dependent variable activity (CO-V) for chromatin accessibility at mouse DHSs (mm10) paired to human using sequence alignment        |
| ENCSR542FJK | FUNCODE conservation scores of context-dependent variable activity (CO-V) for H3K27ac histone modification at human DHSs (hg38) paired to mouse using sequence alignment   |
| ENCSR981BVO | FUNCODE conservation scores of context-dependent variable activity (CO-V) for H3K27ac histone modification at mouse DHSs (mm10) paired to human using sequence alignment   |
| ENCSR825SLF | FUNCODE conservation scores of context-dependent variable activity (CO-V) for H3K4me1 histone modification at human DHSs (hg38) paired to mouse using sequence alignment   |
| ENCSR275KOD | FUNCODE conservation scores of context-dependent variable activity (CO-V) for H3K4me1 histone modification at mouse DHSs (mm10) paired to human using sequence alignment   |
| ENCSR012VAU | FUNCODE conservation scores of context-dependent variable activity (CO-V) for H3K4me3 histone modification at human DHSs (hg38) paired to mouse using sequence alignment   |
| ENCSR811LNW | FUNCODE conservation scores of context-dependent variable activity (CO-V) for H3K4me3 histone modification at mouse DHSs (mm10) paired to human using sequence alignment   |
| ENCSR193MVK | FUNCODE conservation scores of context-dependent variable activity (CO-V) for chromatin accessibility at human DHSs (hg38) paired to mouse using nearby gene homology      |
| ENCSR785AKM | FUNCODE conservation scores of context-dependent variable activity (CO-V) for chromatin accessibility at mouse DHSs (mm10) paired to human using nearby gene homology      |
| ENCSR024SYD | FUNCODE conservation scores of context-dependent variable activity (CO-V) for H3K27ac histone modification at human DHSs (hg38) paired to mouse using nearby gene homology |
| ENCSR389ESN | FUNCODE conservation scores of context-dependent variable activity (CO-V) for H3K27ac histone modification at mouse DHSs (mm10) paired to human using nearby gene homology |
| ENCSR211HIA | FUNCODE conservation scores of context-dependent variable activity (CO-V) for H3K4me1 histone modification at human DHSs (hg38) paired to mouse using nearby gene homology |
| ENCSR786OQS | FUNCODE conservation scores of context-dependent variable activity(CO-V) for H3K4me1 histone modification at mouse DHSs (mm10) paired to human using nearby gene homology  |
| ENCSR760RUZ | FUNCODE conservation scores of context-dependent variable activity (CO-V) for H3K4me3 histone modification at human DHSs (hg38) paired to mouse using nearby gene homology |
| ENCSR168ZNZ | FUNCODE conservation scores of context-dependent variable activity (CO-V) for H3K4me3 histone modification at mouse DHSs (mm10) paired to human using nearby gene homology |


### Details:
To reproduce the FUNCODE scores, use the following steps:

1. Generate pairs of mapped DHS elements based on pairwsie sequence alignment (UCSC Blastz)

2. Download and process sequence read alignment files (bam) from ENCODE data portal

3. Integrate data between human and mouse

4. Compute FUNCODE CO-V and CO-B scores for sequence alignable DHS pairs

5. Compute FUNCODE CO-V scores for sequence DHS pairs with nearby gene homology, and calls significantly conserved pairs

6. Code to produce the score tracks in UCSC genome browser

7. Code to reproduce analyses in the manuscript

For details on each step, see README in the corresponding sub-directories.

Report bugs and provide suggestions by sending email to:

Author and maintainer: Weixiang Fang (wfang9@jh.edu)

Contributers: Chaoran Chen, Boyang Zhang
