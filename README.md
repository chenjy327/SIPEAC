# <p>  <b>SIPEAC</b> </p>

![R](https://img.shields.io/badge/R-4.2-276DC3?logo=r&logoColor=white)
![Python](https://img.shields.io/badge/Python-3-3776AB?logo=python&logoColor=white)
![Shell](https://img.shields.io/badge/Shell-bash-4EAA25?logo=gnubash&logoColor=white)
[![repo size](https://img.shields.io/github/repo-size/chenjy327/SIPEAC)](https://github.com/chenjy327/SIPEAC/)
[![last commit](https://img.shields.io/github/last-commit/chenjy327/SIPEAC)](https://github.com/chenjy327/SIPEAC/commits/main)
[![GitHub stars](https://img.shields.io/github/stars/chenjy327/SIPEAC?style=social)](https://github.com/chenjy327/SIPEAC/)

**SIPEAC (Spatial Identification of Paired Endogenous Antibody Heavy-Light Chains) preferentially targets plasma cells and couples spatial transcriptomics with the structurally paired preference of VH:VL sequences to clone tissue-resident monoclonal antibodies in situ.**

## Code

| File | Description |
|---|---|
| ccs2dedup.sh, run.pacbio.sh | SMRT-seq processing from CCS reads to deduplicated BCR/TCR molecules |
| Step2.ReverseComplement_pTreads.py | Step 2: reverse-complement reads |
| Step3.Barcode.UMI.detection.py | Step 3: spatial barcode and UMI detection |
| Step4.rmpA.py | Step 4: polyA removal |
| Step6.dedup.py | Step 6: UMI deduplication |
| Step7.map2ST.py | Step 7: map molecules to Stereo-seq coordinates |
| Convert2BinX.py | Assign SMRT-seq molecules to BinX spatial bins |
| run.changeo.sh | IgBLAST and Change-O annotation, SHM and clonal grouping; calls add_info.R and rename_column.R |
| add_info.R | Filter Change-O output, compute SHM and the nearest-neighbor distance plot |
| rename_column.R | Keep the clone family ID column after DefineClones |
| run.repair.sh | Run the repair pairing tool and evaluate it against the ground truth (benchmark for SIPEAC) |
| antibody_HL_SIPEAC/ | ESM-based heavy/light chain model (training and inference) |
| Figure1.R | Figure 1 |
| Figure2,E3.R | Figure 2 and Extended Data Fig. 3 |
| Figure3.R | Figure 3 |
| Figure E1.R | Extended Data Fig. 1 |
| Figure E2.R | Extended Data Fig. 2 |
| figure_helpers.R | Shared functions and colors sourced by the Figure scripts (run from the repository root; set `SIPEAC_PYTHON` to change the Python used for spatial plots) |
| Extended_Analyses/ | Additional analyses; see its README |
