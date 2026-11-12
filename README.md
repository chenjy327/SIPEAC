<p align="center">
  <img src="sipeac_title.svg" alt="SIPEAC" width="640">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/R-4.2-276DC3?logo=r&logoColor=white" alt="R">
  <img src="https://img.shields.io/badge/Python-3-3776AB?logo=python&logoColor=white" alt="Python">
  <img src="https://img.shields.io/badge/Shell-bash-4EAA25?logo=gnubash&logoColor=white" alt="Shell">
  <a href="https://github.com/chenjy327/SIPEAC/"><img src="https://img.shields.io/github/repo-size/chenjy327/SIPEAC" alt="repo size"></a>
  <a href="https://github.com/chenjy327/SIPEAC/commits/main"><img src="https://img.shields.io/github/last-commit/chenjy327/SIPEAC" alt="last commit"></a>
  <a href="https://github.com/chenjy327/SIPEAC/"><img src="https://img.shields.io/github/stars/chenjy327/SIPEAC?style=social" alt="GitHub stars"></a>
</p>

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
