# GTalign evaluation scripts

GTalign evaluation scripts for benchmarking [GTalign](https://github.com/minmarg/gtalign_alpha) 
against other structure alignment methods/tools on four different datasets. 

## Description

The following notation is adopted to denote datasets: 

  * scope20840, the SCOPe40 2.08 dataset
  * pdb20, the PDB20 dataset
  * AFSP, the Swiss-Prot dataset
  * HOMSTRAD/homstrad, the HOMSTRAD dataset

The files `db-<dataset>.lst` contain the lists of subject structure files. For the Swiss-Prot dataset, the entire Swiss-Prot archive of AlphaFold2 structures was used. 

The files `queries-<dataset>.lst` contain the lists of query structures. The file `families-homstrad.lst` lists 398 reference protein structure alignment files.

The scripts `commands-run-<dataset>.sh` provide commands for preprocessing structure files and running tools.
The scripts `commands-test-<dataset>.sh` provide commands for evaluating alignments produced by each tool.
The `commands-test-<dataset>-rmsd.sh` and `commands-test-<dataset>-gdtts.sh` scripts provide commands for evaluating alignments by RMSD and GDT\_TS, respectively.

The `commands-test-scope20840-fr.sh` and `commands-test-scope20840-fr2.sh` scripts list commands for SCOPe-based evaluation. The difference between these scripts lies in the GetRecall1 function. The version in `commands-test-scope20840-fr2.sh` additionally outputs the effective number of queries at the family, superfamily, and fold levels.
The `commands-avg-scope20840.sh` script calculates average sensitivity (SCOPe-based ROC1) values and performs statistical tests on sensitivity values between GTalign --speed=0 and other tools.

The commands in the scripts are commented.

The files `plot_TMscores_<dataset>.sh` are scripts to generate the main results figures.
The `plot_TMscores_<dataset>_rmsd.sh` and `plot_TMscores_<dataset>_gdtts.sh` scripts generate figures depicting the results of alignment accuracy evaluation using RMSD and GDT\_TS, respectively.
The `plot_TMscores_scope20840-fr.sh` and `plot_TMscores_scope20840-fr2.sh` scripts display the results obtained from the SCOPe-based evaluation. The difference is in plotting Supplementary Figure S7, where `plot_TMscores_scope20840-fr2.sh` plots sensitivity up to the first false positive against the effective fraction of queries, instead of the total fraction of queries, and is used as the latest version.

The `TMscore_from_alignment.cpp` program is the TM-score (version 2022/02/27) program by Zhang & Skolnick (Proteins 57, 2004) with minimal modifications to normalize GDT\_TS by the number of aligned residue pairs.

The file `PDB_clustered_gtalign_CLC.lst` contains complete-linkage clusters for the entire PDB (08/18/2023; first chains).

The `homstrad.tar.bz2` file is the HOMSTRAD dataset originally obtained from http://yanglab.nankai.edu.cn/mTM-align/benchmark.
