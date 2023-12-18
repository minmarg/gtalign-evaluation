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
The commands in the scripts are commented.
The files `plot_TMscores_<dataset>.sh` are scripts to generate the main results figures.
The file `PDB_clustered_gtalign_CLC.lst` contains complete-linkage clusters for the entire PDB (08/18/2023; first chains).

