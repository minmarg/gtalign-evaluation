export WRKDIR=/home/mindaugas/projects/data

## prepare structures so that all methods interpret them equivalently: no residues w/o " CA [ |A]" atoms;
## leave only residues that have at least N, CA, C, and O atoms, which are valid for Dali;
## also, remove HETATM residues which makes foldseek and TM-align read the same sequence of residues
mkdir ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk
ls -1 ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.org | xargs -i -P 40 perl -e '
  open(F,">$ENV{WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk/{}") or die "ERROR: Failed to open file for writing.";
  $prvchn = "_"; $prvstr="    "; $prvnum=0;
  while(<>) {
    if(/^(?:ATOM|TER)/) {
      if((substr($_,22,4) ne $prvstr) || (substr($_,21,1) ne $prvchn)) {
        $prvchn = substr($_,21,1);
        $prvstr = substr($_,22,4);
        $prvnum++;
        print(F $rec) if($n && $ca && $c && $o);
        $rec = ""; $n = $ca = $c = $o = 0;
      }
##      substr($_,21,1) = "A";
##      substr($_,22,4) = sprintf("%4d", $prvnum);
      $ca = 1 if(substr($_,12,4) eq " CA " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $n = 1  if(substr($_,12,4) eq " N  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $c = 1  if(substr($_,12,4) eq " C  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $o = 1  if(substr($_,12,4) eq " O  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $rec .= $_;
    }
    else {
      print(F $rec) if($n && $ca && $c && $o && $rec);
      $rec = ""; $n = $ca = $c = $o = 0;
      print(F) unless /^HETATM/;
    }
    last if /^ENDMDL/;
  }
  print(F $rec) if($n && $ca && $c && $o && $rec);
  close(F)' ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.org/{}

## remove empty structures and those with <3 residues 
rm $(ls -1S ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk/*|tail -n62)

## verification
ls -1 pdb20_211110_pdb.bmk|xargs -i perl -e 'while(<>){next unless /^(?:ATOM|HETATM)/; if(substr($_,12,4) eq " CA " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A")){ if($rn && $rn eq substr($_,22,5)){print"{}\n$_\n";} $rn = substr($_,22,5); }}' pdb20_211110_pdb.bmk/{}

## files with .pdb extension and their list for fatcat
mkdir ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk.pdb
ls -1 ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk|xargs -i sh -c 'n=$(basename {} .ent); ln -s ../pdb20_211110_pdb.bmk/{} ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk.pdb/${n}.pdb'
(ls -1 ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk|xargs -i sh -c 'echo $(basename {} .ent)') >${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk.lst

## make list of files with extension
(cat ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk.lst|xargs -i echo {}.ent) >${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk.ext.lst


## prepare the query structures the same way
mkdir ${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk
ls -1 ${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets | xargs -i -P 40 perl -e '
  open(F,">$ENV{WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk/{}") or die "ERROR: Failed to open file for writing.";
  $prvchn = "_"; $prvstr="    "; $prvnum=0;
  while(<>) {
    if(/^(?:ATOM|TER)/) {
      if((substr($_,22,4) ne $prvstr) || (substr($_,21,1) ne $prvchn)) {
        $prvchn = substr($_,21,1);
        $prvstr = substr($_,22,4);
        $prvnum++;
        print(F $rec) if($n && $ca && $c && $o);
        $rec = ""; $n = $ca = $c = $o = 0;
      }
##      substr($_,21,1) = "A";
##      substr($_,22,4) = sprintf("%4d", $prvnum);
      $ca = 1 if(substr($_,12,4) eq " CA " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $n = 1  if(substr($_,12,4) eq " N  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $c = 1  if(substr($_,12,4) eq " C  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $o = 1  if(substr($_,12,4) eq " O  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $rec .= $_;
    }
    else {
      print(F $rec) if($n && $ca && $c && $o && $rec);
      $rec = ""; $n = $ca = $c = $o = 0;
      print(F) unless /^HETATM/;
    }
    last if /^ENDMDL/;
  }
  print(F $rec) if($n && $ca && $c && $o && $rec);
  close(F)' ${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets/{}

## verification
ls -1 cameo3d_2021.07.24_2021.10.16_targets.bmk |xargs -i perl -e 'while(<>){next unless /^(?:ATOM|HETATM)/; if(substr($_,12,4) eq " CA " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A")){ if($rn && $rn eq substr($_,22,5)){print"{}\n$_\n";} $rn = substr($_,22,5); }}' cameo3d_2021.07.24_2021.10.16_targets.bmk/{}



## FATCAT
## env var
export FATCAT=/home/mindaugas/install/FATCAT/FATCAT-dist
## run FATCAT (does not read queries with absolute path: use cd)
bash -c '(cd ${WRKDIR}; time (ls -1 cameo3d_2021.07.24_2021.10.16_targets.bmk | xargs -i sh -c "echo {}; time /home/mindaugas/install/FATCAT/FATCAT-dist/FATCATMain/FATCATSearch cameo3d_2021.07.24_2021.10.16_targets.bmk/{} \${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk.lst -i2 \${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk.pdb -r -o \${WRKDIR}/gtalign-benchmark-pdb20/fatcat/fatcat_{}.out -m")) >fatcat.log 2>&1 &'
## check #alignments
ls -1 fatcat/fatcat_*|xargs -i sh -c 'echo {} $(grep -E Align {}|wc -l)'|awk '{print $2}'|uniq -c



## DeepAlign_Search does not produce alignments 
## run DeepAlign for alignments
bash -c '(QRS=($(ls -1 ${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk)); (for((i=0; i<${#QRS[@]}; i++)); do query=${QRS[$i]}; echo -e "\n\n${query}"; time (ls -1 ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk | xargs -i -P 40 sh -c "(/home/mindaugas/install/DeepAlign/DeepAlign/DeepAlign ${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk/${query} ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk/{} || true)") >deepalign/deepalign_${query}.out; done)) >deepalign.log 2>&1 &'
## (rerun the minority of alignments distorted by processes' simultaneous writes)
(grep -E '^ERROR' deepalign__pdb20__processed_1.err | xargs -i sh -c 'rfn=$(echo "{}"|awk "{print \$7}"); qry=$(echo "{}"|awk "{print \$6}"); qryname=$(basename $qry); (/home/mindaugas/install/DeepAlign/DeepAlign/DeepAlign $qry $rfn || true) >>deepalign_rerun/deepalign_${qryname}.out') >deepalign_rerun.log 2>&1
## check #alignments
ls -1 deepalign/deepalign_*|xargs -i sh -c 'echo {} $(grep -E " DeepAlign" {}|wc -l)'|awk '{print $2}'|uniq -c
## total time (40 threads)
grep -E '^real' deepalign.log |perl -e 'while(<>){@a=split(/\s+/);$s+=$1*60+$2 if($a[1]=~/^([^m\s]+)m([^s\s]+)s/);} print($s,"\n")'



## DALI v5
## reformat structure files for Dali format 
## (NOTE: overwriting when a filename includes a multiple-character chain id)
## (hence, use manual ids instead):
(cd ${WRKDIR}; tmpdir="wrk_dir_dali"; [ ! -d ${tmpdir} ] && mkdir ${tmpdir}; cd ${tmpdir}; datdir="../gtalign-benchmark-pdb20/DALIv5_DAT_db"; [ ! -d ${datdir} ] && mkdir ${datdir}; flist=($(ls -1 ../pdb20_211110/pdb20_211110_pdb.bmk)); for((i=0; i<${#flist[@]}; i++)); do f=${flist[$i]}; printf -v id "%04x" $i; echo "$f $id" >>../gtalign-benchmark-pdb20/DALIv5_DAT_db_map.txt; time /home/mindaugas/install/DaliLite_v5/DaliLite.v5/bin/import.pl --pdbfile ../pdb20_211110/pdb20_211110_pdb.bmk/$f --dat ${datdir} --pdbid $id; done; cd ..) >dali5_pdb20_import.log 2>&1
## make list:
(ls -1 DALIv5_DAT_db|xargs -i basename {} .dat) >DALIv5_DAT_db.lst
## reformat query structures:
(cd ${WRKDIR}; tmpdir="wrk_dir_dali"; [ ! -d ${tmpdir} ] && mkdir ${tmpdir}; cd ${tmpdir}; datdir="../gtalign-benchmark-pdb20/DALIv5_DAT_qrs"; [ ! -d ${datdir} ] && mkdir ${datdir}; flist=($(ls -1 ../cameo3d_2021.07.24_2021.10.16_targets.bmk)); for((i=0; i<${#flist[@]}; i++)); do f=${flist[$i]}; printf -v id "%04x" $i; echo "$f $id" >>../gtalign-benchmark-pdb20/DALIv5_DAT_qrs_map.txt; time /home/mindaugas/install/DaliLite_v5/DaliLite.v5/bin/import.pl --pdbfile ../cameo3d_2021.07.24_2021.10.16_targets.bmk/$f --dat ${datdir} --pdbid $id; done; cd ..) >dali5_pdb20_import_qrs.log 2>&1
## queries that failed (no .dat file)
ls -1 DALIv5_DAT_qrs/|cut -b1-4|perl -e 'while(<>){chomp; $nxt=hex($prv)+1; printf("%04x\n",$nxt) if($prv && hex($_)-$nxt); $prv=$_}' | xargs -i grep -E ' {}' DALIv5_DAT_qrs_map.txt
## make list:
(ls -1 DALIv5_DAT_qrs|xargs -i basename {} .dat) >DALIv5_DAT_qrs.lst
## run DaliLitev5
bash -c '(cd dali; time /home/mindaugas/install/DaliLite_v5/DaliLite.v5/bin/dali.pl --np 40 --query ../DALIv5_DAT_qrs.lst --db ../DALIv5_DAT_db.lst --dat1 ../DALIv5_DAT_qrs --dat2 ../DALIv5_DAT_db --outfmt alignments; cd ..) >dali.log 2>&1 &'



## TM-align
## run TM-align
bash -c '(cd ${WRKDIR}; time (ls -1 cameo3d_2021.07.24_2021.10.16_targets.bmk | xargs -i -P 40 sh -c "echo {} \$(time /home/mindaugas/install/TM-align/TMalign cameo3d_2021.07.24_2021.10.16_targets.bmk/{} -dir2 \${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk/ \${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk.ext.lst -het 1 >\${WRKDIR}/gtalign-benchmark-pdb20/tmalign/tmalign_{}.out)")) >tmalign_pdb20.log 2>&1 &'
## check #alignments
ls -1 tmalign/tmalign_*|xargs -i sh -c 'echo {} $(grep -E Aligned {}|wc -l)'|awk '{print $2}'|uniq -c
## run TM-align fast
bash -c '(cd ${WRKDIR}; time (ls -1 cameo3d_2021.07.24_2021.10.16_targets.bmk | xargs -i -P 40 sh -c "echo {} \$(time /home/mindaugas/install/TM-align/TMalign cameo3d_2021.07.24_2021.10.16_targets.bmk/{} -dir2 \${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk/ \${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk.ext.lst -het 1 -fast >\${WRKDIR}/gtalign-benchmark-pdb20/tmalign_fast/tmalign_{}.out)")) >tmalign_fast_pdb20.log 2>&1 &'



## foldseek
## run foldseek using --alignment-type 1
(time ~/install/foldseek/foldseek/bin/foldseek easy-search ${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk/ ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk/ foldseek_tmfast1__pdb20_bmk.res tmp_dir_foldseek --threads 40 --max-seqs 4000 --alignment-type 1 --format-output query,target,evalue,qtmscore,ttmscore,alntmscore,rmsd,qstart,qend,qlen,tstart,tend,tlen,qaln,taln) >foldseek_tmfast1__pdb20_bmk.log 2>&1
## run foldseek using --alignment-type 2 which uses 3d alphabet:
(time ~/install/foldseek/foldseek/bin/foldseek easy-search ${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk/ ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk/ foldseek_alntyp2__pdb20_bmk.res tmp_dir_foldseek --threads 40 --max-seqs 4000 --alignment-type 2 --format-output query,target,evalue,qtmscore,ttmscore,alntmscore,rmsd,qstart,qend,qlen,tstart,tend,tlen,qaln,taln) >foldseek_alntyp2__pdb20_bmk.log 2>&1
## run foldseek using --alignment-type 1 and --tmalign-fast 0
bash -c '(time ~/install/foldseek/foldseek/bin/foldseek easy-search ${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk/ ${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk/ foldseek_tmfast0__pdb20_bmk.res tmp_dir_foldseek --threads 40 --max-seqs 4000 --alignment-type 1 --format-output query,target,evalue,qtmscore,ttmscore,alntmscore,rmsd,qstart,qend,qlen,tstart,tend,tlen,qaln,taln --tmalign-fast 0) >foldseek_tmfast0__pdb20_bmk.log 2>&1 &'



## GTalign speed 0
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk --rfs=${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk -o gtalign_14_speed0_prescore03_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=5000 --speed=0 --pre-score=0.3 -s 0.44 --add-search-by-ss --nhits=4000 --nalns=4000 --dev-N=3) 2>&1 | tee gtalign_14_speed0_prescore03_addss_s044.log
## GTalign speed 9
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk --rfs=${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk -o gtalign_14_speed9_prescore03_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=5000 --speed=9 --pre-score=0.3 -s 0.44 --add-search-by-ss --nhits=4000 --nalns=4000 --dev-N=3) 2>&1 | tee gtalign_14_speed9_prescore03_addss_s044.log
## GTalign speed 9 + --pre-score=0.4
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk --rfs=${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk -o gtalign_14_speed9_prescore04_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=5000 --speed=9 --pre-score=0.4 -s 0.44 --add-search-by-ss --nhits=4000 --nalns=4000 --dev-N=3) 2>&1 | tee gtalign_14_speed9_prescore04_addss_s044.log
## GTalign speed 13
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk --rfs=${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk -o gtalign_14_speed13_prescore03_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=5000 --speed=13 --pre-score=0.3 -s 0.44 --add-search-by-ss --nhits=4000 --nalns=4000 --dev-N=3) 2>&1 | tee gtalign_14_speed13_prescore03_addss_s044.log
## GTalign speed 13 + --pre-score=0.4
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk --rfs=${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk -o gtalign_14_speed13_prescore04_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=5000 --speed=13 --pre-score=0.4 -s 0.44 --add-search-by-ss --nhits=4000 --nalns=4000 --dev-N=3) 2>&1 | tee gtalign_14_speed13_prescore04_addss_s044.log
## GTalign speed 13 + --pre-similarity=15
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/cameo3d_2021.07.24_2021.10.16_targets.bmk --rfs=${WRKDIR}/pdb20_211110/pdb20_211110_pdb.bmk -o gtalign_14_speed13_prescore03_presim15_addss_s03 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=5000 --speed=13 --pre-score=0.3 --pre-similarity=15 -s 0.3 --add-search-by-ss --nhits=4000 --nalns=4000 --dev-N=3 -c cachedir) 2>&1 | tee gtalign_14_speed13_prescore03_presim15_addss_s03.log


