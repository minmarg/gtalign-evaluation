export WRKDIR=/home/mindaugas/projects/data

## select queries at random, one per superfamily
perl -e 'while(<>){next unless /^>(\S+)\s+([a-z]+\.[0-9]+\.[0-9]+)/;push(@{$H{$2}}, $1);} foreach(keys %H){$n=scalar(@{$H{$_}}); $s=int(rand($n)); printf("%-10s %-10s %10d\n",${$H{$_}}[$s],$_,$n);}' astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa >queries.lst

## prepare structures so that all methods interpret them equivalently: single chain, numeration, no residues w/o " CA [ |A]" atoms
ls -1 ${WRKDIR}/scope-2.08/pdbstyle-2.08.org | xargs -i -P 40 perl -e '
  open(F,">$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08/{}") or die "ERROR: Failed to open file for writing.";
  $prvchn = "_"; $prvstr="    "; $prvnum=0;
  while(<>) {
    if(/^(?:ATOM|HETATM|TER)/) {
      if((substr($_,22,4) ne $prvstr) || (substr($_,21,1) ne $prvchn)) {
        $prvchn = substr($_,21,1);
        $prvstr = substr($_,22,4);
        $prvnum++;
        print(F $rec) if $ca;
        $rec = ""; $ca = 0;
      }
      substr($_,21,1) = "A";
      substr($_,22,4) = sprintf("%4d", $prvnum);
      $ca = 1 if(substr($_,12,4) eq " CA " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $rec .= $_;
    }
    else {
      print(F $rec) if($ca && $rec);
      $rec = ""; $ca = 0;
      print(F);
    }
    last if /^ENDMDL/;
  }
  print(F $rec) if($ca && $rec);
  close(F)' ${WRKDIR}/scope-2.08/pdbstyle-2.08.org/{}

## verification
ls -1 pdbstyle-2.08|xargs -i perl -e 'while(<>){next unless /^(?:ATOM|HETATM)/; if(substr($_,12,4) eq " CA " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A")){ if($rn && $rn eq substr($_,22,5)){print"{}\n$_\n";} $rn = substr($_,22,5); }}' pdbstyle-2.08/{}

## for foldseek, additionally remove HETATM residues to get same sequences by TM-align
mkdir ${WRKDIR}/scope-2.08/pdbstyle-2.08.nohet
ls -1 ${WRKDIR}/scope-2.08/pdbstyle-2.08 | xargs -i -P 40 perl -e 'open(F,">$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08.nohet/{}") or die "ERROR: Failed to open file for writing."; while(<>){print(F) unless /^HETATM/;} close(F)' ${WRKDIR}/scope-2.08/pdbstyle-2.08/{}

## for Dali, leave only residues that have at least N, CA, C, and O atoms
ls -1 ${WRKDIR}/scope-2.08/pdbstyle-2.08 | xargs -i -P 40 perl -e '
  open(F,">$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08.dali/{}") or die "ERROR: Failed to open file for writing.";
  $prvchn = "_"; $prvstr="    ";
  while(<>) {
    if(/^(?:ATOM|HETATM|TER)/) {
      if((substr($_,22,4) ne $prvstr) || (substr($_,21,1) ne $prvchn)) {
        $prvchn = substr($_,21,1);
        $prvstr = substr($_,22,4);
        print(F $rec) if($n && $ca && $c && $o);
        $rec = ""; $n = $ca = $c = $o = 0;
      }
      $ca = 1 if(substr($_,12,4) eq " CA " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $n = 1  if(substr($_,12,4) eq " N  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $c = 1  if(substr($_,12,4) eq " C  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $o = 1  if(substr($_,12,4) eq " O  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $rec .= $_;
    }
    else {
      print(F $rec) if($n && $ca && $c && $o && $rec);
      $rec = ""; $n = $ca = $c = $o = 0;
      print(F);
    }
    last if /^ENDMDL/;
  }
  print(F $rec) if($n && $ca && $c && $o && $rec);
  close(F)' ${WRKDIR}/scope-2.08/pdbstyle-2.08/{}

## make list of files with extension
(cat ${WRKDIR}/scope-2.08/pdbstyle-2.08.lst|xargs -i echo {}.ent) >${WRKDIR}/scope-2.08/pdbstyle-2.08.ext.lst

## FATCAT
## env var
export FATCAT=/home/mindaugas/install/FATCAT/FATCAT-dist
## run FATCAT (does not read queries with absolute path: use cd)
bash -c '(cd ${WRKDIR}; time (grep -vE "^#" scope-2.08/queries.lst | awk "{print \$1}" | xargs -i sh -c "echo {}; time /home/mindaugas/install/FATCAT/FATCAT-dist/FATCATMain/FATCATSearch scope-2.08/pdbstyle-2.08/{}.ent \${WRKDIR}/scope-2.08/pdbstyle-2.08.lst -i2 \${WRKDIR}/scope-2.08/pdbstyle-2.08.pdb -r -o \${WRKDIR}/gtalign-benchmark/fatcat/fatcat_{}.out -m")) >fatcat.log 2>&1 &'
## check #alignments
ls -1 fatcat/fatcat_*|xargs -i sh -c 'echo {} $(grep -E Align {}|wc -l)'|awk '{print $2}'|uniq -c
## sort
grep -E 'P\-value' fatcat_d2fp1a_.out |sort -t '-' -n -k2 -k3|less

## DeepAlign_Search does not produce alignments 
## time ~/install/DeepAlign/DeepAlign/DeepAlign_Search.sh -q pdbstyle-2.08/d2fp1a_.ent -l pdbstyle-2.08.lst -d pdbstyle-2.08.pdb -c 40 -O deepalign_d2fp1a_.out -H ~/install/DeepAlign/DeepAlign
## run DeepAlign for alignments
bash -c '(QRS=($(grep -vE "^#" ${WRKDIR}/scope-2.08/queries.lst | awk "{print \$1}")); (for((i=0; i<${#QRS[@]}; i++)); do query=${QRS[$i]}; echo -e "\n\n${query}"; time (ls -1 ${WRKDIR}/scope-2.08/pdbstyle-2.08 | xargs -i -P 40 sh -c "(/home/mindaugas/install/DeepAlign/DeepAlign/DeepAlign ${WRKDIR}/scope-2.08/pdbstyle-2.08/${query}.ent ${WRKDIR}/scope-2.08/pdbstyle-2.08/{} || true)") >deepalign/deepalign_${query}.out; done)) >deepalign.log 2>&1 &'
## check #alignments
ls -1 deepalign/deepalign_*|xargs -i sh -c 'echo {} $(grep -E " DeepAlign" {}|wc -l)'|awk '{print $2}'|uniq -c
## total time (40 threads)
grep -E '^real' deepalign.log |perl -e 'while(<>){@a=split(/\s+/);$s+=$1*60+$2 if($a[1]=~/^([^m\s]+)m([^s\s]+)s/);} print($s,"\n")'

## DALI v5
## reformat structure files for Dali format 
## (NOTE: overwriting when a filename includes a multiple-character chain id)
## (hence, use manual ids instead):
(cd ${WRKDIR}; tmpdir="wrk_dir_dali"; [ ! -d ${tmpdir} ] && mkdir ${tmpdir}; cd ${tmpdir}; datdir="../gtalign-benchmark/DALIv5_DAT_db"; [ ! -d ${datdir} ] && mkdir ${datdir}; flist=($(ls -1 ../scope-2.08/pdbstyle-2.08.dali)); for((i=0; i<${#flist[@]}; i++)); do f=${flist[$i]}; printf -v id "%04x" $i; echo "$f $id" >>../gtalign-benchmark/DALIv5_DAT_db_map.txt; time /home/mindaugas/install/DaliLite_v5/DaliLite.v5/bin/import.pl --pdbfile ../scope-2.08/pdbstyle-2.08.dali/$f --dat ${datdir} --pdbid $id; done; cd ..) >dali5_import.log 2>&1
## make list:
(ls -1 DALIv5_DAT_db|xargs -i basename {} .dat) >DALIv5_DAT_db.lst
## reformat query structures:
(cd ${WRKDIR}; tmpdir="wrk_dir_dali"; [ ! -d ${tmpdir} ] && mkdir ${tmpdir}; cd ${tmpdir}; datdir="../gtalign-benchmark/DALIv5_DAT_qrs"; [ ! -d ${datdir} ] && mkdir ${datdir}; flist=($(grep -vE "^#" ${WRKDIR}/scope-2.08/queries.lst | awk "{print \$1}")); for((i=0; i<${#flist[@]}; i++)); do f="${flist[$i]}.ent"; printf -v id "%04x" $i; echo "$f $id" >>../gtalign-benchmark/DALIv5_DAT_qrs_map.txt; time /home/mindaugas/install/DaliLite_v5/DaliLite.v5/bin/import.pl --pdbfile ../scope-2.08/pdbstyle-2.08.dali/$f --dat ${datdir} --pdbid $id; done; cd ..) >dali5_import_qrs.log 2>&1
## queries that failed (no .dat file)
ls -1 DALIv5_DAT_qrs/|cut -b1-4|perl -e 'while(<>){chomp; $nxt=hex($prv)+1; printf("%04x\n",$nxt) if($prv && hex($_)-$nxt); $prv=$_}' | xargs -i grep -E ' {}' DALIv5_DAT_qrs_map.txt
## make list:
(ls -1 DALIv5_DAT_qrs|xargs -i basename {} .dat) >DALIv5_DAT_qrs.lst
## run DaliLitev5
bash -c '(cd dali; time /home/mindaugas/install/DaliLite_v5/DaliLite.v5/bin/dali.pl --np 40 --query ../DALIv5_DAT_qrs.lst --db ../DALIv5_DAT_db.lst --dat1 ../DALIv5_DAT_qrs --dat2 ../DALIv5_DAT_db --outfmt alignments; cd ..) >dali.log 2>&1 &'

## TM-align
## run TM-align
bash -c '(cd ${WRKDIR}; time (grep -vE "^#" scope-2.08/queries.lst | awk "{print \$1}" | xargs -i -P 40 sh -c "echo {} \$(time /home/mindaugas/install/TM-align/TMalign scope-2.08/pdbstyle-2.08/{}.ent -dir2 \${WRKDIR}/scope-2.08/pdbstyle-2.08/ \${WRKDIR}/scope-2.08/pdbstyle-2.08.ext.lst -het 1 >\${WRKDIR}/gtalign-benchmark/tmalign/tmalign_{}.out)")) >tmalign.log 2>&1 &'
## check #alignments
ls -1 tmalign/tmalign_*|xargs -i sh -c 'echo {} $(grep -E Aligned {}|wc -l)'|awk '{print $2}'|uniq -c
## run TM-align fast
bash -c '(cd ${WRKDIR}; time (grep -vE "^#" scope-2.08/queries.lst | awk "{print \$1}" | xargs -i -P 40 sh -c "echo {} \$(time /home/mindaugas/install/TM-align/TMalign scope-2.08/pdbstyle-2.08/{}.ent -dir2 \${WRKDIR}/scope-2.08/pdbstyle-2.08/ \${WRKDIR}/scope-2.08/pdbstyle-2.08.ext.lst -het 1 -fast >\${WRKDIR}/gtalign-benchmark/tmalign_fast/tmalign_{}.out)")) >tmalign_fast.log 2>&1 &'

## foldseek
## collect queries
mkdir ${WRKDIR}/scope-2.08/pdbstyle-2.08--queries
mkdir ${WRKDIR}/scope-2.08/pdbstyle-2.08--queries.nohet
grep -vE '^#' ${WRKDIR}/scope-2.08/queries.lst | awk '{print $1}' | xargs -i ln -s ../pdbstyle-2.08/{}.ent ${WRKDIR}/scope-2.08/pdbstyle-2.08--queries/
grep -vE '^#' ${WRKDIR}/scope-2.08/queries.lst | awk '{print $1}' | xargs -i ln -s ../pdbstyle-2.08.nohet/{}.ent ${WRKDIR}/scope-2.08/pdbstyle-2.08--queries.nohet/
## run foldseek using --alignment-type 1
(time ~/install/foldseek/foldseek/bin/foldseek easy-search ${WRKDIR}/scope-2.08/pdbstyle-2.08--queries.nohet/ ${WRKDIR}/scope-2.08/pdbstyle-2.08.nohet/ foldseek_tmfast1__scope20840_bmk.res tmp_dir_foldseek --threads 40 --max-seqs 4000 --alignment-type 1 --format-output query,target,evalue,qtmscore,ttmscore,alntmscore,rmsd,qstart,qend,qlen,tstart,tend,tlen,qaln,taln) >foldseek_tmfast1__scope20840_bmk.log 2>&1
## run foldseek using --alignment-type 2 which uses 3d alphabet:
(time ~/install/foldseek/foldseek/bin/foldseek easy-search ${WRKDIR}/scope-2.08/pdbstyle-2.08--queries.nohet/ ${WRKDIR}/scope-2.08/pdbstyle-2.08.nohet/ foldseek_alntyp2__scope20840_bmk.res tmp_dir_foldseek --threads 40 --max-seqs 4000 --alignment-type 2 --format-output query,target,evalue,qtmscore,ttmscore,alntmscore,rmsd,qstart,qend,qlen,tstart,tend,tlen,qaln,taln) >foldseek_alntyp2__scope20840_bmk.log 2>&1
## run foldseek using --alignment-type 1 and --tmalign-fast 0
bash -c '(time ~/install/foldseek/foldseek/bin/foldseek easy-search ${WRKDIR}/scope-2.08/pdbstyle-2.08--queries.nohet/ ${WRKDIR}/scope-2.08/pdbstyle-2.08.nohet/ foldseek_tmfast0__scope20840_bmk.res tmp_dir_foldseek --threads 40 --max-seqs 4000 --alignment-type 1 --format-output query,target,evalue,qtmscore,ttmscore,alntmscore,rmsd,qstart,qend,qlen,tstart,tend,tlen,qaln,taln --tmalign-fast 0) >foldseek_tmfast0__scope20840_bmk.log 2>&1 &'

## GTalign speed 0
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/scope-2.08/pdbstyle-2.08--queries --rfs=${WRKDIR}/scope-2.08/pdbstyle-2.08 -o gtalign_14_speed0_prescore03_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=0 --pre-score=0.3 -s 0.44 --add-search-by-ss --nhits=10000 --nalns=10000 --dev-N=3) 2>&1 | tee gtalign_14_speed0_prescore03_addss_s044.log
## GTalign speed 9
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/scope-2.08/pdbstyle-2.08--queries --rfs=${WRKDIR}/scope-2.08/pdbstyle-2.08 -o gtalign_14_speed9_prescore03_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=9 --pre-score=0.3 -s 0.44 --add-search-by-ss --nhits=10000 --nalns=10000 --dev-N=3) 2>&1 | tee gtalign_14_speed9_prescore03_addss_s044.log
## GTalign speed 9 + --pre-score=0.4
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/scope-2.08/pdbstyle-2.08--queries --rfs=${WRKDIR}/scope-2.08/pdbstyle-2.08 -o gtalign_14_speed9_prescore04_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=9 --pre-score=0.4 -s 0.44 --add-search-by-ss --nhits=10000 --nalns=10000 --dev-N=3) 2>&1 | tee gtalign_14_speed9_prescore04_addss_s044.log
## GTalign speed 13
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/scope-2.08/pdbstyle-2.08--queries --rfs=${WRKDIR}/scope-2.08/pdbstyle-2.08 -o gtalign_14_speed13_prescore03_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=13 --pre-score=0.3 -s 0.44 --add-search-by-ss --nhits=10000 --nalns=10000 --dev-N=3) 2>&1 | tee gtalign_14_speed13_prescore03_addss_s044.log
## GTalign speed 13 + --pre-score=0.4
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/scope-2.08/pdbstyle-2.08--queries --rfs=${WRKDIR}/scope-2.08/pdbstyle-2.08 -o gtalign_14_speed13_prescore04_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=13 --pre-score=0.4 -s 0.44 --add-search-by-ss --nhits=10000 --nalns=10000 --dev-N=3) 2>&1 | tee gtalign_14_speed13_prescore04_addss_s044.log
## GTalign speed 13 + --pre-similarity=15
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/scope-2.08/pdbstyle-2.08--queries --rfs=${WRKDIR}/scope-2.08/pdbstyle-2.08 -o gtalign_14_speed13_prescore03_presim15_addss_s03 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=13 --pre-score=0.3 --pre-similarity=15 -s 0.3 --add-search-by-ss --nhits=10000 --nalns=10000 --dev-N=3 -c cachedir) 2>&1 | tee gtalign_14_speed13_prescore03_presim15_addss_s03.log


