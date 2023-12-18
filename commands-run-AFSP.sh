export WRKDIR=/home/mindaugas/projects/data



## cluster all (CRISPR-)Cas (searched by keyword: can include non-Cas proteins) proteins using a TM-score threshold of 0.4 and coverage of 0.4
/home/mindaugas/projects/share/gtalign/bin/gtalign -v --cls=../wwpdb-CRISPR-Cas -o . --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=1500 --speed=13 --pre-score=0.3 --cpu-threads-reading=10 --add-search-by-ss --cls-coverage=0.4 --cls-threshold=0.4 --dev-N=3 --ter=0 --split=2
## print the length of each protein that constitutes a singleton cluster
awk '{if(NF<2)print}' gtalignclusters.lst | xargs -i perl -e '$pfx=substr("{}",0,4); $C=substr("{}",5,1); open(F,"1k04.pdb__0.out")||die"ERROR: open.";while(<F>){printf("%-s %5d\n","{}",$1) if /^\s+\d+\s+\.\.\S+$pfx\S+\s+Chn:$C\s+\S+\s+\S+\s+\S+\s+\d+\s+\S+\s+\S+\s+(\d+)/} close(F)' >singletons.lst
## select queries evenly distributed in length
## (NOTE: foldseek ignores 6w1s_A)
sort -rnk2 singletons.lst |perl -e 'while(<>){print if $c++%3==0}'|head -n40 >queries.txt
## get structure  chains
mkdir queries.org
awk '{print $1}' queries.txt | xargs -i sh -c 'n=$(echo {}|sed -re "s/_.+\$//"); c=$(echo {}|sed -re "s/^.+_//"); gunzip -c ../wwpdb-CRISPR-Cas/$n.pdb.gz >$n.pdb; python3 /data/installed-software/comer-ws-backend/bin/getchain.py -i $n.pdb -c $c -o queries.org/${n}_${c}.pdb; rm $n.pdb'
## prepare structures so that all methods interpret them equivalently: no residues w/o " CA [ |A]" atoms;
## leave only residues that have at least N, CA, C, and O atoms (valid for Dali);
## also, remove HETATM residues which makes foldseek and TM-align read the same sequence of residues
mkdir queries 
ls -1 queries.org | xargs -i -P 40 perl -e '
  open(F,">queries/{}") or die "ERROR: Failed to open file for writing.";
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
      if(substr($_,12,4) eq " CA " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A")) {
        $ca = 1;
      }
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
  close(F)' queries.org/{}



## gunzip AF structures for TM-align
ls -1 ${WRKDIR}/afsp_db/swissprot_pdb_v4 | xargs -i -P 40 gunzip ${WRKDIR}/afsp_db/swissprot_pdb_v4/{}
## make list of files with extension for TM-align
ls -1 ${WRKDIR}/afsp_db/swissprot_pdb_v4 >${WRKDIR}/afsp_db/swissprot_pdb_v4.ext.lst



## TM-align
## run TM-align
bash -c '(cd ${WRKDIR}; time (ls -1 gtalign-benchmark-afspdb/queries | xargs -i -P 40 sh -c "echo {} \$(time /home/mindaugas/install/TM-align/TMalign gtalign-benchmark-afspdb/queries/{} -dir2 \${WRKDIR}/afsp_db/swissprot_pdb_v4/ \${WRKDIR}/afsp_db/swissprot_pdb_v4.ext.lst -het 1 -outfmt 2 >\${WRKDIR}/gtalign-benchmark-afspdb/tmalign/tmalign_{}.out)")) >tmalign_afspdb.log 2>&1 &'
## check #alignments
ls -1 tmalign/tmalign_*|xargs -i wc -l {}|awk '{print $1}'|uniq -c
## run TM-align fast
bash -c '(cd ${WRKDIR}; time (ls -1 gtalign-benchmark-afspdb/queries | xargs -i -P 40 sh -c "echo {} \$(time /home/mindaugas/install/TM-align/TMalign gtalign-benchmark-afspdb/queries/{} -dir2 \${WRKDIR}/afsp_db/swissprot_pdb_v4/ \${WRKDIR}/afsp_db/swissprot_pdb_v4.ext.lst -het 1 -outfmt 2 -fast >\${WRKDIR}/gtalign-benchmark-afspdb/tmalign_fast/tmalign_{}.out)")) >tmalign_fast_afspdb.log 2>&1 &'



## foldseek
## run foldseek using --alignment-type 2 which uses 3d alphabet:
(time ~/install/foldseek/foldseek/bin/foldseek easy-search ${WRKDIR}/gtalign-benchmark-afspdb/queries/ ${WRKDIR}/afsp_db/archive/swissprot_pdb_v4.tar foldseek_alntyp2__afspdb_bmk.res tmp_dir_foldseek --threads 40 --max-seqs 20000 --alignment-type 2 --format-output query,target,evalue,qtmscore,ttmscore,alntmscore,rmsd,qstart,qend,qlen,tstart,tend,tlen,qaln,taln) >foldseek_alntyp2__afspdb_bmk.log 2>&1
## run foldseek using --alignment-type 1
bash -c '(time ~/install/foldseek/foldseek/bin/foldseek easy-search ${WRKDIR}/gtalign-benchmark-afspdb/queries/ ${WRKDIR}/afsp_db/archive/swissprot_pdb_v4.tar foldseek_tmfast1__afspdb_bmk.res tmp_dir_foldseek --threads 40 --max-seqs 20000 --alignment-type 1 --format-output query,target,evalue,qtmscore,ttmscore,alntmscore,rmsd,qstart,qend,qlen,tstart,tend,tlen,qaln,taln) >foldseek_tmfast1__afspdb_bmk.log 2>&1 &'
## run foldseek using --alignment-type 1 and --tmalign-fast 0
bash -c '(time ~/install/foldseek/foldseek/bin/foldseek easy-search ${WRKDIR}/gtalign-benchmark-afspdb/queries/ ${WRKDIR}/afsp_db/archive/swissprot_pdb_v4.tar foldseek_tmfast0__afspdb_bmk.res tmp_dir_foldseek --threads 40 --max-seqs 20000 --alignment-type 1 --format-output query,target,evalue,qtmscore,ttmscore,alntmscore,rmsd,qstart,qend,qlen,tstart,tend,tlen,qaln,taln --tmalign-fast 0) >foldseek_tmfast0__afspdb_bmk.log 2>&1 &'



## GTalign speed 0
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/gtalign-benchmark-afspdb/queries --rfs=${WRKDIR}/afsp_db/archive/swissprot_pdb_v4.tar -o gtalign_14_speed0_prescore03_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=0 --pre-score=0.3 -s 0.44 --add-search-by-ss --nhits=50000 --nalns=50000 --dev-N=3 -c cachedir) 2>&1 | tee gtalign_14_speed0_prescore03_addss_s044.log
rm cachedir/*
## GTalign speed 9
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/gtalign-benchmark-afspdb/queries --rfs=${WRKDIR}/afsp_db/archive/swissprot_pdb_v4.tar -o gtalign_14_speed9_prescore03_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=9 --pre-score=0.3 -s 0.44 --add-search-by-ss --nhits=50000 --nalns=50000 --dev-N=3 -c cachedir) 2>&1 | tee gtalign_14_speed9_prescore03_addss_s044.log
rm cachedir/*
## GTalign speed 9 + --pre-score=0.4
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/gtalign-benchmark-afspdb/queries --rfs=${WRKDIR}/afsp_db/archive/swissprot_pdb_v4.tar -o gtalign_14_speed9_prescore04_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=9 --pre-score=0.4 -s 0.44 --add-search-by-ss --nhits=50000 --nalns=50000 --dev-N=3 -c cachedir) 2>&1 | tee gtalign_14_speed9_prescore04_addss_s044.log
rm cachedir/*
## GTalign speed 13
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/gtalign-benchmark-afspdb/queries --rfs=${WRKDIR}/afsp_db/archive/swissprot_pdb_v4.tar -o gtalign_14_speed13_prescore03_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=13 --pre-score=0.3 -s 0.44 --add-search-by-ss --nalns=50000 --nhits=50000 --dev-N=3 -c cachedir) 2>&1 | tee gtalign_14_speed13_prescore03_addss_s044.log
rm cachedir/*
## GTalign speed 13 + --pre-score=0.4
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/gtalign-benchmark-afspdb/queries --rfs=${WRKDIR}/afsp_db/archive/swissprot_pdb_v4.tar -o gtalign_14_speed13_prescore04_addss_s044 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=13 --pre-score=0.4 -s 0.44 --add-search-by-ss --nhits=50000 --nalns=50000 --dev-N=3 -c cachedir) 2>&1 | tee gtalign_14_speed13_prescore04_addss_s044.log
rm cachedir/*
## GTalign speed 13 + --pre-similarity=15
(time /home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=${WRKDIR}/gtalign-benchmark-afspdb/queries --rfs=${WRKDIR}/afsp_db/archive/swissprot_pdb_v4.tar -o gtalign_14_speed13_prescore03_presim15_addss_s03 --hetatm --dev-queries-per-chunk=2 --dev-queries-total-length-per-chunk=1500 --dev-min-length=3 --dev-max-length=4000 --speed=13 --pre-score=0.3 --pre-similarity=15 -s 0.3 --add-search-by-ss --nhits=20000 --nalns=20000 --dev-N=3 -c cachedir) 2>&1 | tee gtalign_14_speed13_prescore03_presim15_addss_s03.log


