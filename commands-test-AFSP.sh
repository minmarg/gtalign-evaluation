export WRKDIR=/home/mindaugas/projects/data

##collect TM-align results
export MYPROCDIR=tmalign
export MYPROCDIR=tmalign_fast
(ls -1 ${MYPROCDIR} | grep -E .out | xargs -i -P 40 perl -e '
  $RDIR="$ENV{MYPROCDIR}";
  $QIDfle="{}"; $QID=$QIDfle; $QID=~s/^tmalign_(.+).out$/$1/;
  exit 0 unless $QID;
  open(F, "$RDIR/$QIDfle") || die "ERROR: Failed to open results for: $QIDfle";
  while(<F>) {
    chomp;
    if(/^(\S+)\s+(\S+)\s+([\d\.]+)\s+([\d\.]+)\s+/) {
      $best = ($3<$4)? $4: $3;
      next if $best < 0.4;
      $str = sprintf("%-30s %30s tm1= %8s tm2= %8s best= %8.6f\n", $1, $2, $3, $4, $best);
      syswrite(STDOUT, $str);
    }
  }
  close(F);') >${MYPROCDIR}__afspdb__processed.out 2>${MYPROCDIR}__afspdb__processed.err

##process foldseek results to obtain TM-scores for the alignments
export MYPROCDIR=foldseek_tmfast0
export MYPROCDIR=foldseek_tmfast1
export MYPROCDIR=foldseek_alntyp2
(cat ${MYPROCDIR}__afspdb_bmk.res | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
  $QRYDIR="$ENV{WRKDIR}/gtalign-benchmark-afspdb/queries";
  $RFNDIR="$ENV{WRKDIR}/afsp_db/swissprot_pdb_v4";
  $WDIR="tmp_wrk_dir";
  $_ = "{}";
  if($_) {
    chomp;
    if(/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/) {
      $qid = $1; $rfn = $2; $ftme = $3; $ftm1 = $4; $ftm2 = $5;
      $qstart = $8; $qend = $9; $qlen = $10;  $tstart = $11; $tend = $12; $tlen = $13;
      $qaln = $14; $taln = $15;
      $rfn =~ s/\.gz$//;
      if($qid && $rfn && $qaln && $taln && $qid ne $rfn) {
        $wfile = "$WDIR/${qid}_${rfn}.aln";
        $tmfile = "$WDIR/${qid}_${rfn}.tm";
        $qnnb = $qstart - 1; $qnne = $qlen - $qend;
        $tnnb = $tstart - 1; $tnne = $tlen - $tend;
        $qaln = "X"x$qnnb . $qaln . "X"x$qnne;
        $taln = "-"x$qnnb . $taln . "-"x$qnne;
        $qaln = "-"x$tnnb . $qaln . "-"x$tnne;
        $taln = "X"x$tnnb . $taln . "X"x$tnne;
        open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
        print(W ">$qid\n$qaln\n>$rfn\n$taln\n");
        close(W);
        $prog = "$TMalign $QRYDIR/$qid $RFNDIR/$rfn -I $wfile -het 0 >$tmfile";
        if(system($prog) != 0) {
          print(STDERR "ERROR: Failed to execute: $prog\n");
        } else {
          open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
          while(<W>) {
            $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
            $tm2 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/;
            last if /^\(":" denotes residue pairs/;
          }
          $tmqaln = <W>; $tmtaln = <W>; $tmtaln = <W>; chomp($tmqaln); chomp($tmtaln);
          close(W);
          $strprn = sprintf("%-12s %12s ftme= %10s ftm1= %10s ftm2= %10s  tm1= %8s tm2= %8s best= %8.6f\n",
            $qid,$rfn,$ftme,$ftm1,$ftm2, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1);
          syswrite(STDOUT, $strprn);
          if(length($tmqaln) == length($qaln) && length($tmtaln) == length($taln)) {
            unlink($wfile, $tmfile);
          } else {
            print(STDERR "WARNING: Inconsistent alignments: $qid $rfn\n");
          }
        }
      }
    }
  }
') >${MYPROCDIR}__afspdb__processed.out 2>${MYPROCDIR}__afspdb__processed.err

##process GTalign results to obtain TM-align TM-scores for the alignments
export MYPROCDIR=gtalign_14_speed0_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed9_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed9_prescore04_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore04_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore03_presim15_addss_s03
(ls -1 ${MYPROCDIR} | grep -E .out | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
  $QRYDIR="$ENV{WRKDIR}/gtalign-benchmark-afspdb/queries";
  $RFNDIR="$ENV{WRKDIR}/afsp_db/swissprot_pdb_v4";
  $RDIR="$ENV{MYPROCDIR}";
  $WDIR="tmp_wrk_dir";
  $QIDfle="{}"; $QID=$QIDfle; $QID=~s/__\d+.out$//;
  exit 0 unless $QID;
  open(F, "$RDIR/$QIDfle") || die "ERROR: Failed to open results for: $QIDfle";
  while(<F>) {
    chomp;
    if(eof(F) || /^(\d+)\.\s*$/) {
      $nrtt = "";
      $nrtt = $1 unless eof(F);
      if($rfn && $qaln && $raln) {
        $wfile = "$WDIR/${QID}_${rfn}.aln";
        $tmfile = "$WDIR/${QID}_${rfn}.tm";
        $qnnb = $qstart - 1; $qnne = $qlen - $qend;
        $rnnb = $rstart - 1; $rnne = $rlen - $rend;
        $qaln = "X"x$qnnb . $qaln . "X"x$qnne;
        $raln = "-"x$qnnb . $raln . "-"x$qnne;
        $qaln = "-"x$rnnb . $qaln . "-"x$rnne;
        $raln = "X"x$rnnb . $raln . "X"x$rnne;
        open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
        print(W ">$QID\n$qaln\n>$rfn\n$raln\n");
        close(W);
        $prog = "$TMalign $QRYDIR/${QID}.pdb $RFNDIR/$rfn -I $wfile -het 1 >$tmfile";
        if(system($prog) != 0) {
          print(STDERR "ERROR: Failed to execute: $prog\n");
        } else {
          open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
          while(<W>) {
            $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
            $tm2 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/;
            last if /^\(":" denotes residue pairs/;
          }
          $tmqaln = <W>; $tmraln = <W>; $tmraln = <W>; chomp($tmqaln); chomp($tmraln);
          close(W);
          $strprn = sprintf("%-5s %12s %12s gtm1= %8s gtm2= %8s gbest= %8.6f  tm1= %8s tm2= %8s best= %8.6f\n",
            $nr,"$QID.pdb",$rfn,$gtm1,$gtm2,($gtm1<$gtm2)?$gtm2:$gtm1, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1);
          syswrite(STDOUT, $strprn);
          if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
            unlink($wfile, $tmfile);
          } else {
            print(STDERR "WARNING: Inconsistent alignments: $QID $rfn\n");
          }
        }
      }
      $nr = $nrtt; $gtm1 = $gtm2 = -1; $qaln = $raln = $rfn = ""; $tm1 = $tm2 = -1;
      $qlen = $rlen = $qend = $rend = -1; undef($qstart); undef($rstart);
    }
    do {$rfn = $1; $rfn =~ s/^\S+\:(\S+)\.gz$/$1/;} if /^>(\S+)/;
    do {$rlen = $1; $qlen = $2;} if /^\s+Length:\s+Refn.\s+=\s*(\d+),\s+Query\s+=\s*(\d+)/;
    do {$gtm2 = $1; $gtm1 = $2;} if /^\s+TM-score\s+\(Refn.\/Query\)\s+=\s+([\d\.]+)\s+\/\s+([\d\.]+)/;
    do {$qaln .= $2; $qend = $3; $qstart = $1 unless $qstart;} if /^Query:\s+(\d+)\s+(\S+)\s+(\d+)/;
    do {$raln .= $2; $rend = $3; $rstart = $1 unless $rstart;} if /^Refn.:\s+(\d+)\s+(\S+)\s+(\d+)/;
  }
  close(F);') >${MYPROCDIR}__afspdb__processed.out 2>${MYPROCDIR}__afspdb__processed.err


