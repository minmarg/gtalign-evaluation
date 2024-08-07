export WRKDIR=/home/mindaugas/projects/data

##process FATCAT results to obtain GDT_TS scores for the alignments
(ls -1 fatcat | grep -E .out | xargs -i -P 40 perl -e '
  $TMscore="/home/mindaugas/install/TM-score--2022-2-27/TMscore_from_alignment";
  $QRYDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08--queries";
  $RFNDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08.pdb";
  $RDIR="fatcat";
  $WDIR="tmp_wrk_dir";
  $QIDfle="{}"; $QID=$QIDfle; $QID=~s/^fatcat_(.+).out$/$1/;
  exit 0 unless $QID;
  open(Q, "$QRYDIR/$QID.ent") || die "ERROR: Failed to open query: $QRYDIR/$QID.ent";
  open(F, "$RDIR/$QIDfle") || die "ERROR: Failed to open results for: $QIDfle";
  while(<F>) {
    chomp;
    if(eof(F) || /^Align\s+(\S+)\s+(\d+)\s+with\s+(\S+)\s+(\d+)/) {
      $qidtt = $rfntt = ""; $qlentt = $rlentt = -1;
      unless(eof(F)) {
        print(STDERR "WARNING: Query $QID.ent doesnt match $1.\n") unless "$QID.ent" eq $1;
        $qidtt = $1; $qlentt = $2; $rfntt = $3; $rlentt = $4;
      }
      if($qid && $rfn && $qaln && $raln && substr($qid,0,7) ne substr($rfn,0,7)) {
        $qidStart = $rfnStart = 0;
        seek(Q,0,0);
        while(<Q>) {
          next unless /^ATOM\s+\d+\s+CA [ A]/;
          $resnumc = substr($_,22,5); $resnumc =~ s/\s//g;
          last if $resnumc eq $qstart;
          $qidStart++;
        }
        if(open(R, "$RFNDIR/$rfn")) {
          while(<R>) {
            next unless /^ATOM\s+\d+\s+CA [ A]/;
            $resnumc = substr($_,22,5); $resnumc =~ s/\s//g;
            last if $resnumc eq $rstart;
            $rfnStart++;
          }
          close(R);
        } else {
          print(STDERR "ERROR: Failed to open refn: $RFNDIR/$rfn\n");
        }
        $wfile = "$WDIR/${qid}_${rfn}.aln";
        $tmfile = "$WDIR/${qid}_${rfn}.tm";
        ## FATCAT uses assigned resnum
        $qend = $qidStart + ($qaln =~ tr/a-zA-Z//);
        $rend = $rfnStart + ($raln =~ tr/a-zA-Z//);
        if(length($qaln) > 2 && length($raln) > 2) {
          $qnnb = $qidStart; $qnne = $qlen - $qend;
          $rnnb = $rfnStart; $rnne = $rlen - $rend;
          $qaln = "X"x$qnnb . $qaln . "X"x$qnne;
          $raln = "-"x$qnnb . $raln . "-"x$qnne;
          $qaln = "-"x$rnnb . $qaln . "-"x$rnne;
          $raln = "X"x$rnnb . $raln . "X"x$rnne;
          open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
          print(W ">$qid\n$qaln\n>$rfn\n$raln\n");
          close(W);
          ## FATCAT ignores HETATMs:
          $prog = "$TMscore $QRYDIR/$qid $RFNDIR/$rfn -I $wfile >$tmfile";
          if(system($prog) != 0) {
            print(STDERR "ERROR: Failed to execute: $prog\n");
          } else {
            open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
            while(<W>) {
              $gdt = $1 if /^GDT-TS-score=\s+(\S+)\s+/;
              last if /^\(":" denotes the residue pairs/;
            }
            $tmqaln = <W>; $tmraln = <W>; $tmraln = <W>; chomp($tmqaln); chomp($tmraln);
            close(W);
            $str = sprintf("%-12s %12s Fsc= %8s Pvl= %8s  gdtts= %6s\n",
              $qid,$rfn,$Fsc,$pvl, $gdt);
            syswrite(STDOUT, $str);
            if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
              unlink($wfile, $tmfile);
            } else {
              print(STDERR "WARNING: Inconsistent alignments: $qid $rfn\n");
            }
          }
        }
      }
      $qid = $qidtt; $rfn = $rfntt; $Fsc = $pvl = -1; $qaln = $raln = ""; $gdt = -1;
      $qlen = $qlentt; $rlen = $rlentt; $qend = $rend = 0; undef($qstart); undef($rstart);
    }
    $Fsc = $1 if /^Twists\s+\S+\s+ini-len\s+\S+\s+ini-rmsd\s+\S+\s+opt-equ\s+\S+\s+opt-rmsd\s+\S+\s+chain-rmsd\s+\S+\s+Score\s+(\S+)\s+/;
    $pvl = $1 if /^P-value\s+(\S+)\s+/;
    do {$qaln .= $2; $qstart = $1 unless $qstart;} if /^Chain\s+1:\s+(\S+)\s+(\S+)/;
    do {$raln .= $2; $rstart = $1 unless $rstart;} if /^Chain\s+2:\s+(\S+)\s+(\S+)/;
  }
  close(F);close(Q);') >fatcat__scope20840__processed_gdtts.out 2>fatcat__scope20840__processed_gdtts.err


##process DeepAlign results to obtain TM-scores for the alignments
(ls -1 deepalign | grep -E .out | xargs -i -P 40 perl -e '
  $TMscore="/home/mindaugas/install/TM-score--2022-2-27/TMscore_from_alignment";
  $QRYDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08--queries";
  $RFNDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08";
  $RDIR="deepalign";
  $WDIR="tmp_wrk_dir";
  $QIDfle="{}"; $QID=$QIDfle; $QID=~s/^deepalign_(.+).out$/$1/;
  exit 0 unless $QID;
  open(F, "$RDIR/$QIDfle") || die "ERROR: Failed to open results for: $QIDfle";
  while(<F>) {
    chomp;
    if(eof(F) || /^1st\s+input\s+protein:\s+(\S+)\s+length=\s+(\d+)/) {
      $qidtt = ""; $qlentt = -1;
      unless(eof(F)) {
        print(STDERR "WARNING: Query $QID doesnt match $1.\n") unless $QID eq $1;
        $qidtt = "$1.ent"; $qlentt = $2;
      }
      if($qid && $rfn && $qaln && $raln && $qid ne $rfn) {
        $wfile = "$WDIR/${qid}_${rfn}.aln";
        $tmfile = "$WDIR/${qid}_${rfn}.tm";
        $qnnb = $qstart - 1; $qnne = $qlen - $qend;
        $rnnb = $rstart - 1; $rnne = $rlen - $rend;
        $qaln = "X"x$qnnb . $qaln . "X"x$qnne;
        $raln = "-"x$qnnb . $raln . "-"x$qnne;
        $qaln = "-"x$rnnb . $qaln . "-"x$rnne;
        $raln = "X"x$rnnb . $raln . "X"x$rnne;
        open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
        print(W ">$qid\n$qaln\n>$rfn\n$raln\n");
        close(W);
        $prog = "$TMscore $QRYDIR/$qid $RFNDIR/$rfn -I $wfile -het 1 >$tmfile";
        if(system($prog) != 0) {
          print(STDERR "ERROR: Failed to execute: $prog\n");
        } else {
          open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
          while(<W>) {
            $gdt = $1 if /^GDT-TS-score=\s+(\S+)\s+/;
            last if /^\(":" denotes the residue pairs/;
          }
          $tmqaln = <W>; $tmraln = <W>; $tmraln = <W>; chomp($tmqaln); chomp($tmraln);
          close(W);
          $str = sprintf("%-12s %12s Dsc= %8s dtm= %8s  gdtts= %6s\n",
            $qid,$rfn,$dsc,$dtm, $gdt);
          syswrite(STDOUT, $str);
          if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
            unlink($wfile, $tmfile);
          } else {
            print(STDERR "WARNING: Inconsistent alignments: $qid $rfn\n");
          }
        }
      }
      $qid = $qidtt; $dsc = $dtm = -1; $qaln = $raln = $rfn = ""; $gdt = -1;
      $qlen = $qlentt; $rlen = $qend = $rend = -1; undef($qstart); undef($rstart);
    }
    do {$rfn = "$1.ent"; $rlen = $2;} if /^2nd\s+input\s+protein:\s+(\S+)\s+length=\s+(\d+)/;
    do {$dsc = $1; $dtm = $2;} if /^\s+[\d\-\.]+\s+[\d\-\.]+\s+([\d\-\.]+)\s+[\d\-\.]+\s+[\d\-\.]+\s+[\d\-\.]+\s+([\d\-\.]+)/;
    do {$qaln .= $2; $qend = $3; $qstart = $1 unless $qstart;} if /^T\s+\S+\s+(\d+)\s+(\S+)\s+(\d+)/;
    do {$raln .= $2; $rend = $3; $rstart = $1 unless $rstart;} if /^S\s+\S+\s+(\d+)\s+(\S+)\s+(\d+)/;
  }
  close(F);') >deepalign__scope20840__processed_gdtts.out 2>deepalign__scope20840__processed_gdtts.err


##process DALI results to obtain TM-scores for the alignments
(ls -1 dali | grep -E .txt | xargs -i -P 40 perl -e '
  $TMscore="/home/mindaugas/install/TM-score--2022-2-27/TMscore_from_alignment";
  $QRYDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08.dali";
  $RFNDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08.dali";
  $RDIR="dali";
  $WDIR="tmp_wrk_dir";
  $QIDfle="{}"; $QID=$QIDfle; $QID=~s/.txt$//;
  exit 0 unless $QID;
  open(F, "DALIv5_DAT_qrs_map.txt") || die "ERROR: Failed to open query map.";
  while(<F>) {chomp; @a = split(/\s+/); $hq{$a[1]} = $a[0];}
  close(F);
  open(F, "DALIv5_DAT_db_map.txt") || die "ERROR: Failed to open map of references.";
  while(<F>) {chomp; @a = split(/\s+/); $rf{$a[1]} = $a[0];}
  close(F);
  open(F, "$RDIR/$QIDfle") || die "ERROR: Failed to open results for: $QIDfle";
  die "ERROR: QUERY \"$QID\" not found." unless(exists $hq{substr($QID,0,4)});
  while(<F>) {
    chomp;
    if(eof(F) || /^No\s+(\S+)\s+Query=(\S+)\s+Sbjct=(\S+)\s+Z\-score=(\S+)/) {
      $nrtt = ""; $qidtt = $rfntt = ""; $zsctt = -1;
      do { $nrtt = $1; $qidtt = $2; $rfntt = $3; $zsctt = $4; } unless eof(F);
      if($qid && $rfn && $qseq && $rseq) {
        $wfile = "$WDIR/${qid}_${rfn}.aln";
        $tmfile = "$WDIR/${qid}_${rfn}.tm";
        print(STDERR "WARNING: Query $QID doesnt match $qid.\n") unless $QID eq $qid;
        die "ERROR: Query $qid not found." unless(exists $hq{substr($qid,0,4)});
        die "ERROR: Reference $rfn not found." unless(exists $rf{substr($rfn,0,4)});
        $qstr = $hq{substr($qid,0,4)};
        $rstr = $rf{substr($rfn,0,4)};
        if($qstr ne $rstr) {
          open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
          print(W ">$qstr $qid\n$qseq\n>$rstr $rfn\n$rseq\n");
          close(W);
          $prog = "$TMscore $QRYDIR/$qstr $RFNDIR/$rstr -I $wfile -het 1 >$tmfile";
          if(system($prog) != 0) {
            print(STDERR "ERROR: Failed to execute: $prog\n");
          } else {
            open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
            while(<W>) {
              $gdt = $1 if /^GDT-TS-score=\s+(\S+)\s+/;
              last if /^\(":" denotes the residue pairs/;
            }
            $tmqaln = <W>; $tmraln = <W>; $tmraln = <W>; chomp($tmqaln); chomp($tmraln);
            close(W);
            $strprn = sprintf("No.%-5s %5s %5s %12s %12s Z= %5s  gdtts= %6s\n",
              $nr,$qid,$rfn,$qstr,$rstr,$zsc,$gdt);
            syswrite(STDOUT, $strprn);
            if(length($tmqaln) == length($qseq) && length($tmraln) == length($rseq)) {
              unlink($wfile, $tmfile);
            } else {
              print(STDERR "WARNING: Inconsistent alignments: $qid $rfn: $qstr $rstr\n");
            }
          }
        }
      }
      $nr = $nrtt; $qid = $qidtt; $rfn = $rfntt; $zsc = $zsctt; $qseq = $rseq = ""; $gdt = -1;
    }
    $qseq .= $1 if /^Query\s+(\S+)\s+/;
    $rseq .= $1 if /^Sbjct\s+(\S+)\s+/;
  }
  close(F);') >DALIv5__scope20840__processed_gdtts.out 2>DALIv5__scope20840__processed_gdtts.err


##collect TM-align results
export MYPROCDIR=tmalign
#export MYPROCDIR=tmalign_fast
(ls -1 ${MYPROCDIR} | grep -E .out | xargs -i -P 40 perl -e '
  $TMscore="/home/mindaugas/install/TM-score--2022-2-27/TMscore_from_alignment";
  $QRYDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08--queries";
  $RFNDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08";
  $RDIR="$ENV{MYPROCDIR}";
  $WDIR="tmp_wrk_dir";
  $QIDfle="{}"; $QID=$QIDfle; $QID=~s/^tmalign_(.+).out$/$1/;
  exit 0 unless $QID;
  open(F, "$RDIR/$QIDfle") || die "ERROR: Failed to open results for: $QIDfle";
  while(<F>) {
    chomp;
    if(eof(F) || /^Name\s+of\s+Chain_1:\s+\S+\/(\S+)\s+\(to\s+be/) {
      $qidtt = "";
      unless(eof(F)) {
        print(STDERR "WARNING: Query $QID doesnt match $1.\n") unless "$QID.ent" eq $1;
        $qidtt = $1;
      }
      if($qid && $rfn && $qid ne $rfn) {
        $wfile = "$WDIR/${qid}_${rfn}.aln";
        $tmfile = "$WDIR/${qid}_${rfn}.tm";
        open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
        print(W ">$qid\n$qaln\n>$rfn\n$raln\n");
        close(W);
        $prog = "$TMscore $QRYDIR/$qid $RFNDIR/$rfn -I $wfile -het 1 >$tmfile";
        if(system($prog) != 0) {
          print(STDERR "ERROR: Failed to execute: $prog\n");
        } else {
          open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
          while(<W>) {
            $gdt = $1 if /^GDT-TS-score=\s+(\S+)\s+/;
            last if /^\(":" denotes the residue pairs/;
          }
          $tmqaln = <W>; $tmraln = <W>; $tmraln = <W>; chomp($tmqaln); chomp($tmraln);
          close(W);
          $str = sprintf("%-12s %12s tm1= %8s tm2= %8s best= %8.6f  gdtts= %6s\n",
            $qid,$rfn, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1, $gdt);
          syswrite(STDOUT, $str);
          if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
            unlink($wfile, $tmfile);
          } else {
            print(STDERR "WARNING: Inconsistent alignments: $qid $rfn\n");
          }
        }
      }
      $qid = $qidtt; $qaln = $raln = $rfn = ""; $tm1 = $tm2 = $gdt = -1;
    }
    if(/^\(":" denotes residue pairs/) {
      $qaln = <F>; $raln = <F>; $raln = <F>; chomp($qaln); chomp($raln);
    }
    $rfn = $1 if /^Name\s+of\s+Chain_2:\s+(\S+)/;
    $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
    $tm2 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/;
  }
  close(F);') >${MYPROCDIR}__scope20840__processed_gdtts.out 2>${MYPROCDIR}__scope20840__processed_gdtts.err


##process foldseek results to obtain TM-scores for the alignments
export MYPROCDIR=foldseek_tmfast0
export MYPROCDIR=foldseek_tmfast1
export MYPROCDIR=foldseek_alntyp2
(cat ${MYPROCDIR}__scope20840_bmk.res | xargs -i -P 40 perl -e '
  $TMscore="/home/mindaugas/install/TM-score--2022-2-27/TMscore_from_alignment";
  $QRYDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08--queries";
  $RFNDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08";
  $WDIR="tmp_wrk_dir";
  $_ = "{}";
  if($_) {
    chomp;
    if(/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/) {
      $qid = $1; $rfn = $2; $ftme = $3; $ftm1 = $4; $ftm2 = $5;
      $qstart = $8; $qend = $9; $qlen = $10;  $tstart = $11; $tend = $12; $tlen = $13;
      $qaln = $14; $taln = $15;
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
        $prog = "$TMscore $QRYDIR/$qid $RFNDIR/$rfn -I $wfile -het 0 >$tmfile";
        if(system($prog) != 0) {
          print(STDERR "ERROR: Failed to execute: $prog\n");
        } else {
          open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
          while(<W>) {
            $gdt = $1 if /^GDT-TS-score=\s+(\S+)\s+/;
            last if /^\(":" denotes the residue pairs/;
          }
          $tmqaln = <W>; $tmtaln = <W>; $tmtaln = <W>; chomp($tmqaln); chomp($tmtaln);
          close(W);
          $strprn = sprintf("%-12s %12s ftme= %10s ftm1= %10s ftm2= %10s  gdtts= %6s\n",
            $qid,$rfn,$ftme,$ftm1,$ftm2, $gdt);
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
') >${MYPROCDIR}__scope20840__processed_gdtts.out 2>${MYPROCDIR}__scope20840__processed_gdtts.err


##process GTalign results to obtain TM-align TM-scores for the alignments
export MYPROCDIR=gtalign_14_speed0_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed9_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed9_prescore04_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore04_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore03_presim15_addss_s03
(ls -1 ${MYPROCDIR} | grep -E .out | xargs -i -P 40 perl -e '
  $TMscore="/home/mindaugas/install/TM-score--2022-2-27/TMscore_from_alignment";
  $QRYDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08--queries";
  $RFNDIR="$ENV{WRKDIR}/scope-2.08/pdbstyle-2.08";
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
      if($rfn && $qaln && $raln && "$QID.ent" ne $rfn) {
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
        $prog = "$TMscore $QRYDIR/${QID}.ent $RFNDIR/$rfn -I $wfile -het 1 >$tmfile";
        if(system($prog) != 0) {
          print(STDERR "ERROR: Failed to execute: $prog\n");
        } else {
          open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
          while(<W>) {
            $gdt = $1 if /^GDT-TS-score=\s+(\S+)\s+/;
            last if /^\(":" denotes the residue pairs/;
          }
          $tmqaln = <W>; $tmraln = <W>; $tmraln = <W>; chomp($tmqaln); chomp($tmraln);
          close(W);
          $strprn = sprintf("%-5s %12s %12s gtm1= %8s gtm2= %8s gbest= %8.6f  gdtts= %6s\n",
            $nr,"$QID.ent",$rfn,$gtm1,$gtm2,($gtm1<$gtm2)?$gtm2:$gtm1, $gdt);
          syswrite(STDOUT, $strprn);
          if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
            unlink($wfile, $tmfile);
          } else {
            print(STDERR "WARNING: Inconsistent alignments: $QID $rfn\n");
          }
        }
      }
      $nr = $nrtt; $gtm1 = $gtm2 = -1; $qaln = $raln = $rfn = ""; $gdt = -1;
      $qlen = $rlen = $qend = $rend = -1; undef($qstart); undef($rstart);
    }
    do {$rfn = $1; $rfn =~ s/^\S+\/(\S+\.ent)$/$1/;} if /^>(\S+)/;
    do {$rlen = $1; $qlen = $2;} if /^\s+Length:\s+Refn.\s+=\s*(\d+),\s+Query\s+=\s*(\d+)/;
    do {$gtm2 = $1; $gtm1 = $2;} if /^\s+TM-score\s+\(Refn.\/Query\)\s+=\s+([\d\.]+)\s+\/\s+([\d\.]+)/;
    do {$qaln .= $2; $qend = $3; $qstart = $1 unless $qstart;} if /^Query:\s+(\d+)\s+(\S+)\s+(\d+)/;
    do {$raln .= $2; $rend = $3; $rstart = $1 unless $rstart;} if /^Refn.:\s+(\d+)\s+(\S+)\s+(\d+)/;
  }
  close(F);') >${MYPROCDIR}__scope20840__processed_gdtts.out 2>${MYPROCDIR}__scope20840__processed_gdtts.err

