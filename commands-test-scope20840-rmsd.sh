export WRKDIR=/home/mindaugas/projects/data

##process FATCAT results to obtain TM-scores for the alignments
(ls -1 fatcat | grep -E .out | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
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
          $prog = "$TMalign $QRYDIR/$qid $RFNDIR/$rfn -I $wfile >$tmfile";
          if(system($prog) != 0) {
            print(STDERR "ERROR: Failed to execute: $prog\n");
          } else {
            open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
            while(<W>) {
              $rmsd = $1 if /RMSD=\s+([\d\.]+),/;
              $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
              $tm2 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/;
              last if /^\(":" denotes residue pairs/;
            }
            $tmqaln = <W>; $tmraln = <W>; $tmraln = <W>; chomp($tmqaln); chomp($tmraln);
            close(W);
            $str = sprintf("%-12s %12s Fsc= %8s Pvl= %8s  tm1= %8s tm2= %8s best= %8.6f rmsd= %6s\n",
              $qid,$rfn,$Fsc,$pvl, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1, $rmsd);
            syswrite(STDOUT, $str);
            if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
              unlink($wfile, $tmfile);
            } else {
              print(STDERR "WARNING: Inconsistent alignments: $qid $rfn\n");
            }
          }
        }
      }
      $qid = $qidtt; $rfn = $rfntt; $Fsc = $pvl = -1; $qaln = $raln = ""; $tm1 = $tm2 = $rmsd = -1;
      $qlen = $qlentt; $rlen = $rlentt; $qend = $rend = 0; undef($qstart); undef($rstart);
    }
    $Fsc = $1 if /^Twists\s+\S+\s+ini-len\s+\S+\s+ini-rmsd\s+\S+\s+opt-equ\s+\S+\s+opt-rmsd\s+\S+\s+chain-rmsd\s+\S+\s+Score\s+(\S+)\s+/;
    $pvl = $1 if /^P-value\s+(\S+)\s+/;
    do {$qaln .= $2; $qstart = $1 unless $qstart;} if /^Chain\s+1:\s+(\S+)\s+(\S+)/;
    do {$raln .= $2; $rstart = $1 unless $rstart;} if /^Chain\s+2:\s+(\S+)\s+(\S+)/;
  }
  close(F);close(Q);') >fatcat__scope20840__processed_rmsd.out 2>fatcat__scope20840__processed_rmsd.err

##get the mean #aligned residues for FATCAT
(ls -1 fatcat | grep -E .out | xargs -i cat fatcat/{} | perl -e '
  while(<>) {
    $pvl = $1 if /^P-value\s+(\S+)\s+/;
    do {$naln += length($1); $ngps += ($1=~tr/-//);} if /^Chain\s+1:\s+\S+\s+(\S+)/;
    $ngps += ($1=~tr/-//) if /^Chain\s+2:\s+\S+\s+(\S+)/;
    if(/^Align\s+(\S+)\s+\S+\s+with\s+(\S+)\s+/ || eof) {
      if($naln && ($qry ne $rfn)) {
        $nar = $naln - $ngps;
        push @A, [$nar,$pvl];
      }
      $pvl = 1; $naln = $ngps = 0;
      $qry = $1; $rfn = $2; $qry =~ s/\.(?:ent|pdb)$//; $rfn =~ s/\.(?:ent|pdb)$//;
    }
  }
  $N = 750000; ##max #top alignments considered in the plot
  @A = sort {$a->[1]<=>$b->[1]} @A; ##sort by P-score (ascending order)
  $sum = $sum2 = $nn = 0;
  for($i=0; $i<=$#A && $i<$N; $i++) {
    $sum += $A[$i][0]; $sum2 += $A[$i][0]*$A[$i][0]; $nn++;
  }
  $mn = $sum/$nn; $mn2 = $sum2/$nn; $sd = sqrt(($mn2 - $mn*$mn) * $nn/($nn-1));
  printf("Mean (s.d.) / N: %8.2f (%.2f) / %6d\n", $mn, $sd, $nn);
  ') >fatcat__scope20840__rmsd_meanalnlen.out


##process DeepAlign results to obtain TM-scores for the alignments
(ls -1 deepalign | grep -E .out | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
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
        $prog = "$TMalign $QRYDIR/$qid $RFNDIR/$rfn -I $wfile -het 1 >$tmfile";
        if(system($prog) != 0) {
          print(STDERR "ERROR: Failed to execute: $prog\n");
        } else {
          open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
          while(<W>) {
            $rmsd = $1 if /RMSD=\s+([\d\.]+),/;
            $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
            $tm2 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/;
            last if /^\(":" denotes residue pairs/;
          }
          $tmqaln = <W>; $tmraln = <W>; $tmraln = <W>; chomp($tmqaln); chomp($tmraln);
          close(W);
          $str = sprintf("%-12s %12s Dsc= %8s dtm= %8s  tm1= %8s tm2= %8s best= %8.6f rmsd= %6s\n",
            $qid,$rfn,$dsc,$dtm, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1, $rmsd);
          syswrite(STDOUT, $str);
          if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
            unlink($wfile, $tmfile);
          } else {
            print(STDERR "WARNING: Inconsistent alignments: $qid $rfn\n");
          }
        }
      }
      $qid = $qidtt; $dsc = $dtm = -1; $qaln = $raln = $rfn = ""; $tm1 = $tm2 = $rmsd = -1;
      $qlen = $qlentt; $rlen = $qend = $rend = -1; undef($qstart); undef($rstart);
    }
    do {$rfn = "$1.ent"; $rlen = $2;} if /^2nd\s+input\s+protein:\s+(\S+)\s+length=\s+(\d+)/;
    do {$dsc = $1; $dtm = $2;} if /^\s+[\d\-\.]+\s+[\d\-\.]+\s+([\d\-\.]+)\s+[\d\-\.]+\s+[\d\-\.]+\s+[\d\-\.]+\s+([\d\-\.]+)/;
    do {$qaln .= $2; $qend = $3; $qstart = $1 unless $qstart;} if /^T\s+\S+\s+(\d+)\s+(\S+)\s+(\d+)/;
    do {$raln .= $2; $rend = $3; $rstart = $1 unless $rstart;} if /^S\s+\S+\s+(\d+)\s+(\S+)\s+(\d+)/;
  }
  close(F);') >deepalign__scope20840__processed_rmsd.out 2>deepalign__scope20840__processed_rmsd.err

##get the mean #aligned residues for DeepAlign
(ls -1 deepalign | grep -E .out | xargs -i cat deepalign/{} | perl -e '
  while(<>) {
    $qry = $1 if /^1st\s+input\s+protein:\s+(\S+)/;
    $rfn = $1 if /^2nd\s+input\s+protein:\s+(\S+)/;
    push @A, [$2,$3] if(/^\s+[\d\-\.]+\s+[\d\-\.]+\s+([\d\-\.]+)\s+[\d\-\.]+\s+([\d\-\.]+)\s+[\d\-\.]+\s+([\d\-\.]+)/ && ($qry ne $rfn));
  }
  $N = 750000; ##max #top alignments considered in the plot
  @A = sort {$b->[1]<=>$a->[1]} @A; ##sort by TM-score (descending order)
  $sum = $sum2 = $nn = 0;
  for($i=0; $i<=$#A && $i<$N; $i++) {
    $sum += $A[$i][0]; $sum2 += $A[$i][0]*$A[$i][0]; $nn++;
  }
  $mn = $sum/$nn; $mn2 = $sum2/$nn; $sd = sqrt(($mn2 - $mn*$mn) * $nn/($nn-1));
  printf("Mean (s.d.) / N: %8.2f (%.2f) / %6d\n", $mn, $sd, $nn);
  ') >deepalign__scope20840__rmsd_meanalnlen.out


##process DALI results to obtain TM-scores for the alignments
(ls -1 dali | grep -E .txt | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
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
          $prog = "$TMalign $QRYDIR/$qstr $RFNDIR/$rstr -I $wfile -het 1 >$tmfile";
          if(system($prog) != 0) {
            print(STDERR "ERROR: Failed to execute: $prog\n");
          } else {
            open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
            while(<W>) {
              $rmsd = $1 if /RMSD=\s+([\d\.]+),/;
              $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
              $tm2 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/;
              last if /^\(":" denotes residue pairs/;
            }
            $tmqaln = <W>; $tmraln = <W>; $tmraln = <W>; chomp($tmqaln); chomp($tmraln);
            close(W);
            $strprn = sprintf("No.%-5s %5s %5s %12s %12s Z= %5s  tm1= %8s tm2= %8s best= %8.6f rmsd= %6s\n",
              $nr,$qid,$rfn,$qstr,$rstr,$zsc,$tm1,$tm2,($tm1<$tm2)?$tm2:$tm1, $rmsd);
            syswrite(STDOUT, $strprn);
            if(length($tmqaln) == length($qseq) && length($tmraln) == length($rseq)) {
              unlink($wfile, $tmfile);
            } else {
              print(STDERR "WARNING: Inconsistent alignments: $qid $rfn: $qstr $rstr\n");
            }
          }
        }
      }
      $nr = $nrtt; $qid = $qidtt; $rfn = $rfntt; $zsc = $zsctt; $qseq = $rseq = ""; $tm1 = $tm2 = $rmsd = -1;
    }
    $qseq .= $1 if /^Query\s+(\S+)\s+/;
    $rseq .= $1 if /^Sbjct\s+(\S+)\s+/;
  }
  close(F);') >DALIv5__scope20840__processed_rmsd.out 2>DALIv5__scope20840__processed_rmsd.err

##get the mean #aligned residues for Dali
(ls -1 dali | grep -E .txt | xargs -i cat dali/{} | perl -e '
  while(<>) {
    do {$qs .= $1; $naln += length($1); $ngps += ($1=~tr/-//);} if /^Query\s+(\S+)\s+/;
    do {$ss .= $1; $ngps += ($1=~tr/-//);} if /^Sbjct\s+(\S+)\s+/;
    if(/^No.+Z-score=([\d\.]+)/ || eof) {
      if($zsc && $naln && ($qs ne $ss)) {
        $nar = $naln - $ngps;
        push @A, [$nar,$zsc];
      }
      $zsc = $1; $naln = $ngps = 0; $qs = $ss = "";
    }
  }
  $N = 750000; ##max #top alignments considered in the plot
  @A = sort {$b->[1]<=>$a->[1]} @A; ##sort by Z-score (descending order)
  $sum = $sum2 = $nn = 0;
  for($i=0; $i<=$#A && $i<$N; $i++) {
    $sum += $A[$i][0]; $sum2 += $A[$i][0]*$A[$i][0]; $nn++;
  }
  $mn = $sum/$nn; $mn2 = $sum2/$nn; $sd = sqrt(($mn2 - $mn*$mn) * $nn/($nn-1));
  printf("Mean (s.d.) / N: %8.2f (%.2f) / %6d\n", $mn, $sd, $nn);
  ') >DALIv5__scope20840__rmsd_meanalnlen.out


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
    if(eof(F) || /^Name\s+of\s+Chain_1:\s+\S+\/(\S+)\s+\(to\s+be/) {
      $qidtt = "";
      unless(eof(F)) {
        print(STDERR "WARNING: Query $QID doesnt match $1.\n") unless "$QID.ent" eq $1;
        $qidtt = $1;
      }
      if($qid && $rfn && $qid ne $rfn) {
        $str = sprintf("%-12s %12s  tm1= %8s tm2= %8s best= %8.6f rmsd= %6s\n",
          $qid,$rfn, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1, $rmsd);
        syswrite(STDOUT, $str);
      }
      $qid = $qidtt; $tm1 = $tm2 = $rmsd = -1;
    }
    $rfn = $1 if /^Name\s+of\s+Chain_2:\s+(\S+)/;
    $rmsd = $1 if /RMSD=\s+([\d\.]+),/;
    $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
    $tm2 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/;
  }
  close(F);') >${MYPROCDIR}__scope20840__processed_rmsd.out 2>${MYPROCDIR}__scope20840__processed_rmsd.err

##get the mean #aligned residues for TM-align
export MYPROCDIR=tmalign
export MYPROCDIR=tmalign_fast
(ls -1 ${MYPROCDIR} | grep -E .out | xargs -i cat ${MYPROCDIR}/{} | perl -e '
  while(<>) {
    $nar = $1 if /^Aligned\s+length=\s*(\d+)/;
    $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
    if(/^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/) {
      $tm2 = $1;
      push @A, [$nar,(($tm1<$tm2)?$tm2:$tm1)] unless(($tm1 == $tm2) && (1 <= $tm1)); ##greater TM-score
    }
  }
  $N = 750000; ##max #top alignments considered in the plot
  @A = sort {$b->[1]<=>$a->[1]} @A; ##sort by TM-score (descending order)
  $sum = $sum2 = $nn = 0;
  for($i=0; $i<=$#A && $i<$N; $i++) {
    $sum += $A[$i][0]; $sum2 += $A[$i][0]*$A[$i][0]; $nn++;
  }
  $mn = $sum/$nn; $mn2 = $sum2/$nn; $sd = sqrt(($mn2 - $mn*$mn) * $nn/($nn-1));
  printf("Mean (s.d.) / N: %8.2f (%.2f) / %6d\n", $mn, $sd, $nn);
  ') >${MYPROCDIR}__scope20840__rmsd_meanalnlen.out


##process foldseek results to obtain TM-scores for the alignments
export MYPROCDIR=foldseek_tmfast0
#export MYPROCDIR=foldseek_tmfast1
#export MYPROCDIR=foldseek_alntyp2
(cat ${MYPROCDIR}__scope20840_bmk.res | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
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
        $prog = "$TMalign $QRYDIR/$qid $RFNDIR/$rfn -I $wfile -het 0 >$tmfile";
        if(system($prog) != 0) {
          print(STDERR "ERROR: Failed to execute: $prog\n");
        } else {
          open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
          while(<W>) {
            $rmsd = $1 if /RMSD=\s+([\d\.]+),/;
            $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
            $tm2 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/;
            last if /^\(":" denotes residue pairs/;
          }
          $tmqaln = <W>; $tmtaln = <W>; $tmtaln = <W>; chomp($tmqaln); chomp($tmtaln);
          close(W);
          $strprn = sprintf("%-12s %12s ftme= %10s ftm1= %10s ftm2= %10s  tm1= %8s tm2= %8s best= %8.6f rmsd= %6s\n",
            $qid,$rfn,$ftme,$ftm1,$ftm2, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1, $rmsd);
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
') >${MYPROCDIR}__scope20840__processed_rmsd.out 2>${MYPROCDIR}__scope20840__processed_rmsd.err

##get the mean #aligned residues for Foldseek
export MYPROCDIR=foldseek_tmfast0
export MYPROCDIR=foldseek_tmfast1
export MYPROCDIR=foldseek_alntyp2
(cat ${MYPROCDIR}__scope20840_bmk.res | perl -e '
  while(<>) {
    if(/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/) {
      $qid = $1; $rfn = $2; $ftme = $3; $ftm1 = $4; $ftm2 = $5;
      $qstart = $8; $qend = $9; $qlen = $10;  $tstart = $11; $tend = $12; $tlen = $13;
      $qaln = $14; $taln = $15; $ngps = ($qaln=~tr/-//) + ($taln=~tr/-//); $naln = length($qaln);
      if($qid && $rfn && $qaln && $taln && $qid ne $rfn) {
        $nar = $naln - $ngps;
        push @A, [$nar,(($ftm1<$ftm2)?$ftm2:$ftm1)]; ##greater TM-score
      }
    }
  }
  $N = 750000; ##max #top alignments considered in the plot
  @A = sort {$b->[1]<=>$a->[1]} @A; ##sort by TM-score (descending order)
  $sum = $sum2 = $nn = 0;
  for($i=0; $i<=$#A && $i<$N; $i++) {
    $sum += $A[$i][0]; $sum2 += $A[$i][0]*$A[$i][0]; $nn++;
  }
  $mn = $sum/$nn; $mn2 = $sum2/$nn; $sd = sqrt(($mn2 - $mn*$mn) * $nn/($nn-1));
  printf("Mean (s.d.) / N: %8.2f (%.2f) / %6d\n", $mn, $sd, $nn);
  ') >${MYPROCDIR}__scope20840__rmsd_meanalnlen.out


##process GTalign results to obtain TM-align TM-scores for the alignments
export MYPROCDIR=gtalign_14_speed0_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed9_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed9_prescore04_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore04_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore03_presim15_addss_s03
export MYPROCDIR=gtalign_14_speed13_prescore04_presim15_addss_s03
(ls -1 ${MYPROCDIR} | grep -E .out | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
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
        $prog = "$TMalign $QRYDIR/${QID}.ent $RFNDIR/$rfn -I $wfile -het 1 >$tmfile";
        if(system($prog) != 0) {
          print(STDERR "ERROR: Failed to execute: $prog\n");
        } else {
          open(W, $tmfile) || die "ERROR: Failed to open file: $tmfile";
          while(<W>) {
            $rmsd = $1 if /RMSD=\s+([\d\.]+),/;
            $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
            $tm2 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/;
            last if /^\(":" denotes residue pairs/;
          }
          $tmqaln = <W>; $tmraln = <W>; $tmraln = <W>; chomp($tmqaln); chomp($tmraln);
          close(W);
          $strprn = sprintf("%-5s %12s %12s gtm1= %8s gtm2= %8s gbest= %8.6f  tm1= %8s tm2= %8s best= %8.6f rmsd= %6s\n",
            $nr,"$QID.ent",$rfn,$gtm1,$gtm2,($gtm1<$gtm2)?$gtm2:$gtm1, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1, $rmsd);
          syswrite(STDOUT, $strprn);
          if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
            unlink($wfile, $tmfile);
          } else {
            print(STDERR "WARNING: Inconsistent alignments: $QID $rfn\n");
          }
        }
      }
      $nr = $nrtt; $gtm1 = $gtm2 = -1; $qaln = $raln = $rfn = ""; $tm1 = $tm2 = $rmsd = -1;
      $qlen = $rlen = $qend = $rend = -1; undef($qstart); undef($rstart);
    }
    do {$rfn = $1; $rfn =~ s/^\S+\/(\S+\.ent)$/$1/;} if /^>(\S+)/;
    do {$rlen = $1; $qlen = $2;} if /^\s+Length:\s+Refn.\s+=\s*(\d+),\s+Query\s+=\s*(\d+)/;
    do {$gtm2 = $1; $gtm1 = $2;} if /^\s+TM-score\s+\(Refn.\/Query\)\s+=\s+([\d\.]+)\s+\/\s+([\d\.]+)/;
    do {$qaln .= $2; $qend = $3; $qstart = $1 unless $qstart;} if /^Query:\s+(\d+)\s+(\S+)\s+(\d+)/;
    do {$raln .= $2; $rend = $3; $rstart = $1 unless $rstart;} if /^Refn.:\s+(\d+)\s+(\S+)\s+(\d+)/;
  }
  close(F);') >${MYPROCDIR}__scope20840__processed_rmsd.out 2>${MYPROCDIR}__scope20840__processed_rmsd.err

##get the mean #aligned residues for GTalign
export MYPROCDIR=gtalign_14_speed0_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed9_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed9_prescore04_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore03_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore04_addss_s044
export MYPROCDIR=gtalign_14_speed13_prescore03_presim15_addss_s03
(ls -1 ${MYPROCDIR} | grep -E .out | xargs -i cat ${MYPROCDIR}/{} | perl -e '
  while(<>) {
    do {$gtm2 = $1; $gtm1 = $2;} if /^\s+TM-score\s+\(Refn.\/Query\)\s+=\s+([\d\.]+)\s+\/\s+([\d\.]+)/;
    if(/^\s+(?:Identities|Matched|Gaps)\s+=\s+\d+\/(\d+)/) {
      $naln = $1;
      $ngps = 0;
      $ngps = $1 if(/\s+Gaps\s+=\s+(\d+)\/\d+/);
      $nar = $naln - $ngps;
      push @A, [$nar,(($gtm1<$gtm2)?$gtm2:$gtm1)] unless(($gtm1 == $gtm2) && (1 <= $gtm1)); ##greater TM-score
    }
  }
  $N = 750000; ##max #top alignments considered in the plot
  @A = sort {$b->[1]<=>$a->[1]} @A; ##sort by TM-score (descending order)
  $sum = $sum2 = $nn = 0;
  for($i=0; $i<=$#A && $i<$N; $i++) {
    $sum += $A[$i][0]; $sum2 += $A[$i][0]*$A[$i][0]; $nn++;
  }
  $mn = $sum/$nn; $mn2 = $sum2/$nn; $sd = sqrt(($mn2 - $mn*$mn) * $nn/($nn-1));
  printf("Mean (s.d.) / N: %8.2f (%.2f) / %6d\n", $mn, $sd, $nn);
  ') >${MYPROCDIR}__scope20840__rmsd_meanalnlen.out


