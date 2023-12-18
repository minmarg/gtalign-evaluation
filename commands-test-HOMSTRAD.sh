export WRKDIR=/home/mindaugas/projects/data

##process FATCAT results to obtain TM-scores for the alignments
(ls -1 fatcat | grep -E .aln | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
  $RDIR="fatcat";
  $WDIR="tmp_wrk_dir";
  $QIDfle="{}"; $QID=$QIDfle; $QID=~s/^fatcat__(.+).aln$/$1/;
  exit 0 unless $QID;
  $nonum = "-999999";
  $name = $qry = $rfn = $QID;
  $name =~ s/^(\S+)__\S+__\S+$/$1/; $qry =~ s/^\S+__(\S+)__\S+$/$1/; $rfn =~ s/^\S+__\S+__(\S+)$/$1/;
  $name = lc($name) if $name eq "MMOB";
  $qryfile = "homstrad/$name/$qry.atm";
  $rfnfile = "homstrad/$name/$rfn.atm";
  open(Q, $qryfile) || die "ERROR: Failed to open query: $qryfile";
  open(R, $rfnfile) || die "ERROR: Failed to open refn: $rfnfile";
  open(F, "$RDIR/$QIDfle") || die "ERROR: Failed to open results for: $QIDfle";
  while(<F>) {
    chomp;
    if(eof(F)) {
      $qidStart = $rfnStart = 0;
      seek(Q,0,0);
      while(<Q>) {
        next unless /^ATOM\s+\d+\s+CA [ A]/;
        $resnumc = substr($_,22,5); $resnumc =~ s/\s//g;
        last if $resnumc eq $qstart;
        $qidStart++;
      }
      while(<R>) {
        next unless /^ATOM\s+\d+\s+CA [ A]/;
        $resnumc = substr($_,22,5); $resnumc =~ s/\s//g;
        last if $resnumc eq $rstart;
        $rfnStart++;
      }
      close(R);
      $wfile = "$WDIR/${QID}.aln";
      $tmfile = "$WDIR/${QID}.tm";
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
        print(W ">$qry\n$qaln\n>$rfn\n$raln\n");
        close(W);
        ## FATCAT ignores HETATMs:
        $prog = "$TMalign $qryfile $rfnfile -I $wfile >$tmfile";
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
          $str = sprintf("%-36s Fsc= %8s Pvl= %8s  tm1= %8s tm2= %8s best= %8.6f\n",
            $QID,$Fsc,$pvl, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1);
          syswrite(STDOUT, $str);
          if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
            unlink($wfile, $tmfile);
          } else {
            print(STDERR "WARNING: Inconsistent alignments: $QID\n");
          }
        }
      }
    }
    if(/^Align\s+(\S+)\s+(\d+)\s+with\s+(\S+)\s+(\d+)/) {
      $qid = $1; $qlen = $2; $rfn = $3; $rlen = $4;
      $qaln = $raln = ""; $tm1 = $tm2 = -1;
      $qstart = $rstart = $nonum;
    }
    $Fsc = $1 if /^Twists\s+\S+\s+ini-len\s+\S+\s+ini-rmsd\s+\S+\s+opt-equ\s+\S+\s+opt-rmsd\s+\S+\s+chain-rmsd\s+\S+\s+Score\s+(\S+)\s+/;
    $pvl = $1 if /^P-value\s+(\S+)\s+/;
    do {$qaln .= $2; $qstart = $1 if $qstart eq $nonum;} if /^Chain\s+1:\s+(\S+)\s+(\S+)/;
    do {$raln .= $2; $rstart = $1 if $rstart eq $nonum;} if /^Chain\s+2:\s+(\S+)\s+(\S+)/;
  }
  close(F);close(Q);') >fatcat__homstrad__processed.out 2>fatcat__homstrad__processed.err

##process DeepAlign results to obtain TM-scores for the alignments
(ls -1 deepalign | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
  $RDIR="deepalign";
  $WDIR="tmp_wrk_dir";
  $QIDfle="{}"; $QID=$QIDfle; $QID=~s/^deepalign__(.+)$/$1/;
  exit 0 unless $QID;
  $nonum = -999999;
  $name = $qry = $rfn = $QID;
  $name =~ s/^(\S+)__\S+__\S+$/$1/; $qry =~ s/^\S+__(\S+)__\S+$/$1/; $rfn =~ s/^\S+__\S+__(\S+)$/$1/;
  $name = lc($name) if $name eq "MMOB";
  $qryfile = "homstrad/$name/$qry.atm";
  $rfnfile = "homstrad/$name/$rfn.atm";
  open(F, "$RDIR/$QIDfle") || die "ERROR: Failed to open results for: $QIDfle";
  while(<F>) {
    chomp;
    if(eof(F)) {
      $wfile = "$WDIR/${QID}.aln";
      $tmfile = "$WDIR/${QID}.tm";
      $qnnb = $qstart - 1; $qnne = $qlen - $qend;
      $rnnb = $rstart - 1; $rnne = $rlen - $rend;
      $qaln = "X"x$qnnb . $qaln . "X"x$qnne;
      $raln = "-"x$qnnb . $raln . "-"x$qnne;
      $qaln = "-"x$rnnb . $qaln . "-"x$rnne;
      $raln = "X"x$rnnb . $raln . "X"x$rnne;
      open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
      print(W ">$qry\n$qaln\n>$rfn\n$raln\n");
      close(W);
      $prog = "$TMalign $qryfile $rfnfile -I $wfile -het 1 >$tmfile";
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
        $str = sprintf("%-36s Dsc= %8s dtm= %8s  tm1= %8s tm2= %8s best= %8.6f\n",
          $QID,$dsc,$dtm, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1);
        syswrite(STDOUT, $str);
        if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
          unlink($wfile, $tmfile);
        } else {
          print(STDERR "WARNING: Inconsistent alignments: $QID\n");
        }
      }
    }
    if(/^1st\s+input\s+protein:\s+(\S+)\s+length=\s+(\d+)/) {
      $qid = "$1.atm"; $qlen = $2;
      $dsc = $dtm = -1; $qaln = $raln = ""; $tm1 = $tm2 = -1;
      $rlen = $qend = $rend = -1; $qstart = $rstart = $nonum;
    }
    do {$rfn = "$1.ent"; $rlen = $2;} if /^2nd\s+input\s+protein:\s+(\S+)\s+length=\s+(\d+)/;
    do {$dsc = $1; $dtm = $2;} if /^\s+[\d\-\.]+\s+[\d\-\.]+\s+([\d\-\.]+)\s+[\d\-\.]+\s+[\d\-\.]+\s+[\d\-\.]+\s+([\d\-\.]+)/;
    do {$qaln .= $2; $qend = $3; $qstart = $1 if $qstart == $nonum;} if /^T\s+\S+\s+(\d+)\s+(\S+)\s+(\d+)/;
    do {$raln .= $2; $rend = $3; $rstart = $1 if $rstart == $nonum;} if /^S\s+\S+\s+(\d+)\s+(\S+)\s+(\d+)/;
  }
  close(F);') >deepalign__homstrad__processed.out 2>deepalign__homstrad__processed.err

##process DALI results to obtain TM-scores for the alignments
(ls -1 dali | grep -E .txt | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
  $RDIR="dali";
  $WDIR="tmp_wrk_dir";
  $QIDfle="{}"; $QID=$QIDfle; $QID=~s/^dali_(.+).txt$/$1/;
  exit 0 unless $QID;
  $name = $qry = $rfn = $QID;
  $name =~ s/^(\S+)__\S+__\S+$/$1/; $qry =~ s/^\S+__(\S+)__\S+$/$1/; $rfn =~ s/^\S+__\S+__(\S+)$/$1/;
  $name = lc($name) if $name eq "MMOB";
  $qryfile = "homstrad/$name/$qry.atm";
  $rfnfile = "homstrad/$name/$rfn.atm";
  open(F, "$RDIR/$QIDfle") || die "ERROR: Failed to open results for: $QIDfle";
  while(<F>) {
    chomp;
    if(eof(F)) {
      $wfile = "$WDIR/${QID}.aln";
      $tmfile = "$WDIR/${QID}.tm";
      open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
      print(W ">$qry\n$qseq\n>$rfn\n$rseq\n");
      close(W);
      $prog = "$TMalign $qryfile $rfnfile -I $wfile -het 1 >$tmfile";
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
        $strprn = sprintf("%-36s  Z= %5s  tm1= %8s tm2= %8s best= %8.6f\n",
          $QID,$zsc,$tm1,$tm2,($tm1<$tm2)?$tm2:$tm1);
        syswrite(STDOUT, $strprn);
        if(length($tmqaln) == length($qseq) && length($tmraln) == length($rseq)) {
          unlink($wfile, $tmfile);
        } else {
          print(STDERR "WARNING: Inconsistent alignments: $QID\n");
        }
      }
    }
    if(/^No\s+(\S+)\s+Query=(\S+)\s+Sbjct=(\S+)\s+Z\-score=(\S+)/) {
      $nr = $1; $zsc = $4;
      $qseq = $rseq = ""; $tm1 = $tm2 = -1;
    }
    $qseq .= $1 if /^Query\s+(\S+)\s+/;
    $rseq .= $1 if /^Sbjct\s+(\S+)\s+/;
  }
  close(F);') >DALIv5__homstrad__processed.out 2>DALIv5__homstrad__processed.err

##collect TM-align results
(ls -1 tmalign | grep -E .out | xargs -i -P 40 perl -e '
  $RDIR="tmalign";
  $QIDfle="{}"; $QID=$QIDfle; $QID=~s/^tmalign__(.+).out$/$1/;
  exit 0 unless $QID;
  open(F, "$RDIR/$QIDfle") || die "ERROR: Failed to open results for: $QIDfle";
  while(<F>) {
    chomp;
    if(eof(F)) {
      $str = sprintf("%-36s  tm1= %8s tm2= %8s best= %8.6f\n",
        $QID, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1);
      syswrite(STDOUT, $str);
    }
    if(/^Name\s+of\s+Chain_1:\s+\S+\/(\S+)\s+\(to\s+be/) {
      $qid = $1; $tm1 = $tm2 = -1;
    }
    ##$rid = $1 if /^Name\s+of\s+Chain_2:\s+(\S+)/;
    $tm1 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_1/;
    $tm2 = $1 if /^TM\-score=\s+(\S+)\s+\(if normalized by length of Chain_2/;
  }
  close(F);') >tmalign__homstrad__processed.out 2>tmalign__homstrad__processed.err

##process foldseek results to obtain TM-scores for the alignments
export MYPROCDIR=foldseek
export MYPROCDIR=foldseek_alntyp2
(ls -1 ${MYPROCDIR} | grep -E .res | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
  $RDIR="$ENV{MYPROCDIR}";
  $WDIR="tmp_wrk_dir";
  $QIDfle="{}"; $QID=$QIDfle; $QID=~s/^${RDIR}__(.+).res$/$1/;
  exit 0 unless $QID;
  $name = $qry = $rfn = $QID;
  $name =~ s/^(\S+)__\S+__\S+$/$1/; $qry =~ s/^\S+__(\S+)__\S+$/$1/; $rfn =~ s/^\S+__\S+__(\S+)$/$1/;
  $name = lc($name) if $name eq "MMOB";
  $qryfile = "homstrad/$name/$qry.atm";
  $rfnfile = "homstrad/$name/$rfn.atm";
  open(F, "$RDIR/$QIDfle") || die "ERROR: Failed to open results for: $QIDfle";
  while(<F>) {
    chomp;
    if(/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/) {
      $qid = $1; $rid = $2; $ftme = $3; $ftm1 = $4; $ftm2 = $5;
      $qstart = $8; $qend = $9; $qlen = $10;  $tstart = $11; $tend = $12; $tlen = $13;
      $qaln = $14; $taln = $15;
      $wfile = "$WDIR/${QID}.aln";
      $tmfile = "$WDIR/${QID}.tm";
      $qnnb = $qstart - 1; $qnne = $qlen - $qend;
      $tnnb = $tstart - 1; $tnne = $tlen - $tend;
      $qaln = "X"x$qnnb . $qaln . "X"x$qnne;
      $taln = "-"x$qnnb . $taln . "-"x$qnne;
      $qaln = "-"x$tnnb . $qaln . "-"x$tnne;
      $taln = "X"x$tnnb . $taln . "X"x$tnne;
      open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
      print(W ">$qry\n$qaln\n>$rfn\n$taln\n");
      close(W);
      $prog = "$TMalign $qryfile $rfnfile -I $wfile -het 0 >$tmfile";
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
        $strprn = sprintf("%-36s ftme= %10s ftm1= %10s ftm2= %10s  tm1= %8s tm2= %8s best= %8.6f\n",
          $QID,$ftme,$ftm1,$ftm2, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1);
        syswrite(STDOUT, $strprn);
        if(length($tmqaln) == length($qaln) && length($tmtaln) == length($taln)) {
          unlink($wfile, $tmfile);
        } else {
          print(STDERR "WARNING: Inconsistent alignments: $QID\n");
        }
      }
    }
  }
') >${MYPROCDIR}__homstrad__processed.out 2>${MYPROCDIR}__homstrad__processed.err

##process GTalign results to obtain TM-align TM-scores for the alignments
export speed=0
export MYPROCDIR=gtalign_14_speed${speed}
(ls -1 ${MYPROCDIR} | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
  $RDIR="$ENV{MYPROCDIR}";
  $WDIR="tmp_wrk_dir";
  $QID="{}";
  exit 0 unless $QID;
  $name = $qry = $rfn = $QID;
  $name =~ s/^(\S+)__\S+__\S+$/$1/; $qry =~ s/^\S+__(\S+)__\S+$/$1/; $rfn =~ s/^\S+__\S+__(\S+)$/$1/;
  $name = lc($name) if $name eq "MMOB";
  $qryfile = "homstrad/$name/$qry.atm";
  $rfnfile = "homstrad/$name/$rfn.atm";
  $QIDfle = "${qry}__0.out";
  open(F, "$RDIR/$QID/$QIDfle") || die "ERROR: Failed to open results for: $QID/$QIDfle";
  while(<F>) {
    chomp;
    if(eof(F)) {
      $wfile = "$WDIR/${QID}.aln";
      $tmfile = "$WDIR/${QID}.tm";
      $qnnb = $qstart - 1; $qnne = $qlen - $qend;
      $rnnb = $rstart - 1; $rnne = $rlen - $rend;
      $qaln = "X"x$qnnb . $qaln . "X"x$qnne;
      $raln = "-"x$qnnb . $raln . "-"x$qnne;
      $qaln = "-"x$rnnb . $qaln . "-"x$rnne;
      $raln = "X"x$rnnb . $raln . "X"x$rnne;
      open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
      print(W ">$qry\n$qaln\n>$rfn\n$raln\n");
      close(W);
      $prog = "$TMalign $qryfile $rfnfile -I $wfile -het 1 >$tmfile";
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
        $strprn = sprintf("%-36s gtm1= %8s gtm2= %8s gbest= %8.6f  tm1= %8s tm2= %8s best= %8.6f\n",
          $QID,$gtm1,$gtm2,($gtm1<$gtm2)?$gtm2:$gtm1, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1);
        syswrite(STDOUT, $strprn);
        if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
          unlink($wfile, $tmfile);
        } else {
          print(STDERR "WARNING: Inconsistent alignments: $QID\n");
        }
      }
    }
    if(/^(\d+)\.\s*$/) {
      $nr = $1; $gtm1 = $gtm2 = -1; $qaln = $raln = ""; $tm1 = $tm2 = -1;
      $qlen = $rlen = $qend = $rend = -1; undef($qstart); undef($rstart);
    }
    ##do {$rfn = $1; $rfn =~ s/^\S+\/(\S+\.ent)$/$1/;} if /^>(\S+)/;
    do {$rlen = $1; $qlen = $2;} if /^\s+Length:\s+Refn.\s+=\s*(\d+),\s+Query\s+=\s*(\d+)/;
    do {$gtm2 = $1; $gtm1 = $2;} if /^\s+TM-score\s+\(Refn.\/Query\)\s+=\s+([\d\.]+)\s+\/\s+([\d\.]+)/;
    do {$qaln .= $2; $qend = $3; $qstart = $1 unless $qstart;} if /^Query:\s+(\d+)\s+(\S+)\s+(\d+)/;
    do {$raln .= $2; $rend = $3; $rstart = $1 unless $rstart;} if /^Refn.:\s+(\d+)\s+(\S+)\s+(\d+)/;
  }
  close(F);') >${MYPROCDIR}__homstrad__processed.out 2>${MYPROCDIR}__homstrad__processed.err

## HOMSTRAD reference alignments
(find homstrad -name '*.malf' | xargs -i -P 40 perl -e '
  $TMalign="/home/mindaugas/install/TM-align/TMalign";
  $WDIR="tmp_wrk_dir";
  $name = $dirname = "{}"; $dirname =~ s/^(.+)\/[^\/]+$/$1/; $name =~ s/^.+\/([^\/]+)\.malf$/$1/;
  while(<>) {
    last if /^Percentage/;
    next if (/^$/ || /^[\s\[]+/);
    if(/^(\S+)\s+(\S+)\s*$/) {
      $str = $1; $seq = $2; $seq =~ s/\*+$//;
      $h{$str} .= $seq;
      push @a, $str;
    }
  }
  $qry = $a[0];
  foreach $rfn(keys %h) {
    next if $rfn eq $qry;
    $rfnfile = "${dirname}/$rfn.atm";
    $qryfile = "${dirname}/$qry.atm";
    die "ERROR: File not found: $rfnfile" unless -f $rfnfile;
    die "ERROR: File not found: $qryfile" unless -f $qryfile;
    $QID = "${name}__${qry}__${rfn}";
    $wfile = "$WDIR/${QID}.aln";
    $tmfile = "$WDIR/${QID}.tm";
    $qaln = $h{$qry}; $raln = $h{$rfn};
    die "ERROR: Invalid alignment:\n$qaln\n$raln" if length($qaln) != length($raln);
    for($i = 0; $i < length($qaln);) {
      if((substr($qaln, $i, 1) eq "-") && (substr($raln, $i, 1) eq "-")) {
        substr($qaln, $i, 1) = "";
        substr($raln, $i, 1) = "";
      } else {
        $i++;
      }
    }
    open(W, ">$wfile") || die "ERROR: Failed to open file for writing: $wfile";
    print(W ">$qry\n$qaln\n>$rfn\n$raln\n");
    close(W);
    $prog = "$TMalign $qryfile $rfnfile -I $wfile -het 1 >$tmfile";
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
      $str = sprintf("%-36s  tm1= %8s tm2= %8s best= %8.6f\n",
        $QID, $tm1,$tm2,($tm1<$tm2)?$tm2:$tm1);
      syswrite(STDOUT, $str);
      if(length($tmqaln) == length($qaln) && length($tmraln) == length($raln)) {
        unlink($wfile, $tmfile);
      } else {
        print(STDERR "WARNING: Inconsistent alignments: $QID\n");
      }
    }
  }' {}) >reference__homstrad__processed.out 2>reference__homstrad__processed.err

