#!/bin/bash
export WRKDIR=/home/mindaugas/projects/data

## known cross-fold relationships:
export CFR=("b.66,b.67,b.68,b.69,b.70" "c.2,c.3,c.4,c.5,c.27,c.28,c.30,c.31")
## ignore matches across these fold groups:
export CFRIGNORE=0 #1

## Get recall1 data
function GetRecall1() {
  local ifile=$1 #input file
  local ofile=$2 #output file
  local fqry=$3 #query field number
  local frfn=$4 #reference field number
  local fmsr1=$5 #measure 1 field number (for average)
  local fmsr2=$6 #measure 2 field number (for average)
  local srt=$7 #sort flag
  local hmn=$8 #harmonic mean
  local lstfile="${WRKDIR}/scope-2.08/pdbstyle-2.08.fam.lst";
  local qrslist="${WRKDIR}/scope-2.08/queries.lst" #list of query ids
  FQRY=$fqry FRFN=$frfn FMSR1=$fmsr1 FMSR2=$fmsr2 SRT=$srt HMN=$hmn LSTFILE=$lstfile QRSLIST=$qrslist CFR=${CFR[@]} perl -e '
    @CFRs = ();
    @CFRs = split(/\s+/,$ENV{CFR}) if("$ENV{CFRIGNORE}" eq "1");
    open(F,"$ENV{LSTFILE}") || die "ERROR: Failed to open: $ENV{LSTFILE}";
    while(<F>) {
      chomp; next if /^#/;
      @a = split(/\s+/);
      die "ERROR: Invalid classification." unless($a[1] =~ /^(([a-z]+\.[0-9]+)\.[0-9]+)\.[0-9]+/);
      $fold = $2; $sfam = $1; $fam = $a[1]; $dom = $a[0];
      $Folds{$fold}++; $Sfams{$sfam}++; $Fams{$fam}++;
      $DomFold{$dom} = $fold; $DomSfam{$dom} = $sfam; $DomFam{$dom} = $fam;
    }
    close(F);
    open(F,"$ENV{QRSLIST}") || die "ERROR: Failed to open: $ENV{QRSLIST}";
    while(<F>) {
      chomp; next if /^#/;
      @a = split(/\s+/); $qry = $a[0];
      die "ERROR: Query $qry not found." unless((exists $DomFold{$qry}) && (exists $DomSfam{$qry}) && (exists $DomFam{$qry}));
      $fam = $DomFam{$qry}; $sfam = $DomSfam{$qry}; $fold = $DomFold{$qry};
      $szfam = $Fams{$fam}; $szsfam = $Sfams{$sfam}; $szfold = $Folds{$fold};
      ##effective sizes; -1: self matches excluded:
      $szefffam = $szfam - 1; $szeffsfam = $szsfam - $szfam; $szefffold = $szfold - $szsfam;
      ##effective numbers of queries:
      $QueriesFam++ if($szefffam > 0);
      $QueriesSfam++ if($szeffsfam > 0);
      $QueriesFold++ if($szefffold > 0);
    }
    close(F);
    while(<>) {
      chomp;
      @a = split(/\s+/);
      $qry = $a[$ENV{FQRY}]; $rfn = $a[$ENV{FRFN}]; $msr = $msr1 = $a[$ENV{FMSR1}];
      if($ENV{FMSR1} != $ENV{FMSR2}) {
        $msr2 = $a[$ENV{FMSR2}];
        if($ENV{HMN} eq "H") {
          $msr = (2 * $msr1 * $msr2) / ($msr1 + $msr2) if($msr1+$msr2);
        } else {
          $msr = ($msr1 + $msr2) * 0.5;
        }
      }
      $qry =~ s/.(?:ent|pdb)$//; $rfn =~ s/.(?:ent|pdb)$//;
      push @qrs, $qry; push @rfs, $rfn; push @mss, $msr;
      die "ERROR: Query $qry not found." unless((exists $DomFold{$qry}) && (exists $DomSfam{$qry}) && (exists $DomFam{$qry}));
      die "ERROR: Refn $rfn not found." unless((exists $DomFold{$rfn}) && (exists $DomSfam{$rfn}) && (exists $DomFam{$rfn}));
    }
    if($ENV{SRT} eq "A") {
      @nds=sort{$mss[$a]<=>$mss[$b]} 0..$#mss;
    } else {
      @nds=sort{$mss[$b]<=>$mss[$a]} 0..$#mss;
    }
    for($i = 0; $i <= $#nds; $i++) {
      $ndx = $nds[$i];
      $qry = $qrs[$ndx]; $rfn = $rfs[$ndx];
      $Queries{$qry}++;
      next if $qry eq $rfn;
      if($DomFold{$qry} ne $DomFold{$rfn}) {
        $ignore = 0; $qryfold = $DomFold{$qry}; $rfnfold = $DomFold{$rfn};
        foreach $cfr(@CFRs) {$ignore = 1 if($cfr =~ /(^|,)$qryfold(,|$)/ && $cfr =~ /(^|,)$rfnfold(,|$)/)}
        next if $ignore;
        $DiffFold{$qry}++;
        next
      }
      if($DomFam{$qry} eq $DomFam{$rfn}) {next if $DiffFold{$qry}; $SameFam{$qry}++; next}
      if(($DomFam{$qry} ne $DomFam{$rfn}) && ($DomSfam{$qry} eq $DomSfam{$rfn})) {next if $DiffFold{$qry}; $SameSfam{$qry}++; next}
      if(($DomSfam{$qry} ne $DomSfam{$rfn}) && ($DomFold{$qry} eq $DomFold{$rfn})) {next if $DiffFold{$qry}; $SameFold{$qry}++; next}
    }
    foreach $q(keys %Queries) {
      $fam = $DomFam{$q}; $sfam = $DomSfam{$q}; $fold = $DomFold{$q};
      $szfam = $Fams{$fam}; $szsfam = $Sfams{$sfam}; $szfold = $Folds{$fold};
      $samefam = $SameFam{$q}? $SameFam{$q}: 0;
      $samesfam = $SameSfam{$q}? $SameSfam{$q}: 0;
      $samefold = $SameFold{$q}? $SameFold{$q}: 0;
      $difffold = $DiffFold{$q}? $DiffFold{$q}: 0;
      $rec1fam = $rec1sfam = $rec1fold = 0;
      $rec1fam = $samefam / ($szfam - 1) if($szfam > 1); ##-1: self matches excluded
      $rec1sfam = $samesfam / ($szsfam - $szfam) if($szsfam - $szfam > 0);
      $rec1fold = $samefold / ($szfold - $szsfam) if($szfold - $szsfam > 0);
      $line = sprintf("%-10s %12s %6.4f %6.4f %6.4f  %5d  %5d %5d %5d  %5d %5d %5d", $q,$fam, $rec1fam,$rec1sfam,$rec1fold, $difffold, $szfam,$szsfam,$szfold, $QueriesFam,$QueriesSfam,$QueriesFold);
      $Res{$line} = $rec1fold;
    }
    printf("%s\n",$_) foreach (sort{$Res{$b}<=>$Res{$a}} keys %Res);
  ' $ifile >$ofile
}

## Get precision-recall data
function GetPR() {
  local ifile=$1 #input file
  local ofile=$2 #output file
  local fqry=$3 #query id field number
  local frfn=$4 #reference id field number
  local fmsr1=$5 #measure 1 field number (for average)
  local fmsr2=$6 #measure 2 field number (for average)
  local srt=$7 #sort flag: A, ascending; descending, otherwise
  local hmn=$8 #harmonic mean
  local lstfile="${WRKDIR}/scope-2.08/pdbstyle-2.08.fam.lst" #list of searched scope domain ids with their classification
  local qrslist="${WRKDIR}/scope-2.08/queries.lst" #list of query ids
  FQRY=$fqry FRFN=$frfn FMSR1=$fmsr1 FMSR2=$fmsr2 SRT=$srt HMN=$hmn LSTFILE=$lstfile QRSLIST=$qrslist CFR=${CFR[@]} perl -e '
    @CFRs = ();
    @CFRs = split(/\s+/,$ENV{CFR}) if("$ENV{CFRIGNORE}" eq "1");
    open(F,"$ENV{LSTFILE}") || die "ERROR: Failed to open: $ENV{LSTFILE}";
    while(<F>) {
      chomp; next if /^#/;
      @a = split(/\s+/);
      die "ERROR: Invalid classification." unless($a[1] =~ /^(([a-z]+\.[0-9]+)\.[0-9]+)\.[0-9]+/);
      $fold = $2; $sfam = $1; $fam = $a[1]; $dom = $a[0];
      $Folds{$fold}++; $Sfams{$sfam}++; $Fams{$fam}++;
      $DomFold{$dom} = $fold; $DomSfam{$dom} = $sfam; $DomFam{$dom} = $fam;
    }
    close(F);
    open(F,"$ENV{QRSLIST}") || die "ERROR: Failed to open: $ENV{QRSLIST}";
    while(<F>) {
      chomp; next if /^#/;
      @a = split(/\s+/); $qry = $a[0];
      die "ERROR: Query $qry not found." unless((exists $DomFold{$qry}) && (exists $DomSfam{$qry}) && (exists $DomFam{$qry}));
      $fam = $DomFam{$qry}; $sfam = $DomSfam{$qry}; $fold = $DomFold{$qry};
      $szfam = $Fams{$fam}; $szsfam = $Sfams{$sfam}; $szfold = $Folds{$fold};
      ##effective sizes; -1: self matches excluded:
      $szefffam = $szfam - 1; $szeffsfam = $szsfam - $szfam; $szefffold = $szfold - $szsfam;
      ##max weighted sum for a query equals 1
      $QueriesFam++ if($szefffam > 0);
      $QueriesSfam++ if($szeffsfam > 0);
      $QueriesFold++ if($szefffold > 0);
    }
    close(F);
    while(<>) {
      chomp;
      @a = split(/\s+/);
      $qry = $a[$ENV{FQRY}]; $rfn = $a[$ENV{FRFN}]; $msr = $msr1 = $a[$ENV{FMSR1}];
      if($ENV{FMSR1} != $ENV{FMSR2}) {
        $msr2 = $a[$ENV{FMSR2}];
        if($ENV{HMN} eq "H") {
          $msr = (2 * $msr1 * $msr2) / ($msr1 + $msr2) if($msr1+$msr2);
        } else {
          $msr = ($msr1 + $msr2) * 0.5;
        }
      }
      $qry =~ s/.(?:ent|pdb)$//; $rfn =~ s/.(?:ent|pdb)$//;
      push @qrs, $qry; push @rfs, $rfn; push @mss, $msr;
      die "ERROR: Query $qry not found." unless((exists $DomFold{$qry}) && (exists $DomSfam{$qry}) && (exists $DomFam{$qry}));
      die "ERROR: Refn $rfn not found." unless((exists $DomFold{$rfn}) && (exists $DomSfam{$rfn}) && (exists $DomFam{$rfn}));
    }
    if($ENV{SRT} eq "A") {
      @nds=sort{$mss[$a]<=>$mss[$b]} 0..$#mss;
    } else {
      @nds=sort{$mss[$b]<=>$mss[$a]} 0..$#mss;
    }
    $SameFam = $SameSfam = $SameFold = $DiffFold = 0;
    for($i = 0; $i <= $#nds; $i++) {
      $ndx = $nds[$i];
      $qry = $qrs[$ndx]; $rfn = $rfs[$ndx];
      $fam = $DomFam{$qry}; $sfam = $DomSfam{$qry}; $fold = $DomFold{$qry};
      $szfam = $Fams{$fam}; $szsfam = $Sfams{$sfam}; $szfold = $Folds{$fold};
      ##effective sizes; -1: self matches excluded:
      $szefffam = $szfam - 1; $szeffsfam = $szsfam - $szfam; $szefffold = $szfold - $szsfam;
      if($qry ne $rfn) {
        if($fold ne $DomFold{$rfn}) {
          $ignore = 0; $rfnfold = $DomFold{$rfn};
          foreach $cfr(@CFRs) {$ignore = 1 if($cfr =~ /(^|,)$fold(,|$)/ && $cfr =~ /(^|,)$rfnfold(,|$)/)}
          unless($ignore) {
            $norm = $szefffold; $norm = $szefffam if($szefffold < $szefffam); $norm = 1 if $norm < 1;
            $DiffFold += 1 / $norm;
          }
        }
        elsif($fam eq $DomFam{$rfn}) {
          $SameFam += 1 / $szefffam if($szefffam > 0);
        }
        elsif(($fam ne $DomFam{$rfn}) && ($sfam eq $DomSfam{$rfn})) {
          $SameSfam += 1 / $szeffsfam if($szeffsfam > 0);
        }
        elsif(($sfam ne $DomSfam{$rfn}) && ($fold eq $DomFold{$rfn})) {
          $SameFold += 1 / $szefffold if($szefffold > 0);
        }
      }
      if($i % 1000 == 0 || $i+1 > $#nds) {
        $p1fam = $SameFam + $DiffFold; $prec1fam = 1;
        $p1sfam = $SameSfam + $DiffFold; $prec1sfam = 1;
        $p1fold = $SameFold + $DiffFold; $prec1fold = 1;
        $prec1fam = $SameFam / $p1fam if $p1fam;      $rec1fam = $SameFam / $QueriesFam;
        $prec1sfam = $SameSfam / $p1sfam if $p1sfam;  $rec1sfam = $SameSfam / $QueriesSfam;
        $prec1fold = $SameFold / $p1fold if $p1fold;  $rec1fold = $SameFold / $QueriesFold;
        $line = sprintf("%6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f\n", $prec1fam,$prec1sfam,$prec1fold, $rec1fam,$rec1sfam,$rec1fold);
        $Res{$line} = $rec1fold;
      }
    }
    print foreach (sort{$Res{$a}<=>$Res{$b}} keys %Res);
  ' $ifile >$ofile
}

## Get area under precision-recall curve
function GetAUPRC() {
  local ifile=$1 #input file
  perl -e '
    while(<>) {
      next if /^#/; chomp;
      @a = split(/\s+/);
      push @Fam, [$a[3],$a[0]]; ##recall(x)--precision(y)
      push @Sfam, [$a[4],$a[1]]; ##recall(x)--precision(y)
      push @Fold, [$a[5],$a[2]]; ##recall(x)--precision(y)
    }
    @Fam = sort {$a->[0]<=>$b->[0]} @Fam; ##ascending
    @Sfam = sort {$a->[0]<=>$b->[0]} @Sfam; ##ascending
    @Fold = sort {$a->[0]<=>$b->[0]} @Fold; ##ascending
    $famauc = $sfamauc = $foldauc = 0;
    for($i=0; $i<$#Fam; $i++) { $famauc += $Fam[$i][1] * ($Fam[$i+1][0] - $Fam[$i][0]); }
    for($i=0; $i<$#Sfam; $i++) { $sfamauc += $Sfam[$i][1] * ($Sfam[$i+1][0] - $Sfam[$i][0]); }
    for($i=0; $i<$#Fold; $i++) { $foldauc += $Fold[$i][1] * ($Fold[$i+1][0] - $Fold[$i][0]); }
    printf(" %10s %10s %10s\n %10.6f %10.6f %10.6f\n","Family","S_family","Fold", $famauc,$sfamauc,$foldauc);
  ' $ifile
}


##process FATCAT results for fold recognition analysis
##
NAME=fatcat__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 0 1  5 5  A H
GetPR $INPUT ${OUTPUT}_PR 0 1  5 5  A H


##process DeepAlign results for fold recognition analysis
##
NAME=deepalign__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 0 1  3 3  D H
GetPR $INPUT ${OUTPUT}_PR 0 1  3 3  D H


##process DALI results for fold recognition analysis
##
NAME=DALIv5__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 3 4  6 6  D H
GetPR $INPUT ${OUTPUT}_PR 3 4  6 6  D H


##process TM-align results for fold recognition analysis
##
NAME=tmalign_fast__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 0 1  3 5  D H
GetPR $INPUT ${OUTPUT}_PR 0 1  3 5  D H

NAME=tmalign__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 0 1  3 5  D H
GetPR $INPUT ${OUTPUT}_PR 0 1  3 5  D H


##process foldseek results for fold recognition analysis
##
NAME=foldseek_alntyp2__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 0 1  3 3  A A #H
GetPR $INPUT ${OUTPUT}_PR 0 1  3 3  A A #H

NAME=foldseek_tmfast1__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 0 1  5 7  D A #H
GetPR $INPUT ${OUTPUT}_PR 0 1  5 7  D A #H

NAME=foldseek_tmfast0__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 0 1  5 7  D A #H
GetPR $INPUT ${OUTPUT}_PR 0 1  5 7  D A #H


##process GTalign results for fold recognition analysis (using 2TM-score for ranking; version 0.15.0)
##
VER=15
NAME=gtalign_${VER}_speed0_prescore03_addss_s044__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 1 2  4 6  D H
GetPR $INPUT ${OUTPUT}_PR 1 2  4 6  D H

NAME=gtalign_${VER}_speed9_prescore03_addss_s044__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 1 2  4 6  D H
GetPR $INPUT ${OUTPUT}_PR 1 2  4 6  D H

NAME=gtalign_${VER}_speed9_prescore04_addss_s044__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 1 2  4 6  D H
GetPR $INPUT ${OUTPUT}_PR 1 2  4 6  D H

NAME=gtalign_${VER}_speed13_prescore03_addss_s044__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 1 2  4 6  D H
GetPR $INPUT ${OUTPUT}_PR 1 2  4 6  D H

NAME=gtalign_${VER}_speed13_prescore04_addss_s044__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 1 2  4 6  D H
GetPR $INPUT ${OUTPUT}_PR 1 2  4 6  D H

NAME=gtalign_${VER}_speed13_prescore03_presim15_addss_s03__scope20840__processed; INPUT="${NAME}.out"; OUTPUT="${NAME}__foldrec.out${CFRIGNORE}"
GetRecall1 $INPUT $OUTPUT 1 2  4 6  D H
GetPR $INPUT ${OUTPUT}_PR 1 2  4 6  D H

