export WRKDIR=/home/mindaugas/projects/data

## make sure there are no HETATMs, MODELs, ENDMDLs
find homstrad -name '*.atm'|xargs -i grep -E HETATM {}|less
find homstrad -name '*.atm'|xargs -i grep -E MODEL {}|less
find homstrad -name '*.atm'|xargs -i grep -E ENDMDL {}|less
## make sure the structures are single-chain and N, CA, C, and O 
## atoms are present in the structures for equivalent interpretation;
## result: a single exception, which doesn't have any effect since it's the last residue:
##     homstrad/ghf1/1e73m.atm        333
find homstrad -name '*.atm' | xargs -i perl -e '
  $prvchn = "_"; $prvstr="    "; $prvnum=0;
  while(<>) {
    if(/^(?:ATOM|TER)/) {
      printf("%40s %10d\n", "{}", $prvnum) if (/^ATOM/) && ($prvchn ne "_") && (substr($_,21,1) ne $prvchn);
      if((substr($_,22,4) ne $prvstr) || (substr($_,21,1) ne $prvchn)) {
        $prvchn = substr($_,21,1);
        $prvstr = substr($_,22,4);
        $prvnum++;
        printf("%40s %10d\n", "{}", $prvnum) unless(($n && $ca && $c && $o) || $prvnum == 1);
        $rec = ""; $n = $ca = $c = $o = 0;
      }
      $ca = 1 if(substr($_,12,4) eq " CA " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $n = 1  if(substr($_,12,4) eq " N  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $c = 1  if(substr($_,12,4) eq " C  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $o = 1  if(substr($_,12,4) eq " O  " && (substr($_,16,1) eq " " || substr($_,16,1) eq "A"));
      $rec .= $_;
    }
    else {
      printf("%40s %10d\n", "{}", $prvnum) unless($n && $ca && $c && $o && $rec);
      $rec = ""; $n = $ca = $c = $o = 0;
    }
  }' {}

## FATCAT
## env var
export FATCAT=/home/mindaugas/install/FATCAT/FATCAT-dist
mkdir fatcat
## run FATCAT
(find homstrad -name '*.malf' | xargs -i -P 40 perl -e '
  $name = $dirname = "{}"; $dirname =~ s/^(.+)\/[^\/]+$/$1/; $name =~ s/^.+\/([^\/]+)\.malf$/$1/;
  while(<>) {
    last if ($qry && (/^$/ || /^[\s\[]+/));
    if($qry && /^(\S+)\s+\S+\s*$/) {
      $rfn = $1;
      last if $qry eq $rfn;
      $rfnfile = "${dirname}/$rfn.atm";
      $qryfile = "${dirname}/$qry.atm";
      die "ERROR: File not found: $rfnfile" unless -f $rfnfile;
      die "ERROR: File not found: $qryfile" unless -f $qryfile;
      $outfile = "fatcat/fatcat__${name}__${qry}__${rfn}";
      $cmd = "/home/mindaugas/install/FATCAT/FATCAT-dist/FATCATMain/FATCAT -p1 $qryfile -p2 $rfnfile -r -o $outfile -m";
      $ret = system($cmd);
      print "ERROR: $cmd\n" if($ret != 0 && $ret != 256); ##fatcat returns 256
      next;
    }
    $qry = $1 if(!$qry && /^(\S+)\s+\S+\s*$/);
  }' {}) >fatcat.log 2>&1

## run DeepAlign
mkdir deepalign
(find homstrad -name '*.malf' | xargs -i -P 40 perl -e '
  $name = $dirname = "{}"; $dirname =~ s/^(.+)\/[^\/]+$/$1/; $name =~ s/^.+\/([^\/]+)\.malf$/$1/;
  while(<>) {
    last if ($qry && (/^$/ || /^[\s\[]+/));
    if($qry && /^(\S+)\s+\S+\s*$/) {
      $rfn = $1;
      last if $qry eq $rfn;
      $rfnfile = "${dirname}/$rfn.atm";
      $qryfile = "${dirname}/$qry.atm";
      die "ERROR: File not found: $rfnfile" unless -f $rfnfile;
      die "ERROR: File not found: $qryfile" unless -f $qryfile;
      $outfile = "deepalign/deepalign__${name}__${qry}__${rfn}";
      $cmd = "/home/mindaugas/install/DeepAlign/DeepAlign/DeepAlign $qryfile $rfnfile >$outfile";
      $ret = system($cmd);
      print "ERROR: $cmd\n" if($ret != 0);
      next;
    }
    $qry = $1 if(!$qry && /^(\S+)\s+\S+\s*$/);
  }' {}) >deepalign.log 2>&1

## DALI v5
mkdir dali
(find homstrad -name '*.malf' | xargs -i -P 40 perl -e '
  $name = $dirname = "{}"; $dirname =~ s/^(.+)\/[^\/]+$/$1/; $name =~ s/^.+\/([^\/]+)\.malf$/$1/;
  while(<>) {
    last if ($qry && (/^$/ || /^[\s\[]+/));
    if($qry && /^(\S+)\s+\S+\s*$/) {
      $rfn = $1;
      last if $qry eq $rfn;
      $basename = "${name}__${qry}__${rfn}";
      die "ERROR: mkdir: $basename" unless mkdir("dali/$basename");
      die "ERROR: chdir: $basename" unless chdir("dali/$basename");
      die "ERROR: mkdir: $basename: dat1" unless mkdir("dat1");
      die "ERROR: mkdir: $basename: dat2" unless mkdir("dat2");
      $rfnfile = "../../${dirname}/$rfn.atm";
      $qryfile = "../../${dirname}/$qry.atm";
      die "ERROR: File not found: $rfnfile" unless -f $rfnfile;
      die "ERROR: File not found: $qryfile" unless -f $qryfile;
      $outfile = "dali__${basename}.txt";
      $cmd = "/home/mindaugas/install/DaliLite_v5/DaliLite.v5/bin/dali.pl --pdbfile1 $qryfile --pdbfile2 $rfnfile --dat1 dat1 --dat2 dat2 --outfmt alignments";
      $ret = system($cmd);
      if($ret != 0) {
        print "ERROR: $cmd\n";
      } else {
        $cmd = "mv *.txt ../dali_${basename}.txt";
        $ret = system($cmd);
        if($ret != 0) {
          print "ERROR: $cmd\n";
          chdir("../../");
        } else {
          chdir("../../");
          system("rm -fR dali/$basename");
        }
      }
      next;
    }
    $qry = $1 if(!$qry && /^(\S+)\s+\S+\s*$/);
  }' {}) >dali.log 2>&1

## run TM-align
mkdir tmalign
(find homstrad -name '*.malf' | xargs -i -P 40 perl -e '
  $name = $dirname = "{}"; $dirname =~ s/^(.+)\/[^\/]+$/$1/; $name =~ s/^.+\/([^\/]+)\.malf$/$1/;
  while(<>) {
    last if ($qry && (/^$/ || /^[\s\[]+/));
    if($qry && /^(\S+)\s+\S+\s*$/) {
      $rfn = $1;
      last if $qry eq $rfn;
      $rfnfile = "${dirname}/$rfn.atm";
      $qryfile = "${dirname}/$qry.atm";
      die "ERROR: File not found: $rfnfile" unless -f $rfnfile;
      die "ERROR: File not found: $qryfile" unless -f $qryfile;
      $outfile = "tmalign/tmalign__${name}__${qry}__${rfn}.out";
      $cmd = "/home/mindaugas/install/TM-align/TMalign $qryfile $rfnfile -het 1 >$outfile";
      $ret = system($cmd);
      print "ERROR: $cmd\n" if($ret != 0);
      next;
    }
    $qry = $1 if(!$qry && /^(\S+)\s+\S+\s*$/);
  }' {}) >tmalign.log 2>&1

## run foldseek using --alignment-type 1 and --tmalign-fast 0
mkdir foldseek
(find homstrad -name '*.malf' | xargs -i -P 40 perl -e '
  $name = $dirname = "{}"; $dirname =~ s/^(.+)\/[^\/]+$/$1/; $name =~ s/^.+\/([^\/]+)\.malf$/$1/;
  while(<>) {
    last if ($qry && (/^$/ || /^[\s\[]+/));
    if($qry && /^(\S+)\s+\S+\s*$/) {
      $rfn = $1;
      last if $qry eq $rfn;
      $rfnfile = "${dirname}/$rfn.atm";
      $qryfile = "${dirname}/$qry.atm";
      die "ERROR: File not found: $rfnfile" unless -f $rfnfile;
      die "ERROR: File not found: $qryfile" unless -f $qryfile;
      $outfile = "foldseek/foldseek__${name}__${qry}__${rfn}.res";
      $cmd = "/home/mindaugas/install/foldseek/foldseek/bin/foldseek easy-search $qryfile $rfnfile $outfile tmp_dir_foldseek__${name}__${qry}__${rfn} --threads 1 --alignment-type 1 --format-output query,target,evalue,qtmscore,ttmscore,alntmscore,rmsd,qstart,qend,qlen,tstart,tend,tlen,qaln,taln --tmalign-fast 0 --prefilter-mode 2 -e 1e6";
      $ret = system($cmd);
      print "ERROR: $cmd\n" if($ret != 0);
      next;
    }
    $qry = $1 if(!$qry && /^(\S+)\s+\S+\s*$/);
  }' {}) >foldseek.log 2>&1
## foldseek using --alignment-type 2
mkdir foldseek_alntyp2
(find homstrad -name '*.malf' | xargs -i -P 40 perl -e '
  $name = $dirname = "{}"; $dirname =~ s/^(.+)\/[^\/]+$/$1/; $name =~ s/^.+\/([^\/]+)\.malf$/$1/;
  while(<>) {
    last if ($qry && (/^$/ || /^[\s\[]+/));
    if($qry && /^(\S+)\s+\S+\s*$/) {
      $rfn = $1;
      last if $qry eq $rfn;
      $rfnfile = "${dirname}/$rfn.atm";
      $qryfile = "${dirname}/$qry.atm";
      die "ERROR: File not found: $rfnfile" unless -f $rfnfile;
      die "ERROR: File not found: $qryfile" unless -f $qryfile;
      $outfile = "foldseek_alntyp2/foldseek_alntyp2__${name}__${qry}__${rfn}.res";
      $cmd = "/home/mindaugas/install/foldseek/foldseek/bin/foldseek easy-search $qryfile $rfnfile $outfile tmp_dir_foldseek__${name}__${qry}__${rfn} --threads 1 --alignment-type 2 --format-output query,target,evalue,qtmscore,ttmscore,alntmscore,rmsd,qstart,qend,qlen,tstart,tend,tlen,qaln,taln --prefilter-mode 2 -e 1e6";
      $ret = system($cmd);
      print "ERROR: $cmd\n" if($ret != 0);
      next;
    }
    $qry = $1 if(!$qry && /^(\S+)\s+\S+\s*$/);
  }' {}) >foldseek_alntyp2.log 2>&1

## GTalign speed 0/9/13
export speed=0
mkdir gtalign_14_speed${speed}
(find homstrad -name '*.malf' | xargs -i -P 3 perl -e '
  $name = $dirname = "{}"; $dirname =~ s/^(.+)\/[^\/]+$/$1/; $name =~ s/^.+\/([^\/]+)\.malf$/$1/;
  while(<>) {
    last if ($qry && (/^$/ || /^[\s\[]+/));
    if($qry && /^(\S+)\s+\S+\s*$/) {
      $rfn = $1;
      last if $qry eq $rfn;
      $rfnfile = "${dirname}/$rfn.atm";
      $qryfile = "${dirname}/$qry.atm";
      die "ERROR: File not found: $rfnfile" unless -f $rfnfile;
      die "ERROR: File not found: $qryfile" unless -f $qryfile;
      $outdir = "gtalign_14_speed$ENV{speed}/${name}__${qry}__${rfn}";
      $cmd = "/home/mindaugas/projects/share/gtalign/bin/gtalign --qrs=$qryfile --rfs=$rfnfile -o $outdir --hetatm --dev-min-length=3 --speed=$ENV{speed} --pre-score=0 -s 0 --add-search-by-ss --dev-mem=4096";
      $ret = system($cmd);
      print "ERROR: $cmd\n" if($ret != 0);
      next;
    }
    $qry = $1 if(!$qry && /^(\S+)\s+\S+\s*$/);
  }' {}) >gtalign_14_speed${speed}.log 2>&1

