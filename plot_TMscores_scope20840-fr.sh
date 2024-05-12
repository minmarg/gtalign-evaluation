#!/bin/bash
name=$(basename $0 .sh)
## PARAMS
LEGEND=0 ##legend only
CFRIGNORE=1 ##ignore matches across related fold groups
PR=1 ##P-R plot; Recall1 when 0
##LOGSCALE="+ scale_y_continuous(trans=log10_trans(), breaks=trans_breaks('log10',function(x) 10^x), labels=trans_format('log10',math_format(10^.x)))"
## FILES
WIDTH=7.08; HEIGHT=2.4 #4.8
VER=15
GTA7_speed0_pre03="gtalign_${VER}_speed0_prescore03_addss_s044__scope20840__processed__foldrec.out${CFRIGNORE}"
GTA7_speed9_pre03="gtalign_${VER}_speed9_prescore03_addss_s044__scope20840__processed__foldrec.out${CFRIGNORE}"
GTA7_speed9_pre04="gtalign_${VER}_speed9_prescore04_addss_s044__scope20840__processed__foldrec.out${CFRIGNORE}"
GTA7_speed13_pre03="gtalign_${VER}_speed13_prescore03_addss_s044__scope20840__processed__foldrec.out${CFRIGNORE}"
GTA7_speed13_pre04="gtalign_${VER}_speed13_prescore04_addss_s044__scope20840__processed__foldrec.out${CFRIGNORE}"
GTA7_speed13_pre03_sim15="gtalign_${VER}_speed13_prescore03_presim15_addss_s03__scope20840__processed__foldrec.out${CFRIGNORE}"
TMalign="tmalign__scope20840__processed__foldrec.out${CFRIGNORE}"
TMalignFast="tmalign_fast__scope20840__processed__foldrec.out${CFRIGNORE}"
DeepAlign="deepalign__scope20840__processed__foldrec.out${CFRIGNORE}"
DALI5="DALIv5__scope20840__processed__foldrec.out${CFRIGNORE}"
Fatcat="fatcat__scope20840__processed__foldrec.out${CFRIGNORE}"
FS_="foldseek_alntyp2__scope20840__processed__foldrec.out${CFRIGNORE}"
FS1="foldseek_tmfast1__scope20840__processed__foldrec.out${CFRIGNORE}"
FS0="foldseek_tmfast0__scope20840__processed__foldrec.out${CFRIGNORE}"
## FUNCTIONS
function GetXYRec1() {
  local file=$1
  local fieldv=$2
  local nqrs=2045
  FLDV=$fieldv NQRS=$nqrs perl -e '
    while(<>){
      @a=split(/\s+/);
      push @sv, $a[$ENV{FLDV}];
    } 
    @sndx=sort{$sv[$b]<=>$sv[$a]} 0..$#sv; 
    $y=$x="c("; 
    for($i=0;$i<=$#sndx;$i++){
      do{$y.=",";$x.=","} if $i; 
      $y.=sprintf("%.4f",$sv[$sndx[$i]]); 
      $x.=sprintf("%.4f",($i+1) / $ENV{NQRS});
    } 
    $y.=")"; $x.=")"; 
    print $x, " ", $y' $file
}
function GetXYPR() {
  local file=$1
  local fieldv=$2
  local fieldk=$3
  FLDV=$fieldv FLDK=$fieldk perl -e '
    while(<>){
      @a=split(/\s+/);
      push @sk, $a[$ENV{FLDK}];
      push @sv, $a[$ENV{FLDV}];
    }
    @sndx=sort{$sv[$a]<=>$sv[$b]} 0..$#sv;
    $y=$x="c(";
    for($i=0;$i<=$#sndx;$i++){
      do{$y.=",";$x.=","} if $i;
      $y.=sprintf("%.4f",$sk[$sndx[$i]]);
      $x.=sprintf("%.4f",$sv[$sndx[$i]]);
    }
    $y.=")"; $x.=")";
    print $x, " ", $y' $file
}

## OUTPUT
output="${name}_${VER}_pr${PR}_cfr${CFRIGNORE}"
[ -n "${LOGSCALE}" ] && output="${output}_log"
[ "${LEGEND}" == "1" ] && output="${name}_legend"

GTA7_speed0_pre03_fam_r1=("GTalign --speed=0" $(GetXYRec1 "../$GTA7_speed0_pre03" 2))
GTA7_speed0_pre03_sfam_r1=("GTalign --speed=0" $(GetXYRec1 "../$GTA7_speed0_pre03" 3))
GTA7_speed0_pre03_fold_r1=("GTalign --speed=0" $(GetXYRec1 "../$GTA7_speed0_pre03" 4))

GTA7_speed9_pre03_fam_r1=("GTalign --speed=9 (def)" $(GetXYRec1 "../$GTA7_speed9_pre03" 2))
GTA7_speed9_pre03_sfam_r1=("GTalign --speed=9 (def)" $(GetXYRec1 "../$GTA7_speed9_pre03" 3))
GTA7_speed9_pre03_fold_r1=("GTalign --speed=9 (def)" $(GetXYRec1 "../$GTA7_speed9_pre03" 4))

GTA7_speed9_pre04_fam_r1=("GTalign --speed=9 --pre-score=0.4" $(GetXYRec1 "../$GTA7_speed9_pre04" 2))
GTA7_speed9_pre04_sfam_r1=("GTalign --speed=9 --pre-score=0.4" $(GetXYRec1 "../$GTA7_speed9_pre04" 3))
GTA7_speed9_pre04_fold_r1=("GTalign --speed=9 --pre-score=0.4" $(GetXYRec1 "../$GTA7_speed9_pre04" 4))

GTA7_speed13_pre03_fam_r1=("GTalign --speed=13" $(GetXYRec1 "../$GTA7_speed13_pre03" 2))
GTA7_speed13_pre03_sfam_r1=("GTalign --speed=13" $(GetXYRec1 "../$GTA7_speed13_pre03" 3))
GTA7_speed13_pre03_fold_r1=("GTalign --speed=13" $(GetXYRec1 "../$GTA7_speed13_pre03" 4))

GTA7_speed13_pre04_fam_r1=("GTalign --speed=13 --pre-score=0.4" $(GetXYRec1 "../$GTA7_speed13_pre04" 2))
GTA7_speed13_pre04_sfam_r1=("GTalign --speed=13 --pre-score=0.4" $(GetXYRec1 "../$GTA7_speed13_pre04" 3))
GTA7_speed13_pre04_fold_r1=("GTalign --speed=13 --pre-score=0.4" $(GetXYRec1 "../$GTA7_speed13_pre04" 4))

GTA7_speed13_pre03_sim15_fam_r1=("GTalign --speed=13 --pre-similarity=15" $(GetXYRec1 "../$GTA7_speed13_pre03_sim15" 2))
GTA7_speed13_pre03_sim15_sfam_r1=("GTalign --speed=13 --pre-similarity=15" $(GetXYRec1 "../$GTA7_speed13_pre03_sim15" 3))
GTA7_speed13_pre03_sim15_fold_r1=("GTalign --speed=13 --pre-similarity=15" $(GetXYRec1 "../$GTA7_speed13_pre03_sim15" 4))

TMalign_fam_r1=("TM-align" $(GetXYRec1 "../$TMalign" 2))
TMalign_sfam_r1=("TM-align" $(GetXYRec1 "../$TMalign" 3))
TMalign_fold_r1=("TM-align" $(GetXYRec1 "../$TMalign" 4))

TMalign_fast_fam_r1=("TM-align -fast" $(GetXYRec1 "../$TMalignFast" 2))
TMalign_fast_sfam_r1=("TM-align -fast" $(GetXYRec1 "../$TMalignFast" 3))
TMalign_fast_fold_r1=("TM-align -fast" $(GetXYRec1 "../$TMalignFast" 4))

Deepalign_fam_r1=("DeepAlign" $(GetXYRec1 "../$DeepAlign" 2))
Deepalign_sfam_r1=("DeepAlign" $(GetXYRec1 "../$DeepAlign" 3))
Deepalign_fold_r1=("DeepAlign" $(GetXYRec1 "../$DeepAlign" 4))

DALI5_fam_r1=("DALI" $(GetXYRec1 "../$DALI5" 2))
DALI5_sfam_r1=("DALI" $(GetXYRec1 "../$DALI5" 3))
DALI5_fold_r1=("DALI" $(GetXYRec1 "../$DALI5" 4))

Fatcat_fam_r1=("FATCAT" $(GetXYRec1 "../$Fatcat" 2))
Fatcat_sfam_r1=("FATCAT" $(GetXYRec1 "../$Fatcat" 3))
Fatcat_fold_r1=("FATCAT" $(GetXYRec1 "../$Fatcat" 4))

FS__fam_r1=("Foldseek" $(GetXYRec1 "../$FS_" 2))
FS__sfam_r1=("Foldseek" $(GetXYRec1 "../$FS_" 3))
FS__fold_r1=("Foldseek" $(GetXYRec1 "../$FS_" 4))

FS1_fam_r1=("Foldseek --tmalign-fast 1" $(GetXYRec1 "../$FS1" 2))
FS1_sfam_r1=("Foldseek --tmalign-fast 1" $(GetXYRec1 "../$FS1" 3))
FS1_fold_r1=("Foldseek --tmalign-fast 1" $(GetXYRec1 "../$FS1" 4))

FS0_fam_r1=("Foldseek --tmalign-fast 0" $(GetXYRec1 "../$FS0" 2))
FS0_sfam_r1=("Foldseek --tmalign-fast 0" $(GetXYRec1 "../$FS0" 3))
FS0_fold_r1=("Foldseek --tmalign-fast 0" $(GetXYRec1 "../$FS0" 4))



GTA7_speed0_pre03_fam_pr=("GTalign --speed=0" $(GetXYPR "../$GTA7_speed0_pre03"_PR 3 0))
GTA7_speed0_pre03_sfam_pr=("GTalign --speed=0" $(GetXYPR "../$GTA7_speed0_pre03"_PR 4 1))
GTA7_speed0_pre03_fold_pr=("GTalign --speed=0" $(GetXYPR "../$GTA7_speed0_pre03"_PR 5 2))

GTA7_speed9_pre03_fam_pr=("GTalign --speed=9 (def)" $(GetXYPR "../$GTA7_speed9_pre03"_PR 3 0))
GTA7_speed9_pre03_sfam_pr=("GTalign --speed=9 (def)" $(GetXYPR "../$GTA7_speed9_pre03"_PR 4 1))
GTA7_speed9_pre03_fold_pr=("GTalign --speed=9 (def)" $(GetXYPR "../$GTA7_speed9_pre03"_PR 5 2))

GTA7_speed9_pre04_fam_pr=("GTalign --speed=9 --pre-score=0.4" $(GetXYPR "../$GTA7_speed9_pre04"_PR 3 0))
GTA7_speed9_pre04_sfam_pr=("GTalign --speed=9 --pre-score=0.4" $(GetXYPR "../$GTA7_speed9_pre04"_PR 4 1))
GTA7_speed9_pre04_fold_pr=("GTalign --speed=9 --pre-score=0.4" $(GetXYPR "../$GTA7_speed9_pre04"_PR 5 2))

GTA7_speed13_pre03_fam_pr=("GTalign --speed=13" $(GetXYPR "../$GTA7_speed13_pre03"_PR 3 0))
GTA7_speed13_pre03_sfam_pr=("GTalign --speed=13" $(GetXYPR "../$GTA7_speed13_pre03"_PR 4 1))
GTA7_speed13_pre03_fold_pr=("GTalign --speed=13" $(GetXYPR "../$GTA7_speed13_pre03"_PR 5 2))

GTA7_speed13_pre04_fam_pr=("GTalign --speed=13 --pre-score=0.4" $(GetXYPR "../$GTA7_speed13_pre04"_PR 3 0))
GTA7_speed13_pre04_sfam_pr=("GTalign --speed=13 --pre-score=0.4" $(GetXYPR "../$GTA7_speed13_pre04"_PR 4 1))
GTA7_speed13_pre04_fold_pr=("GTalign --speed=13 --pre-score=0.4" $(GetXYPR "../$GTA7_speed13_pre04"_PR 5 2))

GTA7_speed13_pre03_sim15_fam_pr=("GTalign --speed=13 --pre-similarity=15" $(GetXYPR "../$GTA7_speed13_pre03_sim15"_PR 3 0))
GTA7_speed13_pre03_sim15_sfam_pr=("GTalign --speed=13 --pre-similarity=15" $(GetXYPR "../$GTA7_speed13_pre03_sim15"_PR 4 1))
GTA7_speed13_pre03_sim15_fold_pr=("GTalign --speed=13 --pre-similarity=15" $(GetXYPR "../$GTA7_speed13_pre03_sim15"_PR 5 2))

TMalign_fam_pr=("TM-align" $(GetXYPR "../$TMalign"_PR 3 0))
TMalign_sfam_pr=("TM-align" $(GetXYPR "../$TMalign"_PR 4 1))
TMalign_fold_pr=("TM-align" $(GetXYPR "../$TMalign"_PR 5 2))

TMalign_fast_fam_pr=("TM-align -fast" $(GetXYPR "../$TMalignFast"_PR 3 0))
TMalign_fast_sfam_pr=("TM-align -fast" $(GetXYPR "../$TMalignFast"_PR 4 1))
TMalign_fast_fold_pr=("TM-align -fast" $(GetXYPR "../$TMalignFast"_PR 5 2))

Deepalign_fam_pr=("DeepAlign" $(GetXYPR "../$DeepAlign"_PR 3 0))
Deepalign_sfam_pr=("DeepAlign" $(GetXYPR "../$DeepAlign"_PR 4 1))
Deepalign_fold_pr=("DeepAlign" $(GetXYPR "../$DeepAlign"_PR 5 2))

DALI5_fam_pr=("DALI" $(GetXYPR "../$DALI5"_PR 3 0))
DALI5_sfam_pr=("DALI" $(GetXYPR "../$DALI5"_PR 4 1))
DALI5_fold_pr=("DALI" $(GetXYPR "../$DALI5"_PR 5 2))

Fatcat_fam_pr=("FATCAT" $(GetXYPR "../$Fatcat"_PR 3 0))
Fatcat_sfam_pr=("FATCAT" $(GetXYPR "../$Fatcat"_PR 4 1))
Fatcat_fold_pr=("FATCAT" $(GetXYPR "../$Fatcat"_PR 5 2))

FS__fam_pr=("Foldseek" $(GetXYPR "../$FS_"_PR 3 0))
FS__sfam_pr=("Foldseek" $(GetXYPR "../$FS_"_PR 4 1))
FS__fold_pr=("Foldseek" $(GetXYPR "../$FS_"_PR 5 2))

FS1_fam_pr=("Foldseek --tmalign-fast 1" $(GetXYPR "../$FS1"_PR 3 0))
FS1_sfam_pr=("Foldseek --tmalign-fast 1" $(GetXYPR "../$FS1"_PR 4 1))
FS1_fold_pr=("Foldseek --tmalign-fast 1" $(GetXYPR "../$FS1"_PR 5 2))

FS0_fam_pr=("Foldseek --tmalign-fast 0" $(GetXYPR "../$FS0"_PR 3 0))
FS0_sfam_pr=("Foldseek --tmalign-fast 0" $(GetXYPR "../$FS0"_PR 4 1))
FS0_fold_pr=("Foldseek --tmalign-fast 0" $(GetXYPR "../$FS0"_PR 5 2))



Rtext="
library(ggplot2);
library(scales);
library(gridExtra);
##pdf('${output}.pdf',family='sans',width=3,height=3,pointsize=1);

dat <- list(
  '${GTA7_speed0_pre03_fam_r1[0]}'=cbind(x=${GTA7_speed0_pre03_fam_r1[1]},y=${GTA7_speed0_pre03_fam_r1[2]}), 
  '${GTA7_speed9_pre03_fam_r1[0]}'=cbind(x=${GTA7_speed9_pre03_fam_r1[1]},y=${GTA7_speed9_pre03_fam_r1[2]}), 
  '${GTA7_speed9_pre04_fam_r1[0]}'=cbind(x=${GTA7_speed9_pre04_fam_r1[1]},y=${GTA7_speed9_pre04_fam_r1[2]}), 
  '${GTA7_speed13_pre03_fam_r1[0]}'=cbind(x=${GTA7_speed13_pre03_fam_r1[1]},y=${GTA7_speed13_pre03_fam_r1[2]}), 
  '${GTA7_speed13_pre04_fam_r1[0]}'=cbind(x=${GTA7_speed13_pre04_fam_r1[1]},y=${GTA7_speed13_pre04_fam_r1[2]}), 
  '${GTA7_speed13_pre03_sim15_fam_r1[0]}'=cbind(x=${GTA7_speed13_pre03_sim15_fam_r1[1]},y=${GTA7_speed13_pre03_sim15_fam_r1[2]}), 
  '${TMalign_fam_r1[0]}'=cbind(x=${TMalign_fam_r1[1]},y=${TMalign_fam_r1[2]}), 
  '${TMalign_fast_fam_r1[0]}'=cbind(x=${TMalign_fast_fam_r1[1]},y=${TMalign_fast_fam_r1[2]}), 
  '${Deepalign_fam_r1[0]}'=cbind(x=${Deepalign_fam_r1[1]},y=${Deepalign_fam_r1[2]}), 
  '${DALI5_fam_r1[0]}'=cbind(x=${DALI5_fam_r1[1]},y=${DALI5_fam_r1[2]}), 
  '${Fatcat_fam_r1[0]}'=cbind(x=${Fatcat_fam_r1[1]},y=${Fatcat_fam_r1[2]}), 
  '${FS__fam_r1[0]}'=cbind(x=${FS__fam_r1[1]},y=${FS__fam_r1[2]}), 
  '${FS1_fam_r1[0]}'=cbind(x=${FS1_fam_r1[1]},y=${FS1_fam_r1[2]}), 
  '${FS0_fam_r1[0]}'=cbind(x=${FS0_fam_r1[1]},y=${FS0_fam_r1[2]}));

##credit: Roman Lustrik, stackoverflow
list.names <- names(dat)
lns <- sapply(dat, nrow)
dat <- as.data.frame(do.call('rbind', dat))
dat\$group <- rep(list.names, lns)

mypalette <- c(
  '#191A19', '#0f3113', '#006000', '#008800', '#466f46', '#639763', ##'#56b160', 
  '#a43131', '#ff3535', '#0000be', '#ff7000', '#854100', 
  '#BC6FF1', '#892CDC', '#52057B')

myltypes <- c(
  'solid', 'F1', '224282F2', 'F4448444', 'F282', '2262', ##'73',
  'solid', 'F1',    'solid',   'solid', 'solid', 
  '224282F2', 'F1', 'solid')

p10wlegend <- ggplot(dat, aes(x = x, y = y, colour = group)) +
  theme_bw() +
  theme(legend.text = element_text(size=8), legend.box.background = element_rect(colour = 'black'),
    legend.key.width=unit(0.04,'npc'), legend.key.height=unit(0.02,'npc'),
    legend.position='bottom', legend.box='vertical') +
  geom_line(aes(linetype = group), linewidth = 0.3) +
  scale_linetype_manual(name = '', breaks = list.names, values = myltypes) +
  scale_colour_manual(name = '', breaks = list.names, values = mypalette) +
  guides(linetype = guide_legend(ncol = 3), label.position = 'bottom')

P21 <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  ggtitle('Family') +
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x='Fraction of queries', y='Sensitivity up to 1st FP') +
  ##geom_line(linetype = 1, linewidth = 0.4) +
  geom_line(linewidth = 0.3) +
  scale_linetype_manual(breaks = list.names, values = myltypes) +
  scale_colour_manual(breaks = list.names, values = mypalette) ${LOGSCALE}
  ##scale_colour_discrete(breaks = list.names)


dat <- list(
  '${GTA7_speed0_pre03_sfam_r1[0]}'=cbind(x=${GTA7_speed0_pre03_sfam_r1[1]},y=${GTA7_speed0_pre03_sfam_r1[2]}), 
  '${GTA7_speed9_pre03_sfam_r1[0]}'=cbind(x=${GTA7_speed9_pre03_sfam_r1[1]},y=${GTA7_speed9_pre03_sfam_r1[2]}), 
  '${GTA7_speed9_pre04_sfam_r1[0]}'=cbind(x=${GTA7_speed9_pre04_sfam_r1[1]},y=${GTA7_speed9_pre04_sfam_r1[2]}), 
  '${GTA7_speed13_pre03_sfam_r1[0]}'=cbind(x=${GTA7_speed13_pre03_sfam_r1[1]},y=${GTA7_speed13_pre03_sfam_r1[2]}), 
  '${GTA7_speed13_pre04_sfam_r1[0]}'=cbind(x=${GTA7_speed13_pre04_sfam_r1[1]},y=${GTA7_speed13_pre04_sfam_r1[2]}), 
  '${GTA7_speed13_pre03_sim15_sfam_r1[0]}'=cbind(x=${GTA7_speed13_pre03_sim15_sfam_r1[1]},y=${GTA7_speed13_pre03_sim15_sfam_r1[2]}), 
  '${TMalign_sfam_r1[0]}'=cbind(x=${TMalign_sfam_r1[1]},y=${TMalign_sfam_r1[2]}), 
  '${TMalign_fast_sfam_r1[0]}'=cbind(x=${TMalign_fast_sfam_r1[1]},y=${TMalign_fast_sfam_r1[2]}), 
  '${Deepalign_sfam_r1[0]}'=cbind(x=${Deepalign_sfam_r1[1]},y=${Deepalign_sfam_r1[2]}), 
  '${DALI5_sfam_r1[0]}'=cbind(x=${DALI5_sfam_r1[1]},y=${DALI5_sfam_r1[2]}), 
  '${Fatcat_sfam_r1[0]}'=cbind(x=${Fatcat_sfam_r1[1]},y=${Fatcat_sfam_r1[2]}), 
  '${FS__sfam_r1[0]}'=cbind(x=${FS__sfam_r1[1]},y=${FS__sfam_r1[2]}), 
  '${FS1_sfam_r1[0]}'=cbind(x=${FS1_sfam_r1[1]},y=${FS1_sfam_r1[2]}), 
  '${FS0_sfam_r1[0]}'=cbind(x=${FS0_sfam_r1[1]},y=${FS0_sfam_r1[2]}));

list.names <- names(dat)
lns <- sapply(dat, nrow)
dat <- as.data.frame(do.call('rbind', dat))
dat\$group <- rep(list.names, lns)

P22 <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  ggtitle('Superfamily') +
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x='Fraction of queries', y='') + ##'Sensitivity up to 1st FP') +
  geom_line(linewidth = 0.3) +
  scale_linetype_manual(breaks = list.names, values = myltypes) +
  scale_colour_manual(breaks = list.names, values = mypalette) ${LOGSCALE}


dat <- list(
  '${GTA7_speed0_pre03_fold_r1[0]}'=cbind(x=${GTA7_speed0_pre03_fold_r1[1]},y=${GTA7_speed0_pre03_fold_r1[2]}), 
  '${GTA7_speed9_pre03_fold_r1[0]}'=cbind(x=${GTA7_speed9_pre03_fold_r1[1]},y=${GTA7_speed9_pre03_fold_r1[2]}), 
  '${GTA7_speed9_pre04_fold_r1[0]}'=cbind(x=${GTA7_speed9_pre04_fold_r1[1]},y=${GTA7_speed9_pre04_fold_r1[2]}), 
  '${GTA7_speed13_pre03_fold_r1[0]}'=cbind(x=${GTA7_speed13_pre03_fold_r1[1]},y=${GTA7_speed13_pre03_fold_r1[2]}), 
  '${GTA7_speed13_pre04_fold_r1[0]}'=cbind(x=${GTA7_speed13_pre04_fold_r1[1]},y=${GTA7_speed13_pre04_fold_r1[2]}), 
  '${GTA7_speed13_pre03_sim15_fold_r1[0]}'=cbind(x=${GTA7_speed13_pre03_sim15_fold_r1[1]},y=${GTA7_speed13_pre03_sim15_fold_r1[2]}), 
  '${TMalign_fold_r1[0]}'=cbind(x=${TMalign_fold_r1[1]},y=${TMalign_fold_r1[2]}), 
  '${TMalign_fast_fold_r1[0]}'=cbind(x=${TMalign_fast_fold_r1[1]},y=${TMalign_fast_fold_r1[2]}), 
  '${Deepalign_fold_r1[0]}'=cbind(x=${Deepalign_fold_r1[1]},y=${Deepalign_fold_r1[2]}), 
  '${DALI5_fold_r1[0]}'=cbind(x=${DALI5_fold_r1[1]},y=${DALI5_fold_r1[2]}), 
  '${Fatcat_fold_r1[0]}'=cbind(x=${Fatcat_fold_r1[1]},y=${Fatcat_fold_r1[2]}), 
  '${FS__fold_r1[0]}'=cbind(x=${FS__fold_r1[1]},y=${FS__fold_r1[2]}), 
  '${FS1_fold_r1[0]}'=cbind(x=${FS1_fold_r1[1]},y=${FS1_fold_r1[2]}), 
  '${FS0_fold_r1[0]}'=cbind(x=${FS0_fold_r1[1]},y=${FS0_fold_r1[2]}));

list.names <- names(dat)
lns <- sapply(dat, nrow)
dat <- as.data.frame(do.call('rbind', dat))
dat\$group <- rep(list.names, lns)

P23 <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  ggtitle('Fold') +
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x='Fraction of queries', y='') + ##'Sensitivity up to 1st FP') +
  geom_line(linewidth = 0.3) +
  scale_linetype_manual(breaks = list.names, values = myltypes) +
  scale_colour_manual(breaks = list.names, values = mypalette) ${LOGSCALE}



dat <- list(
  '${GTA7_speed0_pre03_fam_pr[0]}'=cbind(x=${GTA7_speed0_pre03_fam_pr[1]},y=${GTA7_speed0_pre03_fam_pr[2]}), 
  '${GTA7_speed9_pre03_fam_pr[0]}'=cbind(x=${GTA7_speed9_pre03_fam_pr[1]},y=${GTA7_speed9_pre03_fam_pr[2]}), 
  '${GTA7_speed9_pre04_fam_pr[0]}'=cbind(x=${GTA7_speed9_pre04_fam_pr[1]},y=${GTA7_speed9_pre04_fam_pr[2]}), 
  '${GTA7_speed13_pre03_fam_pr[0]}'=cbind(x=${GTA7_speed13_pre03_fam_pr[1]},y=${GTA7_speed13_pre03_fam_pr[2]}), 
  '${GTA7_speed13_pre04_fam_pr[0]}'=cbind(x=${GTA7_speed13_pre04_fam_pr[1]},y=${GTA7_speed13_pre04_fam_pr[2]}), 
  '${GTA7_speed13_pre03_sim15_fam_pr[0]}'=cbind(x=${GTA7_speed13_pre03_sim15_fam_pr[1]},y=${GTA7_speed13_pre03_sim15_fam_pr[2]}), 
  '${TMalign_fam_pr[0]}'=cbind(x=${TMalign_fam_pr[1]},y=${TMalign_fam_pr[2]}), 
  '${TMalign_fast_fam_pr[0]}'=cbind(x=${TMalign_fast_fam_pr[1]},y=${TMalign_fast_fam_pr[2]}), 
  '${Deepalign_fam_pr[0]}'=cbind(x=${Deepalign_fam_pr[1]},y=${Deepalign_fam_pr[2]}), 
  '${DALI5_fam_pr[0]}'=cbind(x=${DALI5_fam_pr[1]},y=${DALI5_fam_pr[2]}), 
  '${Fatcat_fam_pr[0]}'=cbind(x=${Fatcat_fam_pr[1]},y=${Fatcat_fam_pr[2]}), 
  '${FS__fam_pr[0]}'=cbind(x=${FS__fam_pr[1]},y=${FS__fam_pr[2]}), 
  '${FS1_fam_pr[0]}'=cbind(x=${FS1_fam_pr[1]},y=${FS1_fam_pr[2]}), 
  '${FS0_fam_pr[0]}'=cbind(x=${FS0_fam_pr[1]},y=${FS0_fam_pr[2]}));

list.names <- names(dat)
lns <- sapply(dat, nrow)
dat <- as.data.frame(do.call('rbind', dat))
dat\$group <- rep(list.names, lns)

P11 <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  ggtitle('Family') +
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x='Recall', y='Precision') +
  geom_line(linewidth = 0.3) +
  scale_linetype_manual(breaks = list.names, values = myltypes) +
  scale_colour_manual(breaks = list.names, values = mypalette) ${LOGSCALE}


dat <- list(
  '${GTA7_speed0_pre03_sfam_pr[0]}'=cbind(x=${GTA7_speed0_pre03_sfam_pr[1]},y=${GTA7_speed0_pre03_sfam_pr[2]}), 
  '${GTA7_speed9_pre03_sfam_pr[0]}'=cbind(x=${GTA7_speed9_pre03_sfam_pr[1]},y=${GTA7_speed9_pre03_sfam_pr[2]}), 
  '${GTA7_speed9_pre04_sfam_pr[0]}'=cbind(x=${GTA7_speed9_pre04_sfam_pr[1]},y=${GTA7_speed9_pre04_sfam_pr[2]}), 
  '${GTA7_speed13_pre03_sfam_pr[0]}'=cbind(x=${GTA7_speed13_pre03_sfam_pr[1]},y=${GTA7_speed13_pre03_sfam_pr[2]}), 
  '${GTA7_speed13_pre04_sfam_pr[0]}'=cbind(x=${GTA7_speed13_pre04_sfam_pr[1]},y=${GTA7_speed13_pre04_sfam_pr[2]}), 
  '${GTA7_speed13_pre03_sim15_sfam_pr[0]}'=cbind(x=${GTA7_speed13_pre03_sim15_sfam_pr[1]},y=${GTA7_speed13_pre03_sim15_sfam_pr[2]}), 
  '${TMalign_sfam_pr[0]}'=cbind(x=${TMalign_sfam_pr[1]},y=${TMalign_sfam_pr[2]}), 
  '${TMalign_fast_sfam_pr[0]}'=cbind(x=${TMalign_fast_sfam_pr[1]},y=${TMalign_fast_sfam_pr[2]}), 
  '${Deepalign_sfam_pr[0]}'=cbind(x=${Deepalign_sfam_pr[1]},y=${Deepalign_sfam_pr[2]}), 
  '${DALI5_sfam_pr[0]}'=cbind(x=${DALI5_sfam_pr[1]},y=${DALI5_sfam_pr[2]}), 
  '${Fatcat_sfam_pr[0]}'=cbind(x=${Fatcat_sfam_pr[1]},y=${Fatcat_sfam_pr[2]}), 
  '${FS__sfam_pr[0]}'=cbind(x=${FS__sfam_pr[1]},y=${FS__sfam_pr[2]}), 
  '${FS1_sfam_pr[0]}'=cbind(x=${FS1_sfam_pr[1]},y=${FS1_sfam_pr[2]}), 
  '${FS0_sfam_pr[0]}'=cbind(x=${FS0_sfam_pr[1]},y=${FS0_sfam_pr[2]}));

list.names <- names(dat)
lns <- sapply(dat, nrow)
dat <- as.data.frame(do.call('rbind', dat))
dat\$group <- rep(list.names, lns)

P12 <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  ggtitle('Superfamily') +
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x='Recall', y='') +
  geom_line(linewidth = 0.3) +
  scale_linetype_manual(breaks = list.names, values = myltypes) +
  scale_colour_manual(breaks = list.names, values = mypalette) ${LOGSCALE}


dat <- list(
  '${GTA7_speed0_pre03_fold_pr[0]}'=cbind(x=${GTA7_speed0_pre03_fold_pr[1]},y=${GTA7_speed0_pre03_fold_pr[2]}), 
  '${GTA7_speed9_pre03_fold_pr[0]}'=cbind(x=${GTA7_speed9_pre03_fold_pr[1]},y=${GTA7_speed9_pre03_fold_pr[2]}), 
  '${GTA7_speed9_pre04_fold_pr[0]}'=cbind(x=${GTA7_speed9_pre04_fold_pr[1]},y=${GTA7_speed9_pre04_fold_pr[2]}), 
  '${GTA7_speed13_pre03_fold_pr[0]}'=cbind(x=${GTA7_speed13_pre03_fold_pr[1]},y=${GTA7_speed13_pre03_fold_pr[2]}), 
  '${GTA7_speed13_pre04_fold_pr[0]}'=cbind(x=${GTA7_speed13_pre04_fold_pr[1]},y=${GTA7_speed13_pre04_fold_pr[2]}), 
  '${GTA7_speed13_pre03_sim15_fold_pr[0]}'=cbind(x=${GTA7_speed13_pre03_sim15_fold_pr[1]},y=${GTA7_speed13_pre03_sim15_fold_pr[2]}), 
  '${TMalign_fold_pr[0]}'=cbind(x=${TMalign_fold_pr[1]},y=${TMalign_fold_pr[2]}), 
  '${TMalign_fast_fold_pr[0]}'=cbind(x=${TMalign_fast_fold_pr[1]},y=${TMalign_fast_fold_pr[2]}), 
  '${Deepalign_fold_pr[0]}'=cbind(x=${Deepalign_fold_pr[1]},y=${Deepalign_fold_pr[2]}), 
  '${DALI5_fold_pr[0]}'=cbind(x=${DALI5_fold_pr[1]},y=${DALI5_fold_pr[2]}), 
  '${Fatcat_fold_pr[0]}'=cbind(x=${Fatcat_fold_pr[1]},y=${Fatcat_fold_pr[2]}), 
  '${FS__fold_pr[0]}'=cbind(x=${FS__fold_pr[1]},y=${FS__fold_pr[2]}), 
  '${FS1_fold_pr[0]}'=cbind(x=${FS1_fold_pr[1]},y=${FS1_fold_pr[2]}), 
  '${FS0_fold_pr[0]}'=cbind(x=${FS0_fold_pr[1]},y=${FS0_fold_pr[2]}));

list.names <- names(dat)
lns <- sapply(dat, nrow)
dat <- as.data.frame(do.call('rbind', dat))
dat\$group <- rep(list.names, lns)

P13 <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  ggtitle('Fold') +
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x='Recall', y='') +
  geom_line(linewidth = 0.3) +
  scale_linetype_manual(breaks = list.names, values = myltypes) +
  scale_colour_manual(breaks = list.names, values = mypalette) ${LOGSCALE}




# extract legend from ggplot (credit: Joachim Schork)
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1\$grobs, function(x) x\$name) == 'guide-box')
  step3 <- step1\$grobs[[step2]]
  return(step3)
}

legend10 <- extract_legend(p10wlegend);
hgtfct <- 1

if('${LEGEND}'=='1') {
  plot <- grid.arrange(arrangeGrob(legend10, ncol=1), nrow=1); hgtfct <- 1/2;
} else if('${PR}'=='1') {
  plot <- grid.arrange(arrangeGrob(P11,P12,P13, ncol=3), nrow=1);
} else {
  plot <- grid.arrange(arrangeGrob(P21,P22,P23, ncol=3), nrow=1);
}

ggsave(filename='${output}.pdf',plot,width=${WIDTH},height=${HEIGHT}*hgtfct)
#dev.off();
"
echo "$Rtext" | R --vanilla

