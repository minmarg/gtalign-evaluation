#!/bin/bash
name=$(basename $0 .sh)
## PARAMS
##LOGSCALE="+ scale_y_continuous(trans=log10_trans(), limits=c(1,NA), breaks=trans_breaks('log10',function(x) 10^x), labels=trans_format('log10',math_format(10^.x)))"
## FILES
WIDTH=2.36; HEIGHT=2.4
GTA7_speed0_pre03="gtalign_14_speed0_prescore03_addss_s044__pdb20__processed_rmsd.out"
GTA7_speed9_pre03="gtalign_14_speed9_prescore03_addss_s044__pdb20__processed_rmsd.out"
GTA7_speed9_pre04="gtalign_14_speed9_prescore04_addss_s044__pdb20__processed_rmsd.out"
GTA7_speed13_pre03="gtalign_14_speed13_prescore03_addss_s044__pdb20__processed_rmsd.out"
GTA7_speed13_pre04="gtalign_14_speed13_prescore04_addss_s044__pdb20__processed_rmsd.out"
GTA7_speed13_pre03_sim15="gtalign_14_speed13_prescore03_presim15_addss_s03__pdb20__processed_rmsd.out"
TMalign="tmalign__pdb20__processed_rmsd.out"
TMalignFast="tmalign_fast__pdb20__processed_rmsd.out"
DeepAlign="deepalign__pdb20__processed_rmsd.out"
DALI5="DALIv5__pdb20__processed_rmsd.out"
Fatcat="fatcat__pdb20__processed_rmsd.out"
FS_="foldseek_alntyp2__pdb20__processed_rmsd.out"
FS1="foldseek_tmfast1__pdb20__processed_rmsd.out"
FS0="foldseek_tmfast0__pdb20__processed_rmsd.out"
## FUNCTIONS
function GetXYSum() {
  local file=$1
  local cnt=$2
  local fieldv1=$3
  local fieldv2=$4
  local fieldk=$5
  local srt=$6
  NUM=$cnt FLDV1=$fieldv1 FLDV2=$fieldv2 FLDK=$fieldk SRT=$srt perl -e '
    while(<>){
      @a=split(/\s+/); push @sk,$a[$ENV{FLDK}]; 
      push @sv, ($a[$ENV{FLDV1}]<$a[$ENV{FLDV2}])? $a[$ENV{FLDV2}]: $a[$ENV{FLDV1}];
#testing:
#last if 100000 < $c++;
#last if $ENV{NUM} < $c++;
    } 
    if($ENV{SRT} eq "A") {
      @sndx=sort{$sv[$a]<=>$sv[$b]} 0..$#sv; 
    } else {
      @sndx=sort{$sv[$b]<=>$sv[$a]} 0..$#sv; 
    }
    $y=$x="c("; 
    $sum=$sum2=$nn=0; 
    for($i=0;$i<$ENV{NUM} && $i<=$#sndx;$i++){
      $last=($i+1)*0.001 if 0.5<=$sk[$sndx[$i]]; 
      $sum+=$sk[$sndx[$i]]; $sum2+=$sk[$sndx[$i]]*$sk[$sndx[$i]]; $nn++;
      if($i%400==0 || $ENV{NUM}<=$i+1 || $#sndx<$i+1){
        do{$y.=",";$x.=","} if $i; 
        $y.=sprintf("%.4f",$sum*0.001); 
        $x.=sprintf("%.3f",($i+1)*0.001);
      }
    } 
    $y.=")"; $x.=")"; $mn = $sum/$nn; $mn2 = $sum2/$nn; $sd = sqrt(($mn2-$mn*$mn) * $nn/($nn-1));
    print $x, " ", $y, " ", $last, " ", $mn, " ", $sd' $file
}
## OUTPUT
output="${name}"
GTA7F=(8 8 16 D); TMF=(7 7 9 D); DAF=(5 5 13 D); DaF=(6 6 14 D); FCF=(5 5 13 A); FSF=(5 7 15 D)
[ -n "${LOGSCALE}" ] && output="${output}_log"

NN=112000
GTA7_speed0_pre03_dat=("GTalign --speed=0" $(GetXYSum "../$GTA7_speed0_pre03" $NN  ${GTA7F[@]}))
GTA7_speed9_pre03_dat=("GTalign --speed=9 (def)" $(GetXYSum "../$GTA7_speed9_pre03" $NN  ${GTA7F[@]}))
GTA7_speed9_pre04_dat=("GTalign --speed=9 --pre-score=0.4" $(GetXYSum "../$GTA7_speed9_pre04" $NN  ${GTA7F[@]}))
GTA7_speed13_pre03_dat=("GTalign --speed=13" $(GetXYSum "../$GTA7_speed13_pre03" $NN  ${GTA7F[@]}))
GTA7_speed13_pre04_dat=("GTalign --speed=13 --pre-score=0.4" $(GetXYSum "../$GTA7_speed13_pre04" $NN  ${GTA7F[@]}))
GTA7_speed13_pre03_sim15_dat=("GTalign --speed=13 --pre-similarity=15" $(GetXYSum "../$GTA7_speed13_pre03_sim15" $NN  ${GTA7F[@]}))

TMalign_dat=("TM-align" $(GetXYSum "../$TMalign" $NN  ${TMF[@]}))
TMalign_fast_dat=("TM-align -fast" $(GetXYSum "../$TMalignFast" $NN  ${TMF[@]}))
Deepalign_dat=("DeepAlign" $(GetXYSum "../$DeepAlign" $NN  ${DAF[@]})) ##3 3 13 D)) ##:by deepscore
DALI5_dat=("DALI" $(GetXYSum "../$DALI5" $NN  ${DaF[@]}))
Fatcat_dat=("FATCAT" $(GetXYSum "../$Fatcat" $NN  ${FCF[@]}))
FS__dat=("Foldseek" $(GetXYSum "../$FS_" $NN  ${FSF[@]})) ##3 3 15 A)) ##:by Evalue
FS1_dat=("Foldseek --tmalign-fast 1" $(GetXYSum "../$FS1" $NN  ${FSF[@]}))
FS0_dat=("Foldseek --tmalign-fast 0" $(GetXYSum "../$FS0" $NN  ${FSF[@]}))



Rtext="
library(ggplot2);
library(scales);
library(gridExtra);
##pdf('${output}.pdf',family='sans',width=3,height=3,pointsize=1);

dat <- list(
  '${GTA7_speed0_pre03_dat[0]}'=cbind(x=${GTA7_speed0_pre03_dat[1]},y=${GTA7_speed0_pre03_dat[2]}), 
  '${GTA7_speed9_pre03_dat[0]}'=cbind(x=${GTA7_speed9_pre03_dat[1]},y=${GTA7_speed9_pre03_dat[2]}), 
  '${GTA7_speed9_pre04_dat[0]}'=cbind(x=${GTA7_speed9_pre04_dat[1]},y=${GTA7_speed9_pre04_dat[2]}), 
  '${GTA7_speed13_pre03_dat[0]}'=cbind(x=${GTA7_speed13_pre03_dat[1]},y=${GTA7_speed13_pre03_dat[2]}), 
  '${GTA7_speed13_pre04_dat[0]}'=cbind(x=${GTA7_speed13_pre04_dat[1]},y=${GTA7_speed13_pre04_dat[2]}), 
  '${GTA7_speed13_pre03_sim15_dat[0]}'=cbind(x=${GTA7_speed13_pre03_sim15_dat[1]},y=${GTA7_speed13_pre03_sim15_dat[2]}), 
  '${TMalign_dat[0]}'=cbind(x=${TMalign_dat[1]},y=${TMalign_dat[2]}), 
  '${TMalign_fast_dat[0]}'=cbind(x=${TMalign_fast_dat[1]},y=${TMalign_fast_dat[2]}), 
  '${Deepalign_dat[0]}'=cbind(x=${Deepalign_dat[1]},y=${Deepalign_dat[2]}), 
  '${DALI5_dat[0]}'=cbind(x=${DALI5_dat[1]},y=${DALI5_dat[2]}), 
  '${Fatcat_dat[0]}'=cbind(x=${Fatcat_dat[1]},y=${Fatcat_dat[2]}), 
  '${FS__dat[0]}'=cbind(x=${FS__dat[1]},y=${FS__dat[2]}), 
  '${FS1_dat[0]}'=cbind(x=${FS1_dat[1]},y=${FS1_dat[2]}), 
  '${FS0_dat[0]}'=cbind(x=${FS0_dat[1]},y=${FS0_dat[2]}));

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
  theme(legend.text = element_text(size = 8), legend.key.width = unit(0.04,'npc'), legend.key.height = unit(0.02,'npc')) +
  geom_line(aes(linetype = group), linewidth = 0.3) +
  scale_linetype_manual(name = '', breaks = list.names, values = myltypes) +
  scale_colour_manual(name = '', breaks = list.names, values = mypalette)

P1 <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x=expression(paste('# top hits ('%*%10^{3},')')), y='') + ##expression(paste('Cumulative RMSD ('%*%10^{3},')'))) +
  annotate('text', x=-Inf, y=Inf, hjust=-0.05, vjust=1.8, label='PDB', size=4) +
  ##geom_line(linetype = 1, linewidth = 0.4) +
  geom_line(linewidth = 0.3) +
  scale_linetype_manual(breaks = list.names, values = myltypes) +
  scale_colour_manual(breaks = list.names, values = mypalette) ${LOGSCALE}
  ##scale_colour_discrete(breaks = list.names)


# extract legend from ggplot (credit: Joachim Schork)
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1\$grobs, function(x) x\$name) == 'guide-box')
  step3 <- step1\$grobs[[step2]]
  return(step3)
}

legend10 <- extract_legend(p10wlegend);

plot <- grid.arrange(arrangeGrob(P1, ncol=1), nrow=1);

sprintf('%-40s %.2f %.4f','${GTA7_speed0_pre03_dat[0]}',${GTA7_speed0_pre03_dat[4]},${GTA7_speed0_pre03_dat[5]});
sprintf('%-40s %.2f %.4f','${GTA7_speed9_pre03_dat[0]}',${GTA7_speed9_pre03_dat[4]},${GTA7_speed9_pre03_dat[5]});
sprintf('%-40s %.2f %.4f','${GTA7_speed9_pre04_dat[0]}',${GTA7_speed9_pre04_dat[4]},${GTA7_speed9_pre04_dat[5]});
sprintf('%-40s %.2f %.4f','${GTA7_speed13_pre03_dat[0]}',${GTA7_speed13_pre03_dat[4]},${GTA7_speed13_pre03_dat[5]});
sprintf('%-40s %.2f %.4f','${GTA7_speed13_pre04_dat[0]}',${GTA7_speed13_pre04_dat[4]},${GTA7_speed13_pre04_dat[5]});
sprintf('%-40s %.2f %.4f','${GTA7_speed13_pre03_sim15_dat[0]}',${GTA7_speed13_pre03_sim15_dat[4]},${GTA7_speed13_pre03_sim15_dat[5]});
sprintf('%-40s %.2f %.4f','${TMalign_dat[0]}',${TMalign_dat[4]},${TMalign_dat[5]});
sprintf('%-40s %.2f %.4f','${TMalign_fast_dat[0]}',${TMalign_fast_dat[4]},${TMalign_fast_dat[5]});
sprintf('%-40s %.2f %.4f','${Deepalign_dat[0]}',${Deepalign_dat[4]},${Deepalign_dat[5]});
sprintf('%-40s %.2f %.4f','${DALI5_dat[0]}',${DALI5_dat[4]},${DALI5_dat[5]});
sprintf('%-40s %.2f %.4f','${Fatcat_dat[0]}',${Fatcat_dat[4]},${Fatcat_dat[5]});
sprintf('%-40s %.2f %.4f','${FS__dat[0]}',${FS__dat[4]},${FS__dat[5]});
sprintf('%-40s %.2f %.4f','${FS1_dat[0]}',${FS1_dat[4]},${FS1_dat[5]});
sprintf('%-40s %.2f %.4f','${FS0_dat[0]}',${FS0_dat[4]},${FS0_dat[5]});

ggsave(filename='${output}.pdf',plot,width=${WIDTH},height=${HEIGHT})
#dev.off();
"
echo "$Rtext" | R --vanilla

