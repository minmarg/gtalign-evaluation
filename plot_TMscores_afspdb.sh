#!/bin/bash
name=$(basename $0 .sh)
## PARAMS
SORT=besttm
SORT=querytm
## FILES
WIDTH=9.5; HEIGHT=3.2
WIDTH=7.08; HEIGHT=2.4
GTA7_speed0_pre03="gtalign_14_speed0_prescore03_addss_s044__afspdb__processed.out"
GTA7_speed9_pre03="gtalign_14_speed9_prescore03_addss_s044__afspdb__processed.out"
GTA7_speed9_pre04="gtalign_14_speed9_prescore04_addss_s044__afspdb__processed.out"
GTA7_speed13_pre03="gtalign_14_speed13_prescore03_addss_s044__afspdb__processed.out"
GTA7_speed13_pre04="gtalign_14_speed13_prescore04_addss_s044__afspdb__processed.out"
GTA7_speed13_pre03_sim15="gtalign_14_speed13_prescore03_presim15_addss_s03__afspdb__processed.out"
TMalign="tmalign__afspdb__processed.out"
TMalignFast="tmalign_fast__afspdb__processed.out"
FS_="foldseek_alntyp2__afspdb__processed.out"
FS1="foldseek_tmfast1__afspdb__processed.out"
FS0="foldseek_tmfast0__afspdb__processed.out"
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
    } 
    if($ENV{SRT} eq "A") {
      @sndx=sort{$sv[$a]<=>$sv[$b]} 0..$#sv; 
    } else {
      @sndx=sort{$sv[$b]<=>$sv[$a]} 0..$#sv; 
    }
    $y=$x="c("; 
    $sum=0; 
    for($i=0;$i<$ENV{NUM} && $i<=$#sndx;$i++){
      $last=($i+1)*0.001 if 0.5<=$sk[$sndx[$i]]; 
      $sum+=$sk[$sndx[$i]]; 
      if($i%400==0 || $ENV{NUM}<=$i+1 || $#sndx<$i+1){
        do{$y.=",";$x.=","} if $i; 
        $y.=sprintf("%.4f",$sum*0.001); 
        $x.=sprintf("%.3f",($i+1)*0.001);
      }
    } 
    $y.=")"; $x.=")"; 
    print $x, " ", $y, " ", $last, " ", $sum*0.001' $file
}
## OUTPUT
output="${name}_${SORT}"
if [ "$SORT" == "besttm" ]; then
  ## sort by own metric
  NN=312000
  GTA7_speed0_pre03_dat=("GTalign --speed=0" $(GetXYSum "../$GTA7_speed0_pre03" $NN 8 8 14 D))
  GTA7_speed9_pre03_dat=("GTalign --speed=9 (def)" $(GetXYSum "../$GTA7_speed9_pre03" $NN 8 8 14 D))
  GTA7_speed9_pre04_dat=("GTalign --speed=9 --pre-score=0.4" $(GetXYSum "../$GTA7_speed9_pre04" $NN 8 8 14 D))
  GTA7_speed13_pre03_dat=("GTalign --speed=13" $(GetXYSum "../$GTA7_speed13_pre03" $NN 8 8 14 D))
  GTA7_speed13_pre04_dat=("GTalign --speed=13 --pre-score=0.4" $(GetXYSum "../$GTA7_speed13_pre04" $NN 8 8 14 D))
  GTA7_speed13_pre03_sim15_dat=("GTalign --speed=13 --pre-similarity=15" $(GetXYSum "../$GTA7_speed13_pre03_sim15" $NN 8 8 14 D))

  TMalign_dat=("TM-align" $(GetXYSum "../$TMalign" $NN 7 7 7 D))
  TMalign_fast_dat=("TM-align -fast" $(GetXYSum "../$TMalignFast" $NN 7 7 7 D))
  FS__dat=("Foldseek" $(GetXYSum "../$FS_" $NN 5 7 13 D))
  FS1_dat=("Foldseek --tmalign-fast 1" $(GetXYSum "../$FS1" $NN 5 7 13 D))
  FS0_dat=("Foldseek --tmalign-fast 0" $(GetXYSum "../$FS0" $NN 5 7 13 D))


  ## sort by tm-score obtained by tm-align
  GTA7_speed0_pre03_tm_dat=("GTalign --speed=0" $(GetXYSum "../$GTA7_speed0_pre03" $NN 14 14 14 D) 8453.835)
  GTA7_speed9_pre03_tm_dat=("GTalign --speed=9 (def)" $(GetXYSum "../$GTA7_speed9_pre03" $NN 14 14 14 D) 1896.254)
  GTA7_speed9_pre04_tm_dat=("GTalign --speed=9 --pre-score=0.4" $(GetXYSum "../$GTA7_speed9_pre04" $NN 14 14 14 D) 899.154)
  GTA7_speed13_pre03_tm_dat=("GTalign --speed=13" $(GetXYSum "../$GTA7_speed13_pre03" $NN 14 14 14 D) 933.171)
  GTA7_speed13_pre04_tm_dat=("GTalign --speed=13 --pre-score=0.4" $(GetXYSum "../$GTA7_speed13_pre04" $NN 14 14 14 D) 618.006)
  GTA7_speed13_pre03_sim15_tm_dat=("GTalign --speed=13 --pre-similarity=15" $(GetXYSum "../$GTA7_speed13_pre03_sim15" $NN 14 14 14 D) 428.400)

  TMalign_tm_dat=("TM-align" $(GetXYSum "../$TMalign" $NN 7 7 7 D) 879965.485)
  TMalign_fast_tm_dat=("TM-align -fast" $(GetXYSum "../$TMalignFast" $NN 7 7 7 D) 109318.953)
  FS__tm_dat=("Foldseek" $(GetXYSum "../$FS_" $NN 13 13 13 D) 195.579)
  FS1_tm_dat=("Foldseek --tmalign-fast 1" $(GetXYSum "../$FS1" $NN 13 13 13 D) 1138.797)
  FS0_tm_dat=("Foldseek --tmalign-fast 0" $(GetXYSum "../$FS0" $NN 13 13 13 D) 9649.510)


elif [ "$SORT" == "querytm" ]; then
  ## sort by own metric
  NN=100000
  GTA7_speed0_pre03_dat=("GTalign --speed=0" $(GetXYSum "../$GTA7_speed0_pre03" $NN 4 4 10 D))
  GTA7_speed9_pre03_dat=("GTalign --speed=9 (def)" $(GetXYSum "../$GTA7_speed9_pre03" $NN 4 4 10 D))
  GTA7_speed9_pre04_dat=("GTalign --speed=9 --pre-score=0.4" $(GetXYSum "../$GTA7_speed9_pre04" $NN 4 4 10 D))
  GTA7_speed13_pre03_dat=("GTalign --speed=13" $(GetXYSum "../$GTA7_speed13_pre03" $NN 4 4 10 D))
  GTA7_speed13_pre04_dat=("GTalign --speed=13 --pre-score=0.4" $(GetXYSum "../$GTA7_speed13_pre04" $NN 4 4 10 D))
  GTA7_speed13_pre03_sim15_dat=("GTalign --speed=13 --pre-similarity=15" $(GetXYSum "../$GTA7_speed13_pre03_sim15" $NN 4 4 10 D))

  TMalign_dat=("TM-align" $(GetXYSum "../$TMalign" $NN 3 3 3 D))
  TMalign_fast_dat=("TM-align -fast" $(GetXYSum "../$TMalignFast" $NN 3 3 3 D))
  FS__dat=("Foldseek" $(GetXYSum "../$FS_" $NN 3 3 9 A))
  FS1_dat=("Foldseek --tmalign-fast 1" $(GetXYSum "../$FS1" $NN 3 3 9 D))
  FS0_dat=("Foldseek --tmalign-fast 0" $(GetXYSum "../$FS0" $NN 3 3 9 D))


  ## sort by tm-score obtained by tm-align
  GTA7_speed0_pre03_tm_dat=("GTalign --speed=0" $(GetXYSum "../$GTA7_speed0_pre03" $NN 10 10 10 D) 8453.835)
  GTA7_speed9_pre03_tm_dat=("GTalign --speed=9 (def)" $(GetXYSum "../$GTA7_speed9_pre03" $NN 10 10 10 D) 1896.254)
  GTA7_speed9_pre04_tm_dat=("GTalign --speed=9 --pre-score=0.4" $(GetXYSum "../$GTA7_speed9_pre04" $NN 10 10 10 D) 899.154)
  GTA7_speed13_pre03_tm_dat=("GTalign --speed=13" $(GetXYSum "../$GTA7_speed13_pre03" $NN 10 10 10 D) 933.171)
  GTA7_speed13_pre04_tm_dat=("GTalign --speed=13 --pre-score=0.4" $(GetXYSum "../$GTA7_speed13_pre04" $NN 10 10 10 D) 618.006)
  GTA7_speed13_pre03_sim15_tm_dat=("GTalign --speed=13 --pre-similarity=15" $(GetXYSum "../$GTA7_speed13_pre03_sim15" $NN 10 10 10 D) 428.400)

  TMalign_tm_dat=("TM-align" $(GetXYSum "../$TMalign" $NN 3 3 3 D) 879965.485)
  TMalign_fast_tm_dat=("TM-align -fast" $(GetXYSum "../$TMalignFast" $NN 3 3 3 D) 109318.953)
  FS__tm_dat=("Foldseek" $(GetXYSum "../$FS_" $NN 9 9 9 D) 195.579)
  FS1_tm_dat=("Foldseek --tmalign-fast 1" $(GetXYSum "../$FS1" $NN 9 9 9 D) 1138.797)
  FS0_tm_dat=("Foldseek --tmalign-fast 0" $(GetXYSum "../$FS0" $NN 9 9 9 D) 9649.510)


fi

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
  '#a43131', '#ff3535',
  '#BC6FF1', '#892CDC', '#52057B')

myltypes <- c(
  'solid', 'F1', '224282F2', 'F4448444', 'F282', '2262', ##'73',
  'solid', 'F1',
  '224282F2', 'F1', 'solid')

p1wlegend <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size = 7), legend.key.width = unit(0.8,'cm')) +
  geom_line(linewidth = 0.4) +
  geom_point(aes(shape = group)) +
  scale_shape_manual(name = '', breaks = list.names, values = seq(0,length(list.names))) + ##, guide = 'none') +
  scale_linetype_manual(name = '', breaks = list.names, values = myltypes) +
  scale_colour_manual(name = '', breaks = list.names, values = mypalette)

P1 <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  ggtitle(\"sorted by method's measure\") +
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x=expression(paste('# top hits ('%*%10^{3},')')), y=expression(paste('Cumulative TM-score ('%*%10^{3},')'))) +
  annotate('text', x=-Inf, y=Inf, hjust=-0.05, vjust=1.8, label='Swiss-Prot', size=4) +
  ##geom_line(linetype = 1, linewidth = 0.4) +
  geom_line(linewidth = 0.3) +
  scale_linetype_manual(breaks = list.names, values = myltypes) +
  scale_colour_manual(breaks = list.names, values = mypalette)
  ##scale_colour_discrete(breaks = list.names)


dat <- list(
  '${GTA7_speed0_pre03_tm_dat[0]}'=cbind(x=${GTA7_speed0_pre03_tm_dat[1]},y=${GTA7_speed0_pre03_tm_dat[2]}), 
  '${GTA7_speed9_pre03_tm_dat[0]}'=cbind(x=${GTA7_speed9_pre03_tm_dat[1]},y=${GTA7_speed9_pre03_tm_dat[2]}), 
  '${GTA7_speed9_pre04_tm_dat[0]}'=cbind(x=${GTA7_speed9_pre04_tm_dat[1]},y=${GTA7_speed9_pre04_tm_dat[2]}), 
  '${GTA7_speed13_pre03_tm_dat[0]}'=cbind(x=${GTA7_speed13_pre03_tm_dat[1]},y=${GTA7_speed13_pre03_tm_dat[2]}), 
  '${GTA7_speed13_pre04_tm_dat[0]}'=cbind(x=${GTA7_speed13_pre04_tm_dat[1]},y=${GTA7_speed13_pre04_tm_dat[2]}), 
  '${GTA7_speed13_pre03_sim15_tm_dat[0]}'=cbind(x=${GTA7_speed13_pre03_sim15_tm_dat[1]},y=${GTA7_speed13_pre03_sim15_tm_dat[2]}), 
  '${TMalign_tm_dat[0]}'=cbind(x=${TMalign_tm_dat[1]},y=${TMalign_tm_dat[2]}), 
  '${TMalign_fast_tm_dat[0]}'=cbind(x=${TMalign_fast_tm_dat[1]},y=${TMalign_fast_tm_dat[2]}), 
  '${FS__tm_dat[0]}'=cbind(x=${FS__tm_dat[1]},y=${FS__tm_dat[2]}), 
  '${FS1_tm_dat[0]}'=cbind(x=${FS1_tm_dat[1]},y=${FS1_tm_dat[2]}), 
  '${FS0_tm_dat[0]}'=cbind(x=${FS0_tm_dat[1]},y=${FS0_tm_dat[2]}));

list.names <- names(dat)
lns <- sapply(dat, nrow)
dat <- as.data.frame(do.call('rbind', dat))
dat\$group <- rep(list.names, lns)

ylimGT0 <- max(${GTA7_speed0_pre03_tm_dat[2]}) * 1.02;

vlines <- c(
  ${GTA7_speed0_pre03_tm_dat[3]},
  ${GTA7_speed9_pre03_tm_dat[3]},
  ${GTA7_speed9_pre04_tm_dat[3]},
  ${GTA7_speed13_pre03_tm_dat[3]},
  ${GTA7_speed13_pre04_tm_dat[3]},
  ${GTA7_speed13_pre03_sim15_tm_dat[3]},
  ${TMalign_tm_dat[3]},
  ${TMalign_fast_tm_dat[3]},
  ${FS__tm_dat[3]},
  ${FS1_tm_dat[3]},
  ${FS0_tm_dat[3]});

P2 <- ggplot(dat, aes(x = x, y = y, colour = group, linetype = group)) +
  ggtitle('sorted by alignment TM-score') +
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x=expression(paste('# top hits ('%*%10^{3},')')), y='') + ##expression(paste('Cumulative TM-score ('%*%10^{3},')'))) +
  geom_line(linewidth = 0.3) +
  scale_linetype_manual(breaks = list.names, values = myltypes) +
  scale_colour_manual(breaks = list.names, values = mypalette) +
  geom_segment(x = vlines[9], y = ylimGT0, xend = vlines[1], yend = ylimGT0,
    color = 'black', linewidth = 0.2, arrow = arrow(length = unit(0.03,'npc'), ends='both')) +
  geom_vline(xintercept = vlines, linetype = myltypes, color = mypalette, linewidth = 0.2)


dat <- list(
  '${GTA7_speed0_pre03_tm_dat[0]}'=cbind(x=${GTA7_speed0_pre03_tm_dat[5]},y=${GTA7_speed0_pre03_dat[4]}), 
  '${GTA7_speed9_pre03_tm_dat[0]}'=cbind(x=${GTA7_speed9_pre03_tm_dat[5]},y=${GTA7_speed9_pre03_dat[4]}), 
  '${GTA7_speed9_pre04_tm_dat[0]}'=cbind(x=${GTA7_speed9_pre04_tm_dat[5]},y=${GTA7_speed9_pre04_dat[4]}), 
  '${GTA7_speed13_pre03_tm_dat[0]}'=cbind(x=${GTA7_speed13_pre03_tm_dat[5]},y=${GTA7_speed13_pre03_dat[4]}), 
  '${GTA7_speed13_pre04_tm_dat[0]}'=cbind(x=${GTA7_speed13_pre04_tm_dat[5]},y=${GTA7_speed13_pre04_dat[4]}), 
  '${GTA7_speed13_pre03_sim15_tm_dat[0]}'=cbind(x=${GTA7_speed13_pre03_sim15_tm_dat[5]},y=${GTA7_speed13_pre03_sim15_dat[4]}), 
  '${TMalign_tm_dat[0]}'=cbind(x=${TMalign_tm_dat[5]},y=${TMalign_dat[4]}), 
  '${TMalign_fast_tm_dat[0]}'=cbind(x=${TMalign_fast_tm_dat[5]},y=${TMalign_fast_dat[4]}), 
  '${FS__tm_dat[0]}'=cbind(x=${FS__tm_dat[5]},y=${FS__dat[4]}), 
  '${FS1_tm_dat[0]}'=cbind(x=${FS1_tm_dat[5]},y=${FS1_dat[4]}), 
  '${FS0_tm_dat[0]}'=cbind(x=${FS0_tm_dat[5]},y=${FS0_dat[4]}));

list.names <- names(dat)
dat <- as.data.frame(do.call('rbind', dat))
dat\$group <- rep(list.names, 1)

P3 <- ggplot(dat, aes(x = x, y = y, colour = group)) +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'none', panel.grid.minor = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8)) +
  labs(x=expression(paste('Runtime (sec)')), y='') + ##expression(paste('Cumulative TM-score ('%*%10^{3},')'))) +
  ylim(0,NA) +
  geom_point(aes(shape = group), size = 1.4) + #, stroke = 1) +
  annotation_logticks(sides = 'b', short = unit(0.003,'npc'), mid = unit(0.006,'npc'), long = unit(0.01,'npc')) +
  scale_shape_manual(breaks = list.names, values = c(seq(0,length(list.names)-4),seq(length(list.names),length(list.names)+3))) +
  scale_colour_manual(breaks = list.names, values = mypalette) +
  scale_x_continuous(trans = log10_trans(),
    breaks = trans_breaks('log10', function(x) 10^x),
    labels = trans_format('log10', math_format(10^.x)))

# extract legend from ggplot (credit: Joachim Schork)
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1\$grobs, function(x) x\$name) == 'guide-box')
  step3 <- step1\$grobs[[step2]]
  return(step3)
}

legend <- extract_legend(p1wlegend)

plot <- grid.arrange(arrangeGrob(P1, P2, P3, ncol=3), nrow=1);

ggsave(filename='${output}.pdf',plot,width=${WIDTH},height=${HEIGHT})
#dev.off();
"
echo "$Rtext" | R --vanilla

