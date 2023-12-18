#!/bin/bash
name=$(basename $0 .sh)
## PARAMS
WIDTH=9.5; HEIGHT=9.5
output="${name}"
export GTA7_speed0_pre03="gtalign_14_speed0__homstrad__processed.out"
export GTA7_speed9_pre03="gtalign_14_speed9__homstrad__processed.out"
export GTA7_speed13_pre03="gtalign_14_speed13__homstrad__processed.out"
export TMalign="tmalign__homstrad__processed.out"
export DeepAlign="deepalign__homstrad__processed.out"
export DALI5="DALIv5__homstrad__processed.out"
export Fatcat="fatcat__homstrad__processed.out"
export FS_="foldseek_alntyp2__homstrad__processed.out"
export FS0="foldseek__homstrad__processed.out"
export Reference="reference__homstrad__processed.out"
TMPDATFILE=data_compiled.lst
perl -e '
open(F,"../$ENV{GTA7_speed0_pre03}") or die "ERROR: $ENV{GTA7_speed0_pre03}";
while(<F>){@a = split(/\s+/); $H{$a[0]}{GTA0} = $a[12];} close(F);
open(F,"../$ENV{GTA7_speed9_pre03}") or die "ERROR: $ENV{GTA7_speed9_pre03}";
while(<F>){@a = split(/\s+/); $H{$a[0]}{GTA9} = $a[12];} close(F);
open(F,"../$ENV{GTA7_speed13_pre03}") or die "ERROR: $ENV{GTA7_speed13_pre03}";
while(<F>){@a = split(/\s+/); $H{$a[0]}{GTA13} = $a[12];} close(F);
open(F,"../$ENV{TMalign}") or die "ERROR: $ENV{TMalign}";
while(<F>){@a = split(/\s+/); $H{$a[0]}{TMAL} = $a[6];} close(F);
open(F,"../$ENV{DeepAlign}") or die "ERROR: $ENV{DeepAlign}";
while(<F>){@a = split(/\s+/); $H{$a[0]}{DPAL} = $a[10];} close(F);
open(F,"../$ENV{DALI5}") or die "ERROR: $ENV{DALI5}";
while(<F>){@a = split(/\s+/); $H{$a[0]}{DALI} = $a[8];} close(F);
open(F,"../$ENV{Fatcat}") or die "ERROR: $ENV{Fatcat}";
while(<F>){@a = split(/\s+/); $H{$a[0]}{FTCT} = $a[10];} close(F);
open(F,"../$ENV{FS_}") or die "ERROR: $ENV{FS_}";
while(<F>){@a = split(/\s+/); $H{$a[0]}{FS_} = $a[12];} close(F);
open(F,"../$ENV{FS0}") or die "ERROR: $ENV{FS0}";
while(<F>){@a = split(/\s+/); $H{$a[0]}{FS0} = $a[12];} close(F);
open(F,"../$ENV{Reference}") or die "ERROR: $ENV{Reference}";
while(<F>){@a = split(/\s+/); $H{$a[0]}{REFN} = $a[6];} close(F);
foreach $k(keys %H) {
  $H{$k}{GTA0} = "0" unless exists $H{$k}{GTA0};
  $H{$k}{GTA9} = "0" unless exists $H{$k}{GTA9};
  $H{$k}{GTA13} = "0" unless exists $H{$k}{GTA13};
  $H{$k}{TMAL} = "0" unless exists $H{$k}{TMAL};
  $H{$k}{DPAL} = "0" unless exists $H{$k}{DPAL};
  $H{$k}{DALI} = "0" unless exists $H{$k}{DALI};
  $H{$k}{FTCT} = "0" unless exists $H{$k}{FTCT};
  $H{$k}{FS_} = "0" unless exists $H{$k}{FS_};
  $H{$k}{FS0} = "0" unless exists $H{$k}{FS0};
  $H{$k}{REFN} = "0" unless exists $H{$k}{REFN};
  printf("%-36s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",$k,$H{$k}{GTA0},$H{$k}{GTA9},$H{$k}{GTA13},$H{$k}{TMAL}, $H{$k}{DPAL},$H{$k}{DALI},$H{$k}{FTCT}, $H{$k}{FS_},$H{$k}{FS0},$H{$k}{REFN});
}' >${TMPDATFILE}
Rtext="
library(ggplot2);
library(scales);
library(gridExtra);
##pdf('${output}.pdf',family='sans',width=3,height=3,pointsize=1);
Dlst <- scan('${TMPDATFILE}',comment.char='#',what=list(name=character(),gta0=double(),gta9=double(),gta13=double(),tmal=double(), dpal=double(),dali=double(),ftct=double(), fs_=double(),fs0=double(),refn=double()));
N <- length(Dlst\$name);
colors <- rep('steelblue', N);

gta0sum <- sprintf('%.2f',sum(Dlst\$gta0));
gta9sum <- sprintf('%.2f',sum(Dlst\$gta9));
gta13sum <- sprintf('%.2f',sum(Dlst\$gta13));
tmalsum <- sprintf('%.2f',sum(Dlst\$tmal));
dpalsum <- sprintf('%.2f',sum(Dlst\$dpal));
dalisum <- sprintf('%.2f',sum(Dlst\$dali));
ftctsum <- sprintf('%.2f',sum(Dlst\$ftct));
fs_sum <- sprintf('%.2f',sum(Dlst\$fs_));
fs0sum <- sprintf('%.2f',sum(Dlst\$fs0));
refnsum <- sprintf('%.2f',sum(Dlst\$refn));

dalimin0 <- min(Dlst\$dali[Dlst\$dali>0]);

annsize <- 3.5;

df <- data.frame(gta0=Dlst\$gta0,gta9=Dlst\$gta9,gta13=Dlst\$gta13,tmal=Dlst\$tmal, dpal=Dlst\$dpal,dali=Dlst\$dali,ftct=Dlst\$ftct, fs_=Dlst\$fs_,fs0=Dlst\$fs0,refn=Dlst\$refn,name=Dlst\$name,cols=colors);

P1 <- ggplot(data=df, mapping=aes(x=gta9, y=gta0)) +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'none', axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  geom_abline(intercept=0, size=0.1) +
  geom_point(stat='identity', shape=21, size=2, color='black', fill=colors, alpha=0.7) +
  annotate('text', x=-Inf, y=Inf, hjust=-0.05, vjust=1.8, label=paste0('Cumulative TM-score = ',gta0sum), size=annsize) +
  annotate('text', x=Inf, y=-Inf, hjust=1.02, vjust=-0.5, label=paste0('Cumulative TM-score = ',gta9sum), size=annsize) +
##  scale_x_continuous(limits = c(0.2, 1.0)) +
  labs(x = 'GTalign --speed=9', y = 'GTalign --speed=0')

P2 <- ggplot(data=df, mapping=aes(x=gta13, y=gta0)) +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'none', axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  geom_abline(intercept=0, size=0.1) +
  geom_point(stat='identity', shape=21, size=2, color='black', fill=colors, alpha=0.7) +
  annotate('text', x=-Inf, y=Inf, hjust=-0.05, vjust=1.8, label=paste0('Cumulative TM-score = ',gta0sum), size=annsize) +
  annotate('text', x=Inf, y=-Inf, hjust=1.02, vjust=-0.5, label=paste0('Cumulative TM-score = ',gta13sum), size=annsize) +
##  scale_x_continuous(limits = c(0.2, 1.0)) +
  labs(x = 'GTalign --speed=13', y = '') ##'GTalign --speed=0')

P3 <- ggplot(data=df, mapping=aes(x=tmal, y=gta0)) +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'none', axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  geom_abline(intercept=0, size=0.1) +
  geom_point(stat='identity', shape=21, size=2, color='black', fill=colors, alpha=0.7) +
  annotate('text', x=-Inf, y=Inf, hjust=-0.05, vjust=1.8, label=paste0('Cumulative TM-score = ',gta0sum), size=annsize) +
  annotate('text', x=Inf, y=-Inf, hjust=1.02, vjust=-0.5, label=paste0('Cumulative TM-score = ',tmalsum), size=annsize) +
##  scale_x_continuous(limits = c(0.2, 1.0)) +
  labs(x = 'TM-align', y = '') ##'GTalign --speed=0')


P4 <- ggplot(data=df, mapping=aes(x=dpal, y=gta0)) +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'none', axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  geom_abline(intercept=0, size=0.1) +
  geom_point(stat='identity', shape=21, size=2, color='black', fill=colors, alpha=0.7) +
  annotate('text', x=-Inf, y=Inf, hjust=-0.05, vjust=1.8, label=paste0('Cumulative TM-score = ',gta0sum), size=annsize) +
  annotate('text', x=Inf, y=-Inf, hjust=1.02, vjust=-0.5, label=paste0('Cumulative TM-score = ',dpalsum), size=annsize) +
  labs(x = 'DeepAlign', y = 'GTalign --speed=0')

P5 <- ggplot(data=df, mapping=aes(x=dali, y=gta0)) +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'none', axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  geom_abline(intercept=0, size=0.1) +
  geom_point(stat='identity', shape=21, size=2, color='black', fill=colors, alpha=0.7) +
  annotate('text', x=-Inf, y=Inf, hjust=-0.05, vjust=1.8, label=paste0('Cumulative TM-score = ',gta0sum), size=annsize) +
  annotate('text', x=Inf, y=-Inf, hjust=1.02, vjust=-0.5, label=paste0('Cumulative TM-score = ',dalisum), size=annsize) +
  scale_x_continuous(limits = c(dalimin0, 1.0)) +
  labs(x = 'DALI', y = '') ##'GTalign --speed=0')

P6 <- ggplot(data=df, mapping=aes(x=ftct, y=gta0)) +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'none', axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  geom_abline(intercept=0, size=0.1) +
  geom_point(stat='identity', shape=21, size=2, color='black', fill=colors, alpha=0.7) +
  annotate('text', x=-Inf, y=Inf, hjust=-0.05, vjust=1.8, label=paste0('Cumulative TM-score = ',gta0sum), size=annsize) +
  annotate('text', x=Inf, y=-Inf, hjust=1.02, vjust=-0.5, label=paste0('Cumulative TM-score = ',ftctsum), size=annsize) +
  labs(x = 'FATCAT', y = '') ##'GTalign --speed=0')


P7 <- ggplot(data=df, mapping=aes(x=fs0, y=gta0)) +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'none', axis.text=element_text(size=12), axis.title=element_text(size=12), axis.title.x = element_text(hjust=0.0)) +
  geom_abline(intercept=0, size=0.1) +
  geom_point(stat='identity', shape=21, size=2, color='black', fill=colors, alpha=0.7) +
  annotate('text', x=-Inf, y=Inf, hjust=-0.05, vjust=1.8, label=paste0('Cumulative TM-score = ',gta0sum), size=annsize) +
  annotate('text', x=Inf, y=-Inf, hjust=1.02, vjust=-0.5, label=paste0('Cumulative TM-score = ',fs0sum), size=annsize) +
  labs(x = 'Foldseek --tmalign-fast 0', y = 'GTalign --speed=0')

P8 <- ggplot(data=df, mapping=aes(x=fs_, y=gta0)) +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'none', axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  geom_abline(intercept=0, size=0.1) +
  geom_point(stat='identity', shape=21, size=2, color='black', fill=colors, alpha=0.7) +
  annotate('text', x=-Inf, y=Inf, hjust=-0.05, vjust=1.8, label=paste0('Cumulative TM-score = ',gta0sum), size=annsize) +
  annotate('text', x=Inf, y=-Inf, hjust=1.02, vjust=-0.5, label=paste0('Cumulative TM-score = ',fs_sum), size=annsize) +
  labs(x = 'Foldseek', y = '') ##'GTalign --speed=0')

P9 <- ggplot(data=df, mapping=aes(x=refn, y=gta0)) +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'none', axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  geom_abline(intercept=0, size=0.1) +
  geom_point(stat='identity', shape=21, size=2, color='black', fill=colors, alpha=0.7) +
  annotate('text', x=-Inf, y=Inf, hjust=-0.05, vjust=1.8, label=paste0('Cumulative TM-score = ',gta0sum), size=annsize) +
  annotate('text', x=Inf, y=-Inf, hjust=1.02, vjust=-0.5, label=paste0('Cumulative TM-score = ',refnsum), size=annsize) +
  labs(x = 'Reference', y = '') ##'GTalign --speed=0')

plot <- grid.arrange(arrangeGrob(P1, P2, P3, P4, P5, P6, P7, P8, P9, ncol=3), nrow=2, heights = c(10, 0));

ggsave(filename='${output}.pdf',plot,width=${WIDTH},height=${HEIGHT})
#dev.off();
"
echo "$Rtext" | R --vanilla

