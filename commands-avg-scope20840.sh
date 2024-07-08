## average ROC1 scores
ls -1 *__scope20840__processed__foldrec.out0|xargs -i awk '{if($7>1){f+=$3/$10;cf++;} if($8-$7>0){s+=$4/$11;cs++;} if($9-$8>0){fl+=$5/$12;cfl++;}}END{print f,s,fl,cf,cs,cfl,FILENAME}' {}|sort -k7

ls -1 *__scope20840__processed__foldrec.out1|xargs -i awk '{if($7>1){f+=$3/$10;cf++;} if($8-$7>0){s+=$4/$11;cs++;} if($9-$8>0){fl+=$5/$12;cfl++;}}END{print f,s,fl,cf,cs,cfl,FILENAME}' {}|sort -k7

## calculate significance between ROC1 data points
Rtext=""
for CRF in 0 1; do
  Rtext="${Rtext}
    sprintf('%s\\n','Cross-fold relationships ignored= ${CRF}');
  "
  gtafile=gtalign_15_speed0_prescore03_addss_s044__scope20840__processed__foldrec.out${CRF}
  gtafam="$(awk '{if($7>1) print $3 ","}' ${gtafile})"; gtafam=${gtafam%,*}; gtafam="c("$gtafam")"
  gtasfam="$(awk '{if($8-$7>0) print $4 ","}' ${gtafile})"; gtasfam=${gtasfam%,*}; gtasfam="c("$gtasfam")"
  gtafold="$(awk '{if($9-$8>0) print $5 ","}' ${gtafile})"; gtafold=${gtafold%,*}; gtafold="c("$gtafold")"
  for file in gtalign_15_speed9_prescore03_addss_s044__scope20840__processed__foldrec.out${CRF} gtalign_15_speed9_prescore04_addss_s044__scope20840__processed__foldrec.out${CRF} gtalign_15_speed13_prescore03_addss_s044__scope20840__processed__foldrec.out${CRF} gtalign_15_speed13_prescore04_addss_s044__scope20840__processed__foldrec.out${CRF} gtalign_15_speed13_prescore03_presim15_addss_s03__scope20840__processed__foldrec.out${CRF} tmalign__scope20840__processed__foldrec.out${CRF} tmalign_fast__scope20840__processed__foldrec.out${CRF} deepalign__scope20840__processed__foldrec.out${CRF} DALIv5__scope20840__processed__foldrec.out${CRF} fatcat__scope20840__processed__foldrec.out${CRF} foldseek_alntyp2__scope20840__processed__foldrec.out${CRF} foldseek_tmfast1__scope20840__processed__foldrec.out${CRF} foldseek_tmfast0__scope20840__processed__foldrec.out${CRF}; do
    filefam="$(awk '{if($7>1) print $3 ","}' ${file})"; filefam=${filefam%,*}; filefam="c("$filefam")"
    filesfam="$(awk '{if($8-$7>0) print $4 ","}' ${file})"; filesfam=${filesfam%,*}; filesfam="c("$filesfam")"
    filefold="$(awk '{if($9-$8>0) print $5 ","}' ${file})"; filefold=${filefold%,*}; filefold="c("$filefold")"
    Rtext="$Rtext
      xfam<-${gtafam}; yfam<-${filefam};
      xsfam<-${gtasfam}; ysfam<-${filesfam};
      xfold<-${gtafold}; yfold<-${filefold};
      fam<-wilcox.test(xfam, yfam, paired=FALSE);
      sfam<-wilcox.test(xsfam, ysfam, paired=FALSE);
      fold<-wilcox.test(xfold, yfold, paired=FALSE);
      sprintf('%10.4g %10.4g %10.4g  %-80s %-60s', fam\$p.value,sfam\$p.value,fold\$p.value,'${gtafile}','${file}');
  "
  done
done
(echo "$Rtext" | R --vanilla --slave)


