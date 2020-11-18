#! /bin/bash
dihe=/home/jota/000-Tesis/1-FE-NEB/3-FENEB4AMBER/2-Upgrade/4-DisminuirPaso/feneb/getdihe

##########################################################
Start=1
End=2
##########################################################

rm all_bands.xy all_paths.xy
for ((i=Start; i<=End; i++)); # 
#for i in 5400;
#for i in 1 50 100 200 300 400 500 502 504 506 508 510;
	do 
	  cp archivos/getdihe.in TermProd/$i/.	
          awk '{print $1,$2}' TermProd/$i/profile.xy >> all_bands.xy
          echo " " >> all_bands.xy
	  cd TermProd/$i
	  rm phipsi.dat
	  $dihe 
	  cd ../../
          awk '{print $1}' TermProd/$i/phipsi.dat >> all_paths.xy
	  echo " " >> all_paths.xy
        done


