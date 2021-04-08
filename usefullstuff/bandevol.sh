#! /bin/bash
dihe=/home/jota/000-Tesis/1-FE-NEB/3-FENEB4AMBER/2-Upgrade/4-DisminuirPaso/feneb/getdihe

##########################################################
Start=1
End=100
##########################################################

rm all_bands.xy all_paths.xy
for ((i=Start; i<=End; i++)); # 
#for i in  50;
#for i in 270 271 272 273 274 275 276;
	do 
	  cp archivos/getdihe.in TermProd/$i/.	
          cat TermProd/$i/profile.dat >> all_bands.xy
          echo " " >> all_bands.xy
	  cd TermProd/$i
	  rm phipsi.dat
	  $dihe 
	  cd ../../
	  cat TermProd/$i/phipsi.dat >> all_paths.xy
	  echo " " >> all_paths.xy
        done	  

