#! /bin/bash

##########################################################
feneb=/home/jota/000-Tesis/1-FE-NEB/3-FENEB4AMBER/3-NewMaster/feneb/feneb
sander=/home/jota/Programas/amber20/bin/sander
forceintegrator=/home/jota/000-Tesis/1-FE-NEB/2-VerSiRompiAlgo/hybrid/bin/Forceintegrator
dihe=/home/jota/000-Tesis/1-FE-NEB/3-FENEB4AMBER/2-Upgrade/4-DisminuirPaso/feneb/getdihe
nombretop=DALA
name=DALA
ReadExtrema=
Replicas=50
Start=1
MaxSteps=1
BaseStep=0.001d0
deletenc=
##########################################################

converged=F
for ((i=Start; i<=MaxSteps; i++)); # Optimization loop
        do

        if [ $converged == "F" ]
        then

                cd corriendo
        if [ $Replicas == 1 ]
        then
	k=1
	j=$(($i-1|bc))
                cp ../archivos/* .
                cp ../Bands/$j/* .

                $sander -O -i prod.mdin -o prod$k.out -p $nombretop.prmtop -c ${name}_r_$k.rst7 -r ${name}_f_$k.rst7 -x ${name}_f_$k.nc -ref ${name}_r_$k.rst7

	else
	StartRep=1
	EndRep=$Replicas
  if [ $ReadExtrema == 1 ]
    then
      StartRep=2
      EndRep=$(($Replicas-1|bc))
    fi

        for ((k=StartRep; k<=EndRep; k++)); # Band loop
        do
                j=$(($i-1|bc))
                cp ../archivos/* .
                cp ../Bands/$j/* .
	        $sander -O -i prod.mdin -o prod$k.out -p $nombretop.prmtop -c ${name}_r_$k.rst7 -r ${name}_f_$k.rst7 -x ${name}_f_$k.nc -ref ${name}_r_$k.rst7
        done
	fi
        m=$(($i-1|bc))
	      t=$(($Start-1|bc))
        if [ $m == $t ]
         then
           LMFORCE=99999
           STEP=$BaseStep
         else
           if [ $Replicas == 1 ]
           then
             LMFORCE=$(grep "Max force:" ../TermProd/$m/feneb.out|awk '{print $3}')
           else
             LMFORCE=$(grep "Band max force:" ../TermProd/$m/feneb.out|awk '{print $4}')
           fi
#           STEP=$(grep "Step length:" ../TermProd/$m/feneb.out|awk '{print $3}')
# 	   DELTAA=$(grep "DeltaA:" ../TermProd/$m/feneb.out|awk '{print $2}')

        fi
        sed -i s/STEP/$STEP/g feneb.in
        sed -i s/LMFORCE/$LMFORCE/g feneb.in
        $feneb  # Runs FENEB

	converged=$(grep "System" feneb.out|awk '{print $3}')

	if [ $deletenc == "T" ]
	then
        rm *.nc
	fi

#	if [ $converged == "F" ]
#        then


        # for ((k=1; k<=Replicas; k++)); # Copy fles for next movement
        # do
        #         cp ${name}_o_$k.rst7 $i/${name}_r_$k.rst7
        # done

        if [ -d ../Bands/$i ]
          then
          for ((k=StartRep; k<=EndRep; k++)); # Copy fles for next movement
          do
            cp ${name}_o_$k.rst7 ../Bands/$i/${name}_r_$k.rst7
          done
          else
	        mkdir ../Bands/$i


          for ((k=StartRep; k<=EndRep; k++)); # Copy fles for next movement
          do
            cp ${name}_o_$k.rst7 ../Bands/$i/${name}_r_$k.rst7
          done
	fi


#	fi

        if [ $Replicas -gt 1 ]
 	  then
	  FPERP=$(grep "Band max force:" feneb.out|awk '{print $4}')
          FNEB=$(grep "Band max fneb:" feneb.out|awk '{print $4}')
	  FSPRING=$(grep  "Band max fspringLast" feneb.out|awk '{print $4}')
	  RMSFNEB=$(grep "RMS(FNEB):" feneb.out|awk '{print $2}')
          echo $i $FPERP $FNEB $FSPRING $RMSFNEB  >> ../maxstep.dat
	fi

	if [ $Replicas == 1 ]
	  then
          LMFORCE=$(grep "Max force:" feneb.out|awk '{print $3}')
          echo $i $LMFORCE  >> ../maxstep.dat
	fi

#        STEP=$(grep "Step length:" feneb.out|awk '{print $3}')

#        echo $i $LMFORCE $STEP $DELTAA  >> ../maxstep.dat


	if [ $Replicas -gt 1 ]
	  then
	  BARRIER=$(grep Barrier feneb.out|awk '{print $2}')
	  MINPOINT=$(grep Minimum feneb.out|awk '{print $3}')
	  MAXPOINT=$(grep Maximum feneb.out|awk '{print $3}')
	  RMSD=$(cat rmsd.dat)
	  echo $i $BARRIER $MINPOINT $MAXPOINT >> ../barrier.dat
	  echo $i $RMSD >> ../rmsd.dat
	  cat profile.dat >> ../bandevol.dat
	  echo " " >> ../bandevol.dat
	  $dihe
	  cat phipsi.dat >> ../pathevol.dat
	  echo " " >> ../pathevol.dat

#	  $forceintegrator >> salidaintegrator.out
#          awk '{print $1,$7}' atomic_work.dat > temp
#	  first=$(head -n 1 temp|awk '{print $2}')
# 	  awk -v f=$first '{$2 = $2 - f; print}' temp > profile.xy
        fi

       	if [ -d ../TermProd/$i ]
          then
          mv * ../TermProd/$i/.
	else
	  mkdir ../TermProd/$i
          mv * ../TermProd/$i/.
        fi
        echo "Step: "$i
        echo "System converged: "$converged
        cd ..
        fi
done
