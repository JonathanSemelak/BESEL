#! /bin/bash

##########################################################
feneb=/home/jota/000-Tesis/1-FE-NEB/3-FENEB4AMBER/2-Upgrade/4-DisminuirPaso/feneb/feneb
sander=/home/jota/Programas/amber20/bin/sander
forceintegrator=/home/jota/000-Tesis/1-FE-NEB/2-VerSiRompiAlgo/hybrid/bin/Forceintegrator
nombretop=DALA
name=DALA
Replicas=1
Start=1
MaxSteps=1000
BaseStep=0.05
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
        for ((k=2; k<Replicas; k++)); # Band loop
        do
                j=$(($i-1|bc))
                cp ../archivos/* .
                cp ../Bands/$j/* .

	        $sander -O -i prod.mdin -o prod$k.out -p $nombretop.prmtop -c ${name}_r_$k.rst7 -r ${name}_f_$k.rst7 -x ${name}_f_$k.nc -ref ${name}_r_$k.rst7
        done
	fi
        m=$(($i-1|bc))
        if [ $m == 0 ]
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
           STEP=$(grep "Step length:" ../TermProd/$m/feneb.out|awk '{print $3}')
	   DELTAA=$(grep "DeltaA:" ../TermProd/$m/feneb.out|awk '{print $2}')

        fi
        sed -i s/STEP/$STEP/g feneb.in
        sed -i s/LMFORCE/$LMFORCE/g feneb.in
        


        $feneb  # Runs FENEB

	converged=$(grep "System" feneb.out|awk '{print $3}')

        rm *.nc

#	if [ $converged == "F" ]
#        then

        mkdir $i

        for ((k=1; k<=Replicas; k++)); # Copy fles for next movement
        do
                cp ${name}_o_$k.rst7 $i/${name}_r_$k.rst7
        done
	
        mv $i ../Bands/.
        
#	fi

#        LMFORCE=$(grep "Band max force:" feneb.out|awk '{print $4}')
#        STEP=$(grep "Step length:" feneb.out|awk '{print $3}')
	echo $i $LMFORCE $STEP $DELTAA  >> ../maxstep.dat

	if [ $Replicas == 1 ]
	  then
	  $forceintegrator >> salidaintegrator.out
          awk '{print $1,$10}' atomic_work.dat > temp
	  first=$(head -n 1 temp|awk '{print $2}')
 	  awk -v f=$first '{$2 = $2 - f; print}' temp > profile.xy
        fi
         
        mkdir ../TermProd/$i
        mv * ../TermProd/$i/.

        echo "Step: "$i
        echo "System converged: "$converged
        cd ..
        fi
done




