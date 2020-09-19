#! /bin/bash

##########################################################
feneb=/home/jota/000-Tesis/1-FE-NEB/3-FENEB4AMBER/2-Upgrade/2-ArmarTool/feneb/feneb
sander=/home/jota/Programas/amber20/bin/sander
nombretop=DCE
Replicas=25
Start=1
MaxSteps=50

##########################################################

converged=F
for ((i=Start; i<=MaxSteps; i++)); # Optimization loop
        do

        if [ $converged == "F" ]
        then

                cd corriendo
        for ((k=2; k<Replicas; k++)); # Band loop
        do
                j=$(($i-1|bc))
                cp ../archivos/* .
                cp ../Bands/$j/* .

                $sander -O -i prod.mdin -o prod$k.out -p $nombretop.prmtop -c DCE_r_$k.rst7 -r DCE_f_$k.rst7 -x DCE_f_$k.nc -ref DCE_r_$k.rst7
        done

        $feneb  # Runs FENEB

        rm *.nc
        mkdir $i

        for ((k=1; k<=Replicas; k++)); # Copy fles for next movement
        do
                cp DCE_o_$k.rst7 $i/DCE_r_$k.rst7
        done
        mv $i ../Bands/.

        converged=$(grep "System" feneb.out|awk '{print $3}')

/home/jota/000-Tesis/1-FE-NEB/2-VerSiRompiAlgo/hybrid/bin/Forceintegrator >> salidaintegrator.out
        awk '{print $1,$10}' atomic_work.dat > profile.xy

        mkdir ../TermProd/$i
        mv * ../TermProd/$i/.

        echo "Step: "$i
        echo "System converged: "$converged
        cd ..
        fi
done




