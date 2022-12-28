#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# This is an example script to run a free energy optimization
#
# Written by Jonathan A. Semelak
#-------------------------------------------------------------------------------------------------------------
SANDER=sander                                         #Should be replaced by the corresponding path to sander
FENEB=/home/jsemelak/Programs/feneb/feneb             #Should be replaced by the corresponding path to feneb
TOPOLOGY=ALAD.prmtop                                  #Topology file
NAME=ALAD                                             #Prefix for all coordinates files.
STARTSTEP=1                                           #Starting optimization step.
MAXSTEPS=1                                            #Maximum optimization steps to be performed.
DELETENC=1                                            #Delete .nc files after processing ("1" for True)
MDIN=prod.mdin                                        #Sander input
#-------------------------------------------------------------------------------------------------------------
#NOTES:
#[1] The directory where this script is executed must contain $TOPOLOGY, $MDIN, and feneb.in files, as well
#as a directory called called $STARTSTEP-1, with the coordinates corresponding to
#the previous step ("0", in case $STARTSTEP=1).
#This directory is automatically generated, but if this is the first step, must be manually created
#[2] The feneb code requires to "image number" to be 1 for free energy optimizations. This way, the files will
#be named including "_1_".
#-------------------------------------------------------------------------------------------------------------

for ((i=STARTSTEP; i<=MAXSTEPS; i++)); # Optimization loop
     do
        k=1 #feneb code requires to "image number" to be 1 for free energy optimizations
        mkdir -p $i #Creates a directory for the $i-th optimization step
        cd $i
        #Copy input files and topology
        cp ../$MDIN .
        cp ../feneb.in .
        cp ../$TOPOLOGY  .
        #Copy coordinates from previous optimization step ("0" if this is the first optimization step)
        j=$(($i-1|bc)) #Previous optimization step
        if [ $i == "1" ] #If this is the first optimization step, copy _r_ files to _f_
          then
          cp ../$j/${NAME}_r_${k}.rst7 .
          cp ${NAME}_r_${k}.rst7 ${NAME}_fprev_${k}.rst7
        else #Copy _f_ files from previous step and call them _fprev_
          cp ../$j/${NAME}_o_${k}.rst7 ${NAME}_r_${k}.rst7
          cp ../$j/${NAME}_f_${k}.rst7 ${NAME}_fprev_${k}.rst7
        fi
        echo "Running MD"
        $SANDER -O -i prod.mdin \
                   -o prod_${k}.out \
                   -p $TOPOLOGY \
                   -c ${NAME}_fprev_${k}.rst7 \
                   -r ${NAME}_f_${k}.rst7 \
                   -x ${NAME}_f_${k}.nc \
                   -ref ${NAME}_r_${k}.rst7
        #Run feopt
        echo "Running feneb optimization"
        $FENEB
        if [ $DELETENC == "1" ] #Delete .nc files
          then
          rm *.nc
        fi
        echo "Step: "$i " finished"
        cd ..
done
