#!/bin/bash
#AMBER
export AMBERHOME=/home/jota/Programas/amber20_src
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/GCC-4.9.4/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/glibc-2.14/lib
export PATH=$PATH:$AMBERHOME/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$AMBERHOME/lib
#FENEB
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/GCC-4.9.4/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/fpedron/prog/Amber18/lib64/
export PATH=$PATH:/home/jota/Programas/feneb/bin

BANDBUILDERPATH=bandbuilder   #Should be replaced by the corresponding path to bandbuilder
CPPTRAJPATH=cpptraj                             #Should be replaced by the corresponding path to cpptraj
IMAGES=15                                      #Number of images to be generated (according to bandbuilder.in)
TOPOLOGY=ALAD.prmtop              #Topology file
NAME=ALAD                                    #Prefix for all coordinates files.

#Run bandbuilder

$BANDBUILDERPATH

#Generate .nc file for further visualization
for ((i=1; i<=IMAGES; i++));
do
echo "trajin ${NAME}_r_$i.rst7" >> input.cpptraj
done
echo "trajout ${NAME}_BAND_0.nc netcdf" >> input.cpptraj
$CPPTRAJPATH $TOPOLOGY input.cpptraj
rm input.cpptraj
#Move images to the corresponding file
mkdir -p STEP-0
mv *_r_*.rst7 STEP-0/.

