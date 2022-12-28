#!/bin/bash
#--------------------------------------------------------------------------------------------------------------------
# This is an example script to run the bandbuilder tool. A set of structures (images) will be generated, and moved
# to a new folder named 0. A BAND.nc file with all structures will be also generated, using cpptraj.
#
# Written by Jonathan A. Semelak
#--------------------------------------------------------------------------------------------------------------------
BANDBUILDER=/home/jsemelak/Programs/feneb/bandbuilder #Should be replaced by the corresponding path to bandbuilder
CPPTRAJ=cpptraj                                       #Should be replaced by the corresponding path to cpptraj
IMAGES=15                                             #Number of images to be generated (according to bandbuilder.in)
TOPOLOGY=ALAD.prmtop                                  #Topology file
NAME=ALAD                                             #Prefix for all coordinates files.
#--------------------------------------------------------------------------------------------------------------------
#NOTES:
#[1] The directory where this script is executed must contain the corresponding bandbuilder.in, reactants and
#products .rst7 files, and $TOPOLOGY file.
#--------------------------------------------------------------------------------------------------------------------

#Run bandbuilder
$BANDBUILDER
#Generate .nc file for further visualization
for ((i=1; i<=IMAGES; i++));
do
echo "trajin ${NAME}_r_$i.rst7" >> input.cpptraj
done
echo "trajout ${NAME}_BAND_0.nc netcdf" >> input.cpptraj
$CPPTRAJ $TOPOLOGY input.cpptraj
rm input.cpptraj
#Move images to the corresponding file
mkdir -p 0
mv *_r_*.rst7 0/.
