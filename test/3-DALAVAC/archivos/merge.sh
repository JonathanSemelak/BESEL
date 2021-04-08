#!/bin/bash
Start=1
replicas=25
top=DALA.prmtop
name=DALA
cpptraj=/home/jota/Programas/amber20/bin/cpptraj

echo "trajin $name"
for ((i=Start; i<=replicas; i++)); # Optimization loop
do
	echo "trajin ${name}_o_$i.rst7" >> input.cpptraj
done
echo "trajout BAND.nc netcdf" >> input.cpptraj
$cpptraj $top input.cpptraj
rm input.cpptraj
