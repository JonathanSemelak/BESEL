#!/bin/bash
Start=1
replicas=50
inistep=1
steps=20
top=DALA.prmtop
name=DALA
mask=@5,7,9,15,17
cpptraj=/home/jota/Programas/amber20/bin/cpptraj

cp ../archivos/$top .
for ((i=Start; i<=replicas; i++));
do
        cp ../Bands/0/${name}_r_$i.rst7 ${name}_0_$i.rst7
        echo "trajin ${name}_0_$i.rst7" > input.cpptraj

for ((j=inistep; j<=steps; j++));
do
	cp ../Bands/$j/${name}_r_$i.rst7 ${name}_r_${i}_${j}.rst7
	echo "trajin ${name}_r_${i}_${j}.rst7" >> input.cpptraj
done
echo "trajout EVOL_$i.nc netcdf" >> input.cpptraj
echo "rms $mask first out RMSD_$i.agr mass" >> input.cpptraj
$cpptraj $top input.cpptraj
done
rm *.rst7
