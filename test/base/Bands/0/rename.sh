#!/bin/bash
Start=1
replicas=25
top=DCE.prmtop
name=DCE
cpptraj=/home/jota/Programas/amber20/bin/cpptraj

for ((i=Start; i<=replicas; i++)); # Optimization loop
do
	j=$((-180+($i-1)*5))
 cp DCE_r_$i.rst7 toamber/DCE-$j.rst7
 echo $j


done
