#!/bin/bash
export AMBERHOME=/share/apps/amber14
PATH=$PATH:$AMBERHOME/bin
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$AMBERHOME/lib

#Agrego
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/GCC-4.9.4/lib64
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/fpedron/prog/Amber18/lib64/

###############################################
export PATH
export LD_LIBRARY_PATH

feneb=/home/jota/000-Tesis/1-FE-NEB/3-FENEB4AMBER/2-Upgrade/7-OKVAC/feneb/feneb


	$feneb 



