#!/bin/bash
##$ -cwd
##$ -j y
##$ -S /bin/bash
##$ -N TUTORIAL
##$ -pe mpich 1

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

SANDERPATH=sander
FENEBPATH=feneb
STARTSTEP=1
MAXSTEPS=5
DELETENC=T
MDIN=prod.mdin

feopt_wizard -s $SANDERPATH \
             -f $FENEBPATH \
             -x $STARTSTEP \
             -d $MAXSTEPS \
             -m $MDIN \
             -g $DELETENC

