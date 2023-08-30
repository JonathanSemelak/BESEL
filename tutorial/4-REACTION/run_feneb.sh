#!/bin/bash

SANDERPATH=sander
FENEBPATH=feneb
STARTSTEP=1
MAXSTEPS=10
DELETENC=T
NCPU=2
MDIN=prod.mdin

feneb_wizard -s $SANDERPATH \
             -f $FENEBPATH \
             -x $STARTSTEP \
             -d $MAXSTEPS \
             -c $NCPU \
             -m $MDIN \

