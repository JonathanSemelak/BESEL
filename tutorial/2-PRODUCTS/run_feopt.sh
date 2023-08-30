#!/bin/bash

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

