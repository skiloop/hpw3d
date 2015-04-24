#!/bin/bash

# do a series of jobs to calculate 
# high power wave simulations
#

oldpath=`pwd`;
echo $oldpath

mkdir -p common && cd common && cp $HOME/workspace/hpw3d/run.sh run.sh && { sh run.sh; } && cd $oldpath && \
mkdir -p fluidDt && cd fluidDt && { chparam "ftsize" "10" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p mxwllDt && cd mxwllDt && { chparam "mtsize" "80" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p fluidDs && cd fluidDs && { chparam "fgsize" "8" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p pres100 && cd pres100 && { chparam "pressure" "100" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p pres1000 && cd pres1000 && { chparam "pressure" "1000" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p mne1e8 && cd mne1e8 && { chparam "maxne" "1e8" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p amp1bmv && cd amp1bmv && { chparam "amp" "1e8" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p nuType0 && cd nuType0 && { chparam "nuType" "0" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p nuType1 && cd nuType1 && { chparam "nuType" "1" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p nuType3 && cd nuType3 && { chparam "nuType" "3" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p mxwllDs && cd mxwllDs && { chparam "cellsize" "42" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p mxwllDs45 && cd mxwllDs45 && { chparam "cellsize" "45" > run.sh; } && { sh run.sh; } && cd $oldpath && \
mkdir -p fluidDs20 && cd fluidDs20 && { chparam "fgsize" "20" > run.sh; } && { sh run.sh; } && cd $oldpath   # fluid grid size
