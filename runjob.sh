#!/bin/bash


#check param

[ -z $1 ] && exit;

# prepare
num=$1
oldir=`pwd`;

echo $1;
echo $oldir;

# check inputs
[ $num -gt 13 ] && exit;

# run jobs
case $num in
   1)
      mkdir -p common && cd common && cp $HOME/workspace/hpw3d/run.sh run.sh && { sh run.sh; } && cd $oldir;
      ;;
   2)
      mkdir -p fluidDt && cd fluidDt && { chparam "ftsize" "10" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
   3)
      mkdir -p mxwllDt && cd mxwllDt && { chparam "mtsize" "80" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
   4)
      mkdir -p fluidDs && cd fluidDs && { chparam "fgsize" "8" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
   5)
      mkdir -p pres100 && cd pres100 && { chparam "pressure" "100" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
   6)
      mkdir -p pres1000 && cd pres1000 && { chparam "pressure" "1000" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
   7)
      mkdir -p mne1e8 && cd mne1e8 && { chparam "maxne" "1e8" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
   8)
      mkdir -p amp1bmv && cd amp1bmv && { chparam "amp" "1e8" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
   9)
      mkdir -p nuType0 && cd nuType0 && { chparam "nuType" "0" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
   10)
      mkdir -p nuType1 && cd nuType1 && { chparam "nuType" "1" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
   11)
      mkdir -p nuType3 && cd nuType3 && { chparam "nuType" "3" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
   12)
      mkdir -p mxwllDs && cd mxwllDs && { chparam "cellsize" "36" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
   13)
      mkdir -p mxwllDs45 && cd mxwllDs42 && { chparam "cellsize" "42" > run.sh; } && { sh run.sh; } && cd $oldir;
      ;;
esac

