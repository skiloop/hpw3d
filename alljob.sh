#!/bin/bash

# check inputs
[ -z $1 ] && exit -1;

# keep num
num=$1

# do in circle
if [ $1 -lt 14 ]; then
   sh $HOME/workspace/hpw3d/runjob.sh $1;
   sh $HOME/workspace/hpw3d/alljob.sh $((num+1));
fi

