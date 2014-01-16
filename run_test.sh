#!/bin/bash

wavetype=1 # select gaussian wave
pmlwidth=8 # PML layer width
yeecellsize=30
zonesize=1.0 # zone size in wave length
sim_time=30 # simulation time
amp=5e13 # amptitude
nthread=2 # openmp thread count
fgsize=6
iscon=0
maxne=1e18
rei=0 #1e-30
useDensity=1
nuType=3
freq=100e9
echo "../hpw3d --fluid-grid-size=$fgsize --wave-type=$wavetype --pml-width=$pmlwidth --zone-size=$zonesize \
	--amptidute=$amp --yee-cell-size=$yeecellsize  --simulation-time=$sim_time\
   	--openmp-thread=$nthread --is-connecting=$iscon --max-ne=$maxne --rei=$rei --nu-type=$nuType \
   	--use-density=$useDensity --frequency=$freq | tee r.log > /dev/null &"
../hpw3d --fluid-grid-size=$fgsize --wave-type=$wavetype --pml-width=$pmlwidth --zone-size=$zonesize \
	--amptidute=$amp --yee-cell-size=$yeecellsize  --simulation-time=$sim_time \
   	--openmp-thread=$nthread  --is-connecting=$iscon --max-ne=$maxne --rei=$rei --nu-type=$nuType \
   	--use-density=$useDensity --frequency=$freq | tee r.log > /dev/null &


