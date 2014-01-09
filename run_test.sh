#!/bin/bash

wavetype=1 # select gaussian wave
pmlwidth=8 # PML layer width
yeecellsize=30
zonesize=1.0 # zone size in wave length
sim_time=130 # simulation time
amp=1e7 # amptitude
nthread=2 # openmp thread count
fgsize=6
iscon=0
maxne=1e18
rei=0
useDensity=0
nuType=3
echo "../hpw3d --fluid-grid-size=$fgsize --wave-type=$wavetype --pml-width=$pmlwidth --zone-size=$zonesize \
	--amptidute=$amp --yee-cell-size=$yeecellsize  --simulation-time=$sim_time\
   	--openmp-thread=$nthread --is-connecting=$iscon --max-ne=$maxne --rei=$rei --nu-type=$nuType \
   	--use-density=$useDensity | tee r.log > /dev/null &"
../hpw3d --fluid-grid-size=$fgsize --wave-type=$wavetype --pml-width=$pmlwidth --zone-size=$zonesize \
	--amptidute=$amp --yee-cell-size=$yeecellsize  --simulation-time=$sim_time \
   	--openmp-thread=$nthread  --is-connecting=$iscon --max-ne=$maxne --rei=$rei --nu-type=$nuType \
   	--use-density=$useDensity | tee r.log > /dev/null &

