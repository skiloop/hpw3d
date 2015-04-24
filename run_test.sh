#!/bin/bash

wavetype=2 # select gaussian wave
pmlwidth=4 # PML layer width
yeecellsize=50
zonesize=1.2 # zone size in wave length
sim_time=50 # simulation time
amp=5.5e7 # amptitude
nthread=2 # openmp thread count
fgsize=6
iscon=1
rei=1e-13
useDensity=1
nuType=2
echo "../hpw3d --fluid-grid-size=$fgsize --wave-type=$wavetype --pml-width=$pmlwidth --zone-size=$zonesize \
      --rei=$rei --nu-type=$nuType \
	--amptidute=$amp --yee-cell-size=$yeecellsize  --simulation-time=$sim_time\
   	--openmp-thread=$nthread --is-connecting=$iscon --use-density=$useDensity | tee r.log > /dev/null &"
../hpw3d --fluid-grid-size=$fgsize --wave-type=$wavetype --pml-width=$pmlwidth --zone-size=$zonesize \
      --rei=$rei --nu-type=$nuType \
	--amptidute=$amp --yee-cell-size=$yeecellsize  --simulation-time=$sim_time \
   	--openmp-thread=$nthread  --is-connecting=$iscon --use-density=$useDensity | tee r.log > /dev/null &

