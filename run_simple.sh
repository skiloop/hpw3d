#!/bin/bash

wavetype=1 # select gaussian wave
pmlwidth=6 # PML layer width
cellsize=40
zonesize=1.0 # zone size in wave length
sim_time=30 # simulation time
amp=1e10 # amptitude
nthread=2 # openmp thread count
fgsize=10
ftsize=20
mtsize=100
iscon=0
maxne=1e13
rei=1e-13
pressure=760
density=1
nuType=2
freq=100e9
runpath=out
mkdir -p ${runpath}
cd ${runpath}
echo "../hpw3d --fluid-grid-size=$fgsize --wave-type=$wavetype --pml-width=$pmlwidth --zone-size=$zonesize \
	--amptidute=$amp --yee-cell-size=$cellsize  --simulation-time=$sim_time\
	--maxwell-t=$mtsize --fluid-t=$ftsize \
   	--openmp-thread=$nthread --is-connecting=$iscon --max-ne=$maxne --rei=$rei --nu-type=$nuType \
   	--use-density=$density --frequency=$freq --pressure=$pressure > r.log > /dev/null &"
../hpw3d --fluid-grid-size=$fgsize --wave-type=$wavetype --pml-width=$pmlwidth --zone-size=$zonesize \
	--amptidute=$amp --yee-cell-size=$cellsize  --simulation-time=$sim_time \
	--maxwell-t=$mtsize --fluid-t=$ftsize \
   	--openmp-thread=$nthread  --is-connecting=$iscon --max-ne=$maxne --rei=$rei --nu-type=$nuType \
   	--use-density=$density --frequency=$freq --pressure=$pressure > r.log && fetchData r.log r.dat &
cd ..



