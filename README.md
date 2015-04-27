hpw3d
=====
14.01.13 source updating modified but no effect changes

# USAGE
./hpw3d [options]

# options discriptions
	--help,-h	help
	--openmp-thread=n	number of threads for openmp, n must be positive
	--wave-type=[1,2]	source wave form type,1 for gaussian pulse,2 for sine wave
	--fluid-grid-size=	how many fluid grid per Maxwell grid 
	--maxwell-t=	how many Maxwell time steps per T
	--fluid-t=	how many fluid time steps per T
# Density parameters:
	--max-ne=	set Ne max
	--rei=	set recombinition coefficient
	--pressure=	set air pressure
	--nu-type=	set nu type
	--use-density=	whether use density or not
		1	use density
		0	DO NOT use density(default)

# FDTD parameters:
	--pml-width=n	pml size for fdtd
	--yc-size-x=	how many yee cells per wave length of pulse length in x direction
	--yc-size-y=	how many yee cells per wave length of pulse length in y direction
	--yc-size-z=	how many yee cells per wave length of pulse length in z direction
	--yee-cell-size=	how many yee cells per wave length of pulse length,setting yee cell cube
# Gaussian pulse parameters:
	--amptidute=	amptidute for source wave
	--frequency=	wave frequency
	--x-zone-length=	zone length in x-direction,in pulse width
	--y-zone-length=	zone length in y-direction,in pulse width
	--z-zone-length=	zone length in z-direction,in pulse width
	--zone-size=	zone size,in pulse width,set x,y,z zone length together and ,setting yee cell cube
	--simulation-time=	simulation time,in pulse width size
# Sine wave parameters:
	--amptidute=	amptidute for source wave
	--frequency=	wave frequency
	--x-zone-length=	zone length in x-direction,in wave length
	--y-zone-length=	zone length in y-direction,in wave length
	--z-zone-length=	zone length in z-direction,in wave length
	--zone-size=	zone size,in wave length,set x,y,z zone length together,setting yee cell cube
	--simulation-time=	simulation time,in wave length size
	--is-connecting=	1 if use connecting interface 

