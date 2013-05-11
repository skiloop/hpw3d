#ifndef COMMON_H
#define COMMON_H

#include <math.h>
typedef double MyDataF;
typedef MyDataF* pMyDataF;
typedef pMyDataF* ppMyDataF;

#define DEFAULT_THREAD_COUNT (5)
#define DEFAULT_PML_SIZE (10)
#define DEFAULT_FLUID_GRID_SIZE (16)
#define DEFAULT_GRID_SIZE (50)
#define DEFAULT_ZONE_SIZE (1.5)
#define DEFAULT_AMPTIDUTE (1000) // in volt per meter
#define DEFAULT_FREQUENCY (110E9) // in Hz
#define DEFAULT_TIME_ZONE_LENGTH (1.0)

// source type definition
#define GAUSSIAN_WAVE_TYPE (1)
#define SINE_WAVE_TYPE (2)
#define DERIVE_GAUSSIAN_TYPE (3)
#define ZERO_TYPE (4)
#define SINE_PULSE_TYPE (5)
#define DEFAULT_WAVE_TYPE GAUSSIAN_WAVE_TYPE
#define MAX_TYPE_VALUE SINE_PULSE_TYPE

// nu type definition
#define MORROW_AND_LOWKE (1)
#define NIKONOV (2)
#define KANG (3)
#define ALI (4)
#define DEFAULT_NU_FORMAT ALI

// default value for max density
#define DEFAULT_DENSITY_MAX (0)

// some common constants
const MyDataF C = 2.99792458E8; // speed of light in free space
const MyDataF me = 9.110e-31; // electricity mass
const MyDataF e = 1.602e-19; // electricity charge
const MyDataF mu_0 = 4.0 * M_PI * 1.0E-7;
const MyDataF eps_0 = 1.0 / (C * C * mu_0);
const MyDataF M_PI_TWO=M_PI*2;

#endif
