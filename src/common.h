#ifndef COMMON_H
#define COMMON_H
#define _USE_MATH_DEFINES
#include <cmath>
typedef long double MyDataF;
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
#define GAUSSIAN_WAVE (1)
#define SINE_WAVE (2)
#define DERIVATIVE_GAUSSIAN_WAVE (3)
#define COSINE_GAUSSIAN_WAVE (4)
#define ZERO_TYPE (5)
#define ONE_SINE_PULSE (6)
#define SQUARE_PULSE (7)
#define TEST_SOURCE (8)
#define DEFAULT_WAVE_TYPE GAUSSIAN_WAVE
#define MAX_TYPE_VALUE TEST_SOURCE

// nu type definition
#define MORROW_AND_LOWKE (1)
#define NIKONOV (2)
#define KANG (3)
#define ALI (4)
#define DEFAULT_NU_FORMAT ALI

// default value for max density
#define DEFAULT_DENSITY_MAX (1e13)

// some common constants
const MyDataF C = 2.99792458E8; // speed of light in free space
const MyDataF me = 9.110e-31; // electricity mass
const MyDataF e = 1.602e-19; // electricity charge
const MyDataF mu_0 = 4.0 * M_PI * 1.0E-7;
const MyDataF eps_0 = 1.0 / (C * C * mu_0);
const MyDataF M_PI_TWO = M_PI * 2;

// default zone parameters
#define AIR_BUFFER 6
#define NE_BOUND_WIDTH 6

// use density or not
#define USE_DENSITY     1
#define NOT_USE_DENSITY 0
#define IF_USE_DENSITY (NOT_USE_DENSITY)

//#define WITH_DENSITY

//#define WITH_DENSITY
#endif /*  COMMON_H */
