
#include<iostream>
//#define WITH_DENSITY
#include "../src/fdtd.h"
#ifdef _OPENMP
#include <cstdlib>
#include <omp.h>
int thread_count = 1;
#endif
MyDataF epsR;
MyDataF dx, dy, dz;
MyDataF tw;
MyDataF omega;
MyDataF T; // ns
MyDataF Amp;
unsigned pmlw;
void initComData();

int main(int argc, char*argv[]) {
#ifdef _OPENMP
    cout << "OpenMP enabled." << endl;
    if (argc < 2) {
        thread_count = 3;
    } else {
        thread_count = strtol(argv[1], NULL, 10);
    }
    if (thread_count < 0 && thread_count > 100) {
        thread_count = 4;
    }
    cout << "thread count :" << thread_count << endl;
#endif
    initComData();
    fdtd hpw(500, 60, 60, 60, tw, dx, dy, dz, Amp, 10, 12, 4, 1, pmlw);
    hpw.SetSineSource(omega);
#ifdef WITH_DENSITY
    hpw.SetPlasmaVar(0, 760 * 5.3E9, 760, 0);
#endif
    //hpw.initialize();
    hpw.StartUp();
    return 0;
}

void initComData() {
    epsR = 1.0;

    pmlw = 12;

    // sine wave configure
    T = 1 / 110E9;
    omega = 2 * M_PI / T;
    dx = dy = dz = C * T / 100;
    Amp = 1e10;
    tw = 0.3 * T;

    //    // Gaussian Pulse
    //    dx = dy = dz = 1e-3;
    //    Amp = 1e10;
    //    tw = 20e-9;
    //    omega = 2 * M_PI * C / 150 / dx;
}
