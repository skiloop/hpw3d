
#include<iostream>
//#define WITH_DENSITY
#include "../src/fdtd.h"
#ifdef _OPENMP
#include <cstdlib>
#include <omp.h>
int thread_count = 1;
#endif
MyDataF eps_0, epsR;
MyDataF mu_0;
MyDataF dx, dy, dz;
MyDataF pi, C;
MyDataF me, e;
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
    pi = 3.14159265358979;
    C = 2.99792458E8;
    mu_0 = 4.0 * pi * 1.0E-7;
    eps_0 = 1.0 / (C * C * mu_0);
    epsR = 1.0;

    me = 9.110e-31;
    e = 1.602e-19;
    pmlw = 12;

    // sine wave configure
    T = 1 / 110E9;
    omega = 2 * pi / T;
    dx = dy = dz = C * T / 100;
    Amp = 1e10;
    tw = 0.3 * T;

    //    // Gaussian Pulse
    //    dx = dy = dz = 1e-3;
    //    Amp = 1e10;
    //    tw = 20e-9;
    //    omega = 2 * pi * C / 150 / dx;
}
