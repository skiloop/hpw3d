
#include<iostream>
#ifdef _OPENMP
#include <cstdlib>
#include <omp.h>
int thread_count = 1;
#endif
//#define WITH_DENSITY
#include "../src/fdtd.h"
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
        thread_count = 5;
    } else {
        thread_count = strtol(argv[1], NULL, 10);
    }
    if (thread_count < 0 && thread_count > 100) {
        thread_count = 5;
    }
    cout << "thread count :" << thread_count << endl;
#endif
    unsigned xlen, ylen, zlen, tlen;
    unsigned minTimeLen = 500;
    initComData();
    MyDataF zoneLength = tw * 2 * C;
    MyDataF dt = 0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1 / (dz * dz)));
    xlen = zoneLength / dx;
    ylen = zoneLength / dy;
    zlen = zoneLength / dz;
    tlen = tw / dt;
    if (tlen < minTimeLen)tlen = minTimeLen;
    cout << "xlen=" << xlen << endl;
    cout << "ylen=" << ylen << endl;
    cout << "zlen=" << zlen << endl;
    cout << "tlen=" << tlen << endl;
    cout << "dx=" << dx << endl;
    cout << "dt=" << dt << endl;
    fdtd hpw(tlen, xlen, ylen, zlen, tw, dx, dy, dz, Amp, 10, 12, 4, 1, pmlw);
    hpw.setSourceType(fdtd::SOURCE_GAUSSIAN);
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
    //    T = 1 / 110E9;
    //    omega = 2 * pi / T;
    //    dx = dy = dz = C * T / 100;
    //    Amp = 1e10;
    //    tw = 0.3 * T;

    // Gaussian Pulse
    Amp = 1e10;
    tw = 20e-9;
    dx = dy = dz = tw * C / 50;
    omega = 2 * pi * C / 50 / dx;
}
