
#include<iostream>
#ifdef _OPENMP
#include <cstdlib>
#include <omp.h>
int thread_count = 1;
#endif
//#define WITH_DENSITY
#include "fdtd.h"
#include "inputChecker.h"

MyDataF epsR;
MyDataF dx, dy, dz;
MyDataF tw;
MyDataF omega;
MyDataF T; // ns
MyDataF Amp;
unsigned pmlw;

int main(int argc, char*argv[]) {
    inputChecker checker;
    checker.parseInput(argc, argv);
    checker.print();
    //return 0;

    epsR = 1.0;

    // sine wave configure
    T = 1 / checker.frequency;
    omega = 2 * M_PI / T;
    dx = C * T / checker.yeeCellSizeX;
    dy = C * T / checker.yeeCellSizeY;
    dz = C * T / checker.yeeCellSizeZ;
    Amp = checker.amptidute;
    tw = 0.3 * T;
    unsigned xlen, ylen, zlen, tlen;
    unsigned minTimeLen = 500;

    MyDataF dt = 0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1 / (dz * dz)));
    xlen = checker.xZoneLen * C / dx;
    ylen = checker.yZoneLen * C / dy;
    zlen = checker.zZoneLen * C / dz;
    if (checker.waveType == inputChecker::SINE) {
        tlen = T * checker.tZoneLen / dt;
    } else {
        tlen = tw * checker.tZoneLen / dt;
    }

    if (tlen < minTimeLen) {
        tlen = minTimeLen;
    }
    cout << "xlen=" << xlen << endl;
    cout << "ylen=" << ylen << endl;
    cout << "zlen=" << zlen << endl;
    cout << "tlen=" << tlen << endl;
    cout << "dx=" << dx << endl;
    cout << "dt=" << dt << endl;
    fdtd hpw(tlen, xlen, ylen, zlen, tw, dx, dy, dz, Amp, 10, 12, 4, 1, checker.pmlSize);
    hpw.setSourceType(checker.waveType);
#ifdef WITH_DENSITY
    hpw.SetPlasmaVar(0, 760 * 5.3E9, 760, 0);
#endif
    //hpw.initialize();
    hpw.StartUp();
    return 0;
}
