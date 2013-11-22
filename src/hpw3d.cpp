
#include<iostream>
#include <cstdlib>
#ifdef _OPENMP
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

int main(int argc, char*argv[]) {
    inputChecker checker;
    checker.parseInput(argc, argv);
    checker.print();
    //return 0;

    epsR = 1.0;

    T = 1 / checker.frequency;
    switch (checker.waveType) {
        case GAUSSIAN_WAVE:
            tw = 1.5174271293851462339 / M_PI / checker.frequency;
            tw = T;
            break;
        case DERIVATIVE_GAUSSIAN_WAVE:
            break;
        default:
            tw = T;
    }

    omega = 2 * M_PI / T;
    dx = C * T / checker.yeeCellSizeX;
    dy = C * T / checker.yeeCellSizeY;
    dz = C * T / checker.yeeCellSizeZ;
    
#ifdef _OPENMP
    thread_count = checker.threadCount;
#endif
    unsigned xlen, ylen, zlen, tlen;
    unsigned minTimeLen = 2000;

    //    MyDataF dt = 0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1 / (dz * dz)));
    MyDataF dt = dx / 2 / C;
    xlen = (unsigned) (T * checker.xZoneLen * C / dx);
    ylen = (unsigned) (T * checker.yZoneLen * C / dy);
    zlen = (unsigned) (T * checker.zZoneLen * C / dz);
    if (checker.waveType == inputChecker::SINE) {
        tlen = (unsigned) (T * checker.tZoneLen / dt);
    } else {
        tlen = (unsigned) (tw * checker.tZoneLen / dt);
    }

    if (tlen < minTimeLen) {
        tlen = minTimeLen;
    }
    xlen += 2 * checker.pmlSize;
    ylen += 2 * checker.pmlSize;
    zlen += 2 * checker.pmlSize;
    cout << "xlen=" << xlen << endl;
    cout << "ylen=" << ylen << endl;
    cout << "zlen=" << zlen << endl;
    cout << "tlen=" << tlen << endl;
    cout << "dx=" << dx << endl;
    cout << "dt=" << dt << endl;

#ifdef WITH_DENSITY
    int nmaterial = 50;
    cout << "nmaterial=" << nmaterial << endl;
    fdtd hpw(tlen, xlen, ylen, zlen, tw, dx, dy, dz, checker.amptidute, 10, 12, 4, 1, checker.pmlSize, nmaterial, checker.fluidGridSize);
    hpw.SetPlasmaVar(0, 760 * 5.3E9, 760, 0);
#else
    fdtd hpw(tlen, xlen, ylen, zlen, tw, dx, dy, dz, checker.amptidute, 10, 12, 4, 1, checker.pmlSize);
#endif
    hpw.setSourceType(checker.waveType);
    switch (checker.waveType) {
        case GAUSSIAN_WAVE:break;
        case SINE_WAVE:
            hpw.SetSineSource(omega);
            break;
        case DERIVATIVE_GAUSSIAN_WAVE:break;
        case ZERO_TYPE:break;
        case ONE_SINE_PULSE:
            checker.t0 = 0.01 * T;
            checker.omega = omega;
            checker.tUp = 1.01 * T;
            checker.tDown = 0;
            hpw.intSourceSinePulse(checker.t0, checker.omega, checker.tUp, checker.tDown, checker.amptidute);
            break;
        default:
            ;
    }
    //hpw.initialize();
    hpw.StartUp();
    return 0;
}
