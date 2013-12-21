
#include<iostream>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
int thread_count = 1;
#endif
//#define WITH_DENSITY

#include "fdtd.h"
#include "inputChecker.h"
#include "currentSource.h"
#include "SineWaveSource.h"
#include "GaussianWaveSource.h"
#include "CosineGaussianWave.h"

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
    fdtd hpw(tlen, xlen, ylen, zlen, tw, dx, dy, dz, checker.amptidute, 10, 12, 4, 1, checker.pmlSize,
            checker.useConnectingInterface, checker.fluidGridSize);
    hpw.setPlasmaParam(1e-300, 760 * 5.3E9, 760, 0);
#else
    fdtd hpw(tlen, xlen, ylen, zlen, tw, dx, dy, dz, checker.amptidute, 10, 12, 4, 1, checker.pmlSize, 
            checker.useConnectingInterface);
#endif    

    GaussianWaveSource gaussianWave(checker.frequency);
    SineWaveSource sineSource(omega);
    CosineGaussianWave cosGaussian(checker.frequency, 0.5 * checker.frequency);
    Point lower(xlen / 2 - 1, ylen / 2 - 1, zlen / 2 - 1);
    Point upper(xlen / 2 + 1, ylen / 2 + 1, zlen / 2 + 1);
    MyDataF R = 1e-20;
    currentSource cSource(source::Z, R, lower, upper, checker.amptidute * 1e13);

    switch (checker.waveType) {
        case GAUSSIAN_WAVE:
            cSource.setSourceType(&gaussianWave);
            hpw.setSourceType(&gaussianWave);
            break;
        case SINE_WAVE:
            cSource.setSourceType(&sineSource);
            hpw.setSourceType(&sineSource);
            break;
        case DERIVATIVE_GAUSSIAN_WAVE:break;
        case COSINE_GAUSSIAN_WAVE:
            cSource.setSourceType(&cosGaussian);
            hpw.setSourceType(&cosGaussian);
            break;
        case ZERO_TYPE:break;
        case ONE_SINE_PULSE:
            break;
        default:
            cSource.setSourceType(&gaussianWave);
    }
    hpw.setSource(&cSource);
    //hpw.initialize();
    hpw.startUp();
    return 0;
}
