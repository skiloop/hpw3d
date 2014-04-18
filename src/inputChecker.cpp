/* 
 * File:   inputChecker.cpp
 * Author: skiloop
 * 
 * Created on 2013年4月28日, 下午5:01
 */
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <cstring>
#include "inputChecker.h"

using namespace std;

inputChecker::inputChecker()
: waveType(GAUSSIAN)
, threadCount(DEFAULT_THREAD_COUNT)
, pmlSize(DEFAULT_PML_SIZE)
, useConnectingInterface(0)
, yeeCellSizeX(0)
, yeeCellSizeY(0)
, yeeCellSizeZ(0)
, yeeCellSize(DEFAULT_GRID_SIZE)
, fluidGridSize(DEFAULT_FLUID_GRID_SIZE)
, xZoneLen(0)
, yZoneLen(0)
, zZoneLen(0)
, tZoneLen(DEFAULT_TIME_ZONE_LENGTH)
, zoneLen(0)
, frequency(DEFAULT_FREQUENCY)
, amptidute(DEFAULT_AMPTIDUTE)
, maxNe(DEFAULT_DENSITY_MAX)
, rei(0)
, nu_type(DEFAULT_NU_FORMAT)
, useDensity(IF_USE_DENSITY) {
}

inputChecker::inputChecker(const inputChecker& orig) {
}

inputChecker::~inputChecker() {
}

void inputChecker::check() {
    if (waveType < 0 || waveType > MAX_TYPE_VALUE) {
        waveType = GAUSSIAN;
    }
    if (threadCount <= 0) {
        threadCount = DEFAULT_THREAD_COUNT;
    }
    if (pmlSize <= -1) {
        pmlSize = DEFAULT_PML_SIZE;
    }
    if (fluidGridSize <= 0) {
        fluidGridSize = DEFAULT_FLUID_GRID_SIZE;
    }
    if (yeeCellSize > 0) {
        yeeCellSizeX = yeeCellSize;
        yeeCellSizeY = yeeCellSize;
        yeeCellSizeZ = yeeCellSize;
    }
    if (yeeCellSizeX <= 0) {
        yeeCellSizeX = DEFAULT_GRID_SIZE;
    }
    if (yeeCellSizeY <= 0) {
        yeeCellSizeY = DEFAULT_GRID_SIZE;
    }
    if (yeeCellSizeZ <= 0) {
        yeeCellSizeZ = DEFAULT_GRID_SIZE;
    }
    if (zoneLen > 0) {
        xZoneLen = zoneLen;
        yZoneLen = zoneLen;
        zZoneLen = zoneLen;
    }
    if (xZoneLen < 0) {
        xZoneLen = DEFAULT_ZONE_SIZE;
    }
    if (yZoneLen < 0) {
        yZoneLen = DEFAULT_ZONE_SIZE;
    }
    if (zZoneLen < 0) {
        zZoneLen = DEFAULT_ZONE_SIZE;
    }
    if (tZoneLen < 0) {
        tZoneLen = DEFAULT_TIME_ZONE_LENGTH;
    }
    if (frequency < 0) {
        frequency = DEFAULT_FREQUENCY;
    }
    if (useConnectingInterface == 1) {
        useConnectingInterface = 1;
    } else {
        useConnectingInterface = 0;
    }
    if (maxNe < 0) {
        maxNe = DEFAULT_DENSITY_MAX;
    }
    if (USE_DENSITY != useDensity) {
        useDensity = NOT_USE_DENSITY;
    }
}

void inputChecker::help(char *prog) {
    cout << prog << " [options]" << endl;
    const char tab = '\t';

    // HELP
    cout << tab << "--help,-h" << tab << "help" << endl;
    cout << tab << "--openmp-thread=n" << tab << "number of threads for openmp,n is positive" << endl;
    cout << tab << "--wave-type=[1,2]" << tab << "source wave form type,1 for gaussian pulse,2 for sine wave" << endl;
    cout << tab << "--fluid-grid-size=" << tab << "how many fluid grid size per Maxwell grid " << endl;
    //=======================================
    // fdtd parameters
    //=======================================
    cout << "FDTD parameters:" << endl;
    cout << tab << "--pml-width=n" << tab << "pml size for fdtd" << endl;
    cout << tab << "--yc-size-x=" << tab << "how many yee cells per wave length of pulse length in x direction" << endl;
    cout << tab << "--yc-size-y=" << tab << "how many yee cells per wave length of pulse length in y direction" << endl;
    cout << tab << "--yc-size-z=" << tab << "how many yee cells per wave length of pulse length in z direction" << endl;
    cout << tab << "--yee-cell-size=" << tab << "how many yee cells per wave length of pulse length,setting yee cell cube" << endl;

    //=======================================
    // Gaussian pulse parameters
    //=======================================
    cout << "Gaussian pulse parameters:" << endl;
    //cout << tab << "--pulse-width=pw" << tab << "Gaussian pulse width,in naseconds" << endl;
    cout << tab << "--amptidute=" << tab << "amptidute for source wave" << endl;
    cout << tab << "--frequency=" << tab << "wave frequency" << endl;
    cout << tab << "--x-zone-length=" << tab << "zone length in x-direction,in pulse width" << endl;
    cout << tab << "--y-zone-length=" << tab << "zone length in y-direction,in pulse width" << endl;
    cout << tab << "--z-zone-length=" << tab << "zone length in z-direction,in pulse width" << endl;
    cout << tab << "--zone-size=" << tab << "zone size,in pulse width,set x,y,z zone length together and ,setting yee cell cube" << endl;
    cout << tab << "--simulation-time=" << tab << "simulation time,in pulse width size" << endl;

    //=======================================
    // Sine wave parameters
    //=======================================
    cout << "Sine wave parameters:" << endl;
    cout << tab << "--amptidute=" << tab << "amptidute for source wave" << endl;
    cout << tab << "--frequency=" << tab << "wave frequency" << endl;
    cout << tab << "--x-zone-length=" << tab << "zone length in x-direction,in wave length" << endl;
    cout << tab << "--y-zone-length=" << tab << "zone length in y-direction,in wave length" << endl;
    cout << tab << "--z-zone-length=" << tab << "zone length in z-direction,in wave length" << endl;
    cout << tab << "--zone-size=" << tab << "zone size,in wave length,set x,y,z zone length together,setting yee cell cube" << endl;
    cout << tab << "--simulation-time=" << tab << "simulation time,in wave length size" << endl;
    cout << tab << "--is-connecting=" << tab << "1 if use connecting interface " << endl;
    cout << tab << "--max-ne=" << tab << "set Ne max" << endl;
    cout << tab << "--rei=" << tab << "set R_ei" << endl;
    cout << tab << "--nu-type=" << tab << "set nu type" << endl;
    cout << tab << "--use-density=" << tab << "whether use density or not" << endl;
    cout << tab << tab << USE_DENSITY << tab << "use density";
    cout << (IF_USE_DENSITY == USE_DENSITY ? "(default)" : "") << endl;
    cout << tab << tab << NOT_USE_DENSITY << tab << "DO NOT use density";
    cout << (IF_USE_DENSITY == NOT_USE_DENSITY ? "(default)" : "") << endl;
    exit(0);
}

void inputChecker::parseInput(int argc, char *argv[]) {
    if (argc <= 1) {
        help(argv[0]);
    }
    for (int i = 1; i < argc; i++) {
        if (0 == strcmp("-h", argv[i]) || 0 == strcmp(argv[i], "--help")) {
            help(argv[0]);
        } else
            if (0 == strncmp(argv[i], "--openmp-thread=", 16)) {
            threadCount = strtol(argv[i] + 16, NULL, 10);
        } else if (0 == strncmp(argv[i], "--wave-type=", 12)) {
            waveType = strtol(argv[i] + 12, NULL, 10);
        } else if (0 == strncmp(argv[i], "--pml-width=", 12)) {
            pmlSize = strtol(argv[i] + 12, NULL, 10);
        } else if (0 == strncmp(argv[i], "--x-zone-length=", 16)) {
            xZoneLen = atof(argv[i] + 16);
        } else if (0 == strncmp(argv[i], "--y-zone-length=", 16)) {
            yZoneLen = atof(argv[i] + 16);
        } else if (0 == strncmp(argv[i], "--z-zone-length=", 16)) {
            zZoneLen = atof(argv[i] + 16);
        } else if (0 == strncmp(argv[i], "--zone-size=", 12)) {
            zoneLen = atof(argv[i] + 12);
        } else if (0 == strncmp(argv[i], "--simulation-time=", 18)) {
            tZoneLen = atof(argv[i] + 18);
        } else if (0 == strncmp(argv[i], "--amptidute=", 12)) {
            amptidute = atof(argv[i] + 12);
        } else if (0 == strncmp(argv[i], "--frequency=", 12)) {
            frequency = atof(argv[i] + 12);
        } else if (0 == strncmp(argv[i], "--yc-size-x=", 12)) {
            yeeCellSizeX = strtol(argv[i] + 12, NULL, 10);
        } else if (0 == strncmp(argv[i], "--yc-size-y=", 12)) {
            yeeCellSizeY = strtol(argv[i] + 12, NULL, 10);
        } else if (0 == strncmp(argv[i], "--yc-size-z=", 12)) {
            yeeCellSizeZ = strtol(argv[i] + 12, NULL, 10);
        } else if (0 == strncmp(argv[i], "--yee-cell-size=", 16)) {
            yeeCellSize = strtol(argv[i] + 16, NULL, 10);
        } else if (0 == strncmp(argv[i], "--fluid-grid-size=", 18)) {
            fluidGridSize = strtol(argv[i] + 18, NULL, 10);
        } else if (0 == strncmp(argv[i], "--is-connecting=", 16)) {
            useConnectingInterface = strtol(argv[i] + 16, NULL, 10);
        } else if (0 == strncmp(argv[i], "--max-ne=", 9)) {
            maxNe = atof(argv[i] + 9);
        } else if (0 == strncmp(argv[i], "--rei=", 6)) {
            rei = atof(argv[i] + 6);
        } else if (0 == strncmp(argv[i], "--nu-type=", 10)) {
            nu_type = strtol(argv[i] + 10, NULL, 10);
        } else if (0 == strncmp(argv[i], "--use-density=", 14)) {
            useDensity = strtol(argv[i] + 14, NULL, 10);
        }
    }
    check();
}

void inputChecker::print() {
    cout << "===========< input parameters >===============" << endl;
    cout << "waveType=" << waveType << endl;
    cout << "threadCount=" << threadCount << endl;
    cout << "pmlSize=" << pmlSize << endl;
    cout << "yeeCellSizeX=" << yeeCellSizeX << endl;
    cout << "yeeCellSizeY=" << yeeCellSizeY << endl;
    cout << "yeeCellSizeZ=" << yeeCellSizeZ << endl;
    cout << "yeeCellSize=" << yeeCellSize << endl;
    cout << "fluidGridSize=" << fluidGridSize << endl;
    cout << "xZoneLen=" << xZoneLen << " lambda" << endl;
    cout << "yZoneLen=" << yZoneLen << " lambda" << endl;
    cout << "zZoneLen=" << zZoneLen << " lambda" << endl;
    cout << "tZoneLen=" << tZoneLen << " T" << endl;
    cout << "zoneLen=" << zoneLen << " lambda" << endl;
    cout << "frequency=" << frequency << " Hz" << endl;
    cout << "amptidute=" << amptidute << " V/m" << endl;
    cout << "useConnecting=" << useConnectingInterface << endl;
    cout << "max-ne=" << maxNe << endl;
    cout << "rei=" << rei << endl;
    cout << "nu-type=" << nu_type << endl;
    cout << "if-use-density=" << useDensity << endl;
    cout << "============================================" << endl;
}
