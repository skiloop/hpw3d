/* 
 * File:   inputChecker.h
 * Author: skiloop
 *
 * Created on 2013年4月28日, 下午5:01
 */

#ifndef INPUTCHECKER_H
#define	INPUTCHECKER_H

#include "common.h"

class inputChecker {
public:
    inputChecker();
    inputChecker(const inputChecker& orig);
    virtual ~inputChecker();
    void help(char *argv);
    void print();
    void parseInput(int argc, char *argv[]);
public:
    int waveType;
    int threadCount;
    int pmlSize;
    int yeeCellSizeX;
    int yeeCellSizeY;
    int yeeCellSizeZ;
    int yeeCellSize;
    int fluidGridSize;
    double xZoneLen;
    double yZoneLen;
    double zZoneLen;
    double tZoneLen;
    double zoneLen;
    double frequency;
    double amptidute;
    //parameter for sine pulse
    double t0;
    double omega;
    double tUp;
    double tDown;
    
    const static int GAUSSIAN = GAUSSIAN_WAVE_TYPE;
    const static int SINE = SINE_WAVE_TYPE;
private:
    void checkInput();
    void check();
};

#endif	/* INPUTCHECKER_H */

