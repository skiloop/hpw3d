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
    float xZoneLen;
    float yZoneLen;
    float zZoneLen;
    float tZoneLen;
    float zoneLen;
    double frequency;
    double amptidute;
    const static int GAUSSIAN = GAUSSIAN_WAVE_TYPE;
    const static int SINE = SINE_WAVE_TYPE;
private:
    void checkInput();
    void check();
};

#endif	/* INPUTCHECKER_H */

