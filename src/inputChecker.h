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

    /**
     * constructor for inputChecker
     */
    inputChecker();

    /**
     * copy constructor
     * @param orig
     */
    inputChecker(const inputChecker& orig);

    /**
     * deconstructor
     */
    virtual ~inputChecker();

    /**
     * print help for program
     * @param argv
     */
    void help(char *argv);

    /**
     * print parramters from input
     */
    void print();

    /**
     * parse input variables
     * @param argc
     * @param argv
     */
    void parseInput(int argc, char *argv[]);

public:
    // wave form type
    int waveType;

    // thread count for openmp
    int threadCount;

    // pml width
    int pmlSize;
    
    // if use connecting interface
    int useConnectingInterface;

    // 
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
    double amplidute;
    double maxNe;
    //parameter for sine pulse
    double t0;
    double omega;
    double tUp;
    double tDown;
    double rei;
    double pressure;
    
    int nu_type;
    
    int useDensity;
    
    int MaxwellTimeNumber;
    int FluidTimeNumber;

    const static int GAUSSIAN = GAUSSIAN_WAVE;
    const static int SINE = SINE_WAVE;
private:
    void checkInput();
    void check();
};

#endif	/* INPUTCHECKER_H */

