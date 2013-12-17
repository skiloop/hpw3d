/* 
 * File:   Capacitor.h
 * Author: skiloop
 *
 * Created on 2013年12月15日, 下午10:09
 */

#ifndef CAPACITOR_H
#define	CAPACITOR_H
#include "Point.h"

class Capacitor {
public:
    Capacitor();
    Capacitor(const Capacitor& orig);
    virtual ~Capacitor();
private:
    int mDirection;
    double mCapasitance;
    Point mLowerIndex;
    Point mUpperIndex;
};

#endif	/* CAPACITOR_H */

