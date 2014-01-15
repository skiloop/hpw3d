/* 
 * File:   SquareWave.cpp
 * Author: skiloop
 * 
 * Created on 2014年1月15日, 上午9:42
 */

#include <math.h>

#include "SquarePulse.h"

SquarePulse::SquarePulse(MyDataF lower_t, MyDataF upper_t)
: mLowerTime(lower_t)
, mUpperTime(upper_t) {
}

SquarePulse::SquarePulse(const SquarePulse& orig)
: mLowerTime(orig.mLowerTime)
, mUpperTime(orig.mUpperTime) {
}

SquarePulse::~SquarePulse() {
}

MyDataF SquarePulse::valueAtTime(MyDataF t) {
    if (t >= mLowerTime && t <= mUpperTime) {
        return 1.0;
    } else {
        return 0.0;
    }
}

