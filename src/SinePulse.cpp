/* 
 * File:   SinePulse.cpp
 * Author: skiloop
 * 
 * Created on 2014年1月14日, 下午9:44
 */

#include <cmath>

#include "SinePulse.h"

SinePulse::SinePulse(MyDataF T, MyDataF t0)
: mT(fabs(T))
, mTimeDelay(t0) {
    mFinishingTime = mT + t0;
    mOmega = M_PI_TWO / T;
}

SinePulse::SinePulse(const SinePulse& orig)
: mT(orig.mT)
, mTimeDelay(orig.mTimeDelay)
, mOmega(orig.mOmega)
, mFinishingTime(orig.mFinishingTime) {
}

SinePulse::~SinePulse() {
}

MyDataF SinePulse::valueAtTime(MyDataF t) {
    if (t > mTimeDelay && t < mFinishingTime) {
        return sin(M_PI_TWO * (t - mTimeDelay) * mOmega);
    } else {
        return 0.0;
    }
}
