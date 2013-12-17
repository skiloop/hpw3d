/* 
 * File:   CosineGaussianWave.cpp
 * Author: skiloop
 * 
 * Created on 2013年12月17日, 上午10:36
 */

#include "CosineGaussianWave.h"

CosineGaussianWave::CosineGaussianWave(MyDataF modulation_frequency, MyDataF bandWidth)
: mOmega(M_PI_TWO *modulation_frequency)
, mBandWidth(bandWidth) {
    MyDataF tau = 0.966 / mBandWidth;
    mTimeDelay = 4.5 * tau;
    mTauSquare=tau*tau;
}

CosineGaussianWave::CosineGaussianWave(const CosineGaussianWave& orig)
: mOmega(orig.mOmega)
, mBandWidth(orig.mBandWidth)
, mTauSquare(orig.mTauSquare)
, mTimeDelay(orig.mTimeDelay) {
}

CosineGaussianWave::~CosineGaussianWave() {
}

MyDataF CosineGaussianWave::valueAtTime(MyDataF t) {
    return cos(mOmega * (t - mTimeDelay))
            * exp(-(mTimeDelay - t)*(mTimeDelay - t) / mTauSquare);
}
