/* 
 * File:   GaussianWaveSource.cpp
 * Author: skiloop
 * 
 * Created on 2013年12月16日, 下午3:19
 */

#include "common.h"
#include "GaussianWaveSource.h"

GaussianWaveSource::GaussianWaveSource(MyDataF maxf)
: mMaxFrequency(maxf) {
    mTau = 1 / mMaxFrequency / 2;
    mTimeDelay = 4.5 * mTau;
}

GaussianWaveSource::GaussianWaveSource(const GaussianWaveSource& orig)
: mMaxFrequency(orig.mMaxFrequency)
, mTimeDelay(orig.mTimeDelay)
, mTau(orig.mTimeDelay) {

}

GaussianWaveSource::~GaussianWaveSource() {
}

MyDataF GaussianWaveSource::valueAtTime(MyDataF t) {
    return exp(-(mTimeDelay - t)*(mTimeDelay - t) / mTau / mTau);
}
