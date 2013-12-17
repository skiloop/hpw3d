/* 
 * File:   GaussianWaveSource.h
 * Author: skiloop
 *
 * Created on 2013年12月16日, 下午3:19
 */

#ifndef GAUSSIANWAVESOURCE_H
#define	GAUSSIANWAVESOURCE_H
#include "common.h"
#include "sourceType.h"

class GaussianWaveSource : public sourceType {
public:
    GaussianWaveSource(MyDataF maxf);
    GaussianWaveSource(const GaussianWaveSource& orig);
    virtual ~GaussianWaveSource();
    MyDataF valueAtTime(MyDataF t);
private:
    MyDataF mMaxFrequency; // number of cells per wavelength
    MyDataF mTimeDelay; // delay t0
    MyDataF mTau; // tau
};

#endif	/* GAUSSIANWAVESOURCE_H */

