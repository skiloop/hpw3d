/* 
 * File:   CosineGaussianWave.h
 * Author: skiloop
 *
 * Created on 2013年12月17日, 上午10:36
 */

#ifndef COSINEGAUSSIANWAVE_H
#define	COSINEGAUSSIANWAVE_H

#include "sourceType.h"

class CosineGaussianWave : public sourceType {
public:
    CosineGaussianWave(MyDataF modulation_frequency, MyDataF bandWidth);
    CosineGaussianWave(const CosineGaussianWave& orig);
    virtual ~CosineGaussianWave();
    MyDataF valueAtTime(MyDataF t);
private:
    MyDataF mOmega;
    MyDataF mBandWidth;
    MyDataF mTauSquare;
    MyDataF mTimeDelay;

};

#endif	/* COSINEGAUSSIANWAVE_H */

