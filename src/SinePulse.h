/* 
 * File:   SinePulse.h
 * Author: skiloop
 *
 * Created on 2014年1月14日, 下午9:44
 */

#ifndef SINEPULSE_H
#define	SINEPULSE_H
#include "sourceType.h"

class SinePulse : public sourceType {
public:
    SinePulse(MyDataF T, MyDataF t0);
    SinePulse(const SinePulse& orig);
    virtual ~SinePulse();
    MyDataF valueAtTime(MyDataF t);
private:
    MyDataF mT;
    MyDataF mTimeDelay;
    MyDataF mOmega;
    MyDataF mFinishingTime;
};

#endif	/* SINEPULSE_H */

