/* 
 * File:   SquareWave.h
 * Author: skiloop
 *
 * Created on 2014年1月15日, 上午9:41
 */

#ifndef SQUAREWAVE_H
#define	SQUAREWAVE_H

#include "sourceType.h"

class SquarePulse:public sourceType {
public:
    SquarePulse(MyDataF lower_t,MyDataF upper_t);
    SquarePulse(const SquarePulse& orig);
    virtual ~SquarePulse();
    MyDataF valueAtTime(MyDataF t);
private:
    MyDataF mLowerTime;
    MyDataF mUpperTime;    
};

#endif	/* SQUAREWAVE_H */

