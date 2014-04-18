/* 
 * File:   PointArrays.h
 * Author: skiloop
 *
 * Created on 2014年1月19日, 下午7:47
 */

#ifndef POINTARRAYS_H
#define	POINTARRAYS_H

#include "Point.h"
#include "common.h"

class PointArrays {
public:
    PointArrays();
    PointArrays(const PointArrays& orig);
    virtual ~PointArrays();
    int n;
    MyDataF *pVal;
    Point *pArray;
private:

};

#endif	/* POINTARRAYS_H */

