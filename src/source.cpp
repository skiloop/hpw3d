/* 
 * File:   source.cpp
 * Author: skiloop
 * 
 * Created on 2013年12月16日, 下午12:11
 */

#include "source.h"

source::source(int direction, MyDataF R, Point lower, Point upper, MyDataF amp)
: mDirection(direction)
, mResistancePerComponent(R)
, mAmplitude(amp)
, mLowerIndex(lower)
, mUpperIndex(upper)
, mPSource(NULL) {
    Cese.create3DArray(upper.x - lower.x + 1, upper.y - lower.y + 1, upper.z - lower.z + 1);
}

source::source(const source& orig) {
    if (this != &orig) {
        mDirection = orig.mDirection;
        mResistancePerComponent = orig.mResistancePerComponent;
        mAmplitude = orig.mAmplitude;
        mLowerIndex = orig.mLowerIndex;
        mUpperIndex = orig.mUpperIndex;
        Cese.freeArray();
        Cese.create3DArray(orig.Cese);
        mPSource = orig.mPSource;
    }
}

source::~source() {
}

