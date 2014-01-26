/* 
 * File:   source.cpp
 * Author: skiloop
 * 
 * Created on 2013年12月16日, 下午12:11
 */

#include "source.h"

source::source(int direction, MyDataF R, MyDataF amp)
: mDirection(direction)
, mResistancePerComponent(R)
, mAmplitude(amp)
, mPSource(NULL) {

}

source::source(const source& orig) {
    if (this != &orig) {
        mDirection = orig.mDirection;
        mResistancePerComponent = orig.mResistancePerComponent;
        mAmplitude = orig.mAmplitude;
        mPSource = orig.mPSource;
    }
}

source::~source() {
        
}

void source::add(Point point) {
    mSrcIndexes.push_back(point);
}
