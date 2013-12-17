/* 
 * File:   SineWaveSource.cpp
 * Author: skiloop
 * 
 * Created on 2013年12月16日, 下午3:53
 */

#include <math.h>

#include "SineWaveSource.h"

SineWaveSource::SineWaveSource(MyDataF omega) : mOmega(omega) {
}

SineWaveSource::SineWaveSource(const SineWaveSource& orig) : mOmega(orig.mOmega) {
}

SineWaveSource::~SineWaveSource() {
}

MyDataF SineWaveSource::valueAtTime(MyDataF t) {
    return sin(M_PI_TWO * t * mOmega);
}
