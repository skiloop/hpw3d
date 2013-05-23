/* 
 * File:   source.cpp
 * Author: skiloop
 * 
 * Created on 2013年5月10日, 下午9:58
 */
#include <math.h>

#include "common.h"
#include "source.h"

Source::Source() {
}

Source::Source(const Source& orig) {
}

Source::~Source() {
}

MyDataF Source::SinePulse(MyDataF t, MyDataF omega, MyDataF t_max, MyDataF t_min) {
    if (t >= t_min && t <= t_max) {
        return cos(M_PI_TWO * omega * t);
    } else {
        return 0.0;
    }
}

MyDataF Source::SineWave(MyDataF t, MyDataF omega) {
    return cos(M_PI_TWO * omega * t);
}