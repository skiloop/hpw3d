/* 
 * File:   source.cpp
 * Author: skiloop
 * 
 * Created on 2013年5月10日, 下午9:58
 */
#include <cmath>

#include "common.h"
#include "sourceType.h"

sourceType::sourceType() {
}

sourceType::sourceType(const sourceType& orig) {
}

sourceType::~sourceType() {
}

MyDataF sourceType::SinePulse(MyDataF t, MyDataF omega, MyDataF t_max, MyDataF t_min) {
    if (t >= t_min && t <= t_max) {
        return cos(M_PI_TWO * omega * t);
    } else {
        return 0.0;
    }
}

MyDataF sourceType::SineWave(MyDataF t, MyDataF omega) {
    return cos(M_PI_TWO * omega * t);
}

MyDataF sourceType::valueAtTime(MyDataF t) {
    return 0.0;
}
