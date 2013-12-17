/* 
 * File:   currentSource.cpp
 * Author: skiloop
 * 
 * Created on 2013年12月16日, 下午12:05
 */

#include "currentSource.h"
#include "sourceType.h"

currentSource::currentSource(int direction, MyDataF R, Point lower, Point upper, MyDataF amp)
: source(direction, R, lower, upper, amp) {

}

currentSource::currentSource(const currentSource& orig) : source(orig) {
}

currentSource::~currentSource() {
}

void currentSource::initUpdateCoefficients(data3d<MyDataF>& Cere, data3d<MyDataF>& Cerhw,
        data3d<MyDataF>& Cerhv, MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) {
    MyDataF dr, dw, dv;
    switch (mDirection) {
        case X:
            dr = dx;
            dw = dz;
            dv = dy;
            break;
        case Y:
            dr = dy;
            dw = dx;
            dv = dz;
            break;
        case Z:
            dr = dz;
            dw = dy;
            dv = dx;
            break;
    }
    MyDataF C = dt * dr / mResistancePerComponent / dw / dv;
    MyDataF eps2 = eps_0 * 2;
    MyDataF dt2 = dt * 2;
    for (unsigned i = 0, in = mLowerIndex.x; in <= mUpperIndex.x; i++, in++) {
        for (unsigned j = 0, jn = mLowerIndex.y; jn <= mUpperIndex.y; j++, jn++) {
            for (unsigned k = 0, kn = mLowerIndex.z; kn <= mUpperIndex.z; k++, kn++) {
                Cere.p[in][jn][kn] = (eps2 - C) / (eps2 + C);
                Cerhw.p[in][jn][kn] = dt2 / (eps2 + C) / dv;
                Cerhv.p[in][jn][kn] = -dt2 / (eps2 + C) / dw;
                Cese.p[i][j][k] = -dt2 / (eps2 + C) / dw / dv;
            }
        }
    }

}

void currentSource::updateSource(data3d<MyDataF>& Ex, data3d<MyDataF>& Ey,
        data3d<MyDataF>& Ez, MyDataF t) {
    switch (mDirection) {
        case X:
            for (unsigned i = mLowerIndex.x, is = 0; i <= mUpperIndex.x; i++, is++) {
                for (unsigned j = mLowerIndex.y, js = 0; j <= mUpperIndex.y; j++, js++) {
                    for (unsigned k = mLowerIndex.z, ks = 0; k <= mUpperIndex.z; k++, ks++) {
                        Ex.p[i][j][k] += Cese.p[is][js][ks] * mAmplitude * mPSource->valueAtTime(t);
                    }
                }
            }
            break;
        case Y:
            for (unsigned i = mLowerIndex.x, is = 0; i <= mUpperIndex.x; i++, is++) {
                for (unsigned j = mLowerIndex.y, js = 0; j <= mUpperIndex.y; j++, js++) {
                    for (unsigned k = mLowerIndex.z, ks = 0; k <= mUpperIndex.z; k++, ks++) {
                        Ey.p[i][j][k] += Cese.p[is][js][ks] * mAmplitude * mPSource->valueAtTime(t);
                    }
                }
            }
            break;
        case Z:
            for (unsigned i = mLowerIndex.x, is = 0; i <= mUpperIndex.x; i++, is++) {
                for (unsigned j = mLowerIndex.y, js = 0; j <= mUpperIndex.y; j++, js++) {
                    for (unsigned k = mLowerIndex.z, ks = 0; k <= mUpperIndex.z; k++, ks++) {
                        Ez.p[i][j][k] += Cese.p[is][js][ks] * mAmplitude * mPSource->valueAtTime(t);
                    }
                }
            }
            break;
    }
}