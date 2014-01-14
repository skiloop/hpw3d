/* 
 * File:   currentSource.cpp
 * Author: skiloop
 * 
 * Created on 2013年12月16日, 下午12:05
 */

#include "sourceType.h"
#include "currentSource.h"

currentSource::currentSource(int direction, MyDataF R, Point lower, Point upper, MyDataF amp)
: source(direction, R, lower, upper, amp) {

}

currentSource::currentSource(const currentSource& orig) : source(orig) {
}

currentSource::~currentSource() {
}

void currentSource::desideXYZ(MyDataF& dr, MyDataF& dw, MyDataF& dv,
        const MyDataF dx, const MyDataF dy, const MyDataF dz) {
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
}

void currentSource::initUpdateCoefficients(data3d<MyDataF>& Cere, data3d<MyDataF>& Cerhw,
        data3d<MyDataF>& Cerhv, MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) {
    MyDataF dr, dw, dv;
    desideXYZ(dr, dw, dv, dx, dy, dz);
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

void currentSource::initUpdateCoefficients(data3d<MyDataF>& Cere, data3d<MyDataF>& Cerhw,
        data3d<MyDataF>& Cerhv, data3d<MyDataF>& Cervr, const data3d<MyDataF>& Beta,
        const data3d<MyDataF>& Ne, const MyDataF vm, const Point& nes, const Point& bes,
        unsigned m, MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) {
    MyDataF dr, dw, dv;
    desideXYZ(dr, dw, dv, dx, dy, dz);
    MyDataF psi = dr / mResistancePerComponent / 2 / dw / dv; //psi
    MyDataF alpha, a, bb, b, theta, phi;
    MyDataF ne;
    alpha = vm * dt / 2;
    alpha = (1.0 - alpha) / (1.0 + alpha);
    a = eps_0 / dt;
    bb = e / 2;
    theta = 1.0 / dw / dv;
    for (unsigned i = 0, in = mLowerIndex.x; in <= mUpperIndex.x; i++, in++) {
        for (unsigned j = 0, jn = mLowerIndex.y; jn <= mUpperIndex.y; j++, jn++) {
            for (unsigned k = 0, kn = mLowerIndex.z; kn <= mUpperIndex.z; k++, kn++) {
                ne = Ne.p[i * m + nes.x][j * m + nes.y][k * m + nes.z];
                b = bb * Beta.p[i * m + bes.x][j * m + bes.y][k * m + bes.z] * ne + psi;
                Cere.p[in][jn][kn] = (a - b) / (a + b);
                Cervr.p[in][jn][kn] = bb * ne * (1 + alpha) / (a + b);
                Cerhw.p[in][jn][kn] = 1.0 / dv / eps_0 / (a + b);
                Cerhv.p[in][jn][kn] = -1.0 / dw / eps_0 / (a + b);
                Cese.p[i][j][k] = -theta / (a + b);
            }
        }
    }
}

void currentSource::initUpdateCoefficients(data3d<MyDataF>& Cere, data3d<MyDataF>& Cerhw,
        data3d<MyDataF>& Cerhv, data3d<MyDataF>& Cervr, const data3d<MyDataF>& Beta,
        const data3d<MyDataF>& Ne, const data3d<MyDataF>& Nu_c, const Point& nes, const Point& bes,
        unsigned m, MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) {
    MyDataF dr, dw, dv;
    desideXYZ(dr, dw, dv, dx, dy, dz);
    MyDataF psi = dr / mResistancePerComponent / 2 / dw / dv; //psi
    MyDataF dt2 = dt / 2;
    MyDataF alpha, a, bb, b, theta, phi;
    MyDataF ne;

    a = eps_0 / dt;
    bb = e / 2;
    theta = 1.0 / dw / dv;
    for (unsigned i = 0, in = mLowerIndex.x; in <= mUpperIndex.x; i++, in++) {
        for (unsigned j = 0, jn = mLowerIndex.y; jn <= mUpperIndex.y; j++, jn++) {
            for (unsigned k = 0, kn = mLowerIndex.z; kn <= mUpperIndex.z; k++, kn++) {
                ne = Ne.p[i * m + nes.x][j * m + nes.y][k * m + nes.z];
                b = bb * Beta.p[i * m + bes.x][j * m + bes.y][k * m + bes.z] * ne + psi;
                alpha = Nu_c.p[i * m + bes.x][j * m + bes.y][k * m + bes.z] * dt2;
                alpha = (1.0 - alpha) / (1.0 + alpha);
                Cere.p[in][jn][kn] = (a - b) / (a + b);
                Cervr.p[in][jn][kn] = bb * ne * (1 + alpha) / (a + b);
                Cerhw.p[in][jn][kn] = 1.0 / dv / eps_0 / (a + b);
                Cerhv.p[in][jn][kn] = -1.0 / dw / eps_0 / (a + b);
                Cese.p[i][j][k] = -theta / (a + b);
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