/* 
 * File:   currentSource.cpp
 * Author: skiloop
 * 
 * Created on 2013年12月16日, 下午12:05
 */

#include "sourceType.h"
#include "currentSource.h"

currentSource::currentSource(int direction, MyDataF R, MyDataF amp)
: source(direction, R, amp) {

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

void currentSource::initCoefficients(data3d<MyDataF>& Cere, data3d<MyDataF>& Cerhw,
        data3d<MyDataF>& Cerhv, MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) {
    MyDataF dr, dw, dv;
    desideXYZ(dr, dw, dv, dx, dy, dz);
    MyDataF C = dt * dr / mResistancePerComponent / dw / dv;
    MyDataF eps2 = eps_0 * 2;
    MyDataF dt2 = dt * 2;
    //    for (unsigned i = 0, in = mLowerIndex.x; in <= mUpperIndex.x; i++, in++) {
    //        for (unsigned j = 0, jn = mLowerIndex.y; jn <= mUpperIndex.y; j++, jn++) {
    //            for (unsigned k = 0, kn = mLowerIndex.z; kn <= mUpperIndex.z; k++, kn++) {
    //                Cere.p[in][jn][kn] = (eps2 - C) / (eps2 + C);
    //                Cerhw.p[in][jn][kn] = dt2 / (eps2 + C) / dv;
    //                Cerhv.p[in][jn][kn] = -dt2 / (eps2 + C) / dw;
    //                Cese.p[i][j][k] = -dt2 / (eps2 + C) / dw / dv;
    //            }
    //        }
    //    }
    int i = 0;
    while (i < mSrcIndexes.size()) {
        Cere.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = (eps2 - C) / (eps2 + C);
        Cerhw.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = dt2 / (eps2 + C) / dv;
        Cerhv.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = -dt2 / (eps2 + C) / dw;
        Cese.push_back(-dt2 / (eps2 + C) / dw / dv);
        i++;
    }

}

void currentSource::initCoefficients(data3d<MyDataF>& Cere, data3d<MyDataF>& Cerhw,
        data3d<MyDataF>& Cerhv, data3d<MyDataF>& Cervr, const data3d<MyDataF>& Beta,
        const data3d<MyDataF>& Ne, const MyDataF vm, const Point& nes, const Point& bes,
        unsigned m, MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) {
    MyDataF dr, dw, dv;
    desideXYZ(dr, dw, dv, dx, dy, dz);
    MyDataF psi = dr / mResistancePerComponent / 2 / dw / dv; //psi
    MyDataF alpha, a, bb, b, theta;
    MyDataF ne;
    alpha = vm * dt / 2;
    alpha = (1.0 - alpha) / (1.0 + alpha);
    a = eps_0 / dt;
    bb = e / 2;
    theta = 1.0 / dw / dv;

    unsigned i = 0;
    while (i < mSrcIndexes.size()) {

        ne = Ne.p[mSrcIndexes[i].x * m - nes.x][mSrcIndexes[i].y * m - nes.y][mSrcIndexes[i].z * m - nes.z];
        b = bb * Beta.p[mSrcIndexes[i].x * m - bes.x][mSrcIndexes[i].y * m + bes.y][mSrcIndexes[i].z * m - bes.z] * ne + psi;
        Cere.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = (a - b) / (a + b);
        Cervr.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = bb * ne * (1 + alpha) / (a + b);
        Cerhw.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = 1.0 / dv / eps_0 / (a + b);
        Cerhv.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = -1.0 / dw / eps_0 / (a + b);
        Cese.push_back(-theta / (a + b));

        i++;
    }
}

void currentSource::initCoefficients(data3d<MyDataF>& Cere, data3d<MyDataF>& Cerhw,
        data3d<MyDataF>& Cerhv, data3d<MyDataF>& Cervr, const data3d<MyDataF>& Beta,
        const data3d<MyDataF>& Ne, const data3d<MyDataF>& Nu_c, const Point& nes, const Point& bes,
        unsigned m, MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) {
    MyDataF dr, dw, dv;
    desideXYZ(dr, dw, dv, dx, dy, dz);
    MyDataF psi = dr / mResistancePerComponent / 2 / dw / dv; //psi
    MyDataF dt2 = dt / 2;
    MyDataF alpha, a, bb, b, theta;
    MyDataF ne;

    a = eps_0 / dt;
    bb = e / 2;
    theta = 1.0 / dw / dv;

    unsigned i = 0;
    while (i < mSrcIndexes.size()) {
        ne = Ne.p[mSrcIndexes[i].x * m - nes.x][mSrcIndexes[i].y * m - nes.y][mSrcIndexes[i].z * m - nes.z];
        b = bb * Beta.p[mSrcIndexes[i].x * m - bes.x][mSrcIndexes[i].y * m - bes.y][mSrcIndexes[i].z * m - bes.z] * ne + psi;
        alpha = Nu_c.p[mSrcIndexes[i].x * m - bes.x][mSrcIndexes[i].y * m - bes.y][mSrcIndexes[i].z * m - bes.z] * dt2;
        alpha = (1.0 - alpha) / (1.0 + alpha);
        Cere.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = (a - b) / (a + b);
        Cervr.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = bb * ne * (1 + alpha) / (a + b);
        Cerhw.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = 1.0 / dv / eps_0 / (a + b);
        Cerhv.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = -1.0 / dw / eps_0 / (a + b);
        Cese.push_back(-theta / (a + b));

        i++;
    }

}

///////////////////////////////////////////
// UPDAE COEFF
//////////////////////////////////////////

void currentSource::updateCoefficients(data3d<MyDataF>& Cere, data3d<MyDataF>& Cerhw,
        data3d<MyDataF>& Cerhv, MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) {
    MyDataF dr, dw, dv;
    desideXYZ(dr, dw, dv, dx, dy, dz);
    MyDataF C = dt * dr / mResistancePerComponent / dw / dv;
    MyDataF eps2 = eps_0 * 2;
    MyDataF dt2 = dt * 2;
    //    for (unsigned i = 0, in = mLowerIndex.x; in <= mUpperIndex.x; i++, in++) {
    //        for (unsigned j = 0, jn = mLowerIndex.y; jn <= mUpperIndex.y; j++, jn++) {
    //            for (unsigned k = 0, kn = mLowerIndex.z; kn <= mUpperIndex.z; k++, kn++) {
    //                Cere.p[in][jn][kn] = (eps2 - C) / (eps2 + C);
    //                Cerhw.p[in][jn][kn] = dt2 / (eps2 + C) / dv;
    //                Cerhv.p[in][jn][kn] = -dt2 / (eps2 + C) / dw;
    //                Cese.p[i][j][k] = -dt2 / (eps2 + C) / dw / dv;
    //            }
    //        }
    //    }
    int i = 0;
    while (i < mSrcIndexes.size()) {
        Cere.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = (eps2 - C) / (eps2 + C);
        Cerhw.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = dt2 / (eps2 + C) / dv;
        Cerhv.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = -dt2 / (eps2 + C) / dw;
        Cese[i] = (-dt2 / (eps2 + C) / dw / dv);
        i++;
    }

}

void currentSource::updateCoefficients(data3d<MyDataF>& Cere, data3d<MyDataF>& Cerhw,
        data3d<MyDataF>& Cerhv, data3d<MyDataF>& Cervr, const data3d<MyDataF>& Beta,
        const data3d<MyDataF>& Ne, const MyDataF vm, const Point& nes, const Point& bes,
        unsigned m, MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) {
    MyDataF dr, dw, dv;
    desideXYZ(dr, dw, dv, dx, dy, dz);
    MyDataF psi = dr / mResistancePerComponent / 2 / dw / dv; //psi
    MyDataF alpha, a, bb, b, theta;
    MyDataF ne;
    alpha = vm * dt / 2;
    alpha = (1.0 - alpha) / (1.0 + alpha);
    a = eps_0 / dt;
    bb = e / 2;
    theta = 1.0 / dw / dv;

    unsigned i = 0;
    while (i < mSrcIndexes.size()) {

        ne = Ne.p[mSrcIndexes[i].x * m - nes.x][mSrcIndexes[i].y * m - nes.y][mSrcIndexes[i].z * m - nes.z];
        b = bb * Beta.p[mSrcIndexes[i].x * m - bes.x][mSrcIndexes[i].y * m + bes.y][mSrcIndexes[i].z * m - bes.z] * ne + psi;
        Cere.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = (a - b) / (a + b);
        Cervr.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = bb * ne * (1 + alpha) / (a + b);
        Cerhw.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = 1.0 / dv / eps_0 / (a + b);
        Cerhv.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = -1.0 / dw / eps_0 / (a + b);
        Cese[i] = (-theta / (a + b));

        i++;
    }
}

void currentSource::updateCoefficients(data3d<MyDataF>& Cere, data3d<MyDataF>& Cerhw,
        data3d<MyDataF>& Cerhv, data3d<MyDataF>& Cervr, const data3d<MyDataF>& Beta,
        const data3d<MyDataF>& Ne, const data3d<MyDataF>& Nu_c, const Point& nes, const Point& bes,
        unsigned m, MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) {
    MyDataF dr, dw, dv;
    desideXYZ(dr, dw, dv, dx, dy, dz);
    MyDataF psi = dr / mResistancePerComponent / 2 / dw / dv; //psi
    MyDataF dt2 = dt / 2;
    MyDataF alpha, a, bb, b, theta;
    MyDataF ne;

    a = eps_0 / dt;
    bb = e / 2;
    theta = 1.0 / dw / dv;

    unsigned i = 0;
    while (i < mSrcIndexes.size()) {
        ne = Ne.p[mSrcIndexes[i].x * m - nes.x][mSrcIndexes[i].y * m - nes.y][mSrcIndexes[i].z * m - nes.z];
        b = bb * Beta.p[mSrcIndexes[i].x * m - bes.x][mSrcIndexes[i].y * m - bes.y][mSrcIndexes[i].z * m - bes.z] * ne + psi;
        alpha = Nu_c.p[mSrcIndexes[i].x * m - bes.x][mSrcIndexes[i].y * m - bes.y][mSrcIndexes[i].z * m - bes.z] * dt2;
        alpha = (1.0 - alpha) / (1.0 + alpha);
        Cere.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = (a - b) / (a + b);
        Cervr.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = bb * ne * (1 + alpha) / (a + b);
        Cerhw.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = 1.0 / dv / eps_0 / (a + b);
        Cerhv.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = -1.0 / dw / eps_0 / (a + b);
        Cese[i] = (-theta / (a + b));

        i++;
    }

}

///////////////////////////////////////////
// update source
///////////////////////////////////////////

void currentSource::updateSource(data3d<MyDataF>& Ex, data3d<MyDataF>& Ey,
        data3d<MyDataF>& Ez, MyDataF t) {
    unsigned i = 0;
    switch (mDirection) {
        case X:
            while (i < mSrcIndexes.size()) {
                Ex.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] += Cese[i] * mAmplitude * mPSource->valueAtTime(t);
                i++;
            }
            break;
        case Y:
            while (i < mSrcIndexes.size()) {
                Ey.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] += Cese[i] * mAmplitude * mPSource->valueAtTime(t);
                i++;
            }
            break;
        case Z:
            while (i < mSrcIndexes.size()) {
                Ez.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] += Cese[i] * mAmplitude * mPSource->valueAtTime(t);
                i++;
            }
            break;
    }
}

void currentSource::updateHardSource(data3d<MyDataF>& Ex, data3d<MyDataF>& Ey,
        data3d<MyDataF>& Ez, MyDataF t) {
    unsigned i = 0;
    switch (mDirection) {
        case X:
            while (i < mSrcIndexes.size()) {
                Ex.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = Cese[i] * mAmplitude * mPSource->valueAtTime(t);
                i++;
            }
            break;
        case Y:
            while (i < mSrcIndexes.size()) {
                Ey.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = Cese[i] * mAmplitude * mPSource->valueAtTime(t);
                i++;
            }
            break;
        case Z:
            while (i < mSrcIndexes.size()) {
                Ez.p[mSrcIndexes[i].x][mSrcIndexes[i].y][mSrcIndexes[i].z] = Cese[i] * mAmplitude * mPSource->valueAtTime(t);
                i++;
            }
            break;
    }
}