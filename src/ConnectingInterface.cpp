/* 
 * File:   ConnectingInterface.cpp
 * Author: skiloop
 * 
 * Created on 2013年12月17日, 下午3:39
 */

#include <sstream>

#include "ConnectingInterface.h"

ConnectingInterface::ConnectingInterface() {
}

ConnectingInterface::ConnectingInterface(const ConnectingInterface& orig) {
}

ConnectingInterface::~ConnectingInterface() {
}

void ConnectingInterface::setLowerIndex(Point &p) {
    mLowerIndex = p;
}

void ConnectingInterface::setUpperIndex(Point &p) {
    mUpperIndex = p;
}

void ConnectingInterface::setIncidentAngle(MyDataF theta, MyDataF phi, MyDataF psi, MyDataF lamda, MyDataF c, MyDataF dt, MyDataF ds) {
    mTheta = theta;
    mPhi = phi;
    mEx = -sin(psi) * sin(phi) + cos(psi) * cos(theta) * cos(phi);
    mEy = sin(psi) * cos(phi) + cos(psi) * cos(theta) * sin(phi);
    mEz = -cos(psi) * sin(theta);
    mHx = -cos(psi) * sin(phi) - sin(psi) * cos(theta) * cos(phi);
    mHy = cos(psi) * cos(phi) - sin(psi) * cos(theta) * sin(phi);
    mHz = sin(psi) * sin(phi);
    mDisC1 = sin(mTheta) * cos(mPhi);
    mDisC2 = sin(mTheta) * sin(mPhi);
    mDisC3 = cos(mTheta);
    if (fabs(mDisC1) < 1e-10)mDisC1 = 0;
    if (fabs(mDisC2) < 1e-10)mDisC2 = 0;
    if (fabs(mDisC3) < 1e-10)mDisC3 = 0;
    mPhaseVelocityRatio = PhaseVelAngle(0.0, 0.0, lamda, dt, ds, c) / PhaseVelAngle(phi, theta, lamda, dt, ds, c);
}

void ConnectingInterface::initCoefficients(MyDataF ds, MyDataF dt) {
    MyDataF delta = ds * mPhaseVelocityRatio;
    mChe = dt / mu_0 / delta;
    mCeh = dt / eps_0 / delta;
    mCMur = (C * dt - delta) / (C * dt + delta);
}

void ConnectingInterface::updateEConnect(data3d<MyDataF>& Ex, const data3d<MyDataF> &Cexhy, const data3d<MyDataF> &Cexhz,
        data3d<MyDataF>& Ey, const data3d<MyDataF> &Ceyhx, const data3d<MyDataF> &Ceyhz,
        data3d<MyDataF>& Ez, const data3d<MyDataF> &Cezhy, const data3d<MyDataF> &Cezhx) {
    //update left side
    MyDataF delay, df;
    int d;
    Point lower(mLowerIndex);
    Point upper(mLowerIndex);
    upper.y = mUpperIndex.y;
    for (unsigned i = 0; i < mHxDelayFrontBack.nx; i++, lower.x++, upper.x++) {
        lower.z = mLowerIndex.z;
        upper.z = mLowerIndex.z;
        for (unsigned k = 0; k < mHxDelayFrontBack.nz; k++, lower.z++, upper.z++) {
            delay = mHxDelayFrontBack.p[i][0][k];
            d = (int) floor(delay);
            df = delay - d;
            Ez.p[lower.x][lower.y][lower.z] -= Cezhx[lower] * mHx * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
            delay = mHxDelayFrontBack.p[i][1][k];
            d = (int) floor(delay);
            df = delay - d;
            Ez.p[upper.x][upper.y][upper.z] += Cezhx[upper] * mHx * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
        }
    }
    lower = mLowerIndex;
    upper = mLowerIndex;
    upper.y = mUpperIndex.y;
    for (unsigned i = 0; i < mHzDelayFrontBack.nx; i++, lower.x++, upper.x++) {
        lower.z = mLowerIndex.z;
        upper.z = mLowerIndex.z;
        for (unsigned k = 0; k < mHzDelayFrontBack.nz; k++, lower.z++, upper.z++) {
            delay = mHzDelayFrontBack.p[i][0][k];
            d = (int) floor(delay);
            df = delay - d;
            Ex.p[lower.x][lower.y][lower.z] -= Cexhz[lower] * mHz * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
            delay = mHzDelayFrontBack.p[i][1][k];
            d = (int) floor(delay);
            df = delay - d;
            Ex.p[upper.x][upper.y][upper.z] += Cexhz[upper] * mHz * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
        }
    }
    // front and back
    lower = mLowerIndex;
    upper = mLowerIndex;
    upper.x = mUpperIndex.x;
    for (unsigned j = 0; j < mHzDelayLeftRight.ny; j++, lower.y++, upper.y++) {
        lower.z = mLowerIndex.z;
        upper.z = mLowerIndex.z;
        for (unsigned k = 0; k < mHzDelayLeftRight.nz; k++, lower.z++, upper.z++) {
            delay = mHzDelayLeftRight.p[0][j][k];
            d = (int) floor(delay);
            df = delay - d;
            Ey.p[lower.x][lower.y][lower.z] -= Ceyhz[lower] * mHz * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
            delay = mHzDelayLeftRight.p[1][j][k];
            d = (int) floor(delay);
            df = delay - d;
            Ey.p[upper.x][upper.y][upper.z] += Ceyhz[upper] * mHz * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
        }
    }
    lower = mLowerIndex;
    upper = mLowerIndex;
    upper.x = mUpperIndex.x;
    for (unsigned j = 0; j < mHyDelayLeftRight.ny; j++, lower.y++, upper.y++) {
        lower.z = mLowerIndex.z;
        upper.z = mLowerIndex.z;
        for (unsigned k = 0; k < mHyDelayLeftRight.nz; k++, lower.z++, upper.z++) {
            delay = mHyDelayLeftRight.p[0][j][k];
            d = (int) floor(delay);
            df = delay - d;
            Ez.p[lower.x][lower.y][lower.z] -= Cezhy[lower] * mHy * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
            delay = mHyDelayLeftRight.p[1][j][k];
            d = (int) floor(delay);
            df = delay - d;
            Ez.p[upper.x][upper.y][upper.z] += Cezhy[upper] * mHy * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
        }
    }

    // top and back
    lower = mLowerIndex;
    upper = mLowerIndex;
    upper.z = mUpperIndex.z;
    for (unsigned j = 0; j < mHyDelayTopBottom.ny; j++, lower.y++, upper.y++) {
        lower.x = mLowerIndex.x;
        upper.x = mLowerIndex.x;
        for (unsigned i = 0; i < mHyDelayTopBottom.nx; i++, lower.x++, upper.x++) {
            delay = mHyDelayTopBottom.p[i][j][0];
            d = (int) floor(delay);
            df = delay - d;
            Ex.p[lower.x][lower.y][lower.z] -= Cexhy[lower] * mHy * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
            delay = mHyDelayTopBottom.p[i][j][1];
            d = (int) floor(delay);
            df = delay - d;
            Ex.p[upper.x][upper.y][upper.z] += Cexhy[upper] * mHy * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
        }
    }
    lower = mLowerIndex;
    upper = mLowerIndex;
    upper.z = mUpperIndex.z;
    for (unsigned j = 0; j < mHxDelayTopBottom.ny; j++, lower.y++, upper.y++) {
        lower.x = mLowerIndex.x;
        upper.x = mLowerIndex.x;
        for (unsigned i = 0; i < mHxDelayTopBottom.nx; i++, lower.x++, upper.x++) {
            delay = mHxDelayTopBottom.p[i][j][0];
            d = (int) floor(delay);
            df = delay - d;
            Ey.p[lower.x][lower.y][lower.z] -= Ceyhx[lower] * mHx * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
            delay = mHxDelayTopBottom.p[i][j][1];
            d = (int) floor(delay);
            df = delay - d;
            Ey.p[upper.x][upper.y][upper.z] += Ceyhx[upper] * mHx * (mHinc.p[d] + df * (mHinc.p[d + 1] - mHinc.p[d]));
        }
    }
}

void ConnectingInterface::updateMConnect(data3d<MyDataF> &Hx, data3d<MyDataF> &Chxey, const data3d<MyDataF> &Chxez,
        data3d<MyDataF>& Hy, const data3d<MyDataF> &Chyex, const data3d<MyDataF> &Chyez,
        data3d<MyDataF>& Hz, const data3d<MyDataF> &Chzey, const data3d<MyDataF> &Chzex) {
    //update left side
    MyDataF df;
    int d;
    Point lower(mLowerIndex);
    Point upper(mLowerIndex);
    lower.y--;
    upper.y = mUpperIndex.y;
    // left right
    for (unsigned i = 0; i < mEzDelayFrontBack.nx; i++, lower.x++, upper.x++) {
        lower.z = mLowerIndex.z;
        upper.z = mLowerIndex.z;
        for (unsigned k = 0; k < mEzDelayFrontBack.nz; k++, lower.z++, upper.z++) {
            d = (int) floor(mEzDelayFrontBack.p[i][0][k]);
            df = mEzDelayFrontBack.p[i][0][k] - d;
            Hx.p[lower.x][lower.y][lower.z] += Chxez[lower] * mEz * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
            d = (int) floor(mEzDelayFrontBack.p[i][1][k]);
            df = mEzDelayFrontBack.p[i][1][k] - d;
            Hx.p[upper.x][upper.y][upper.z] -= Chxez[upper] * mEz * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
        }
    }
    lower = mLowerIndex;
    upper = mLowerIndex;
    lower.y--;
    upper.y = mUpperIndex.y;
    for (unsigned i = 0; i < mExDelayFrontBack.nx; i++, lower.x++, upper.x++) {
        lower.z = mLowerIndex.z;
        upper.z = mLowerIndex.z;
        for (unsigned k = 0; k < mExDelayFrontBack.nz; k++, lower.z++, upper.z++) {
            d = (int) floor(mExDelayFrontBack.p[i][0][k]);
            df = mExDelayFrontBack.p[i][0][k] - d;
            Hz.p[lower.x][lower.y][lower.z] += Chzex[lower] * mEx * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
            d = (int) floor(mExDelayFrontBack.p[i][1][k]);
            df = mExDelayFrontBack.p[i][1][k] - d;
            Hz.p[upper.x][upper.y][upper.z] -= Chzex[upper] * mEx * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
        }
    }
    // front and back
    lower = mLowerIndex;
    upper = mLowerIndex;
    lower.x--;
    upper.x = mUpperIndex.x;
    for (unsigned j = 0; j < mEyDelayLeftRight.ny; j++, lower.y++, upper.y++) {
        lower.z = mLowerIndex.z;
        upper.z = mLowerIndex.z;
        for (unsigned k = 0; k < mEyDelayLeftRight.nz; k++, lower.z++, upper.z++) {
            d = (int) floor(mEyDelayLeftRight.p[0][j][k]);
            df = mEyDelayLeftRight.p[0][j][k] - d;
            Hz.p[lower.x][lower.y][lower.z] += Chzey[lower] * mEy * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
            d = (int) floor(mEyDelayLeftRight.p[1][j][k]);
            df = mEyDelayLeftRight.p[1][j][k] - d;
            Hz.p[upper.x][upper.y][upper.z] -= Chzey[upper] * mEy * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
        }
    }
    lower = mLowerIndex;
    upper = mLowerIndex;
    lower.x--;
    upper.x = mUpperIndex.x;
    for (unsigned j = 0; j < mEzDelayLeftRight.ny; j++, lower.y++, upper.y++) {
        lower.z = mLowerIndex.z;
        upper.z = mLowerIndex.z;
        for (unsigned k = 0; k < mEzDelayLeftRight.nz; k++, lower.z++, upper.z++) {
            d = (int) floor(mEzDelayLeftRight.p[0][j][k]);
            df = mEzDelayLeftRight.p[0][j][k] - d;
            Hy.p[lower.x][lower.y][lower.z] += Chyez[lower] * mEz * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
            d = (int) floor(mEzDelayLeftRight.p[1][j][k]);
            df = mEzDelayLeftRight.p[1][j][k] - d;
            Hy.p[upper.x][upper.y][upper.z] -= Chyez[upper] * mEz * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
        }
    }

    // top and bottom
    lower = mLowerIndex;
    upper = mLowerIndex;
    lower.z--;
    upper.z = mUpperIndex.z;
    for (unsigned j = 0; j < mExDelayTopBottom.ny; j++, lower.y++, upper.y++) {
        lower.x = mLowerIndex.x;
        upper.x = mLowerIndex.x;
        for (unsigned i = 0; i < mExDelayTopBottom.nx; i++, lower.x++, upper.x++) {
            d = (int) floor(mExDelayTopBottom.p[i][j][0]);
            df = mExDelayTopBottom.p[i][j][0] - d;
            Hy.p[lower.x][lower.y][lower.z] += Chyex[lower] * mEx * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
            d = (int) floor(mExDelayTopBottom.p[i][j][1]);
            df = mExDelayTopBottom.p[i][j][1] - d;
            Hy.p[upper.x][upper.y][upper.z] -= Chyex[upper] * mEx * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
        }
    }
    lower = mLowerIndex;
    upper = mLowerIndex;
    lower.z--;
    upper.z = mUpperIndex.z;
    for (unsigned j = 0; j < mEyDelayTopBottom.ny; j++, lower.y++, upper.y++) {
        lower.x = mLowerIndex.x;
        upper.x = mLowerIndex.x;
        for (unsigned i = 0; i < mEyDelayTopBottom.nx; i++, lower.x++, upper.x++) {
            d = (int) floor(mEyDelayTopBottom.p[i][j][0]);
            df = mEyDelayTopBottom.p[i][j][0] - d;
            Hx.p[lower.x][lower.y][lower.z] += Chxey[lower] * mEy * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
            d = (int) floor(mEyDelayTopBottom.p[i][j][1]);
            df = mEyDelayTopBottom.p[i][j][1] - d;
            Hx.p[upper.x][upper.y][upper.z] -= Chxey[upper] * mEy * (mEinc.p[d] + df * (mEinc.p[d + 1] - mEinc.p[d]));
        }
    }
}

void ConnectingInterface::updateESource(MyDataF source) {
    MyDataF elast2 = mEinc.p[mEinc.n - 2];
    mEinc.p[0] = source;
    for (unsigned i = 1; i <= mEinc.n - 2; i++) {
        mEinc.p[i] += mCeh * (mHinc.p[i] - mHinc.p[i - 1]);
    }
    mEinc.p[mEinc.n - 1] = elast2 + mCMur * (mEinc.p[mEinc.n - 2] - mEinc.p[mEinc.n - 1]);
}

void ConnectingInterface::updateMSource() {
    for (unsigned i = 0; i < mHinc.n; i++) {
        mHinc.p[i] += mChe * (mEinc.p[i + 1] - mEinc.p[i]);
    }
}

void ConnectingInterface::setLowerAndUpper(Point& lower, Point& upper) {
    setLowerIndex(lower);
    setUpperIndex(upper);
}

void ConnectingInterface::initData() {
    unsigned si = (mUpperIndex.x - mLowerIndex.x);
    unsigned sj = (mUpperIndex.y - mLowerIndex.y);
    unsigned sk = (mUpperIndex.z - mLowerIndex.z);
    unsigned sx = si + 1;
    unsigned sy = sj + 1;
    unsigned sz = sk + 1;
    MyDataF tmp;

    unsigned len = ceil(sqrt(si * si + sj * sj + sk * sk)) + 3;

    mEinc.createArray(len, 0.0);
    mHinc.createArray(len - 1, 0.0);
    //    // create coeffcients array
    //    mCezhy.create3DArray(2, sy, sk);
    //    mCezhx.create3DArray(sx, 2, sk);
    //    mCexhy.create3DArray(si, 2, sz);
    //    mCexhz.create3DArray(si, sy, 2);
    //    mCeyhx.create3DArray(2, sj, sz);
    //    mCeyhz.create3DArray(sx, sj, 2);
    //    mChzey.create3DArray(2, sy, sk);
    //    mChzex.create3DArray(sx, 2, sk);
    //    mChxey.create3DArray(si, 2, sz);
    //    mChxez.create3DArray(si, sy, 2);
    //    mChyex.create3DArray(2, sj, sz);
    //    mChyez.create3DArray(sx, sj, 2);
    // create E array
    mEzDelayFrontBack.create3DArray(sx, 2, sk);
    mEzDelayLeftRight.create3DArray(2, sy, sk);
    mExDelayFrontBack.create3DArray(si, 2, sz);
    mExDelayTopBottom.create3DArray(si, sy, 2);
    mEyDelayLeftRight.create3DArray(2, sj, sz);
    mEyDelayTopBottom.create3DArray(sx, sj, 2);
    // create H array
    mHzDelayFrontBack.create3DArray(si, 2, sz);
    mHzDelayLeftRight.create3DArray(2, sj, sz);
    mHxDelayFrontBack.create3DArray(sx, 2, sk);
    mHxDelayTopBottom.create3DArray(sx, sj, 2);
    mHyDelayLeftRight.create3DArray(2, sy, sk);
    mHyDelayTopBottom.create3DArray(si, sy, 2);
    // inittial delays    
    for (unsigned k = 0; k < sk; k++) {
        tmp = mLowerIndex.z + k + 0.5;
        for (unsigned i = 0, ik = mLowerIndex.x; i < sx; i++, ik++) {
            mEzDelayFrontBack.p[i][0][k] = 1.0 + distance(ik, mLowerIndex.y, tmp);
            mEzDelayFrontBack.p[i][1][k] = 1.0 + distance(ik, mUpperIndex.y, tmp);
            mHxDelayFrontBack.p[i][0][k] = 0.5 + distance(ik, mLowerIndex.y - 0.5, tmp);
            mHxDelayFrontBack.p[i][1][k] = 0.5 + distance(ik, mUpperIndex.y + 0.5, tmp);
        }
        for (unsigned j = 0, ik = mLowerIndex.y; j < sy; j++, ik++) {
            mEzDelayLeftRight.p[0][j][k] = 1.0 + distance(mLowerIndex.x, ik, tmp);
            mEzDelayLeftRight.p[1][j][k] = 1.0 + distance(mUpperIndex.x, ik, tmp);
            mHyDelayLeftRight.p[0][j][k] = 0.5 + distance(mLowerIndex.x - 0.5, ik, tmp);
            mHyDelayLeftRight.p[1][j][k] = 0.5 + distance(mUpperIndex.x + 0.5, ik, tmp);
        }
    }
    for (unsigned i = 0; i < si; i++) {
        tmp = mLowerIndex.x + i + 0.5;
        for (unsigned k = 0, ik = mLowerIndex.z; k < sz; k++, ik++) {
            mExDelayFrontBack.p[i][0][k] = 1.0 + distance(tmp, mLowerIndex.y, ik);
            mExDelayFrontBack.p[i][1][k] = 1.0 + distance(tmp, mUpperIndex.y, ik);
            mHzDelayFrontBack.p[i][0][k] = 0.5 + distance(tmp, mLowerIndex.y - 0.5, ik);
            mHzDelayFrontBack.p[i][1][k] = 0.5 + distance(tmp, mUpperIndex.y + 0.5, ik);
        }
        for (unsigned j = 0, ik = mLowerIndex.y; j < sy; j++, ik++) {
            mExDelayTopBottom.p[i][j][0] = 1.0 + distance(tmp, ik, mLowerIndex.z);
            mExDelayTopBottom.p[i][j][1] = 1.0 + distance(tmp, ik, mUpperIndex.z);
            mHyDelayTopBottom.p[i][j][0] = 0.5 + distance(tmp, ik, mLowerIndex.z - 0.5);
            mHyDelayTopBottom.p[i][j][1] = 0.5 + distance(tmp, ik, mUpperIndex.z + 0.5);
        }
    }

    for (unsigned j = 0; j < sj; j++) {
        tmp = mLowerIndex.y + j + 0.5;
        for (unsigned i = 0, ik = mLowerIndex.x; i < sx; i++, ik++) {
            mEyDelayTopBottom.p[i][j][0] = 1.0 + distance(ik, tmp, mLowerIndex.z);
            mEyDelayTopBottom.p[i][j][1] = 1.0 + distance(ik, tmp, mUpperIndex.z);
            mHxDelayTopBottom.p[i][j][0] = 0.5 + distance(ik, tmp, mLowerIndex.z - 0.5);
            mHxDelayTopBottom.p[i][j][1] = 0.5 + distance(ik, tmp, mUpperIndex.z + 0.5);
        }
        for (unsigned k = 0, ik = mLowerIndex.z; k < sz; k++, ik++) {
            mEyDelayLeftRight.p[0][j][k] = 1.0 + distance(mLowerIndex.x, tmp, ik);
            mEyDelayLeftRight.p[1][j][k] = 1.0 + distance(mUpperIndex.x, tmp, ik);
            mHzDelayLeftRight.p[0][j][k] = 0.5 + distance(mLowerIndex.x - 0.5, tmp, ik);
            mHzDelayLeftRight.p[1][j][k] = 0.5 + distance(mUpperIndex.x + 0.5, tmp, ik);
        }
    }
}

MyDataF ConnectingInterface::distance(MyDataF i, MyDataF j, MyDataF k) {
#if(DEBUG>=1)
    if ((i < mOrigPoint.x) || (j < mOrigPoint.y) || (k < mOrigPoint.z)) {
        i = i;
    }
#endif
    return mDisC1 * (i - mOrigPoint.x) + mDisC2 * (j - mOrigPoint.y) + mDisC3 * (k - mOrigPoint.z);
}

void ConnectingInterface::desideOrigPoint() {
    mPhi = mPhi - M_PI * 2 * floor(mPhi / (2 * M_PI));
    mTheta = mTheta - M_PI * floor(mTheta / M_PI);
    switch ((int) floor(mTheta * 2 / M_PI)) {
        case 0:
            switch ((int) floor(mPhi * 2 / M_PI)) {
                case 0:
                    mOrigPoint = mLowerIndex;
                    break;
                case 1:
                    mOrigPoint.x = mUpperIndex.x;
                    mOrigPoint.y = mLowerIndex.y;
                    mOrigPoint.z = mLowerIndex.z;
                    break;
                case 2:
                    mOrigPoint.x = mUpperIndex.x;
                    mOrigPoint.y = mUpperIndex.y;
                    mOrigPoint.z = mLowerIndex.z;
                    break;
                case 3:
                    mOrigPoint.x = mLowerIndex.x;
                    mOrigPoint.y = mUpperIndex.y;
                    mOrigPoint.z = mLowerIndex.z;
                    break;
            }
            break;
        case 1:
            switch ((int) floor(mPhi * 2 / M_PI)) {
                case 0:
                    mOrigPoint.x = mLowerIndex.x;
                    mOrigPoint.y = mLowerIndex.y;
                    mOrigPoint.z = mUpperIndex.z;
                    break;
                case 1:
                    mOrigPoint.x = mUpperIndex.x;
                    mOrigPoint.y = mLowerIndex.y;
                    mOrigPoint.z = mUpperIndex.z;
                    break;
                case 2:
                    mOrigPoint = mUpperIndex;
                    break;
                case 3:
                    mOrigPoint.x = mLowerIndex.x;
                    mOrigPoint.y = mUpperIndex.y;
                    mOrigPoint.z = mUpperIndex.z;
                    break;
            }
            break;
    }

}

void ConnectingInterface::invalidate() {
    desideOrigPoint();
    initData();
    mEzDelayFrontBack.setName("ezdl");
    mEzDelayFrontBack.save(1);
    mEzDelayLeftRight.setName("ezdf");
    mEzDelayLeftRight.save(1);
    mExDelayFrontBack.setName("exdl");
    mExDelayFrontBack.save(1);
    mExDelayTopBottom.setName("exdt");
    mExDelayTopBottom.save(1);
    mEyDelayLeftRight.setName("eydf");
    mEyDelayLeftRight.save(1);
    mHyDelayLeftRight.setName("hydf");
    mHyDelayLeftRight.save(1);
    mHyDelayTopBottom.setName("hydt");
    mHyDelayTopBottom.save(1);
    mHzDelayLeftRight.setName("hzdf");
    mHzDelayLeftRight.save(1);
    mHzDelayFrontBack.setName("hzdl");
    mHzDelayFrontBack.save(1);
    mEyDelayLeftRight.setName("eydf");
    mEyDelayLeftRight.save(1);
}

MyDataF ConnectingInterface::PhaseVelAngle(MyDataF phi, MyDataF theta, MyDataF lamda, MyDataF dt, MyDataF dx, MyDataF c) {

    MyDataF A, B, C, D, S, eta, N, eta_k;

    N = lamda / dx;
    S = c * dt / dx;
    A = M_PI * cos(phi) * sin(theta) / N;
    B = M_PI * sin(phi) * sin(theta) / N;
    C = M_PI * cos(theta) / N;
    MyDataF tmp = M_PI * S / N;
    D = sin(tmp) * sin(tmp) / S / S;

    eta = 0.5;
    eta_k = 1;
    while (fabs(eta - eta_k) > 1e-4) {
        eta_k = eta;
        eta = eta_k - (sin(A * eta_k) * sin(A * eta_k) + sin(B * eta_k) * sin(B * eta_k) + sin(C * eta_k) * sin(C * eta_k) - D)
                / (A * sin(2.0 * A * eta_k) + B * sin(2.0 * B * eta_k) + C * sin(2.0 * C * eta_k));
    }

    return lamda / c / eta_k;
}

void ConnectingInterface::saveEMInc(int step) {
    stringstream ss;
    ss << "einc_" << step << ".dat";
    mEinc.save(ss.str());
    stringstream sm;
    sm << "minc_" << step << ".dat";
    mHinc.save(sm.str());
}