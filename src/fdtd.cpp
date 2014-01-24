

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>     
#include <assert.h>
#include <fstream>
#include <sstream>
#include <exception>

#ifdef _OPENMP
#include <omp.h>
extern int thread_count;
#endif

#include "cpml.h"
#include "fdtd.h"
#include "data1d.h"
#include "data3d.h"
#include "sourceType.h"
#include "InonizationFormula.h"
#include "Point.h"
#include "ConnectingInterface.h"

extern MyDataF epsR;
//extern MyDataF dt, dx, dy, dz;
extern MyDataF T;

using namespace std;

void checkmax(unsigned &u_2check, unsigned min, unsigned max) {
    if (u_2check >= max || u_2check < min)u_2check = (min + max) / 2;
}

//#ifdef WITH_DENSITY

fdtd::fdtd(int useDensity, unsigned _totalTimeSteps, unsigned xzoneSize, unsigned yzoneSize, unsigned zzoneSize,
        MyDataF _tw, MyDataF _dx, MyDataF _dy, MyDataF _dz,
        MyDataF _amp, unsigned _savemodulus, unsigned _ksource,
        unsigned _m, unsigned _ma, unsigned pmlw, int useConnect, unsigned _neGrid, MyDataF maxNe)
: mIsUseDensity(useDensity)
, mTotalTimeSteps(_totalTimeSteps)
, tw(_tw), mDx(_dx), mDy(_dy), mDz(_dz)
, mAmplitude(_amp), mSaveModulus(_savemodulus), mKSource(_ksource)
, mPMLOrder(_m), mAlphaOrder(_ma), mPMLWidth(pmlw)
, mAirBufferWidth(AIR_BUFFER)
, mNeGridSize(_neGrid)
, Ne0(maxNe)
, mSrcType(SOURCE_SINE)
, mUseConnectingInterface(!(useConnect == 0))
, mEpsilon(NULL), mSigma(NULL), mMu(NULL)
, pSource(NULL), pSourceType(NULL), CA(NULL), CB(NULL)
, mNeBoundWidth(NE_BOUND_WIDTH) {
    mMaxIndex.setValue(xzoneSize + 2 * (pmlw + mAirBufferWidth),
            yzoneSize + 2 * (pmlw + mAirBufferWidth),
            zzoneSize + 2 * (pmlw + mAirBufferWidth));
    desideDomainZone();
}
//#else
//
//fdtd::fdtd(unsigned _totalTimeSteps, unsigned _imax, unsigned _jmax, unsigned _kmax,
//        MyDataF _tw, MyDataF _dx, MyDataF _dy, MyDataF _dz,
//        MyDataF _amp, unsigned _savemodulus, unsigned _ksource,
//        unsigned _m, unsigned _ma, unsigned pmlw, int useConnect)
//: mTotalTimeSteps(_totalTimeSteps), mMaxIndex(_imax, _jmax, _kmax)
//, tw(_tw), mDx(_dx), mDy(_dy), mDz(_dz)
//, mAmplitude(_amp), mSaveModulus(_savemodulus), mKSource(_ksource)
//, mPMLOrder(_m), mAlphaOrder(_ma), mPMLWidth(pmlw)
//, mAirBufferWidth(AIR_BUFFER)
//, mSrcType(SOURCE_SINE)
//, mUseConnectingInterface(~(useConnect == 0))
//, mEpsilon(NULL), mSigma(NULL), mMu(NULL)
//, pSource(NULL), pSourceType(NULL)
//, CA(NULL), CB(NULL) {
//    mMaxIndex.setValue(_imax + 2 * (pmlw + mAirBufferWidth),
//            _jmax + 2 * (pmlw + mAirBufferWidth),
//            _kmax + 2 * (pmlw + mAirBufferWidth));
//    desideDomainZone();
//}
//#endif

fdtd::~fdtd(void) {
    if (mEpsilon != NULL)delete []mEpsilon;
    if (mSigma != NULL)delete []mSigma;
    if (mMu != NULL)delete[]mMu;
    if (CA != NULL)delete[]CA;
    if (CB != NULL)delete[]CB;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//#ifdef WITH_DENSITY
const MyDataF fdtd::mNeutralGasDensity = 2.44e19;

void fdtd::setPlasmaParam(MyDataF _rei, MyDataF _vm, MyDataF _p, int _ftype) {
    mRei = _rei;
    mNu_m = _vm;
    mAirPressure = _p;
    mNiuType = _ftype;
}

void fdtd::captureEFieldForEeff(void) {
    unsigned i, j, k;
    unsigned io, jo, ko;
    MyDataF sxIJK, syIJK, szIJK;

    Point start(0, 0, 0);
    switch (mSrcType) {
        case SOURCE_SINE:
            for (i = mDomainStartIndex.x, io = start.x; i <= mDomainEndIndex.x; i++, io += mNeGridSize) {
                for (j = mDomainStartIndex.y, jo = start.y; j <= mDomainEndIndex.y; j++, jo += mNeGridSize) {
                    for (k = mDomainStartIndex.z, ko = start.z; k <= mDomainEndIndex.z; k++, ko += mNeGridSize) {
                        sxIJK = (Ex.p[i][j][k] + Ex.p[i + 1][j][k]) / 2;
                        syIJK = (Ey.p[i][j][k] + Ey.p[i][j + 1][k]) / 2;
                        szIJK = (Ez.p[i][j][k] + Ez.p[i][j][k + 1]) / 2;
                        MyDataF tmp = sqrt(szIJK * szIJK + sxIJK * sxIJK + syIJK * syIJK);
                        if (tmp > Eeff.p[io][jo][ko]) {
                            Eeff.p[io][jo][ko] = tmp;
#if DEBUG>=4
                            if (isnan(Eeff.p[io][jo][ko]) || isinf(Eeff.p[io][jo][ko])) {
                                tmp = 0.0;
                            }
#endif
                        }
                    }
                }
            }
            break;
        case SOURCE_GAUSSIAN:
        default:
            for (i = mDomainStartIndex.x, io = start.x; i <= mDomainEndIndex.x; i++, io += mNeGridSize) {
                for (j = mDomainStartIndex.y, jo = start.y; j <= mDomainEndIndex.y; j++, jo += mNeGridSize) {
                    for (k = mDomainStartIndex.z, ko = start.z; k <= mDomainEndIndex.z; k++, ko += mNeGridSize) {
                        sxIJK = (Ex.p[i][j][k] + Ex.p[i + 1][j][k]) / 2;
                        syIJK = (Ey.p[i][j][k] + Ey.p[i][j + 1][k]) / 2;
                        szIJK = (Ez.p[i][j][k] + Ez.p[i][j][k + 1]) / 2;
                        Eeff.p[io][jo][ko] += (szIJK * szIJK + sxIJK * sxIJK + syIJK * syIJK);
#if DEBUG>=4
                        if (isnan(Eeff.p[io][jo][ko]) || isinf(Eeff.p[io][jo][ko]) /*|| fabs(Eeff.p[io][jo][ko]) > 1e30*/) {
                            szIJK = 0.0;
                        }
#endif
                    }
                }
            }
    }
}

void fdtd::updateCollisionFrequency() {
    if (fdtd::SOURCE_GAUSSIAN == mSrcType) {
        int i, j, k;
        MyDataF EeffDivP;
        MyDataF DivParam = 100 * mAirPressure * 133.3;
        MyDataF C1 = 5.20e8 * mAirPressure;
        MyDataF C2 = 2.93e8 * mAirPressure;
        MyDataF C3 = 3.24e8 * mAirPressure;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k,EeffDivP) //shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif
        for (i = 0; i < Nu_c.nx; i++) {
            for (j = 0; j < Nu_c.ny; j++) {
                for (k = 0; k < Nu_c.nz; k++) {
                    EeffDivP = Eeff.p[i][j][k] / DivParam;
                    if (EeffDivP >= 120) {
                        Nu_c.p[i][j][k] = C1 * sqrt(EeffDivP);
                    } else if (EeffDivP >= 54) {
                        Nu_c.p[i][j][k] = C2 * EeffDivP / (1 + 0.041 * EeffDivP);
                    } else {
                        Nu_c.p[i][j][k] = C3 * EeffDivP / (1 + 0.04 * EeffDivP);
                    }
#if DEBUG>=4
                    if (isnan(Nu_c.p[i][j][k])) {
                        EeffDivP = 0;
                    }
#endif
                }
            }
        }
    }
}

void fdtd::vec2Eeff() {
    MyDataF tmp = mDt / mDtFluid;
    for (unsigned i = 0; i < Eeff.nx; i += mNeGridSize) {
        for (unsigned j = 0; j < Eeff.ny; j += mNeGridSize) {
            for (unsigned k = 0; k < Eeff.nz; k += mNeGridSize) {
#if DEBUG>=4
                MyDataF smp=Eeff.p[i][j][k];
#endif
                Eeff.p[i][j][k] = sqrt(tmp * Eeff.p[i][j][k]);
#if DEBUG>=4
                if (isnan(Eeff.p[i][j][k]) || isinf(Eeff.p[i][j][k])) {
                    smp += mDt / mDtFluid;
                }
#endif
            }
        }
    }
}

void fdtd::updateEeff() {
    unsigned is, js, ks;
    unsigned in, jn, kn;
    unsigned i, j, k;
    unsigned im, jm, km;
    unsigned iu, ju, ku;
    unsigned ngred = mNeGridSize * mNeGridSize*mNeGridSize;
    Point start(0, 0, 0);
    Point end(Eeff.nx, Eeff.ny, Eeff.nz);
    if (fdtd::SOURCE_SINE != mSrcType) {
        vec2Eeff();
    }
    for (is = start.x, in = is + mNeGridSize; in < end.x; is = in, in += mNeGridSize) {
        for (js = start.y, jn = js + mNeGridSize; jn < end.y; js = jn, jn += mNeGridSize) {
            for (ks = start.z, kn = ks + mNeGridSize; kn < end.z; ks = kn, kn += mNeGridSize) {
                // integrate Erms
                for (i = is, im = 0, iu = mNeGridSize; i <= in; i++, im++, iu--) {
                    for (j = js, jm = 0, ju = mNeGridSize; j <= jn; j++, jm++, ju--) {
                        for (k = ks, km = 0, ku = mNeGridSize; k <= kn; k++, km++, ku--) {
                            Eeff.p[i][j][k] = (iu * ju * ku * Eeff.p[is][js][ks] + im * ju * ku * Eeff.p[in][js][ks] +
                                    im * jm * ku * Eeff.p[in][jn][ks] + im * jm * km * Eeff.p[in][jn][kn] +
                                    iu * jm * ku * Eeff.p[is][jn][ks] + iu * jm * km * Eeff.p[is][jn][kn] +
                                    iu * ju * km * Eeff.p[is][js][kn] + im * ju * km * Eeff.p[in][js][kn]) / ngred;
#if DEBUG>=4
                            if (isnan(Eeff.p[i][j][k]) || isinf(Eeff.p[i][j][k])) {
                                Eeff.p[i][j][k] += 1;
                            }
#endif
                        }
                    }
                }
            }
        }
    }
}

void fdtd::calIonizationParam(int i, int j, int k, MyDataF &va, MyDataF &vi, MyDataF &Deff) {
    MyDataF EeffVPerCM;
    switch (mSrcType) {
        case SOURCE_SINE:
            EeffVPerCM = Eeff.p[i - mNeBoundWidth][j - mNeBoundWidth][k - mNeBoundWidth] / 100 * pow(1 / (1 + mOmega * mOmega / mNu_m / mNu_m), 0.5);
            break;
        default:
            EeffVPerCM = Eeff.p[i - mNeBoundWidth][j - mNeBoundWidth][k - mNeBoundWidth] / 100; //convert to V/cm
    }

    switch (mNiuType) {
        case MORROW_AND_LOWKE:
            Niu_MorrowAndLowke(&vi, &va, EeffVPerCM, mNeutralGasDensity);
            break;
        case NIKONOV:
            Niu_Nikonov(&vi, &va, EeffVPerCM, mAirPressure);
            break;
        case KANG:
            Niu_Kang(&vi, &va, EeffVPerCM);
            break;
        case ALI:
        default:
            Niu_Ali(&vi, &va, EeffVPerCM, mAirPressure);
    }
    if (Ne.p[i][j][k] < 1) {
        Deff = mDe;
    } else {
        // tau_m = eps_0 / (e * Ne.p[i][j][k] * (mu_e + mu_i));
        MyDataF kasi = vi * eps_0 / (e * Ne.p[i][j][k] * (mMu_e + mMu_i));
        Deff = (kasi * mDe + mDa) / (kasi + 1);
    }
}

void fdtd::updateDensity(void) {

    int i, j, k;
    Point ms(mNeBoundWidth, mNeBoundWidth, mNeBoundWidth);
    Point me(Ne.nx - mNeBoundWidth, Ne.ny - mNeBoundWidth, Ne.nz - mNeBoundWidth);
    MyDataF Ne_ijk, Neip1, Neim1, Nejm1, Nejp1, Nekp1, Nekm1;
    MyDataF Deff;
    MyDataF maxvi = 0, minvi = 0;
    MyDataF vi, va;
    MyDataF dtfDivDsfSquare = mDtFluid / mDsFluid / mDsFluid;

    Ne_pre = Ne;

#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) \
        schedule(dynamic) private(i,j,k,Ne_ijk, Neip1, Neim1, Nejm1, Nejp1, Nekp1, Nekm1,vi,va,Deff)
#endif
    for (i = ms.x; i < me.x; i++) {
        for (j = ms.y; j < me.y; j++) {
            for (k = ms.z; k < me.z; k++) {

                Ne_ijk = Ne_pre.p[i][j][k];
                Neip1 = Ne_pre.p[i + 1][j][k];
                Neim1 = Ne_pre.p[i - 1][j][k];
                Nejp1 = Ne_pre.p[i][j + 1][k];
                Nejm1 = Ne_pre.p[i][j - 1][k];
                Nekp1 = Ne_pre.p[i][j][k + 1];
                Nekm1 = Ne_pre.p[i][j][k - 1];

                calIonizationParam(i, j, k, va, vi, Deff);

                Ne.p[i][j][k] = (Ne_ijk * (1.0 + mDtFluid * vi) + Deff * dtfDivDsfSquare *
                        (Neip1 + Neim1 + Nejp1 + Nejm1 + Nekp1 + Nekm1 - 6 * Ne_ijk))
                        / (1.0 + mDtFluid * (va + mRei * Ne_ijk));
                if (Ne.p[i][j][k] < 0) {
                    Ne.p[i][j][k] = 0;
                }
                if (vi > maxvi) {
                    maxvi = vi;
                    //                    ci = i;
                    //                    cj = j;
                    //                    ck = k;
                }
                if (vi < minvi) minvi = vi;
#ifdef DEBUG
                if (mOfNeCheck.is_open() && i == mNeSrcPos.x && j == mNeSrcPos.y && k == mNeSrcPos.z) {
                    mOfNeCheck << Eeff.p[i][j][k] << '\t' << Ne_ijk * mDtFluid * vi << '\t' << Deff * dtfDivDsfSquare * \
                            (Neip1 + Neim1 + Nejp1 + Nejm1 + Nekp1 + Nekm1 - 6 * Ne_ijk) << \
                            '\t' << mDtFluid * (va + mRei * Ne_ijk) << endl;
                }
                if(isnan(Ne.p[i][j][k])){
                    vi=0;
                }
#endif
            }
        }
    }
    wallCircleBound(Ne);
    //    cout << Ne.p[Ne.nx / 2][Ne.ny / 2][Ne.nz / 2] << '\t';
    //    cout << maxvi << '\t' << minvi << '\t' << Ne.p[ci][cj][ck] << '\t' << Erms.p[ci][cj][ck] << '\t';
}

void fdtd::updateVelocity(void) {
    updateVx();
    updateVy();
    updateVz();
}

void fdtd::updateElectricFields() {
    updateEx();
    updateEy();
    updateEz();
}

void fdtd::wallCircleBound(data3d<MyDataF> &stru) {
    unsigned i, j, k;

    for (int widthLeft = mNeBoundWidth - 1; widthLeft >= 0; widthLeft--) {

        unsigned start = widthLeft;
        unsigned end = stru.nz - widthLeft - 1;
        unsigned start1 = widthLeft + 1;
        unsigned start2 = widthLeft + 2;
        unsigned end1 = end - 1;
        unsigned end2 = end - 2;
        //bottom and top
        for (i = start; i <= end; i++) {
            for (j = start; j <= end; j++) {
                stru.p[i][j][start] = 2 * stru.p[i][j][start1] - stru.p[i][j][start2];
                stru.p[i][j][end] = 2 * stru.p[i][j][end1] - stru.p[i][j][end2];
            }
        }

        //left and right
        for (j = start; j <= end; j++) {
            for (k = start; k <= end; k++) {
                stru.p[start][j][k] = 2 * stru.p[start1][j][k] - stru.p[start2][j][k];
                stru.p[end][j][k] = 2 * stru.p[end1][j][k] - stru.p[end2][j][k];
            }
        }

        //front and back
        for (i = start; i <= end; i++) {
            for (k = start; k <= end; k++) {
                stru.p[i][start][k] = 2 * stru.p[i][start1][k] - stru.p[i][start2][k];
                stru.p[i][end][k] = 2 * stru.p[i][end1][k] - stru.p[i][end2][k];
            }
        }
    }
}

void fdtd::createDensityRelatedArrays() {
    // velocity coefficients
    //Cvxex.create3DArray(Vx,0.0);
    //Cvyey.create3DArray(Vy,0.0);
    //Cvzez.create3DArray(Vz,0.0);
    // electricity coefficients
    //Cexe.create3DArray(Ex, 0.0);
    //Ceye.create3DArray(Ey, 0.0);
    //Ceze.create3DArray(Ez, 0.0);
    Cexvx.create3DArray(Ex, 0.0);
    Ceyvy.create3DArray(Ey, 0.0);
    Cezvz.create3DArray(Ez, 0.0);
    // beta
    Beta.create3DArray(Eeff.nx, Eeff.ny, Eeff.nz, 0.0);
    // velocity coefficients
    if (fdtd::SOURCE_SINE != mSrcType) {
        Cvxvx.create3DArray(Vx, mAlpha);
        Cvyvy.create3DArray(Vy, mAlpha);
        Cvzvz.create3DArray(Vz, mAlpha);
        Cvxex.create3DArray(Vx, 0.0);
        Cvyey.create3DArray(Vy, 0.0);
        Cvzez.create3DArray(Vz, 0.0);
        Nu_c.create3DArray(Eeff.nx, Eeff.ny, Eeff.nz, 0.0);
        Nu_c.setName("nuc");
    }


}

void fdtd::initCoeffForDensity() {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Velocity Coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    mA = mNu_m * mDt / 2;
    mGamma = 1 + mA;
    mAlpha = (1 - mA) / mGamma;
    mCvxex = mCvyey = mCvzez = e * mDt / 2 / me / mGamma;
    mCoeffVelocity = mHalf_e * mDt / me;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // initials collision frequency
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    updateCollisionFrequency();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //electricity coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    updateCoeffWithDensity();

}

void fdtd::updateCoeffWithDensity() {

    updateBeta();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //electricity coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int i, j, k;
    unsigned im, jm, km;
    MyDataF tmp = eMDtDiv2DivEps0 * (1 + mAlpha);
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k,im,jm,km)//shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif    
    for (j = mDomainStartIndex.y; j <= mDomainEndIndex.y; j++) {
        jm = (j - mDomainStartIndex.y) * mNeGridSize;
        for (i = mDomainStartIndex.x, im = mHalfNeGridSize; i < mDomainEndIndex.x; i++, im += mNeGridSize) {
            for (k = mDomainStartIndex.z, km = 0; k <= mDomainEndIndex.z; k++, km += mNeGridSize) {
                MyDataF kappa = (1 + Beta.p[im][jm][km]);
                Cexe.p[i][j][k] = (1 - Beta.p[im][jm][km]) / kappa;
                Cexhy.p[i][j][k] = -dtDivEps0DivDz / kappa;
                Cexhz.p[i][j][k] = dtDivEps0DivDy / kappa;
                if (fdtd::SOURCE_SINE != mSrcType) {
                    MyDataF a = mHalfDelta_t * Nu_c.p[im][jm][km];
                    MyDataF gamma = 1 + a;
                    Cvxex.p[i][j][k] = mCoeffVelocity / gamma;
                    Cvxvx.p[i][j][k] = (1 - a) / gamma;
                    Cexvx.p[i][j][k] = eMDtDiv2DivEps0 * (1 + (1 - a) / gamma)
                            * Ne.p[im][jm][km] / kappa;
#if DEBUG>=4
                    if (isnan(Cvxvx.p[i][j][k]) || isnan(Cvxex.p[i][j][k]) || isnan(Cexvx.p[i][j][k])) {
                        a = 0;
                    }
#endif
                } else {
                    Cexvx.p[i][j][k] = tmp * Ne.p[im][jm][km] / kappa;
                }
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k,im,jm,km)//shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif    
    for (j = mDomainStartIndex.y; j < mDomainEndIndex.y; j++) {
        jm = (j - mDomainStartIndex.y) * mNeGridSize + mHalfNeGridSize;
        for (i = mDomainStartIndex.x, im = 0; i <= mDomainEndIndex.x; i++, im += mNeGridSize) {
            for (k = mDomainStartIndex.z, km = 0; k <= mDomainEndIndex.z; k++, km += mNeGridSize) {
                MyDataF kappa = (1 + Beta.p[im][jm][km]);
                Ceye.p[i][j][k] = (1 - Beta.p[im][jm][km]) / kappa;
                Ceyhx.p[i][j][k] = dtDivEps0DivDz / kappa;
                Ceyhz.p[i][j][k] = -dtDivEps0DivDx / kappa;
                if (fdtd::SOURCE_SINE != mSrcType) {
                    MyDataF a = mHalfDelta_t * Nu_c.p[im][jm][km];
                    MyDataF gamma = 1 + a;
                    Cvyvy.p[i][j][k] = (1 - a) / gamma;
                    Cvyey.p[i][j][k] = mCoeffVelocity / gamma;
                    Ceyvy.p[i][j][k] = eMDtDiv2DivEps0 * (1 + (1 - a) / gamma)
                            * Ne.p[im][jm][km] / kappa;
#if DEBUG>=4
                    if (isnan(Cvyvy.p[i][j][k]) || isnan(Cvyey.p[i][j][k]) || isnan(Ceyvy.p[i][j][k])) {
                        a = 0;
                    }
#endif
                } else {
                    Ceyvy.p[i][j][k] = tmp * Ne.p[im][jm][km] / kappa;
                }
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k,im,jm,km)//shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif   
    for (j = mDomainStartIndex.y; j <= mDomainEndIndex.y; j++) {
        jm = (j - mDomainStartIndex.y) * mNeGridSize;
        for (i = mDomainStartIndex.x, im = 0; i <= mDomainEndIndex.x; i++, im += mNeGridSize) {
            for (k = mDomainStartIndex.z, km = mHalfNeGridSize; k < mDomainEndIndex.z; k++, km += mNeGridSize) {
                MyDataF kappa = (1 + Beta.p[im][jm][km]);
                Ceze.p[i][j][k] = (1 - Beta.p[im][jm][km]) / kappa;
                Cezhy.p[i][j][k] = dtDivEps0DivDx / kappa;
                Cezhx.p[i][j][k] = -dtDivEps0DivDy / kappa;
                if (fdtd::SOURCE_SINE != mSrcType) {
                    MyDataF a = mHalfDelta_t * Nu_c.p[im][jm][km];
                    MyDataF gamma = 1 + a;
                    Cvzvz.p[i][j][k] = (1 - a) / gamma;
                    Cvzez.p[i][j][k] = mCoeffVelocity / gamma;
                    Cezvz.p[i][j][k] = eMDtDiv2DivEps0 * (1 + (1 - a) / gamma)
                            * Ne.p[im][jm][km] / kappa;
#if DEBUG>=4
                    if (isnan(Cvzvz.p[i][j][k]) || isnan(Cvzez.p[i][j][k]) || isnan(Cezvz.p[i][j][k])) {
                        a = 0;
                    }
#endif
                } else {
                    Cezvz.p[i][j][k] = tmp * Ne.p[im][jm][km] / kappa;
                }
            }
        }
    }
}

void fdtd::updateBeta() {
    Point start(mNeBoundWidth, mNeBoundWidth, mNeBoundWidth);
    Point end(Ne.nx - mNeBoundWidth, Ne.nz - mNeBoundWidth, Ne.nz - mNeBoundWidth);
    MyDataF temp = e2Dt2Div4DivEps0DivMe;
    if (fdtd::SOURCE_SINE == mSrcType) {
        temp = temp / mGamma;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)  //shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif
        for (int i = 0; i < Eeff.nx; i++) {
            int im = i + start.x;
            for (unsigned j = 0, jm = start.y; j < Eeff.ny; j++, jm++) {
                for (unsigned k = 0, km = start.z; k < Eeff.nz; k++, km++) {
                    Beta.p[i][j][k] = temp * Ne.p[im][jm][km];
                }
            }
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif
        for (int i = 0; i < Eeff.nx; i++) {
            int im = i + start.x;
            for (unsigned j = 0, jm = start.y; j < Eeff.ny; j++, jm++) {
                for (unsigned k = 0, km = start.z; k < Eeff.nz; k++, km++) {
                    Beta.p[i][j][k] = temp / (1 + mHalfDelta_t * Nu_c.p[i][j][k]) * Ne.p[im][jm][km];
                }
            }
        }
    }
}

void fdtd::initDensity() {
    MyDataF tmp = 2 * pow(4 * mDsFluid, 2); // 4 Maxwell grid size width    
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)
#endif   
    for (int i = 0; i < Ne.nx; i++) {
        for (int j = 0; j < Ne.ny; j++) {
            for (int k = 0; k < Ne.nz; k++) {
#ifdef DEBUG
                MyDataF sx, sy, sz;
                MyDataF px, py, pz;
                MyDataF ea;
                sx = (i - (int) mNeSrcPos.x) * mDsFluid;
                sy = (j - (int) mNeSrcPos.y) * mDsFluid;
                sz = (k - (int) mNeSrcPos.z) * mDsFluid;
                px = sx*sx;
                py = sy*sy;
                pz = sz*sz;
                ea = exp(-(px + py + pz) / tmp);
                Ne.p[i][j][k] = Ne0*ea;
#else
                Ne.p[i][j][k] = Ne0 * exp(-(pow((i - (int) mNeSrcPos.x) * mDsFluid, 2)
                        + pow((j - (int) mNeSrcPos.y) * mDsFluid, 2)
                        + pow((k - (int) mNeSrcPos.z) * mDsFluid, 2)) / tmp);
#endif
            }
        }
    }
}

void fdtd::createDensityArrays() {

    Vz.setName("Vz");
    Vx.setName("Vx");
    Vy.setName("Vy");
    Vx.create3DArray(Ex, 0.0);
    Vy.create3DArray(Ey, 0.0);
    Vz.create3DArray(Ez, 0.0);

    Exn.create3DArray(Ex, 0.0);
    Eyn.create3DArray(Ey, 0.0);
    Ezn.create3DArray(Ez, 0.0);

    unsigned nx = (mDomainEndIndex.x - mDomainStartIndex.x + 1) * mNeGridSize + 1;
    unsigned ny = (mDomainEndIndex.y - mDomainStartIndex.y + 1) * mNeGridSize + 1;
    unsigned nz = (mDomainEndIndex.z - mDomainStartIndex.z + 1) * mNeGridSize + 1;
    Ne.create3DArray(nx + mNeBoundWidth, ny + mNeBoundWidth, nz + mNeBoundWidth, 0.0);
    Eeff.create3DArray(nx, ny, nz, 0.0);
    Ne_pre.create3DArray(Ne, 0.0);

    createDensityRelatedArrays();
    Eeff.setName("erms");
    Ne.setName("Ne");
}
//#endif

void fdtd::createFieldArray() {
    // initial PML
    //    pml.InitialMuEps();
    //    pml.Initial(mMaxIndex.x, mMaxIndex.y, mMaxIndex.z, pmlWidth);
#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << "mMaxIndex.x=" << mMaxIndex.x << endl;
    cout << "mMaxIndex.y=" << mMaxIndex.y << endl;
    cout << "mMaxIndex.z=" << mMaxIndex.z << endl;
#endif
    Ez.create3DArray(mMaxIndex.x + 1, mMaxIndex.y + 1, mMaxIndex.z, 0);
    Ey.create3DArray(mMaxIndex.x + 1, mMaxIndex.y, mMaxIndex.z + 1, 0);
    Ex.create3DArray(mMaxIndex.x, mMaxIndex.y + 1, mMaxIndex.z + 1, 0);

    Hx.create3DArray(mMaxIndex.x + 1, mMaxIndex.y, mMaxIndex.z, 0);
    Hy.create3DArray(mMaxIndex.x, mMaxIndex.y + 1, mMaxIndex.z, 0);
    Hz.create3DArray(mMaxIndex.x, mMaxIndex.y, mMaxIndex.z + 1, 0);

    //coefficients
    Cexe.create3DArray(Ex, 0.0);
    Ceye.create3DArray(Ey, 0.0);
    Ceze.create3DArray(Ez, 0.0);
    Chxh.create3DArray(Hx, 0.0);
    Chyh.create3DArray(Hy, 0.0);
    Chzh.create3DArray(Hz, 0.0);
    Cexhy.create3DArray(Ex, 0.0);
    Cexhz.create3DArray(Ex, 0.0);
    Chxey.create3DArray(Hx, 0.0);
    Chxez.create3DArray(Hx, 0.0);
    Ceyhx.create3DArray(Ey, 0.0);
    Ceyhz.create3DArray(Ey, 0.0);
    Chyex.create3DArray(Hy, 0.0);
    Chyez.create3DArray(Hy, 0.0);
    Cezhy.create3DArray(Ez, 0.0);
    Cezhx.create3DArray(Ez, 0.0);
    Chzey.create3DArray(Hz, 0.0);
    Chzex.create3DArray(Hz, 0.0);

    Ez.setName("Ez");
    Ex.setName("Ex");
    Ey.setName("Ey");
    Hz.setName("Hz");
    Hx.setName("Hx");
    Hy.setName("Hy");
#ifdef DEBUG
    Cezhx.setName("cezhx");
    Cezhy.setName("cezhy");
    Ceze.setName("ceze");
#endif
}

void fdtd::setUp() {
    //Time step
    //    dt = 0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) +
    //            1.0 / (dz * dz)));
    mDt = mDx / 2 / C;

    //delay
    if (fdtd::SOURCE_SINE != mSrcType) {
        t0 = 4.5 * tw;
    } else {
        t0 = tw;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //temporary variables that often used
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    dtDivEps0DivDxyz = mDt / eps_0 / mDx / mDy / mDz;
    if (USE_DENSITY == mIsUseDensity) {
        mDsFluid = mDx / mNeGridSize;
        mHalfDelta_t = mDt / 2;
        mHalf_e = e / 2;
        dtDivEps0DivDx = mDt / eps_0 / mDx;
        dtDivEps0DivDy = mDt / eps_0 / mDy;
        dtDivEps0DivDz = mDt / eps_0 / mDz;
        e2Dt2Div4DivEps0DivMe = 0.25 * e * e * mDt * mDt / me / eps_0;
        eMDtDiv2DivEps0 = mHalf_e * mDt / eps_0;
        mHalfNeGridSize = mNeGridSize / 2;


        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // fluid variables
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        mMu_e = e / me / mNu_m; //3.7e-2;
        mMu_i = mMu_e / 100.0; //mu_e/mu_i ranges from 100 to 200
        mDe = mMu_e * 2 * 1.602e-19 / e; //
        mDa = mMu_i * 2 * 1.602e-19 / e; //
        MyDataF Dmax = mDe > mDa ? mDe : mDa;
        //Fine Time Step Size
        mDtFluid = 10 * mDt; //0.01 * mDsFluid * mDsFluid / 2 / Dmax;
        mNeSkipStep = mDtFluid / mDt;

        cout << "neSkipStep=" << mNeSkipStep << endl;
        cout << tw / mDt / mNeSkipStep << endl;
        createDensityArrays();
        //exit(0);
    }


    // source position    
#ifdef DEBUG
    mSourceIndex.setValue(mMaxIndex.x / 2, mMaxIndex.y / 2, mMaxIndex.z / 2);
#else
    mSourceIndex.setValue(mDomainStartIndex.x + (unsigned) (((float) (mDomainEndIndex.x - mDomainStartIndex.x)*2.25) / 3.0),
            mMaxIndex.y / 2, mMaxIndex.z / 2);
#endif    

    checkmax(mSourceIndex.x, 1, mMaxIndex.x);
    checkmax(mSourceIndex.y, 1, mMaxIndex.y);
    checkmax(mSourceIndex.z, 1, mMaxIndex.z);

    if (mNonPMLEndIndex.x < mNonPMLStartIndex.x)mNonPMLEndIndex.x = mNonPMLStartIndex.x + 1;
    if (mNonPMLEndIndex.y < mNonPMLStartIndex.y)mNonPMLEndIndex.y = mNonPMLStartIndex.y + 1;
    if (mNonPMLEndIndex.z < mNonPMLStartIndex.z)mNonPMLEndIndex.z = mNonPMLStartIndex.z + 1;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  COMPUTING FIELD UPDATE EQUATION COEFFICIENTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    initCoeficients();

    if (USE_DENSITY == mIsUseDensity) {
        mNeSrcPos.setValue((mSourceIndex.x - mDomainStartIndex.x) * mNeGridSize + mNeBoundWidth,
                (mSourceIndex.y - mDomainStartIndex.y - 3) * mNeGridSize + mNeBoundWidth,
                (mSourceIndex.z - mDomainStartIndex.z) * mNeGridSize + mNeBoundWidth);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // initial density
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        initDensity();

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Initial Coefficients for Density
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        initCoeffForDensity();
#ifdef DEBUG
        //Ne.save();
        Cvxex.setName("cv");
        Cezvz.setName("cezvz");
        Cvzez.setName("cvzzg");
#endif

    }

    if (mUseConnectingInterface) {
        mConnectingInterface.setLowerAndUpper(mDomainStartIndex, mDomainEndIndex);
        mConnectingInterface.setIncidentAngle(0.5 * M_PI, 0.0 * M_PI, 1.0 * M_PI, tw*C, C, mDt, mDx);
        mConnectingInterface.initCoefficients(mDx, mDt);
        mConnectingInterface.invalidate();
    } else {
        // initial coefficients at source position
        if (USE_DENSITY == mIsUseDensity) {
            Point nes((mSourceIndex.x - mDomainStartIndex.x) * mNeGridSize + mNeBoundWidth,
                    (mSourceIndex.y - mDomainStartIndex.y) * mNeGridSize + mNeBoundWidth,
                    (mSourceIndex.z - mDomainStartIndex.z) * mNeGridSize + mNeBoundWidth);
            Point bes((mSourceIndex.x - mDomainStartIndex.x) * mNeGridSize,
                    (mSourceIndex.y - mDomainStartIndex.y) * mNeGridSize,
                    (mSourceIndex.z - mDomainStartIndex.z) * mNeGridSize);
            initSourceCoeff(nes, bes);
        } else {
            initSourceCoeff();
        }
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  PML parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MyDataF sigmaMax = 0.601;
    MyDataF kappaMax = 0.5;
    MyDataF alphaMax = 0.1;
    int pmlOrder = 4;
    MyDataF alphaOrder = 1;

    cout << "alphaMax=" << alphaMax << endl;
    cout << "kappaMax=" << kappaMax << endl;
    cout << "sigmaMax=" << sigmaMax << endl;
    cout << "alphaOrder=" << alphaOrder << endl;
    cout << "pmlOrder=" << pmlOrder << endl;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  PML initials
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mPML.setCPMLRegion(mPMLWidth);
    mPML.createCPMLArrays(mMaxIndex.x, mMaxIndex.y, mMaxIndex.z);
    mPML.initCoefficientArrays(pmlOrder, alphaOrder, sigmaMax, kappaMax, alphaMax, epsR, mDt, mDx, mDy, mDz,
            Ceyhz, Cezhy, Chyez, Chzey,
            Cexhz, Cezhx, Chxez, Chzex,
            Ceyhx, Cexhy, Chyex, Chxey);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // print parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    printParam();
#ifdef DEBUG
    mOfNeCheck.open("nec.dat");
#endif
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SET CPML PARAMETERS IN EACH DIRECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void fdtd::compute() {

    unsigned n;
    Point capturePosition(mSourceIndex.x + 10, mSourceIndex.y + 10, mSourceIndex.z + 10);

    Point bes((mDomainStartIndex.x) * mNeGridSize,
            (mDomainStartIndex.y) * mNeGridSize,
            (mDomainStartIndex.z) * mNeGridSize);
    Point nes(bes.x - mNeBoundWidth, bes.y - mNeBoundWidth, bes.z - mNeBoundWidth);

#if DEBUG>=3
    if (mSourceIndex.y + 30 < mMaxIndex.y) {
        capturePosition.y = mSourceIndex.y + 30;
    }
#endif

    if (!capturePosition.checkMax(mMaxIndex.x, mMaxIndex.y, mMaxIndex.z)) {
        if (capturePosition.y >= mMaxIndex.y) {
            capturePosition.y = mMaxIndex.y / 2 + 3;
        }
        if (capturePosition.x >= mMaxIndex.x) {
            capturePosition.x = mMaxIndex.x / 2 + 3;
        }
        if (capturePosition.z >= mMaxIndex.z) {
            capturePosition.z = mMaxIndex.z / 2 + 3;
        }
    }
    cout << "capturePosition.x=" << capturePosition.x << endl;
    cout << "capturePosition.y=" << capturePosition.y << endl;
    cout << "capturePosition.z=" << capturePosition.z << endl;
    cout << "mSourceIndex.x=" << mSourceIndex.x << endl;
    cout << "mSourceIndex.y=" << mSourceIndex.y << endl;
    cout << "mSourceIndex.z=" << mSourceIndex.z << endl;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  BEGIN TIME STEP
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cout << endl;
    cout << "Begin time-stepping..." << endl;
#ifdef MATLAB_SIMULATION
    if (initMatlabSimulation() < 0) {
        return;
    }
#endif
    for (n = 1; n <= mTotalTimeSteps; ++n) {

        cout << "Ez at time step " << n << " at (" << capturePosition.x << ", " << capturePosition.y << ", " << capturePosition.z;
        cout << ") :  ";
        cout << Ez.p[capturePosition.x][capturePosition.y][capturePosition.z] << '\t';
        cout << Ez.p[mSourceIndex.x][capturePosition.y][mSourceIndex.z] << '\t';
        cout << Ez.p[capturePosition.x][mSourceIndex.y][mSourceIndex.z] << '\t';
        cout << Ez.p[mSourceIndex.x][mSourceIndex.y][capturePosition.z] << '\t';
        cout << Ez.p[mSourceIndex.x][mSourceIndex.y][mSourceIndex.z] << '\t';
        cout << n * mDt / 1e-9 << " ns" << "\t";
        cout << endl;
        //        cout << "Source Value:";
        //        cout << Ez.p[mSourceIndex.x][mSourceIndex.y][mSourceIndex.z] << '\t';
        //        cout << endl;
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+10][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+15][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+20][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+25][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+30][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+35][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+40][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+45][mSourceIndex.z] << endl;
        if (mUseConnectingInterface) {

            updateMagneitcFields();
            mConnectingInterface.updateMSource();
            mConnectingInterface.updateMConnect(Hx, Chxey, Chxez, Hy, Chyex, Chyez, Hz, Chzey, Chzex);
            mPML.updateCPML_M_Fields(Hx, Hy, Hz, Ex, Ey, Ez);

            updateElectricFields();
            if (USE_DENSITY == mIsUseDensity)updateVelocity();
            mConnectingInterface.updateESource(pSourceType->valueAtTime(n * mDt) * mAmplitude);
            mConnectingInterface.updateEConnect(Ex, Cexhy, Cexhz, Ey, Ceyhx, Ceyhz, Ez, Cezhy, Cezhx);
            mPML.updateCPML_E_Fields(Ex, Ey, Ez, Hx, Hy, Hz);
        } else {
            updateMagneitcFields();
            mPML.updateCPML_M_Fields(Hx, Hy, Hz, Ex, Ey, Ez);

            updateElectricFields();
            pSource->updateHardSource(Ex, Ey, Ez, n * mDt);
            if (USE_DENSITY == mIsUseDensity)updateVelocity();
            mPML.updateCPML_E_Fields(Ex, Ey, Ez, Hx, Hy, Hz);
        }

        if (USE_DENSITY == mIsUseDensity) {
            captureEFieldForEeff();

            if (n % mNeSkipStep == 0) {
                updateEeff();
#if DEBUG>=4
                Eeff.save(Eeff.nx / 2, 1, n, 1);
                Eeff.save(Eeff.ny / 2, 1, n, 2);
                Eeff.save(Eeff.nz / 2, 1, n, 3);
#endif
                if (!mUseConnectingInterface) {
                    updateSourceCoeff(nes, bes);
                }
                updateCollisionFrequency();
                updateDensity();
                updateCoeffWithDensity();
#if DEBUG>=4
                Ne.save(Ne.nz / 2, 1, n, 3);
                Ne.save(Ne.nz / 2, 1, n, 2);
                Ne.save(Ne.nz / 2, 1, n, 1);
                if (NULL != Nu_c.p) {
                    Nu_c.save(Nu_c.nz / 2, 1, n, 3);
                }
#else
                Ne.save(Ne.nx / 2, 1, n, 1);
                Ne.save(Ne.nz / 2, 1, n, 3);
                Eeff.save(Eeff.nx / 2, 1, n, 1);
                Ne.save(mSourceIndex.y*mNeGridSize, 1, n, 2);
#endif                
                Eeff.resetArray();
            }
        }
        if ((n % mSaveModulus) == 0) {

            //writeField(n);
            //Ez.save(mSourceIndex.x + 10, 1, n, 1);
            //Ez.save(mSourceIndex.y + 10, 1, n, 2);
            Ez.save(Ez.nz / 2, 1, n, 3);
            Ex.save(Ex.nz / 2, 1, n, 3);
            Ey.save(Ey.nz / 2, 1, n, 3);
#ifdef DEBUG 
            if (mUseConnectingInterface) {
                mConnectingInterface.saveEMInc(n);
                Ez.save(mMaxIndex.z / 2, 1, n, 3);
                Ez.save(mMaxIndex.y / 2, 1, n, 2);
                Ez.save(mMaxIndex.x / 2, 1, n, 1);
                Ex.save(mMaxIndex.z / 2, 1, n, 3);
                Ex.save(mMaxIndex.y / 2, 1, n, 2);
                Ex.save(mMaxIndex.x / 2, 1, n, 1);
                Ey.save(mMaxIndex.z / 2, 1, n, 3);
                Ey.save(mMaxIndex.y / 2, 1, n, 2);
                Ey.save(mMaxIndex.x / 2, 1, n, 1);
            }

            //            pml.Psi_exz_zp.setName("psi");
            //            pml.Psi_exz_zp.save(0, 1, n, 3);
            //            pml.Psi_exz_zp.save(4, 1, n, 3);
#ifdef DEBUG
            Cezhx.save(Ez.nz / 2, 1, n, 3);
            Cezhy.save(Ez.ny / 2, 1, n, 3);
            Ceze.save(Ez.nz / 2, 1, n, 3);
#endif
            //            if (USE_DENSITY == mIsUseDensity) {
            //                Cezvz.save(Ez.nz / 2, 1, n, 3);
            //                Cvzez.save(Ez.nz / 2, 1, n, 3);
            //            }

#else /* not define DEBUG*/
            Ez.save(mSourceIndex.x, 1, n, 1);
            Ez.save(mSourceIndex.y, 1, n, 2);
            Ez.save(mSourceIndex.z, 1, n, 3);
#endif   /*DEBUG*/

        }
#ifdef MATLAB_SIMULATION
        doMatlabSimulation();
#endif

    }
#ifdef MATLAB_SIMULATION
    finishMatlabSimulation();
#endif

    if (USE_DENSITY == mIsUseDensity) {
#ifdef DEBUG
        Ne.save();
#endif
    }
#ifdef DEBUG
    if (mOfNeCheck.is_open()) {

        mOfNeCheck.close();
    }
#endif
    //  END TIME STEP
    cout << "Done time-stepping..." << endl;

}

void fdtd::updateSource(unsigned n) {
    MyDataF source;
    switch (mSrcType) {
        case SOURCE_GAUSSIAN:
            source = mAmplitude * -2.0 * ((n * mDt - t0) / tw / tw)
                    * exp(-pow(((n * mDt - t0) / tw), 2)); //Differentiated Gaussian pulse
            //cout << "Gaussian source:" << source << endl;
            break;
        case SOURCE_SINE:
            // sine wave
            source = M_PI_TWO * mOmega * mAmplitude * cos((n * mDt - t0) * M_PI_TWO * mOmega);
            break;
        case ONE_SINE_PULSE:
            source = M_PI_TWO * mOmega * mAmplitude * sourceType::SinePulse(n * mDt - t0, mOmega, t_up, t_down);

            break;
        default:
            source = 0;
    }
    Ez.p[mSourceIndex.x][mSourceIndex.y][mSourceIndex.z] = Ez[mSourceIndex] + dtDivEps0DivDxyz * source;
    //cout<<"source="<<source<<"\t"<<amp<<"\t"<<n<<"\t"<<dt<<"\t"<<
    //        amp * -2.0 * ((n * dt - t0) / tw / tw) * exp(-pow(((n * dt - t0) / tw), 2))<<endl;

}
//Builds an object

void fdtd::buildObject() {

    buildSphere();
    //buildBrick();
}

//Builds a sphere (Sample code - NOt used in this program)

void fdtd::buildSphere() {

    MyDataF dist; //distance
    MyDataF rad = 2; //(MyDataF)mMaxIndex.x / 5.0; // sphere radius
    MyDataF sc = (MyDataF) mMaxIndex.x / 2.0; //sphere centre
    //MyDataF rad2 = 0.3; //(MyDataF)mMaxIndex.x / 5.0 - 3.0; // sphere radius

    unsigned i, j, k;
    Point *p;
    for (i = mDomainStartIndex.x; i <= mDomainEndIndex.x; ++i) {
        for (j = mDomainStartIndex.y; j <= mDomainEndIndex.y; ++j) {
            for (k = mDomainStartIndex.z; k <= mDomainEndIndex.z; ++k) {
                //compute distance form centre to the point i, j, k
                dist = sqrt((i + 0.5 - sc) * (i + 0.5 - sc) +
                        (j + 0.5 - sc) * (j + 0.5 - sc) +
                        (k + 0.5 - sc) * (k + 0.5 - sc));

                if (dist <= rad) {
                    p = new Point(i, j, k);
                    pSource->add(*p);
                }
            }
        }
    }

}

//Builds a dipole

void fdtd::buildBrick() {
    unsigned i, j, k;
    Point *p;
    const unsigned w = 2;
    Point lower(mMaxIndex.x / 2 - w, mMaxIndex.y / 2 - w, mMaxIndex.z / 2 - w);
    Point upper(mMaxIndex.x / 2 + w, mMaxIndex.y / 2 + w, mMaxIndex.z / 2 + w);

    for (i = lower.x; i <= upper.x; ++i) {

        for (j = lower.y; j <= upper.y; ++j) {

            for (k = lower.z; k <= upper.z; ++k) {

                p = new Point(i, j, k);
                pSource->add(*p);

            }
        }
    }
}

//creates a dielectric cube (yee cell) made up of the selected material

void fdtd::yeeCube(unsigned I, unsigned J, unsigned K, unsigned mType) {

    //set face 1 (for EX)

    ID1.p[I][J][K] = mType;
    ID1.p[I][J][K + 1] = mType;
    ID1.p[I][J + 1][K + 1] = mType;
    ID1.p[I][J + 1][K] = mType;

    //set face 2 (for EY)
    ID2.p[I][J][K] = mType;
    ID2.p[I + 1][J][K] = mType;
    ID2.p[I + 1][J][K + 1] = mType;
    ID2.p[I][J][K + 1] = mType;

    //set face 3 (for EZ)
    ID3.p[I][J][K] = mType;
    ID3.p[I + 1][J][K] = mType;
    ID3.p[I + 1][J + 1][K] = mType;
    ID3.p[I][J + 1][K] = mType;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Saving Output Data to files
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void fdtd::writeField(unsigned iteration) {
    unsigned i, j;
    ofstream out;
    stringstream ss;
    // form file name
    ss << "E_Field_ks_" << iteration << ".dat";
    string fileBaseName(ss.str());
    // open file
    out.open(fileBaseName.c_str());
    if (out.is_open()) {

        for (i = 0; i < mMaxIndex.x - 1; i++) {
            for (j = 0; j < mMaxIndex.y - 1; j++) { // |E|
                out << sqrt(pow(Ex.p[i][j][mKSource], 2) +
                        pow(Ey.p[i][j][mKSource], 2) + pow(Ez.p[i][j][mKSource], 2)) << '\t';
            }
            out << endl;
        }
        out.close();
    }
    stringstream sc;
    sc << "E_Field_j_" << iteration << ".dat";
    fileBaseName = sc.str();
    // open file
    out.open(fileBaseName.c_str());
    if (out.is_open()) {
        unsigned jsource = mKSource + 10;
        checkmax(jsource, 1, mMaxIndex.y);
        for (i = 0; i < mMaxIndex.x - 1; i++) {
            for (j = 0; j < mMaxIndex.z - 1; j++) { // |E|

                out << sqrt(pow(Ex.p[i][jsource][j], 2) +
                        pow(Ey.p[i][jsource][j], 2) + pow(Ez.p[i][jsource][j], 2)) << '\t';
            }
            out << endl;
        }
        out.close();
    }
}

//start up

void fdtd::startUp() {

    cout << "initializing(in Startup)..." << endl;
    createFieldArray();
    cout << "buildObject (in Startup)" << endl;
    buildObject();
    cout << "setUp (in Startup)" << endl;
    setUp();
    cout << "computing (in Startup)" << endl;
    compute();
    cout << "exit Startup" << endl;
}

void fdtd::printParam() {

    cout << "dx = " << mDx << endl;
    cout << "dy = " << mDy << endl;
    cout << "dz = " << mDz << endl;
    cout << "(mMaxIndex.x,mMaxIndex.y,mMaxIndex.z) = (" << mMaxIndex.x << "," << mMaxIndex.y << "," << mMaxIndex.z << ")" << endl;
    // time step increment
    cout << "dt = " << mDt << endl;

    //  Specify the Impulsive Source (Differentiated Gaussian) parameters
    cout << "tw = " << tw << endl; //pulse width
    cout << "t0 = " << t0 << endl; //delay    
    cout << "Amplitude = " << mAmplitude << endl; // Amplitude
    cout << "omega = " << mOmega << endl; // angle speed for sine wave

    //Specify the Time Step at which the data has to be saved for Visualization
    cout << "save_modulus = " << mSaveModulus << endl;

    //  Specify the dipole Boundaries(A cuboidal rode- NOT as a cylinder)
    cout << "(mStartIndex.x,mStartIndex.y, mStartIndex.z) = (" << mNonPMLStartIndex.x << ',' << mNonPMLStartIndex.y << ',' << mNonPMLStartIndex.z << ')' << endl;
    cout << "(mEndIndex.x,  mEndIndex.y, mEndIndex.z) = (" << mNonPMLEndIndex.x << ',' << mNonPMLEndIndex.y << ',' << mNonPMLEndIndex.z << ')' << endl;

    //Output recording point
    cout << "ksource = " << mKSource << endl;

    //  Specify the CPML Order and Other Parameters:
    cout << " PML Order = " << mPMLOrder << endl;
    cout << " Alpha Order = " << mAlphaOrder << endl;
    //#ifdef WITH_DENSITY 
    cout << "neGridSize=" << mNeGridSize << endl;
    cout << "neSkipStep=" << mNeSkipStep << endl;
    cout << "DtFluid=" << mDtFluid << endl;
    cout << "DsFluid=" << mDsFluid << endl;
    cout << "mu_i=" << mMu_i << endl;
    cout << "mu_e=" << mMu_e << endl;
    cout << "Da=" << mDa << endl;
    cout << "De=" << mDe << endl;
    //#endif
    cout << endl << "Time step = " << mDt << endl;
    cout << endl << "Number of steps = " << mTotalTimeSteps << endl;
    cout << endl << "Total Simulation time = " << mTotalTimeSteps * mDt << " Seconds" << endl;
}

void fdtd::setSourceType(sourceType* pSrcType) {

    pSourceType = pSrcType;
}

void fdtd::defineSineSource(MyDataF omega_) {

    mSrcType = SOURCE_SINE;
    mOmega = omega_;
}

void fdtd::updateHx() {
    int i, j, k;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) //shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif
    for (k = 0; k < Hx.nz; ++k) {
        for (j = 0; j < Hx.ny; ++j) {
            for (i = 0; i < Hx.nx; ++i) {

                Hx.p[i][j][k] = Chxh.p[i][j][k] * Hx.p[i][j][k] +
                        Chxez.p[i][j][k]*(Ez.p[i][j + 1][k] - Ez.p[i][j][k]) +
                        Chxey.p[i][j][k]*(Ey.p[i][j][k + 1] - Ey.p[i][j][k]);

                //#if (DEBUG>=4&&!_OPENMP)
                //                Hx.isValid(i, j, k);
                //                Hx.whenLargerThan(i, j, k, 1e30 / 188, NULL);
                //#endif

            }
        }
    }
}

void fdtd::updateHy() {
    int i, j, k;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) //shared(Hy,Ez,Ex,pml,DA,DB,dx,dz)
#endif
    for (k = 0; k < Hy.nz; ++k) {
        for (i = 0; i < Hy.nx; ++i) {
            for (j = 0; j < Hy.ny; ++j) {

                Hy.p[i][j][k] = Chyh.p[i][j][k] * Hy.p[i][j][k] +
                        Chyez.p[i][j][k]*(Ez.p[i + 1][j][k] - Ez.p[i][j][k]) +
                        Chyex.p[i][j][k]*(Ex.p[i][j][k + 1] - Ex.p[i][j][k]);

                //#if (DEBUG>=4&&!_OPENMP)
                //                Hy.isValid(i, j, k);
                //                Hy.whenLargerThan(i, j, k, 1e30 / 188, NULL);
                //#endif

            }
        }
    }
}

void fdtd::updateHz() {
    int i, j, k;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  UPDATE Hz
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (k = 0; k < Hz.nz; ++k) {
        for (i = 0; i < Hz.nx; ++i) {
            for (j = 0; j < Hz.ny; ++j) {

                Hz.p[i][j][k] = Chzh.p[i][j][k] * Hz.p[i][j][k] +
                        (Ey.p[i + 1][j][k] - Ey.p[i][j][k]) * Chzey.p[i][j][k] +
                        (Ex.p[i][j + 1][k] - Ex.p[i][j][k]) * Chzex.p[i][j][k];
                //#if (DEBUG>=4&&!_OPENMP)
                //                Hz.isValid(i, j, k);
                //                Hz.whenLargerThan(i, j, k, 1e30 / 188, NULL);
                //#endif

            }
        }
    }
}

void fdtd::updateEx() {
    int i, j, k;
    if (USE_DENSITY == mIsUseDensity)Exn = Ex;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)  private(i,j,k)
#endif // _OPENMP
    for (k = 1; k < mMaxIndex.z; ++k) {
        for (j = 1; j < mMaxIndex.y; ++j) {
            for (i = 0; i < mMaxIndex.x; ++i) {
                Ex.p[i][j][k] = Cexe.p[i][j][k] * Ex.p[i][j][k] +
                        (Hz.p[i][j][k] - Hz.p[i][j - 1][k]) * Cexhz.p[i][j][k] +
                        (Hy.p[i][j][k] - Hy.p[i][j][k - 1]) * Cexhy.p[i][j][k];
            }
        }
    }
    ///////////////////////////////////////////
    // add velocity
    ///////////////////////////////////////////
    if (USE_DENSITY == mIsUseDensity) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)  private(i,j,k)
#endif // _OPENMP
        for (k = 1; k < mMaxIndex.z; ++k) {
            for (j = 1; j < mMaxIndex.y; ++j) {
                for (i = 0; i < mMaxIndex.x; ++i) {
                    Ex.p[i][j][k] += Cexvx.p[i][j][k] * Vx.p[i][j][k];
#if DEBUG>=4
                    //                Ex.whenLargerThan(i, j, k, 1e30, NULL);
                    if (isnan(Ex.p[i][j][k])) {
                        Ex.p[i][j][k] += 0.0;
                    }
#endif

                }
            }
        }
    }
}

void fdtd::updateVx() {
    int i, j, k;
    if (fdtd::SOURCE_SINE != mSrcType) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) 
#endif
        for (i = 0; i < Vx.nx; i++) {
            for (j = 0; j < Vx.ny; j++) {
                for (k = 0; k < Vx.nz; k++) {
                    Vx.p[i][j][k] = Cvxvx.p[i][j][k] * Vx.p[i][j][k] - Cvxex.p[i][j][k] * (Exn.p[i][j][k] + Ex.p[i][j][k]);
                }
            }
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) 
#endif
        for (i = 0; i < Vx.nx; i++) {
            for (j = 0; j < Vx.ny; j++) {
                for (k = 0; k < Vx.nz; k++) {
                    Vx.p[i][j][k] = mAlpha * Vx.p[i][j][k] - mCvxex * (Exn.p[i][j][k] + Ex.p[i][j][k]);
                }
            }
        }
    }
}

void fdtd::updateEy() {
    int i, j, k;
    if (USE_DENSITY == mIsUseDensity)Eyn = Ey;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) 
#endif
    for (k = 1; k < mMaxIndex.z; ++k) {
        for (i = 1; i < mMaxIndex.x; ++i) {
            for (j = 0; j < mMaxIndex.y; ++j) {
                Ey.p[i][j][k] = Ceye.p[i][j][k] * Ey.p[i][j][k] +
                        (Hz.p[i][j][k] - Hz.p[i - 1][j][k]) * Ceyhz.p[i][j][k] +
                        (Hx.p[i][j][k] - Hx.p[i][j][k - 1]) * Ceyhx.p[i][j][k];
            }
        }
    }

    ////////////////////////////////////////
    // add velocity
    ///////////////////////////////////////
    if (USE_DENSITY == mIsUseDensity) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) 
#endif
        for (k = 1; k < mMaxIndex.z; ++k) {
            for (i = 1; i < mMaxIndex.x; ++i) {
                for (j = 0; j < mMaxIndex.y; ++j) {
                    Ey.p[i][j][k] += Ceyvy.p[i][j][k] * Vy.p[i][j][k];
#if (DEBUG>=4&&!_OPENMP)
                    //                Ey.whenLargerThan(i, j, k, 1e30, NULL);
                    //                Ey.isValid(i, j, k);
                    if (isnan(Ey.p[i][j][k])) {

                        Ey.p[i][j][k] += 0.0;
                    }
#endif
                }
            }
        }
    }
}

void fdtd::updateVy() {
    int i, j, k;
    if (fdtd::SOURCE_SINE != mSrcType) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) 
#endif
        for (i = 0; i < Vy.nx; i++) {
            for (j = 0; j < Vy.ny; j++) {
                for (k = 0; k < Vy.nz; k++) {
                    Vy.p[i][j][k] = Cvyvy.p[i][j][k] * Vy.p[i][j][k] - Cvyey.p[i][j][k] * (Eyn.p[i][j][k] + Ey.p[i][j][k]);
                }
            }
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) 
#endif
        for (i = 0; i < Vy.nx; i++) {
            for (j = 0; j < Vy.ny; j++) {
                for (k = 0; k < Vy.nz; k++) {
                    Vy.p[i][j][k] = mAlpha * Vy.p[i][j][k] - mCvyey * (Eyn.p[i][j][k] + Ey.p[i][j][k]);
                }
            }
        }
    }
}

void fdtd::updateEz() {
    int i, j, k;
    if (USE_DENSITY == mIsUseDensity)Ezn = Ez;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) 
#endif

    for (i = 1; i < mMaxIndex.x; ++i) {
        for (k = 0; k < mMaxIndex.z; ++k) {
            for (j = 1; j < mMaxIndex.y; ++j) {
                Ez.p[i][j][k] = Ceze.p[i][j][k] * Ez.p[i][j][k] +
                        (Hy.p[i][j][k] - Hy.p[i - 1][j][k]) * Cezhy.p[i][j][k]+
                        (Hx.p[i][j][k] - Hx.p[i][j - 1][k]) * Cezhx.p[i][j][k];
            }
        }
    }
    //===========================
    // ADD VELOCITY
    //===========================
    if (USE_DENSITY == mIsUseDensity) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) 
#endif
        for (i = 1; i < mMaxIndex.x; ++i) {
            for (k = 0; k < mMaxIndex.z; ++k) {
                for (j = 1; j < mMaxIndex.y; ++j) {
                    Ez.p[i][j][k] += Cezvz.p[i][j][k] * Vz.p[i][j][k];
#if DEBUG>=4
                    //Ez.whenLargerThan(i, j, k, 1e30, NULL);
                    if (isnan(Ez.p[i][j][k])) {
                        Ez.p[i][j][k] += 0.0;
                    }
#endif
                }
            }
        }
    }
}

void fdtd::updateVz() {
    int i, j, k;
    if (fdtd::SOURCE_SINE != mSrcType) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) 
#endif
        for (i = 0; i < Vz.nx; i++) {
            for (j = 0; j < Vz.ny; j++) {
                for (k = 0; k < Vz.nz; k++) {
                    Vz.p[i][j][k] = Cvzvz.p[i][j][k] * Vz.p[i][j][k] - Cvzez.p[i][j][k] * (Ezn.p[i][j][k] + Ez.p[i][j][k]);
                }
            }
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) 
#endif
        for (i = 0; i < Vz.nx; i++) {
            for (j = 0; j < Vz.ny; j++) {
                for (k = 0; k < Vz.nz; k++) {
                    Vz.p[i][j][k] = mAlpha * Vz.p[i][j][k] - mCvzez * (Ezn.p[i][j][k] + Ez.p[i][j][k]);
                }
            }
        }
    }
}

void fdtd::initCoeficients() {
    //////////////////////////
    // E Field coefficients
    //////////////////////////
    MyDataF dtDivEps0DivDz, dtDivEps0DivDy, dtDivEps0DivDx;
    dtDivEps0DivDz = mDt / eps_0 / mDz;
    dtDivEps0DivDy = mDt / eps_0 / mDy;
    dtDivEps0DivDx = mDt / eps_0 / mDx;
    for (unsigned i = 0; i < Ex.nx; i++) {
        for (unsigned j = 0; j < Ex.ny; j++) {
            for (unsigned k = 0; k < Ex.nz; k++) {
                Cexhy.p[i][j][k] = -dtDivEps0DivDz;
                Cexhz.p[i][j][k] = dtDivEps0DivDy;
                Cexe.p[i][j][k] = 1;
            }
        }
    }

    for (unsigned i = 0; i < Ey.nx; i++) {
        for (unsigned j = 0; j < Ey.ny; j++) {
            for (unsigned k = 0; k < Ey.nz; k++) {
                Ceyhx.p[i][j][k] = dtDivEps0DivDz;
                Ceyhz.p[i][j][k] = -dtDivEps0DivDx;
                Ceye.p[i][j][k] = 1;
            }
        }
    }
    for (unsigned i = 0; i < Ez.nx; i++) {
        for (unsigned j = 0; j < Ez.ny; j++) {
            for (unsigned k = 0; k < Ez.nz; k++) {
                Cezhx.p[i][j][k] = -dtDivEps0DivDy;
                Cezhy.p[i][j][k] = dtDivEps0DivDx;
                Ceze.p[i][j][k] = 1;
            }
        }
    }

    //////////////////////////
    // M Field coefficients
    //////////////////////////
    MyDataF dtDivMu0DivDx, dtDivMu0DivDy, dtDivMu0DivDz;
    dtDivMu0DivDx = mDt / mu_0 / mDx;
    dtDivMu0DivDy = mDt / mu_0 / mDy;
    dtDivMu0DivDz = mDt / mu_0 / mDz;
    for (unsigned i = 0; i < Hx.nx; i++) {
        for (unsigned j = 0; j < Hx.ny; j++) {
            for (unsigned k = 0; k < Hx.nz; k++) {
                Chxey.p[i][j][k] = dtDivMu0DivDz;
                Chxez.p[i][j][k] = -dtDivMu0DivDy;
                Chxh.p[i][j][k] = 1;
            }
        }
    }
    for (unsigned i = 0; i < Hy.nx; i++) {
        for (unsigned j = 0; j < Hy.ny; j++) {
            for (unsigned k = 0; k < Hy.nz; k++) {
                Chyex.p[i][j][k] = -dtDivMu0DivDz;
                Chyez.p[i][j][k] = dtDivMu0DivDx;
                Chyh.p[i][j][k] = 1;
            }
        }
    }
    for (unsigned i = 0; i < Hz.nx; i++) {
        for (unsigned j = 0; j < Hz.ny; j++) {
            for (unsigned k = 0; k < Hz.nz; k++) {

                Chzex.p[i][j][k] = dtDivMu0DivDy;
                Chzey.p[i][j][k] = -dtDivMu0DivDx;
                Chzh.p[i][j][k] = 1;
            }
        }
    }
}

void fdtd::updateMagneitcFields() {

    updateHx();
    updateHy();
    updateHz();
}

void fdtd::intSourceSinePulse(MyDataF t_0, MyDataF omega_, MyDataF tUp, MyDataF tDown, MyDataF amplitude) {

    t0 = t_0;
    mOmega = omega_;
    t_up = tUp;
    t_down = tDown;
    mAmplitude = amplitude;
}

// =================================================================
// MATLAB SIMULATION
// =================================================================
#ifdef MATLAB_SIMULATION

int fdtd::initMatlabSimulation() {
    if (!Ez.isMatlabEngineStarted()) {
        Ez.initMatlabEngine();
    }
    if (Ez.isMatlabEngineStarted()) {
        Ez.preparePlotting();

        return 0;
    }
    return -1;

}

void fdtd::doMatlabSimulation() {
    if (Ez.isMatlabEngineStarted()) {

        Ez.plotArrays();
    }
}

void fdtd::finishMatlabSimulation() {
    if (Ez.isMatlabEngineStarted()) {

        Ez.clearMatlabEngineArray();
        Ez.closeMatlabEngine();
    }
}
#endif
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// END OF PROGRAM CPMLFDTD3D
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void fdtd::setSrcType(int srcType) {

    mSrcType = srcType;
}

void fdtd::desideDomainZone() {
    //  Specify the dipole size 

    mNonPMLStartIndex.setValue(mPMLWidth, mPMLWidth, mPMLWidth);
    mNonPMLEndIndex.setValue(mMaxIndex.x - mPMLWidth, mMaxIndex.y - mPMLWidth, mMaxIndex.z - mPMLWidth);
    mDomainStartIndex.setValue(mNonPMLStartIndex.x + mAirBufferWidth,
            mNonPMLStartIndex.y + mAirBufferWidth,
            mNonPMLStartIndex.z + mAirBufferWidth);
    mDomainEndIndex.setValue(mNonPMLEndIndex.x - mAirBufferWidth,
            mNonPMLEndIndex.y - mAirBufferWidth,
            mNonPMLEndIndex.z - mAirBufferWidth);
}

void fdtd::initSourceCoeff() {
    // initial coefficients at source position
    switch (pSource->getDirection()) {
        case source::Z:
            pSource->initCoefficients(Ceze, Cezhy, Cezhx, mDz, mDy, mDx, mDt);
            break;
        case source::X:
            pSource->initCoefficients(Cexe, Cexhz, Cexhy, mDx, mDz, mDy, mDt);
            break;
        case source::Y:
            pSource->initCoefficients(Ceye, Ceyhx, Ceyhz, mDy, mDx, mDz, mDt);
            break;
    }
}

void fdtd::updateSourceCoeff() {
    // initial coefficients at source position
    switch (pSource->getDirection()) {
        case source::Z:
            pSource->updateCoefficients(Ceze, Cezhy, Cezhx, mDz, mDy, mDx, mDt);
            break;
        case source::X:
            pSource->updateCoefficients(Cexe, Cexhz, Cexhy, mDx, mDz, mDy, mDt);
            break;
        case source::Y:
            pSource->updateCoefficients(Ceye, Ceyhx, Ceyhz, mDy, mDx, mDz, mDt);
            break;
    }
}

void fdtd::initSourceCoeff(const Point& nes, const Point & bes) {
    // initial coefficients at source position
    if (fdtd::SOURCE_SINE != mSrcType) {
        switch (pSource->getDirection()) {
            case source::Z:
                pSource->initCoefficients(Ceze, Cezhy, Cezhx, Cexvx, Beta, Ne, Nu_c, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
            case source::X:
                pSource->initCoefficients(Cexe, Cexhz, Cexhy, Ceyvy, Beta, Ne, Nu_c, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
            case source::Y:
                pSource->initCoefficients(Ceye, Ceyhx, Ceyhz, Cezvz, Beta, Ne, Nu_c, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
        }
    } else {
        switch (pSource->getDirection()) {
            case source::Z:
                pSource->initCoefficients(Ceze, Cezhy, Cezhx, Cexvx, Beta, Ne, mNu_m, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
            case source::X:
                pSource->initCoefficients(Cexe, Cexhz, Cexhy, Ceyvy, Beta, Ne, mNu_m, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
            case source::Y:
                pSource->initCoefficients(Ceye, Ceyhx, Ceyhz, Cezvz, Beta, Ne, mNu_m, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
        }
    }
}

void fdtd::updateSourceCoeff(const Point& nes, const Point & bes) {
    // initial coefficients at source position
    if (fdtd::SOURCE_SINE != mSrcType) {
        switch (pSource->getDirection()) {
            case source::Z:
                pSource->updateCoefficients(Ceze, Cezhy, Cezhx, Cexvx, Beta, Ne, Nu_c, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
            case source::X:
                pSource->updateCoefficients(Cexe, Cexhz, Cexhy, Ceyvy, Beta, Ne, Nu_c, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
            case source::Y:
                pSource->updateCoefficients(Ceye, Ceyhx, Ceyhz, Cezvz, Beta, Ne, Nu_c, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
        }
    } else {
        switch (pSource->getDirection()) {
            case source::Z:
                pSource->updateCoefficients(Ceze, Cezhy, Cezhx, Cexvx, Beta, Ne, mNu_m, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
            case source::X:
                pSource->updateCoefficients(Cexe, Cexhz, Cexhy, Ceyvy, Beta, Ne, mNu_m, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
            case source::Y:
                pSource->updateCoefficients(Ceye, Ceyhx, Ceyhz, Cezvz, Beta, Ne, mNu_m, nes, bes, mNeGridSize, mDx, mDy, mDz, mDt);
                break;
        }
    }
}
