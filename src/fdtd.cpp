

#include <math.h>
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
#include "source.h"
#include "InonizationFormula.h"
#include "Point.h"

extern MyDataF epsR;
//extern MyDataF dt, dx, dy, dz;
extern MyDataF T;

using namespace std;

void checkmax(unsigned &u_2check, unsigned max, unsigned min) {
    if (u_2check >= max || u_2check < min)u_2check = (min + max) / 2;
}

#ifdef WITH_DENSITY

fdtd::fdtd(unsigned _totalTimeSteps, unsigned _imax, unsigned _jmax, unsigned _kmax,
        MyDataF _tw, MyDataF _dx, MyDataF _dy, MyDataF _dz,
        MyDataF _amp, unsigned _savemodulus, unsigned _ksource,
        unsigned _m, unsigned _ma, unsigned pmlw, unsigned _neGrid)
: mTotalTimeSteps(_totalTimeSteps), mMaxIndex(_imax, _jmax, _kmax)
, tw(_tw), mDx(_dx), mDy(_dy), mDz(_dz)
, mAmplitude(_amp), save_modulus(_savemodulus), mKSource(_ksource)
, mPMLOrder(_m), mAlphaOrder(_ma), mPMLWidth(pmlw)
, mNeGridSize(_neGrid)
, Ne0(DEFAULT_DENSITY_MAX)
, mSrcType(SOURCE_GAUSSIAN)
, mEpsilon(NULL), mSigma(NULL), mMu(NULL), CA(NULL), CB(NULL) {
}
#else

fdtd::fdtd(unsigned _totalTimeSteps, unsigned _imax, unsigned _jmax, unsigned _kmax,
        MyDataF _tw, MyDataF _dx, MyDataF _dy, MyDataF _dz,
        MyDataF _amp, unsigned _savemodulus, unsigned _ksource,
        unsigned _m, unsigned _ma, unsigned pmlw)
: mTotalTimeSteps(_totalTimeSteps), mMaxIndex(_imax, _jmax, _kmax)
, tw(_tw), mDx(_dx), mDy(_dy), mDz(_dz)
, mAmplitude(_amp), save_modulus(_savemodulus), mKSource(_ksource)
, mPMLOrder(_m), mAlphaOrder(_ma), mPMLWidth(pmlw)
, mSrcType(SOURCE_GAUSSIAN)
, mEpsilon(NULL), mSigma(NULL), mMu(NULL), CA(NULL), CB(NULL) {
}
#endif

fdtd::~fdtd(void) {
    if (mEpsilon != NULL)delete []mEpsilon;
    if (mSigma != NULL)delete []mSigma;
    if (mMu != NULL)delete[]mMu;
    if (CA != NULL)delete[]CA;
    if (CB != NULL)delete[]CB;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef WITH_DENSITY
const MyDataF fdtd::mNeutralGasDensity = 2.44e19;

void fdtd::setPlasmaParam(MyDataF _rei, MyDataF _vm, MyDataF _p, int _ftype) {
    mRei = _rei;
    mNiu_m = _vm;
    mAirPressure = _p;
    mNiuType = _ftype;
}

void fdtd::integerEeff() {
    unsigned i, j, k;
    unsigned io, jo, ko;
    //MyDataF vxIJK,vyIJK;
    for (i = mStartIndex.x, io = mStartIndex.x * mNeGridSize; i <= mEndIndex.x; i++, io += mNeGridSize) {
        for (j = mStartIndex.y, jo = mStartIndex.y * mNeGridSize; j <= mEndIndex.y; j++, jo += mNeGridSize) {
            for (k = mStartIndex.z, ko = mStartIndex.z * mNeGridSize; k <= mEndIndex.z; k++, ko += mNeGridSize) {
                Erms.p[io][jo][ko] = mPMLOrder / e * sqrt(Erms.p[io][jo][ko] / mDeltaTimeFluid * Nu_c.p[i][j][k] / 2);
            }
        }
    }
}

void fdtd::updateErms(void) {
    unsigned i, j, k;
    unsigned io, jo, ko;
    switch (mSrcType) {
        case SOURCE_GAUSSIAN:
            MyDataF vxIJK, vyIJK;
            for (i = mStartIndex.x, io = mStartIndex.x * mNeGridSize; i <= mEndIndex.x; i++, io += mNeGridSize) {
                for (j = mStartIndex.y, jo = mStartIndex.y * mNeGridSize; j <= mEndIndex.y; j++, jo += mNeGridSize) {
                    for (k = mStartIndex.z, ko = mStartIndex.z * mNeGridSize; k <= mEndIndex.z; k++, ko += mNeGridSize) {
                        vxIJK = (Vx.p[i - 1][j][k - 1] + Vx.p[i + 1][j][k - 1] + Vx.p[i - 1][j][k + 1] + Vx.p[i + 1][j][k + 1]) / 4;
                        vyIJK = (Vy.p[i][j - 1][k - 1] + Vy.p[i][j + 1][k - 1] + Vy.p[i][j - 1][k + 1] + Vy.p[i][j + 1][k + 1]) / 4;
                        Erms.p[io][jo][ko] += (Vz.p[i][j][k] * Vz.p[i][j][k] + vxIJK * vxIJK + vyIJK * vyIJK) * mDt;
                    }
                }
            }
            break;
        case SOURCE_SINE:
        default:
            MyDataF exIJK, eyIJK;
            for (i = mStartIndex.x, io = mStartIndex.x * mNeGridSize; i <= mEndIndex.x; i++, io += mNeGridSize) {
                for (j = mStartIndex.y, jo = mStartIndex.y * mNeGridSize; j <= mEndIndex.y; j++, jo += mNeGridSize) {
                    for (k = mStartIndex.z, ko = mStartIndex.z * mNeGridSize; k <= mEndIndex.z; k++, ko += mNeGridSize) {
                        exIJK = (Ex.p[i - 1][j][k - 1] + Ex.p[i + 1][j][k - 1] + Ex.p[i - 1][j][k + 1] + Ex.p[i + 1][j][k + 1]) / 4;
                        eyIJK = (Ey.p[i][j - 1][k - 1] + Ey.p[i][j + 1][k - 1] + Ey.p[i][j - 1][k + 1] + Ey.p[i][j + 1][k + 1]) / 4;
                        Erms.p[io][jo][ko] = sqrt(Ez.p[i][j][k] * Ez.p[i][j][k] + exIJK * exIJK + eyIJK * eyIJK);
                    }
                }
            }
    }
}

void fdtd::updateCollisionFrequency() {
    int i, j, k;
    //unsigned io, jo, ko;
    MyDataF EeffDivP;
    MyDataF DivParam = 100 * mAirPressure * 133.3;
    MyDataF C1 = 5.20e8 * mAirPressure;
    MyDataF C2 = 2.93e8 * mAirPressure;
    MyDataF C3 = 3.24e8 * mAirPressure;
    //    for (i = mStartIndex.x, io = mStartIndex.x * neGrid; i <= mEndIndex.x; i++, io += neGrid) {
    //        for (j = mStartIndex.y, jo = mStartIndex.y * neGrid; j <= mEndIndex.y; j++, jo += neGrid) {
    //            for (k = mStartIndex.z, ko = mStartIndex.z * neGrid; k <= mEndIndex.z; k++, ko += neGrid) {
    //                EeffDivP = Erms.p[io][jo][ko] / DivParam;
    //                if (EeffDivP >= 120) {
    //                    Nu_c.p[i][j][k] = C1 * sqrt(EeffDivP);
    //                } else if (EeffDivP >= 54) {
    //                    Nu_c.p[i][j][k] = C2 * EeffDivP / (1 + 0.041 * EeffDivP);
    //                } else {
    //                    Nu_c.p[i][j][k] = C3 * EeffDivP / (1 + 0.04 * EeffDivP);
    //                }
    //            }
    //        }
    //    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k,EeffDivP) //shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif
    for (i = 0; i < Nu_c.nx; i++) {
        for (j = 0; j < Nu_c.ny; j++) {
            for (k = 0; k < Nu_c.nz; k++) {
                EeffDivP = Erms.p[i][j][k] / DivParam;
                if (EeffDivP >= 120) {
                    Nu_c.p[i][j][k] = C1 * sqrt(EeffDivP);
                } else if (EeffDivP >= 54) {
                    Nu_c.p[i][j][k] = C2 * EeffDivP / (1 + 0.041 * EeffDivP);
                } else {
                    Nu_c.p[i][j][k] = C3 * EeffDivP / (1 + 0.04 * EeffDivP);
                }
            }
        }
    }
}

void fdtd::interpErms() {
    unsigned is, js, ks;
    unsigned in, jn, kn;
    unsigned i, j, k;
    unsigned im, jm, km;
    unsigned iu, ju, ku;
    unsigned ngred = mNeGridSize * mNeGridSize*mNeGridSize;
    for (is = mStartIndex.x, in = mStartIndex.x + mNeGridSize; is < mEndIndex.x; is = in, in += mNeGridSize) {
        for (js = mStartIndex.y, jn = mStartIndex.y + mNeGridSize; js < mEndIndex.y; js = jn, jn += mNeGridSize) {
            for (ks = mStartIndex.z, kn = mStartIndex.z + mNeGridSize; ks < mEndIndex.z; ks = kn, kn += mNeGridSize) {
                // integrate Erms
                for (i = is, im = 0, iu = mNeGridSize; i < in; i++, im++, iu--) {
                    for (j = js, jm = 0, ju = mNeGridSize; j < jn; j++, jm++, ju--) {
                        for (k = ks, km = 0, ku = mNeGridSize; k < kn; k++, km++, ku--) {
                            //if (im == 0 && jm == 0 && km == 0)continue;
                            Erms.p[i][j][k] = (iu * ju * ku * Erms.p[is][js][ks] + im * ju * ku * Erms.p[in][js][ks] +
                                    im * jm * ku * Erms.p[in][jn][ks] + im * jm * km * Erms.p[in][jn][kn] +
                                    iu * jm * ku * Erms.p[is][jn][ks] + iu * jm * km * Erms.p[is][jn][kn] +
                                    iu * ju * km * Erms.p[is][js][kn] + im * ju * km * Erms.p[in][js][kn]) / ngred;
                        }
                    }
                }
            }
        }
    }
}

void fdtd::calculateIonizationParameters(int i, int j, int k, MyDataF &va, MyDataF &vi, MyDataF &Deff) {
    MyDataF Eeff;
    switch (mSrcType) {
        case SOURCE_GAUSSIAN:
            Eeff = Erms.p[i][j][k] / 100; //convert to V/cm
            break;
        default:
            Eeff = Erms.p[i][j][k] / 100 * pow(1 / (1 + mOmega * mOmega / mNiu_m / mNiu_m), 0.5);
    }

    switch (mNiuType) {
        case MORROW_AND_LOWKE:
            Niu_MorrowAndLowke(&vi, &va, Eeff, mNeutralGasDensity);
            break;
        case NIKONOV:
            Niu_Nikonov(&vi, &va, Eeff, mAirPressure);
            break;
        case KANG:
            Niu_Kang(&vi, &va, Eeff);
            break;
        case ALI:
        default:
            Niu_Ali(&vi, &va, Eeff, mAirPressure);
    }
    if (Ne.p[i][j][k] < 1) {
        Deff = mDe;
    } else {
        // tau_m = eps_0 / (e * Ne.p[i][j][k] * (mu_e + mu_i));
        MyDataF kasi = vi * eps_0 / (e * Ne.p[i][j][k] * (mMu_e + mMu_i));
        Deff = (kasi * mDe + mDa) / (kasi + 1);
    }
}

/************************************************************************/
/* update density                                                                     */

/************************************************************************/
void fdtd::updateDensity(void) {

    int i, j, k, mt = 1;
    MyDataF Ne_ijk, Neip1, Neim1, Nejm1, Nejp1, Nekp1, Nekm1;
    MyDataF Deff;
    MyDataF maxvi = 0, minvi = 0;
    MyDataF vi, va;

    unsigned ci = 0, cj = 0, ck = 0;

    Ne_pre = Ne;

#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) \
        schedule(dynamic) private(i,j,k,Ne_ijk, Neip1, Neim1, Nejm1, Nejp1, Nekp1, Nekm1,vi,va,Deff)
#endif
    for (i = mt; i < Ne.nx - mt; i++) {
        for (j = mt; j < Ne.ny - mt; j++) {
            for (k = mt; k < Ne.nz - mt; k++) {

                Ne_ijk = Ne_pre.p[i][j][k];
                Neip1 = Ne_pre.p[i + 1][j][k];
                Neim1 = Ne_pre.p[i - 1][j][k];
                Nejp1 = Ne_pre.p[i][j + 1][k];
                Nejm1 = Ne_pre.p[i][j - 1][k];
                Nekp1 = Ne_pre.p[i][j][k + 1];
                Nekm1 = Ne_pre.p[i][j][k - 1];

                calculateIonizationParameters(i, j, k, va, vi, Deff);

                Ne.p[i][j][k] = (Ne_ijk * (1 + mDeltaTimeFluid * vi) + Deff * mDeltaTimeFluid *
                        (Neip1 + Neim1 + Nejp1 + Nejm1 + Nekp1 + Nekm1 - 6 * Ne_ijk)
                        / mDeltaSizeFluid / mDeltaSizeFluid / mDeltaSizeFluid) / (1 + mDeltaTimeFluid * (va + mRei * Ne_ijk));
                if (vi > maxvi) {
                    maxvi = vi;
                    ci = i;
                    cj = j;
                    ck = k;
                }
                if (vi < minvi) minvi = vi;
            }
        }
    }
    wallCircleBound(Ne);
    cout << Ne.p[Ne.nx / 2][Ne.ny / 2][Ne.nz / 2] << '\t';
    cout << maxvi << '\t' << minvi << '\t' << Ne.p[ci][cj][ck] << '\t' << Erms.p[ci][cj][ck] << '\t';
}

void fdtd::updateVeloity(void) {
    //return 0;
}

void fdtd::wallCircleBound(data3d<MyDataF> &stru) {
    unsigned i, j, k;
    unsigned endx, endy, endz;

    endx = stru.nx - 1;
    endy = stru.ny - 1;
    endz = stru.nz - 1;

    //bottom and top
    unsigned endz1 = endz - 1;
    unsigned endz2 = endz - 2;
    for (i = 1; i < endx; i++) {
        for (j = 1; j < endy; j++) {
            stru.p[i][j][0] = 2 * stru.p[i][j][1] - stru.p[i][j][2];
            stru.p[i][j][endz] = 2 * stru.p[i][j][endz1] - stru.p[i][j][endz2];
        }
    }

    //left and right
    unsigned endx1 = endx - 1;
    unsigned endx2 = endx - 2;
    for (j = 1; j < endy; j++) {
        for (k = 1; k < endz; k++) {
            stru.p[0][j][k] = 2 * stru.p[1][j][k] - stru.p[2][j][k];
            stru.p[endx][j][k] = 2 * stru.p[endx1][j][k] - stru.p[endx2][j][k];
        }
    }

    //front and back
    unsigned endy1 = endy - 1;
    unsigned endy2 = endy - 2;
    for (i = 1; i < endx; i++) {
        for (k = 1; k < endz; k++) {
            stru.p[i][0][k] = 2 * stru.p[i][1][k] - stru.p[i][2][k];
            stru.p[i][endy][k] = 2 * stru.p[i][endy1][k] - stru.p[i][endy2][k];
        }
    }
}

void fdtd::createCoeff() {
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
    Beta.create3DArray(Ne, 0.0);
    // velocity coefficients
    if (mSrcType == fdtd::SOURCE_GAUSSIAN) {
        Cvxex_guassian.create3DArray(Vx, 0.0);
        Cvyey_guassian.create3DArray(Vy, 0.0);
        Cvzez_guassian.create3DArray(Vz, 0.0);
    }
}

void fdtd::initCoeff() {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Velocity Coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    mA = mNiu_m * mDt / 2;
    mGamma = 1 + mA;
    mAlpha = (1 - mA) / mGamma;
    Cvxex = Cvyey = Cvzez = e * mDt / 2 / me / mGamma;
    Coeff_velocity = mHalf_e * mDt / me;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update collision frequency
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (mSrcType == fdtd::SOURCE_GAUSSIAN) {
        Nu_c.resetArray(mNiu_m);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //electricity coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    updateCoeff();

}

void fdtd::updateCoeff() {

    updateBeta();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //electricity coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int i, j, k;
    unsigned im, jm = 0, km;
    MyDataF tmp = eMDtDiv2DivEps0 * (1 + mAlpha);
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k,im,jm,km)//shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif    
    for (j = 0; j < Ex.ny; j++, jm += mNeGridSize) {
        for (i = 0, im = mHalfNeGridSize; i < Ex.nx; i++, im += mNeGridSize) {
            for (k = 0, km = mHalfNeGridSize; k < Ex.nz; k++, km += mNeGridSize) {
                MyDataF kappa = (1 + Beta.p[im][jm][km]);
                Cexe.p[i][j][k] = (1 - Beta.p[im][jm][km]) / kappa;
                Cexhy.p[i][j][k] = -dtDivEps0DivDz / kappa;
                Cexhz.p[i][j][k] = dtDivEps0DivDy / kappa;
                if (mSrcType == fdtd::SOURCE_GAUSSIAN) {
                    MyDataF a = mHalfDelta_t * Nu_c.p[im][jm][km];
                    MyDataF gamma = 1 + a;
                    MyDataF tmpc = eMDtDiv2DivEps0 * (1 + (1 - a) / gamma);
                    Cvxex_guassian.p[i][j][k] = Coeff_velocity / gamma;
                    Cexvx.p[i][j][k] = tmpc * Ne.p[im][jm][km] / kappa;
                } else {
                    Cexvx.p[i][j][k] = tmp * Ne.p[im][jm][km] / kappa;
                }
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k,im,jm,km)//shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif    
    for (i = 0; i < Ey.nx; i++) {
        im = i*mNeGridSize;
        for (j = 0, jm = mHalfNeGridSize; j < Ey.ny; j++, jm += mNeGridSize) {
            for (k = 0, km = mHalfNeGridSize; k < Ey.nz; k++, km += mNeGridSize) {
                MyDataF kappa = (1 + Beta.p[im][jm][km]);
                Ceye.p[i][j][k] = (1 - Beta.p[im][jm][km]) / kappa;
                Ceyhx.p[i][j][k] = dtDivEps0DivDz / kappa;
                Ceyhz.p[i][j][k] = -dtDivEps0DivDx / kappa;
                if (mSrcType == fdtd::SOURCE_GAUSSIAN) {
                    MyDataF a = mHalfDelta_t * Nu_c.p[im][jm][km];
                    MyDataF gamma_t = 1 + a;
                    MyDataF tmpc = eMDtDiv2DivEps0 * (1 + (1 - a) / gamma_t);
                    Cvyey_guassian.p[i][j][k] = Coeff_velocity / gamma_t;
                    Ceyvy.p[i][j][k] = tmpc * Ne.p[im][jm][km] / kappa;
                } else {
                    Ceyvy.p[i][j][k] = tmp * Ne.p[im][jm][km] / kappa;
                }
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k,im,jm,km)//shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif   
    for (i = 0; i < Ez.nx; i++) {
        im = i*mNeGridSize;
        for (j = 0, jm = 0; j < Ez.ny; j++, jm += mNeGridSize) {
            for (k = 0, km = 0; k < Ez.nz; k++, km += mNeGridSize) {
                MyDataF kappa = (1 + Beta.p[im][jm][km]);
                Ceze.p[i][j][k] = (1 - Beta.p[im][jm][km]) / kappa;
                Cezhy.p[i][j][k] = dtDivEps0DivDx / kappa;
                Cezhx.p[i][j][k] = -dtDivEps0DivDy / kappa;
                if (mSrcType == fdtd::SOURCE_GAUSSIAN) {
                    MyDataF a = mHalfDelta_t * Nu_c.p[im][jm][km];
                    MyDataF gamma_t = 1 + a;
                    MyDataF tmpc = eMDtDiv2DivEps0 * (1 + (1 - a) / gamma_t);
                    Cvzez_guassian.p[i][j][k] = Coeff_velocity / gamma_t;
                    Cezvz.p[i][j][k] = tmpc * Ne.p[im][jm][km] / kappa;
                } else {
                    Cezvz.p[i][j][k] = tmp * Ne.p[im][jm][km] / kappa;
                }
            }
        }
    }
}

void fdtd::updateBeta() {
    Point start(mStartIndex.x*mNeGridSize, mStartIndex.y*mNeGridSize, mStartIndex.z * mNeGridSize);
    Point end(mEndIndex.x*mNeGridSize, mEndIndex.y*mNeGridSize, mEndIndex.z * mNeGridSize);
    MyDataF temp = e2Dt2Div4DivEps0DivMe;
    if (mSrcType != fdtd::SOURCE_GAUSSIAN) {
        temp = temp / mGamma;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)  //shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif
        for (int i = start.x; i < end.x; i++) {
            for (unsigned j = start.y; j < end.y; j++) {
                for (unsigned k = start.z; k < end.z; k++) {
                    Beta.p[i][j][k] = temp * Ne.p[i][j][k];
                }
            }
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif
        for (int i = start.x; i < end.x; i++) {
            for (unsigned j = start.y; j < end.y; j++) {
                for (unsigned k = start.z; k < end.z; k++) {
                    Beta.p[i][j][k] = temp / (1 + mHalfDelta_t * Nu_c.p[i][j][k]) * Ne.p[i][j][k];
                }
            }
        }
    }
}

void fdtd::initDensity() {
    MyDataF tmp = 2 * pow(4 * mDx, 2);
    Point srcPos(mSourceIndex.x*mNeGridSize, mSourceIndex.y*mNeGridSize, mSourceIndex.z * mNeGridSize);

    for (int i = 0; i < Ne.nx; i++) {
        for (int j = 0; j < Ne.ny; j++) {
            for (int k = 0; k < Ne.nz; k++) {
#ifdef DEBUG
                MyDataF sx, sy, sz;
                MyDataF px, py, pz;
                MyDataF ea;
                sx = (i - (int) srcPos.x) * mDx;
                sy = (j - (int) srcPos.y) * mDy;
                sz = (k - (int) srcPos.z) * mDz;
                px = sx*sx;
                py = sy*sy;
                pz = sz*sz;
                ea = exp(-(px + py + pz) / tmp);
                Ne.p[i][j][k] = Ne0*ea;
#else
                Ne.p[i][j][k] = Ne0 * exp(-(pow((i - (int) srcPos.x) * mDx, 2) + pow((j - (int) srcPos.y) * mDy, 2) + pow((k - (int) srcPos.z) * mDz, 2)) / tmp);
#endif
            }
        }
    }
}
#endif

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

#ifdef WITH_DENSITY
    Vz.setName("Vz");
    Vx.setName("Vx");
    Vy.setName("Vy");
    Vx.create3DArray(Ex, 0.0);
    Vy.create3DArray(Ey, 0.0);
    Vz.create3DArray(Ez, 0.0);

#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << " neGrid = " << mNeGridSize << endl;
#endif
    Ne.create3DArray((mMaxIndex.x + 1) * mNeGridSize, (mMaxIndex.y + 1) * mNeGridSize, (mMaxIndex.z + 1) * mNeGridSize, Ne0);
    Erms.create3DArray(Ne, 0.0);
    Ne_pre.create3DArray(Ne, 0.0);

    // Gaussian need niu_c from previous step
    if (mSrcType == fdtd::SOURCE_GAUSSIAN) {
        Nu_c.create3DArray(Ne, 0.0);
        Nu_c.setName("nu_c");
    }

    createCoeff();
    Ne.setName("Ne");

#endif
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void fdtd::setUp() {
    //Time step
    //    dt = 0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) +
    //            1.0 / (dz * dz)));
    mDt = mDx / 2 / C;

    //delay
    if (mSrcType == fdtd::SOURCE_GAUSSIAN) {
        t0 = 4.5 * tw;
    } else {
        t0 = tw;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //temporary variables that often used
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    dtDivEps0DivDxyz = mDt / eps_0 / mDx / mDy / mDz;
#ifdef WITH_DENSITY
    mDeltaSizeFluid = mDx / mNeGridSize;
    mHalfDelta_t = mDt / 2;
    mHalf_e = e / 2;
    dtDivEps0DivDx = mDt / eps_0 / mDx;
    dtDivEps0DivDy = mDt / eps_0 / mDy;
    dtDivEps0DivDz = mDt / eps_0 / mDz;
    e2Dt2Div4DivEps0DivMe = 0.25 * e * e * mDt * mDt / me / eps_0;
    eMDtDiv2DivEps0 = mHalf_e * mDt / eps_0;
    mHalfNeGridSize = mNeGridSize / 2;
#endif

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // fluid variables
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef WITH_DENSITY	
    mMu_e = e / me / mNiu_m; //3.7e-2;
    mMu_i = mMu_e / 100.0; //mu_e/mu_i ranges from 100 to 200
    mDe = mMu_e * 2 * 1.602e-19 / e; //
    mDa = mMu_i * 2 * 1.602e-19 / e; //
    MyDataF Dmax = mDe > mDa ? mDe : mDa;
    //Fine Time Step Size
    mDeltaTimeFluid = 0.05 * mDeltaSizeFluid * mDeltaSizeFluid / 2 / Dmax;
    mNeSkipStep = mDeltaTimeFluid / mDt;


    cout << "neSkipStep=" << mNeSkipStep << endl;
    cout << tw / mDt / mNeSkipStep << endl;
    //exit(0);
#endif
    //  Specify the dipole size 
    mStartIndex.setValue(mPMLWidth, mPMLWidth, mPMLWidth);
    mEndIndex.setValue(mMaxIndex.x - mPMLWidth, mMaxIndex.y - mPMLWidth, mMaxIndex.z - mPMLWidth);

    // source position    
#ifdef DEBUG
    mSourceIndex.setValue(mMaxIndex.x / 2, mMaxIndex.y / 2, mMaxIndex.z / 2);
#else
    mSourceIndex.setValue(mMaxIndex.x / 2, mMaxIndex.y - mPMLWidth - 35, mMaxIndex.z / 2);
#endif    

    checkmax(mSourceIndex.x, 1, mMaxIndex.x);
    checkmax(mSourceIndex.y, 1, mMaxIndex.y);
    checkmax(mSourceIndex.z, 1, mMaxIndex.z);

    if (mEndIndex.x < mStartIndex.x)mEndIndex.x = mStartIndex.x + 1;
    if (mEndIndex.y < mStartIndex.y)mEndIndex.y = mStartIndex.y + 1;
    if (mEndIndex.z < mStartIndex.z)mEndIndex.z = mStartIndex.z + 1;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  COMPUTING FIELD UPDATE EQUATION COEFFICIENTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    initCoeficients();

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
    pml.setCPMLRegion(mPMLWidth);
    pml.createCPMLArrays(mMaxIndex.x, mMaxIndex.y, mMaxIndex.z);
    pml.initCoefficientArrays(pmlOrder, alphaOrder, sigmaMax, kappaMax, alphaMax, epsR, mDt, mDx, mDy, mDz,
            Ceyhz, Cezhy, Chyez, Chzey,
            Cexhz, Cezhx, Chxez, Chzex,
            Ceyhx, Cexhy, Chyex, Chxey);

#ifdef WITH_DENSITY
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // initial density
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    initDensity();
#ifdef DEBUG
    Ne.save(Ne.ny / 2, 0, 0, 3);
    Ne.save(Ne.ny / 2, 0, 0, 1);
    Ne.save(Ne.ny / 2, 0, 0, 2);
    Ne.save();
#endif
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initial Coefficients for Density
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    initCoeff();
#endif
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // print parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    printParameters();
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SET CPML PARAMETERS IN EACH DIRECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void fdtd::compute() {

    unsigned n;
    Point capturePosition(mSourceIndex.x, mSourceIndex.y + 30, mSourceIndex.z);

#if DEBUG>=3
    if (mSourceIndex.y + 30 < mMaxIndex.y) {
        capturePosition.y = mSourceIndex.y + 30;
    }
#endif

    if (!capturePosition.checkMax(mMaxIndex.x, mMaxIndex.y, mMaxIndex.z)) {
        if (capturePosition.y >= mMaxIndex.y) {
            capturePosition.y = mMaxIndex.y - 1;
        }
        if (capturePosition.x >= mMaxIndex.x) {
            capturePosition.x = mMaxIndex.x - 1;
        }
        if (capturePosition.z >= mMaxIndex.z) {
            capturePosition.z = mMaxIndex.z - 1;
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
        cout << endl;
        cout << "Source Value:";
        cout << Ez.p[mSourceIndex.x][mSourceIndex.y][mSourceIndex.z] << '\t';
        cout << endl;
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+10][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+15][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+20][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+25][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+30][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+35][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+40][mSourceIndex.z] << '\t';
        //cout	<< Ez.p[mSourceIndex.x][mSourceIndex.y+45][mSourceIndex.z] << endl;

        updateMagneitcFields();
        pml.updateCPML_M_Fields(Hx, Hy, Hz, Ex, Ey, Ez);
        updateElectricAndVeloityFields();
        pml.updateCPML_E_Fields(Ex, Ey, Ez, Hx, Hy, Hz);
        //====================================
        // update Source
        //====================================
        updateSource(n);

#ifdef WITH_DENSITY
        updateErms();
        if (n % mNeSkipStep == 0) {
            interpErms();
            if (mSrcType == fdtd::SOURCE_GAUSSIAN)
                updateCollisionFrequency();
            updateDensity();
            updateCoeff();
        }
#endif
        if ((n % save_modulus) == 0) {
            writeField(n);
            //Ez.save(mSourceIndex.x + 10, 1, n, 1);
            //Ez.save(mSourceIndex.y + 10, 1, n, 2);
            Ez.save(mMaxIndex.z - mPMLWidth - 5, 1, n, 3);
            Ex.save(mMaxIndex.z - mPMLWidth - 5, 1, n, 3);
            Ey.save(mMaxIndex.z - mPMLWidth - 5, 1, n, 3);
            /*
            pml.Psi_exz_zp.setName("psi");
            pml.Psi_exz_zp.save(0, 1, n, 3);
            pml.Psi_exz_zp.save(4, 1, n, 3);
             */
#ifdef WITH_DENSITY
            Ne.save(Ne.nz / 2, 1, n, 2);
#endif
        }
#ifdef MATLAB_SIMULATION
        doMatlabSimulation();
#endif

    }
#ifdef MATLAB_SIMULATION
    finishMatlabSimulation();
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
            cout << "Gaussian source:" << source << endl;
            break;
        case SOURCE_SINE:
            // sine wave
            source = M_PI_TWO * mOmega * mAmplitude * cos((n * mDt - t0) * M_PI_TWO * mOmega);
            break;
        case ONE_SINE_PULSE:
            source = M_PI_TWO * mOmega * mAmplitude * Source::SinePulse(n * mDt - t0, mOmega, t_up, t_down);
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

    //buildSphere();
    //buildDipole();
}

//Builds a sphere (Sample code - NOt used in this program)

void fdtd::buildSphere() {

    MyDataF dist; //distance
    MyDataF rad = 8; //(MyDataF)mMaxIndex.x / 5.0; // sphere radius
    MyDataF sc = (MyDataF) mMaxIndex.x / 2.0; //sphere centre
    //MyDataF rad2 = 0.3; //(MyDataF)mMaxIndex.x / 5.0 - 3.0; // sphere radius

    unsigned i, j, k;

    for (i = 0; i < mMaxIndex.x; ++i) {
        for (j = 0; j < mMaxIndex.y; ++j) {
            for (k = 0; k < mMaxIndex.z; ++k) {
                //compute distance form centre to the point i, j, k
                dist = sqrt((i + 0.5 - sc) * (i + 0.5 - sc) +
                        (j + 0.5 - sc) * (j + 0.5 - sc) +
                        (k + 0.5 - sc) * (k + 0.5 - sc));

                //if point is within the sphere
                if (dist <= rad) {
                    //set the material at that point
                    yeeCube(i, j, k, 6);

                }
            }
        }
    }

}

//Builds a dipole

void fdtd::buildDipole() {
    unsigned i, j, k;
    unsigned centre = (mStartIndex.y + mEndIndex.y) / 2;

    for (i = mStartIndex.x; i <= mEndIndex.x; ++i) {

        for (j = mStartIndex.y; j <= mEndIndex.y; ++j) {

            for (k = mStartIndex.z; k <= mEndIndex.z; ++k) {

                if (j != centre) {

                    yeeCube(i, j, k, 1); //PEC material
                }
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
    //    cout << "initial pml (in Statup)" << endl;
    //    pml.Initial(mMaxIndex.x, mMaxIndex.y, mMaxIndex.z, 11);
    cout << "setUp (in Startup)" << endl;
    setUp();
    cout << "buildObject (in Startup)" << endl;
    buildObject();
    cout << "computing (in Startup)" << endl;
    compute();
    cout << "exit Startup" << endl;
}

void fdtd::printParameters() {
    cout << "dx = " << mDx << endl;
    cout << "dy = " << mDy << endl;
    cout << "dz = " << mDz << endl;
    cout << "(mMaxIndex.x,mMaxIndex.y,mMaxIndex.z) = (" << mMaxIndex.x << "," << mMaxIndex.y << "," << mMaxIndex.z << ")" << endl;
    // time step increment
    cout << "dt = " << mDt << endl;

    //  Specify the Impulsive Source (Differentiated Gaussian) parameters
    cout << "tw = " << tw << endl; //pulse width
    cout << "t0 = " << t0 << endl; //delay
    cout << "source = " << source << endl; //Differentiated Gaussian source
    cout << "amp = " << mAmplitude << endl; // Amplitude
    cout << "omega = " << mOmega << endl; // angle speed for sine wave

    //Specify the Time Step at which the data has to be saved for Visualization
    cout << "save_modulus = " << save_modulus << endl;

    //  Specify the dipole Boundaries(A cuboidal rode- NOT as a cylinder)
    cout << "(mStartIndex.x,mStartIndex.y, mStartIndex.z) = (" << mStartIndex.x << ',' << mStartIndex.y << ',' << mStartIndex.z << ')' << endl;
    cout << "(mEndIndex.x,  mEndIndex.y, mEndIndex.z) = (" << mEndIndex.x << ',' << mEndIndex.y << ',' << mEndIndex.z << ')' << endl;

    //Output recording point
    cout << "ksource = " << mKSource << endl;

    //  Specify the CPML Order and Other Parameters:
    cout << " m = " << mPMLOrder << endl;
    cout << " ma = " << mAlphaOrder << endl;
#ifdef WITH_DENSITY 
    cout << "neGrid=" << mNeGridSize << endl;
    cout << "neSkipStep=" << mNeSkipStep << endl;
    cout << "dtf=" << mDeltaTimeFluid << endl;
    cout << "dsf=" << mDeltaSizeFluid << endl;
    cout << "mu_i=" << mMu_i << endl;
    cout << "mu_e=" << mMu_e << endl;
    cout << "Da=" << mDa << endl;
    cout << "De=" << mDe << endl;
#endif
    cout << endl << "Time step = " << mDt << endl;
    cout << endl << "Number of steps = " << mTotalTimeSteps << endl;
    cout << endl << "Total Simulation time = " << mTotalTimeSteps * mDt << " Seconds" << endl;
}

void fdtd::setSourceType(int sourceType) {
    mSrcType = sourceType;
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
    for (k = 0; k < mMaxIndex.z; ++k) {
        for (j = 0; j < mMaxIndex.y; ++j) {
            for (i = 0; i <= mMaxIndex.x; ++i) {
                Hx.p[i][j][k] = Chxh.p[i][j][k] * Hx.p[i][j][k] +
                        Chxez.p[i][j][k]*(Ez.p[i][j + 1][k] - Ez.p[i][j][k]) +
                        Chxey.p[i][j][k]*(Ey.p[i][j][k + 1] - Ey.p[i][j][k]);
#ifdef WITH_DENSITY
#if (DEBUG>=4&&!_OPENMP)
                Hx.isValid(i, j, k);
#endif
#endif
            }
        }
    }
}

void fdtd::updateHy() {
    int i, j, k;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) //shared(Hy,Ez,Ex,pml,DA,DB,dx,dz)
#endif
    for (k = 0; k < mMaxIndex.z; ++k) {
        for (i = 0; i < mMaxIndex.x; ++i) {
            for (j = 0; j <= mMaxIndex.y; ++j) {
                Hy.p[i][j][k] = Chyh.p[i][j][k] * Hy.p[i][j][k] +
                        Chyez.p[i][j][k]*(Ez.p[i + 1][j][k] - Ez.p[i][j][k]) +
                        Chyex.p[i][j][k]*(Ex.p[i][j][k + 1] - Ex.p[i][j][k]);
#ifdef WITH_DENSITY
#if (DEBUG>=4&&!_OPENMP)
                Hy.isValid(i, j, k);
#endif
#endif
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
    for (k = 0; k <= mMaxIndex.z; ++k) {
        for (i = 0; i < mMaxIndex.x; ++i) {
            for (j = 0; j < mMaxIndex.y; ++j) {
                Hz.p[i][j][k] = Chzh.p[i][j][k] * Hz.p[i][j][k] +
                        (Ey.p[i + 1][j][k] - Ey.p[i][j][k]) * Chzey.p[i][j][k] +
                        (Ex.p[i][j + 1][k] - Ex.p[i][j][k]) * Chzex.p[i][j][k];
#ifdef WITH_DENSITY
#if (DEBUG>=4&&!_OPENMP)
                Hz.isValid(i, j, k);
#endif
#endif
            }
        }
    }
}

void fdtd::updateEx() {
    int i, j, k;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  UPDATE Ex
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k)//shared(Ex,Hz,Hy,pml,CA,CB,ID1,dy,dz)
#endif
    for (k = 1; k < mMaxIndex.z; ++k) {
        for (j = 1; j < mMaxIndex.y; ++j) {
            for (i = 0; i < mMaxIndex.x; ++i) {
#ifdef WITH_DENSITY
                MyDataF Exp = Ex.p[i][j][k];
#endif
                Ex.p[i][j][k] = Cexe.p[i][j][k] * Ex.p[i][j][k] +
                        (Hz.p[i][j][k] - Hz.p[i][j - 1][k]) * Cexhz.p[i][j][k] +
                        (Hy.p[i][j][k] - Hy.p[i][j][k - 1]) * Cexhy.p[i][j][k];
#ifdef WITH_DENSITY
                Ex.p[i][j][k] += Cexvx.p[i][j][k] * Vx.p[i][j][k];

                if (mSrcType == fdtd::SOURCE_GAUSSIAN) {
                    MyDataF a = Nu_c.p[i * mNeGridSize][j * mNeGridSize + mHalfNeGridSize][k * mNeGridSize];
                    Vx.p[i][j][k] = (1 - a) / (1 + a) * Vx.p[i][j][k] - Cvxex_guassian.p[i][j][k] * (Exp + Ex.p[i][j][k]);
                } else {
                    Vx.p[i][j][k] = mAlpha * Vx.p[i][j][k] - Cvxex * (Exp + Ex.p[i][j][k]);
                }
#endif
            }
        }
    }
}

void fdtd::updateEy() {
    int i, j, k;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  UPDATE Ey
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) //shared(Ex,Hz,Hy,pml,CA,CB,ID1,dy,dz)
#endif
    for (k = 1; k < mMaxIndex.z; ++k) {
        for (i = 1; i < mMaxIndex.x; ++i) {
            for (j = 0; j < mMaxIndex.y; ++j) {
#ifdef WITH_DENSITY
                MyDataF Eyp = Ey.p[i][j][k];
#endif  /* WITH_DENSITY */

                Ey.p[i][j][k] = Ceye.p[i][j][k] * Ey.p[i][j][k] +
                        (Hz.p[i][j][k] - Hz.p[i - 1][j][k]) * Ceyhz.p[i][j][k] +
                        (Hx.p[i][j][k] - Hx.p[i][j][k - 1]) * Ceyhx.p[i][j][k];
#ifdef WITH_DENSITY
                Ey.p[i][j][k] += Ceyvy.p[i][j][k] * Vy.p[i][j][k];

#if (DEBUG>=4&&!_OPENMP)
                Ey.isValid(i, j, k);
#endif
                // Vy.p[i][j][k] = alpha * Vy.p[i][j][k] - Cvyey * (Eyp + Ey.p[i][j][k]);
                if (mSrcType == fdtd::SOURCE_GAUSSIAN) {
                    MyDataF a = Nu_c.p[i * mNeGridSize + mHalfNeGridSize][j * mNeGridSize][k * mNeGridSize];
                    Vy.p[i][j][k] = (1 - a) / (1 + a) * Vy.p[i][j][k] - Cvyey_guassian.p[i][j][k] * (Eyp + Ey.p[i][j][k]);
                } else {
                    Vy.p[i][j][k] = mAlpha * Vy.p[i][j][k] - Cvyey * (Eyp + Ey.p[i][j][k]);
                }
#endif /* WITH_DENSITY */
            }
        }
    }
}

void fdtd::updateEz() {
    int i, j, k;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  UPDATE Ez
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k) //shared(Ex,Hz,Hy,pml,CA,CB,ID1,dy,dz)
#endif

    for (i = 1; i < mMaxIndex.x; ++i) {
        for (k = 0; k < mMaxIndex.z; ++k) {
            for (j = 1; j < mMaxIndex.y; ++j) {
#ifdef WITH_DENSITY
                MyDataF Ezp = Ez.p[i][j][k];
#endif
                Ez.p[i][j][k] = Ceze.p[i][j][k] * Ez.p[i][j][k] +
                        (Hy.p[i][j][k] - Hy.p[i - 1][j][k]) * Cezhy.p[i][j][k]+
                        (Hx.p[i][j][k] - Hx.p[i][j - 1][k]) * Cezhx.p[i][j][k];
#ifdef WITH_DENSITY
                Ez.p[i][j][k] += Cezvz.p[i][j][k] * Vz.p[i][j][k];

                // Vz.p[i][j][k] = alpha * Vz.p[i][j][k] - Cvzez * (Ezp + Ez.p[i][j][k]);
                if (mSrcType == fdtd::SOURCE_GAUSSIAN) {
                    MyDataF a = Nu_c.p[i * mNeGridSize ][j * mNeGridSize][k * mNeGridSize];
                    Vz.p[i][j][k] = (1 - a) / (1 + a) * Vz.p[i][j][k] - Cvzez_guassian.p[i][j][k] * (Ezp + Ez.p[i][j][k]);
                } else {
                    Vz.p[i][j][k] = mAlpha * Vz.p[i][j][k] - Cvzez * (Ezp + Ez.p[i][j][k]);
                }
#endif
#if DEBUG>=4
                if (isnan(Ez.p[i][j][k])) {
                    cout << "Ez is nan at " << i << "," << j << "," << k << endl;
                }
#endif
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

void fdtd::updateElectricAndVeloityFields() {
    updateEx();
    updateEy();
    updateEz();
}

void fdtd::intSourceSinePulse(MyDataF t_0, MyDataF omega_, MyDataF tUp, MyDataF tDown, MyDataF amptidute) {
    t0 = t_0;
    mOmega = omega_;
    t_up = tUp;
    t_down = tDown;
    mAmplitude = amptidute;
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
