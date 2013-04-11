

#include <math.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <fstream>
#include <strstream>

#include "fdtd.h"
#include "InonizationFormula.h"

extern MyDataF eps_0, epsR;
extern MyDataF mu_0;
extern MyDataF dt, dx, dy, dz;
extern MyDataF pi, C, me, e, T;

using namespace std;
#ifdef WITH_DENSITY

fdtd::fdtd(unsigned _nmax, unsigned _imax, unsigned _jmax, unsigned _kmax,
        MyDataF _tw, MyDataF _dx, MyDataF _dy, MyDataF _dz,
        unsigned _amp, unsigned _savemodulus, unsigned _ksource,
        unsigned _m, unsigned _ma, unsigned pmlw, unsigned _nmaterial, unsigned _neGrid)
: nMax(_nmax), Imax(_imax), Jmax(_jmax), Kmax(_kmax)
, tw(_tw), dx(_dx), dy(_dy), dz(_dz)
, amp(_amp), save_modulus(_savemodulus), ksource(_ksource)
, m(_m), ma(_ma), pmlWith(pmlw)
, neGrid(_neGrid)
, numMaterials(_nmaterial)
, Ne0(1e7)
, epsilon(NULL), sigma(NULL), mu(NULL), CA(NULL), CB(NULL) {
}
#else

fdtd::fdtd(unsigned _nmax, unsigned _imax, unsigned _jmax, unsigned _kmax,
        MyDataF _tw, MyDataF _dx, MyDataF _dy, MyDataF _dz,
        unsigned _amp, unsigned _savemodulus, unsigned _ksource,
        unsigned _m, unsigned _ma, unsigned pmlw, unsigned _nmaterial)
: nMax(_nmax), Imax(_imax), Jmax(_jmax), Kmax(_kmax)
, tw(_tw), dx(_dx), dy(_dy), dz(_dz)
, amp(_amp), save_modulus(_savemodulus), ksource(_ksource)
, m(_m), ma(_ma), pmlWith(pmlw)
, numMaterials(_nmaterial), epsilon(NULL), sigma(NULL), mu(NULL), CA(NULL), CB(NULL) {
}
#endif

fdtd::~fdtd(void) {
    if (epsilon != NULL)delete []epsilon;
    if (sigma != NULL)delete []sigma;
    if (mu != NULL)delete[]mu;
    if (CA != NULL)delete[]CA;
    if (CB != NULL)delete[]CB;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef WITH_DENSITY

void fdtd::SetPlasmaVar(MyDataF _rei, MyDataF _vm, MyDataF _p, int _ftype) {
    rei = _rei;
    vm = _vm;
    p = _p;
    niutype = _ftype;
}

int fdtd::UpdateErms(void) {
    unsigned i, j, k;
    unsigned io, jo, ko;
    for (i = istart, io = istart * neGrid; i <= iend; i++, io += neGrid) {
        for (j = jstart, jo = jstart * neGrid; j <= jend; j++, jo += neGrid) {
            for (k = kstart, ko = kstart * neGrid; k <= kend; k++, ko += neGrid) {
                MyDataF exIJK = (Ex.p[i - 1][j][k - 1] + Ex.p[i + 1][j][k - 1] + Ex.p[i - 1][j][k + 1] + Ex.p[i + 1][j][k + 1]) / 4;
                MyDataF eyIJK = (Ey.p[i][j - 1][k - 1] + Ey.p[i][j + 1][k - 1] + Ey.p[i][j - 1][k + 1] + Ey.p[i][j + 1][k + 1]) / 4;
                Erms.p[io][jo][ko] = sqrt(Ez.p[i][j][k] * Ez.p[i][j][k] + exIJK * exIJK + eyIJK * eyIJK);
            }
        }
    }
    return 0;
}

int fdtd::InterpErms() {
    unsigned is, js, ks;
    unsigned in, jn, kn;
    unsigned i, j, k;
    unsigned im, jm, km;
    unsigned iu, ju, ku;
    unsigned ngred = neGrid * neGrid*neGrid;
    for (is = istart, in = istart + neGrid; is < iend; is = in, in += neGrid) {
        for (js = jstart, jn = jstart + neGrid; js < jend; js = jn, jn += neGrid) {
            for (ks = kstart, kn = kstart + neGrid; ks < kend; ks = kn, kn += neGrid) {
                // integrate Erms
                for (i = is, im = 0, iu = neGrid; i < in; i++, im++, iu--) {
                    for (j = js, jm = 0, ju = neGrid; j < jn; j++, jm++, ju--) {
                        for (k = ks, km = 0, ku = neGrid; k < kn; k++, km++, ku--) {
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
    return 0;
}

int fdtd::UpdateDensity(void) {

    unsigned i, j, k, mt = 1;

    MyDataF Eeff, alpha_t, tau_m, kasi;
    MyDataF Ne_ijk, Neip1, Neim1, Nejm1, Nejp1, Nekp1, Nekm1;
    MyDataF Deff;
    MyDataF maxvi = 0, minvi = 0;
    MyDataF vi, va;

    unsigned ci = 0, cj = 0, ck = 0;
    Ne_pre = Ne;

    for (i = mt; i < Ne.nx - mt; i++) {
        for (j = mt; j < Ne.ny - mt; j++) {
            for (k = mt; k < Ne.nz - mt; k++) {

                Eeff = Erms.p[i][j][k] / 100 * pow(1 / (1 + omega * omega / vm / vm), 0.5);

                Ne_ijk = Ne_pre.p[i][j][k];
                Neip1 = Ne_pre.p[i + 1][j][k];
                Neim1 = Ne_pre.p[i - 1][j][k];
                Nejp1 = Ne_pre.p[i][j + 1][k];
                Nejm1 = Ne_pre.p[i][j - 1][k];
                Nekp1 = Ne_pre.p[i][j][k + 1];
                Nekm1 = Ne_pre.p[i][j][k - 1];

                switch (niutype) {
                    case 1:
                        Niu_MorrowAndLowke(&vi, &va, Eeff, Ne_ijk * 1e6);
                        break;
                    case 2:
                        Niu_Nikonov(&vi, &va, Eeff, p);
                        break;
                    case 3:
                        Niu_Kang(&vi, &va, Eeff);
                        break;
                    default:
                        alpha_t = Eeff / p;
                        if (alpha_t < 30) {
                            if (alpha_t < 1e-12) {
                                vi = 0;
                            } else if (alpha_t >= 1) {
                                vi = (1.45 + 0.01 * pow(alpha_t, 1.5))*2.5e7 * exp(-208 / alpha_t) * p;
                            } else {
                                vi = 5.14e11 * exp(-73 * pow(alpha_t, -0.44)) * p;
                            }
                        } else if (alpha_t > 120) {
                            if (alpha_t <= 3000) {
                                vi = 54.08e6 * pow(alpha_t, 0.5) * exp(-359 / alpha_t) * p;
                            } else {
                                vi = 5.14e11 * exp(-73 * pow(alpha_t, -0.44)) * p;
                            }
                        } else if (alpha_t > 54) {
                            vi = (1.32 + 0.054 * alpha_t)*1e7 * exp(-208 / alpha_t) * p;
                        } else {
                            vi = (5.0 + 0.19 * alpha_t)*1e7 * exp(-273.8 / alpha_t) * p;
                        }
                        va = 7.6e-4 * pow(alpha_t / (alpha_t + 218), 2) / p;
                }
                if (Ne_ijk < 1) {
                    Deff = De;
                } else {
                    tau_m = eps_0 / (e * Ne_ijk * (mu_e + mu_i));
                    kasi = vi * tau_m;
                    Deff = (kasi * De + Da) / (kasi + 1);
                }
                Ne.p[i][j][k] =
                        (
                        Ne_ijk * (1 + dtf * vi)
                        + Deff * dtf * (Neip1 + Neim1 + Nejp1 + Nejm1 + Nekp1 + Nekm1 - 6 * Ne_ijk) / dtf / dtf
                        ) / (1 + dtf * (va + rei * Ne_ijk));
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
    WallCircleBound(Ne);
    cout << Ne.p[Ne.nx / 2][Ne.ny / 2][Ne.nz / 2] << '\t';
    cout << maxvi << '\t' << minvi << '\t' << Ne.p[ci][cj][ck] << '\t' << Erms.p[ci][cj][ck] << '\t';
    return 0;
}

int fdtd::UpdateVeloity(void) {
    return 0;
}

void fdtd::WallCircleBound(data3d<MyDataF> &stru) {
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
    //    Cvxex.CreateStruct(Vx,0.0);
    //    Cvyey.CreateStruct(Vy,0.0);
    //    Cvzez.CreateStruct(Vz,0.0);
    // electricity coefficients
    Cexex.CreateStruct(Ex, 0.0);
    Ceyey.CreateStruct(Ey, 0.0);
    Cezez.CreateStruct(Ez, 0.0);
    Cexh.CreateStruct(Ex, 0.0);
    Ceyh.CreateStruct(Ey, 0.0);
    Cezh.CreateStruct(Ez, 0.0);
    Cexvx.CreateStruct(Ex, 0.0);
    Ceyvy.CreateStruct(Ey, 0.0);
    Cezvz.CreateStruct(Ez, 0.0);
    // beta
    beta.CreateStruct(Ne, 0.0);
}

void fdtd::initCoeff() {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Velocity Coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    a = vm * dt / 2;
    gamma = 1 + a;
    alpha = (1 - a) / gamma;
    Cvxex = Cvyey = Cvzez = e * dt / 2 / me / gamma;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update beta
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //updateBeta();

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
    unsigned i, j, k;
    unsigned im, jm, km;
    unsigned halfNeGrid = neGrid / 2;
    MyDataF tmp = 0.5 * e * (1 + alpha);
    for (i = 0, im = halfNeGrid; i < Ex.nx; i++, im += neGrid) {
        for (j = 0, jm = 0; j < Ex.ny; j++, jm += neGrid) {
            for (k = 0, km = halfNeGrid; k < Ex.nz; k++, km += neGrid) {
                MyDataF kappa = (1 + beta.p[im][jm][km]);
                Cexex.p[i][j][k] = (1 - beta.p[im][jm][km]) / kappa;
                Cexh.p[i][j][k] = 1 / kappa;
                Cexvx.p[i][j][k] = tmp * Ne.p[im][jm][km] / kappa;
            }
        }
    }
    for (i = 0, im = 0; i < Ey.nx; i++, im += neGrid) {
        for (j = 0, jm = halfNeGrid; j < Ey.ny; j++, jm += neGrid) {
            for (k = 0, km = halfNeGrid; k < Ey.nz; k++, km += neGrid) {
                MyDataF kappa = (1 + beta.p[im][jm][km]);
                Ceyey.p[i][j][k] = (1 - beta.p[im][jm][km]) / kappa;
                Ceyh.p[i][j][k] = 1 / kappa;
                Ceyvy.p[i][j][k] = tmp * Ne.p[im][jm][km] / kappa;
            }
        }
    }
    for (i = 0, im = 0; i < Ez.nx; i++, im += neGrid) {
        for (j = 0, jm = 0; j < Ez.ny; j++, jm += neGrid) {
            for (k = 0, km = 0; k < Ez.nz; k++, km += neGrid) {
                MyDataF kappa = (1 + beta.p[im][jm][km]);
                Cezez.p[i][j][k] = (1 - beta.p[im][jm][km]) / kappa;
                Cezh.p[i][j][k] = 1 / kappa;
                Cezvz.p[i][j][k] = tmp * Ne.p[im][jm][km] / kappa;
            }
        }
    }
}

void fdtd::updateBeta() {
    unsigned is = istart*neGrid;
    unsigned ie = iend*neGrid;
    unsigned js = jstart*neGrid;
    unsigned je = jend*neGrid;
    unsigned ks = kstart*neGrid;
    unsigned ke = kend*neGrid;
    MyDataF temp = 0.25 * e * e * dt * dt / me / eps_0 / gamma;
    for (unsigned i = is; i < ie; i++) {
        for (unsigned j = js; j < je; j++) {
            for (unsigned k = ks; k < ke; k++) {
                beta.p[i][j][k] = temp * Ne.p[i][j][k];
            }
        }
    }

}

void fdtd::initDensity() {

}
#endif

void fdtd::initialize() {

    unsigned i;

    // initial PML
    pml.InitialMuEps();
    pml.Initial(Imax, Jmax, Kmax, pmlWith);
#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << "numMaterials = " << numMaterials << endl;
#endif
    //Dynamic memory allocation
    epsilon = new double[numMaterials];

    for (i = 0; i < numMaterials; i++) {

        epsilon[i] = eps_0;
    }

    mu = new double[numMaterials];

    for (i = 0; i < numMaterials; i++) {

        mu[i] = mu_0;
    }

    sigma = new double[numMaterials];

    for (i = 0; i < numMaterials; i++) {

        sigma[i] = 0.0;
    }

    CA = new double[numMaterials];

    for (i = 0; i < numMaterials; i++) {

        CA[i] = 0.0;
    }

    CB = new double[numMaterials];

    for (i = 0; i < numMaterials; i++) {

        CB[i] = 0.0;
    }
#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << "Imax=" << Imax << endl;
    cout << "Jmax=" << Jmax << endl;
    cout << "Kmax=" << Kmax << endl;
#endif
    Ez.CreateStruct(Imax, Jmax, Kmax, 0);
    Ey.CreateStruct(Imax, Jmax - 1, Kmax - 1, 0);
    Ex.CreateStruct(Imax - 1, Jmax, Kmax - 1, 0);

    Hx.CreateStruct(Imax, Jmax - 1, Kmax, 0);
    Hy.CreateStruct(Imax - 1, Jmax, Kmax, 0);
    Hz.CreateStruct((Imax - 1), (Jmax - 1), (Kmax - 1), 0);

#ifdef WITH_DENSITY
    Vx.CreateStruct(Ex, 0.0);
    Vy.CreateStruct(Ey, 0.0);
    Vz.CreateStruct(Ez, 0.0);
#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << " neGrid = " << neGrid << endl;
#endif
    Ne.CreateStruct(Imax*neGrid, Jmax*neGrid, Kmax*neGrid, Ne0);
    Erms.CreateStruct(Ne, 0.0);
    Ne_pre.CreateStruct(Ne, 0.0);
    createCoeff();
#endif

#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << "creating ID1..." << endl;
#endif
    ID1.CreateStruct(Imax, Jmax, Kmax, 0);

#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << "creating ID2..." << endl;
#endif
    ID2.CreateStruct(Imax, Jmax, Kmax, 0);

#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << "creating ID3..." << endl;
#endif
    ID3.CreateStruct(Imax, Jmax, Kmax, 0);
    pml.createCPMLArray();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void fdtd::setUp() {
    unsigned i;
    //Time step
    dt = 0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) +
            1.0 / (dz * dz)));
    //delay
    t0 = 4.0 * tw;
#ifdef WITH_DENSITY
    //Fine Grid size
    dsf = dx / neGrid;

    mu_e = e / me / vm; //3.7e-2;
    mu_i = mu_e / 100.0; //mu_e/mu_i ranges from 100 to 200
    De = mu_e * 2 * 1.602e-19 / e; //
    Da = De * mu_i / mu_e;
    MyDataF Dmax = (De > Da ? De : Da);
    //Fine Time Step Size
    dtf = 0.6 * dsf * dsf / 2 / Dmax;
    neTotalStep = dtf / dt;
#endif
    //  Specify the dipole size 
    istart = pmlWith;
    iend = Imax - pmlWith;
    jstart = pmlWith;
    jend = Jmax - pmlWith;
    kstart = pmlWith;
    kend = Kmax - pmlWith;

    // source position
    isp = Imax / 2;
    jsp = Jmax / 2;
    ksp = Kmax / 2;

    if (iend < istart)iend = istart + 1;
    if (jend < jstart)jend = jstart + 1;
    if (kend < kstart)kend = kstart + 1;
    //Material properties
    //Location '0' is for free space and '1' is for PEC
    epsilon[2] = 4.0 * eps_0;
    sigma[2] = 0.005;
    epsilon[3] = 8.0 * eps_0;
    sigma[3] = 3.96E7; // aluminum
    epsilon[4] = 9.5 * eps_0;
    sigma[4] = 5.76E7; //copper
    epsilon[5] = 9.0 * eps_0;
    sigma[5] = 2e6; //steel
    epsilon[6] = 2.1 * eps_0;
    sigma[6] = 7.8e-4; //teflon
    epsilon[7] = 81 * eps_0;
    sigma[7] = 1e-2; //water

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  COMPUTING FIELD UPDATE EQUATION COEFFICIENTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DA = 1.0;
    DB = dt / mu_0;

    for (i = 0; i < numMaterials; ++i) {

        CA[i] = (1.0 - sigma[i] * dt / (2.0 * epsilon[i])) /
                (1.0 + sigma[i] * dt / (2.0 * epsilon[i]));
        CB[i] = (dt / (epsilon[i])) /
                (1.0 + sigma[i] * dt / (2.0 * epsilon[i]));

    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  PML parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pml.initParmeters(dx, dy, dz, m, ma);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initial Coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    initCoeff();

    cout << endl << "TIme step = " << dt << endl;
    cout << endl << "Number of steps = " << nMax << endl;
    cout << endl << "Total Simulation time = " << nMax * dt << " Seconds" << endl;

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SET CPML PARAMETERS IN EACH DIRECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void fdtd::compute() {

    unsigned id, n;
    unsigned i, j, k;
    unsigned ic, jc, kc; //capture field
    ic = isp + 1;
    jc = jsp + 1;
    kc = ksp + 2;

    assert(ic < Imax && jc < Jmax && kc < Kmax);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  BEGIN TIME STEP
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cout << endl;
    cout << "Begin time-stepping..." << endl;
#ifdef MATLAB_SIMULATION
    initMatlabSimulation();
#endif
    for (n = 1; n <= nMax; ++n) {

        cout << "Ez at time step " << n << " at (" << ic << ", " << jc << ", " << kc;
        cout << ") :  " << Ez.p[ic][jc][kc] << endl;

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Hx
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 1; k < Kmax - 1; ++k) {

            for (i = 0; i < Imax - 1; ++i) {

                for (j = 0; j < Jmax - 1; ++j) {

                    Hx.p[i][j][k] = DA * Hx.p[i][j][k] + DB *
                            ((Ez.p[i][j][k] - Ez.p[i][j + 1][k]) * pml.den_hy.p[j] +
                            (Ey.p[i][j][k] - Ey.p[i][j][k - 1]) * pml.den_hz.p[k]);
                }
            }
            pml.updateHxIn(k, Hx, Ez, DB, dy);
        }
        pml.updateHxOut(Hx, Ey, DB, dz);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Hy
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 1; k < Kmax - 1; ++k) {

            for (i = 0; i < Imax - 1; ++i) {

                for (j = 0; j < Jmax - 1; ++j) {

                    Hy.p[i][j][k] = DA * Hy.p[i][j][k] + DB *
                            ((Ez.p[i + 1][j][k] - Ez.p[i][j][k]) * pml.den_hx.p[i] +
                            (Ex.p[i][j][k - 1] - Ex.p[i][j][k]) * pml.den_hz.p[k]);
                }
            }
            pml.updateHyIn(k, Hy, Ez, DB, dx);
        }
        pml.updateHyOut(Hy, Ex, DB, dz);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Hz
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 0; k < Kmax - 1; ++k) {

            for (i = 0; i < Imax - 1; ++i) {

                for (j = 0; j < Jmax - 1; ++j) {

                    Hz.p[i][j][k] = DA * Hz.p[i][j][k] + DB
                            * ((Ey.p[i][j][k] - Ey.p[i + 1][j][k]) * pml.den_hx.p[i] +
                            (Ex.p[i][j + 1][k] - Ex.p[i][j][k]) * pml.den_hy.p[j]);
                }
            }
            pml.updateHz(k, Hz, Ex, Ey, DB, dx, dy);
        }

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Ex
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 0; k < Kmax - 1; ++k) {

            for (i = 0; i < Imax - 1; ++i) {

                for (j = 1; j < Jmax - 1; ++j) {
#ifdef WITH_DENSITY
                    MyDataF Exp = Ex.p[i][j][k];
#endif
                    id = ID1.p[i][j][k];
                    if (id == 1) { // PEC

                        Ex.p[i][j][k] = 0;

                    } else {
#ifdef WITH_DENSITY
                        Ex.p[i][j][k] = CA[id] * Ex.p[i][j][k] * Cexex.p[i][j][k] + CB[id] * Cexh.p[i][j][k]*
                                ((Hz.p[i][j][k] - Hz.p[i][j - 1][k]) * pml.den_ey.p[j] +
                                (Hy.p[i][j][k] - Hy.p[i][j][k + 1]) * pml.den_ez.p[k]) +
                                Cexvx.p[i][j][k] * Vx.p[i][j][k];
#else
                        Ex.p[i][j][k] = CA[id] * Ex.p[i][j][k] + CB[id] *
                                ((Hz.p[i][j][k] - Hz.p[i][j - 1][k]) * pml.den_ey.p[j] +
                                (Hy.p[i][j][k] - Hy.p[i][j][k + 1]) * pml.den_ez.p[k]);
#endif                                  
                    }
#ifdef WITH_DENSITY
                    Vx.p[i][j][k] = alpha * Vx.p[i][j][k] - Cvxex * (Exp + Ex.p[i][j][k]);
#endif
                }
            }
            pml.updateExIn(k, Ex, Hz, ID1, CB, dy);
        }
        pml.updateExOut(Ex, Hy, ID1, CB, dz);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Ey
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 0; k < Kmax - 1; ++k) {

            for (i = 1; i < Imax - 1; ++i) {

                for (j = 0; j < Jmax - 1; ++j) {
#ifdef WITH_DENSITY
                    MyDataF Eyp = Ey.p[i][j][k];
#endif
                    id = ID2.p[i][j][k];
                    if (id == 1) { // PEC

                        Ey.p[i][j][k] = 0;

                    } else {
#ifdef WITH_DENSITY
                        Ey.p[i][j][k] = CA[id] * Ey.p[i][j][k] * Ceyey.p[i][j][k] + CB[id] * Ceyh.p[i][j][k]*
                                ((Hz.p[i - 1][j][k] - Hz.p[i][j][k]) * pml.den_ex.p[i] +
                                (Hx.p[i][j][k + 1] - Hx.p[i][j][k]) * pml.den_ez.p[k]) +
                                Ceyvy.p[i][j][k] * Vy.p[i][j][k];
#else
                        Ey.p[i][j][k] = CA[id] * Ey.p[i][j][k] + CB[id] *
                                ((Hz.p[i - 1][j][k] - Hz.p[i][j][k]) * pml.den_ex.p[i] +
                                (Hx.p[i][j][k + 1] - Hx.p[i][j][k]) * pml.den_ez.p[k]);
#endif
                    }
#ifdef WITH_DENSITY
                    if (isnan(Ey.p[i][j][k])) {
                        cout << "(" << i << "," << j << "," << k << ")" << endl;
                    }
                    Vy.p[i][j][k] = alpha * Vy.p[i][j][k] - Cvyey * (Eyp + Ey.p[i][j][k]);
#endif
                }
            }
            pml.updateEyIn(k, Ey, Hz, ID2, CB, dx);
        }
        pml.updateEyOut(Ey, Hx, ID2, CB, dy);

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Ez
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 1; k < Kmax - 1; ++k) {

            for (i = 1; i < Imax - 1; ++i) {

                for (j = 1; j < Jmax - 1; ++j) {
#ifdef WITH_DENSITY
                    MyDataF Ezp = Ez.p[i][j][k];
#endif
                    id = ID3.p[i][j][k];
                    if (id == 1) { // PEC

                        Ez.p[i][j][k] = 0;

                    } else {
#ifdef WITH_DENSITY
                        Ez.p[i][j][k] = CA[id] * Ez.p[i][j][k] * Cezez.p[i][j][k] + CB[id] * Cezh.p[i][j][k]
                                * ((Hy.p[i][j][k] - Hy.p[i - 1][j][k]) * pml.den_ex.p[i] +
                                (Hx.p[i][j - 1][k] - Hx.p[i][j][k]) * pml.den_ey.p[j]) +
                                Cezvz.p[i][j][k] * Vz.p[i][j][k];
#else
                        Ez.p[i][j][k] = CA[id] * Ez.p[i][j][k] + CB[id]
                                * ((Hy.p[i][j][k] - Hy.p[i - 1][j][k]) * pml.den_ex.p[i] +
                                (Hx.p[i][j - 1][k] - Hx.p[i][j][k]) * pml.den_ey.p[j]);
#endif
                    }
#ifdef WITH_DENSITY
                    Vz.p[i][j][k] = alpha * Vz.p[i][j][k] - Cvzez * (Ezp + Ez.p[i][j][k]);
#endif
                }
            }
            pml.updateEz(k, Ez, Hx, Hy, ID3, CB, dx, dy);
        }

        //-----------------------------------------------------------
        //   Apply a point source (Soft)
        //-----------------------------------------------------------

        source = amp * -2.0 * ((n * dt - t0) / tw)
                * exp(-pow(((n * dt - t0) / tw), 2)); //Differentiated Gaussian pulse

        Ez.p[isp][jsp][ksp] = Ez.p[isp][jsp][ksp] - CB[ID3.p[isp][jsp][ksp]] * source;


        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  WRITE TO OUTPUT FILES
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef WITH_DENSITY
        UpdateErms();
        if (n % neTotalStep == 0) {
            InterpErms();
            UpdateDensity();
            updateCoeff();
        }
#endif
        if ((n % save_modulus) == 0) {
            writeField(n);
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

//Builds an object

void fdtd::buildObject() {

    //buildSphere();
    //buildDipole();
}

//Builds a sphere (Sample code - NOt used in this program)

void fdtd::buildSphere() {

    double dist; //distance
    double rad = 8; //(double)Imax / 5.0; // sphere radius
    double sc = (double) Imax / 2.0; //sphere centre
    //double rad2 = 0.3; //(double)Imax / 5.0 - 3.0; // sphere radius

    unsigned i, j, k;

    for (i = 0; i < Imax; ++i) {

        for (j = 0; j < Jmax; ++j) {

            for (k = 0; k < Kmax; ++k) {

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
    unsigned centre = (jstart + jend) / 2;

    for (i = istart; i <= iend; ++i) {

        for (j = jstart; j <= jend; ++j) {

            for (k = kstart; k <= kend; ++k) {

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
    strstream ss;
    // form file name
    ss << "E_Field_" << iteration << ".txt";
    //string fileBaseName(ss.str());
    // open file
    out.open(ss.str());
    if (out.is_open()) {

        for (i = 0; i < Imax - 1; i++) {
            for (j = 0; j < Jmax - 1; j++) { // |E|
                out<<sqrt(pow(Ex.p[i][j][ksource], 2) +
                        pow(Ey.p[i][j][ksource], 2) + pow(Ez.p[i][j][ksource], 2))<<'\t';
            }
            out<<endl;
        }
        out.close();
    }
}

//start up

void fdtd::StartUp() {
    cout << "initializing(in Startup)..." << endl;
    initialize();
    //    cout << "initial pml (in Statup)" << endl;
    //    pml.Initial(Imax, Jmax, Kmax, 11);
    cout << "setUp (in Startup)" << endl;
    setUp();
    cout << "buildObject (in Startup)" << endl;
    buildObject();
    cout << "initial CPML (in Startup)" << endl;
    pml.initCPML(dt, dx, dy, dz);
    cout << "computing (in Startup)" << endl;
    compute();
    cout << "exit Startup" << endl;
}

void fdtd::putvars() {
    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;
    cout << "dz = " << dz << endl;
    cout << "(Imax,Jmax,Kmax) = (" << Imax << "," << Jmax << "," << Kmax << ")" << endl;
    // time step increment
    cout << "dt = " << dt << endl;

    //  Specify the Impulsive Source (Differentiated Gaussian) parameters
    cout << "tw = " << tw << endl; //pulse width
    cout << "t0 = " << t0 << endl; //delay
    cout << "source = " << source << endl; //Differentiated Gaussian source
    cout << "amp = " << amp << endl; // Amplitude

    //Specify the Time Step at which the data has to be saved for Visualization
    cout << "save_modulus = " << save_modulus << endl;

    //  Specify the dipole Boundaries(A cuboidal rode- NOT as a cylinder)
    cout << "(istart, iend, jstart) = (" << istart << ',' << iend << ',' << jstart << ')' << endl;
    cout << "(jend, kstart, kend) = (" << jend << ',' << kstart << ',' << kend << ')' << endl;

    //Output recording point
    cout << "ksource = " << ksource << endl;

    //  Specify the CPML Order and Other Parameters:
    cout << " m = " << m << endl;
    cout << " ma = " << ma << endl;
}
// =================================================================
// MATLAB SIMULATION
// =================================================================
#ifdef MATLAB_SIMULATION

void fdtd::initMatlabSimulation() {
    Ez.InitMatlabEngine();
    Ez.InitPlot();

}

void fdtd::doMatlabSimulation() {
    Ez.PlotArrays();
}

void fdtd::finishMatlabSimulation() {
    Ez.ClearSim();
    Ez.CloseEngine();
}
#endif
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// END OF PROGRAM CPMLFDTD3D
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
