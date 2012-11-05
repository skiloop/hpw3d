

#include <math.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>

#include "fdtd.h"
#include "InonizationFormula.h"

extern MyDataF eps_0, epsR;
extern MyDataF mu_0;
extern MyDataF dt, dx, dy, dz;
extern MyDataF pi, C, me, e;

using namespace std;
#ifdef WITH_DENSITY

fdtd::fdtd(unsigned _nmax, unsigned _imax, unsigned _jmax, unsigned _kmax,
        MyDataF _tw, MyDataF _dx, MyDataF _dy, MyDataF _dz,
        unsigned _amp, unsigned _savemodulus, unsigned _ksource,
        unsigned _m, unsigned _ma, unsigned _nmaterial, unsigned _neGrid)
: nMax(_nmax), Imax(_imax), Jmax(_jmax), Kmax(_kmax)
, tw(_tw), dx(_dx), dy(_dy), dz(_dz)
, amp(_amp), save_modulus(_savemodulus), ksource(_ksource)
, m(_m), ma(_ma)
, neGrid(_neGrid)
, numMaterials(_nmaterial)
, Ne0(1e7)
, pml(_m, _ma) {
}
#else

fdtd::fdtd(unsigned _nmax, unsigned _imax, unsigned _jmax, unsigned _kmax,
        MyDataF _tw, MyDataF _dx, MyDataF _dy, MyDataF _dz,
        unsigned _amp, unsigned _savemodulus, unsigned _ksource,
        unsigned _m, unsigned _ma, unsigned _nmaterial)
: nMax(_nmax), Imax(_imax), Jmax(_jmax), Kmax(_kmax)
, tw(_tw), dx(_dx), dy(_dy), dz(_dz)
, amp(_amp), save_modulus(_savemodulus), ksource(_ksource)
, m(_m), ma(_ma)
, numMaterials(_nmaterial)
, pml(_m, _ma) {
}
#endif

fdtd::~fdtd(void) {
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
    return 0;
}

int fdtd::UpdateDensity(void) {

    unsigned i, j, k, mt = 1;

    MyDataF Eeff, alpha_t, tau_m, kasi;
    MyDataF Ne_ijk, Neip1, Neim1, Nejm1, Nejp1, Nekp1, Nekm1;
    MyDataF Deff;
    MyDataF maxvi = 0, minvi = 0;
    //MyDataF tmp=0;
    MyDataF vi, va;

    int ci = 0, cj = 0, ck = 0;
    Ne_pre = Ne;
    //printf("\n%4.3e\t",Ne.p[0][0][0]);
    for (i = mt; i < Ne.nx - mt; i++)
        for (j = mt; j < Ne.ny - mt; j++)
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
                    case 1:Niu_MorrowAndLowke(&vi, &va, Eeff, Ne_ijk * 1e6);
                        break;
                    case 2:Niu_Nikonov(&vi, &va, Eeff, p);
                        break;
                    case 3:Niu_Kang(&vi, &va, Eeff);
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
                //				if((i==spx*FiNeGridRatio)&&(j==spy*FiNeGridRatio)&&(k==spz*FiNeGridRatio))
                //				{
                //					printf("Deff = %6.5e\n",Deff);
                //					printf("vi\t=\t%6.5e\n",vi);
                //					printf("in\t=\t%7.6e\n",Deff*dtf*(Neip1+Neim1+Nejp1+Nejm1+Nekp1+Nekm1-6*Ne_ijk)/dtf/dtf);
                //					printf("Erms\t=\t%7.6e\n",Erms.p[i][j][k]);
                //					printf("va\t=\t%6.5e\n",va);
                //					printf("alpha_t\t=\t%6.5e\n",alpha_t);
                //				}
                //printf("%6.4e\t%6.4e\t",Ne.p[spx*FiNeGridRatio][spy*FiNeGridRatio][spz*FiNeGridRatio],Ne.p[nx][ny][nz]);
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
    //printf("%4.3e\n",Ne.p[0][0][0]);
    WallCircleBound(Ne);
    //ApplyDensityBound(Ne,mt);
    //printf("%4.3e\n",Ne.p[0][0][0]);
    printf("%6.5e\t", Ne.p[Ne.nx / 2][Ne.ny / 2][Ne.nz / 2]);
    printf("\t%6.4e\t%6.4e\t%6.4e\t%6.4e\n", maxvi, minvi, Ne.p[ci][cj][ck], Erms.p[ci][cj][ck]);
    return 0;
}

int fdtd::UpdateVeloity(void) {
    return 0;
}

void fdtd::WallCircleBound(data3d &stru) {
    int i, j, k;
    int endx, endy, endz;

    endx = stru.nx - 1;
    endy = stru.ny - 1;
    endz = stru.nz - 1;

    //Buttom and top

    for (i = 1; i < endx; i++)
        for (j = 1; j < endy; j++) {
            stru.p[i][j][0] = 2 * stru.p[i][j][1] - stru.p[i][j][2];
            stru.p[i][j][endz] = 2 * stru.p[i][j][endz - 1] - stru.p[i][j][endz - 2];
        }
    //left and right
    for (j = 1; j < endy; j++)
        for (k = 1; k < endz; k++) {
            stru.p[0][j][k] = 2 * stru.p[1][j][k] - stru.p[2][j][k];
            stru.p[endx][j][k] = 2 * stru.p[endx - 1][j][k] - stru.p[endx - 2][j][k];
        }
    //front and back
    for (i = 1; i < endx; i++)
        for (k = 1; k < endz; k++) {
            stru.p[i][0][k] = 2 * stru.p[i][1][k] - stru.p[i][2][k];
            stru.p[i][endy][k] = 2 * stru.p[i][endy - 1][k] - stru.p[i][endy - 2][k];
        }
}

#endif

void fdtd::initialize() {

    unsigned i, j, k;
#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << "numMaterials = " << numMaterials << endl;
#endif
    //Dynamic memory allocation
    epsilon = (double *) malloc((numMaterials) * sizeof (double));

    for (i = 0; i < numMaterials; i++) {

        epsilon[i] = eps_0;
    }

    mu = (double *) malloc((numMaterials) * sizeof (double));

    for (i = 0; i < numMaterials; i++) {

        mu[i] = mu_0;
    }

    sigma = (double *) malloc((numMaterials) * sizeof (double));

    for (i = 0; i < numMaterials; i++) {

        sigma[i] = 0.0;
    }

    CA = (double *) malloc((numMaterials) * sizeof (double));

    for (i = 0; i < numMaterials; i++) {

        CA[i] = 0.0;
    }

    CB = (double *) malloc((numMaterials) * sizeof (double));

    for (i = 0; i < numMaterials; i++) {

        CB[i] = 0.0;
    }
#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << "Imax=" << Imax << endl;
    cout << "Jmax=" << Jmax << endl;
    cout << "Kmax=" << Kmax << endl;
#endif
    Ez.CreateStruct(Imax, Jmax, Kmax,0);   
    Ey.CreateStruct(Imax, Jmax - 1, Kmax - 1,0);
    Ex.CreateStruct(Imax - 1, Jmax, Kmax - 1,0);

    Hx.CreateStruct(Imax, Jmax - 1, Kmax,0);
    Hy.CreateStruct(Imax - 1, Jmax, Kmax,0);
    Hz.CreateStruct((Imax - 1), (Jmax - 1), (Kmax - 1),0);

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
#endif

#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << "creating ID1..." << endl;
#endif
    ID1.CreateStruct(Imax,Jmax,Kmax,0);

#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << "creating ID2..." << endl;
#endif
    ID2.CreateStruct(Imax,Jmax,Kmax,0);

#if(DEBUG>=3)
    cout << __FILE__ << ":" << __LINE__ << endl;
    cout << "creating ID3..." << endl;
#endif
    ID3.CreateStruct(Imax,Jmax,Kmax,0);
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

    //Fine Time Step Size
    dtf = dt / neGrid;

    mu_e = e / me / vm; //3.7e-2;
    mu_i = mu_e / 100.0; //mu_e/mu_i ranges from 100 to 200
    De = mu_e * 2 * 1.602e-19 / e; //
    Da = De * mu_i / mu_e;
#endif
    //  Specify the dipole size 
    istart = 24;
    iend = 26;
    jstart = 55;
    jend = 71;
    kstart = 11;
    kend = 13;

    // source position
    isp = 25;
    jsp = 63;
    ksp = 12;

    if (iend > Imax)iend = Imax - 1;
    if (jend > Jmax)jend = Jmax - 1;
    if (kend > Kmax)kend = Kmax - 1;
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

    printf("\nTIme step = %e", dt);
    printf("\n Number of steps = %d", nMax);
    printf("\n Total Simulation time = %e Seconds", nMax * dt);

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SET CPML PARAMETERS IN EACH DIRECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void fdtd::compute() {

    unsigned id, n;
    unsigned i, j, k;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  BEGIN TIME STEP
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    printf("\nBegin time-stepping...\n");

    for (n = 1; n <= nMax; ++n) {

        printf("Ez at time step %d at (%d, %d, %d) :  %f\n", n, isp, 40, ksp, Ez.p[isp][40][ksp]);

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
        }
        pml.UpdatePMLForHx(Hx, Ey, Ez, DB);
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
        }
        pml.UpdatePMLForHy(Hy, Ez, Ex, DB);
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
        }
        pml.UpdatePMLForHz(Hz, Ex, Ey, DB);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Ex
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 0; k < Kmax - 1; ++k) {

            for (i = 0; i < Imax - 1; ++i) {

                for (j = 1; j < Jmax - 1; ++j) {

                    id = ID1.p[i][j][k];
                    if (id == 1) { // PEC

                        Ex.p[i][j][k] = 0;

                    } else {

                        Ex.p[i][j][k] = CA[id] * Ex.p[i][j][k] + CB[id] *
                                ((Hz.p[i][j][k] - Hz.p[i][j - 1][k]) * pml.den_ey.p[j] +
                                (Hy.p[i][j][k] - Hy.p[i][j][k + 1]) * pml.den_ez.p[k]);
                    }
                }
            }
        }
        pml.UpdatePMLForEx(Ex, Hy, Hz, CB, ID1);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Ey
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 0; k < Kmax - 1; ++k) {

            for (i = 1; i < Imax - 1; ++i) {

                for (j = 0; j < Jmax - 1; ++j) {

                    id = ID2.p[i][j][k];
                    if (id == 1) { // PEC

                        Ey.p[i][j][k] = 0;

                    } else {

                        Ey.p[i][j][k] = CA[id] * Ey.p[i][j][k] + CB[id] *
                                ((Hz.p[i - 1][j][k] - Hz.p[i][j][k]) * pml.den_ex.p[i] +
                                (Hx.p[i][j][k + 1] - Hx.p[i][j][k]) * pml.den_ez.p[k]);
                    }
                }
            }
        }
        pml.UpdatePMLForEy(Ey, Hz, Hx, CB, ID2);

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Ez
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 1; k < Kmax - 1; ++k) {

            for (i = 1; i < Imax - 1; ++i) {

                for (j = 1; j < Jmax - 1; ++j) {

                    id = ID3.p[i][j][k];
                    if (id == 1) { // PEC

                        Ez.p[i][j][k] = 0;

                    } else {

                        Ez.p[i][j][k] = CA[id] * Ez.p[i][j][k] + CB[id]
                                * ((Hy.p[i][j][k] - Hy.p[i - 1][j][k]) * pml.den_ex.p[i] +
                                (Hx.p[i][j - 1][k] - Hx.p[i][j][k]) * pml.den_ey.p[j]);
                    }
                }
            }
        }
        pml.UpdatePMLForEz(Ez, Hx, Hy, CB, ID3);

        //-----------------------------------------------------------
        //   Apply a point source (Soft)
        //-----------------------------------------------------------

        source = amp * -2.0 * ((n * dt - t0) / tw)
                * exp(-pow(((n * dt - t0) / tw), 2)); //Differentiated Gaussian pulse

        Ez.p[isp][jsp][ksp] = Ez.p[isp][jsp][ksp] - CB[ID3.p[isp][jsp][ksp]] * source;


        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  WRITE TO OUTPUT FILES
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if ((n % save_modulus) == 0) {

            writeField(n);
        }

    }
    //  END TIME STEP
    printf("Done time-stepping...\n");

}

//Builds an object

void fdtd::buildObject() {

    //buildSphere();
    buildDipole();
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

//creates a dielctric cube (yee cell) made up of the selected material

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
    FILE *ptr;
    char step[10];
    char fileBaseName[100] = "E_Field_";
    sprintf(step, "%d", iteration);
    strcat(fileBaseName, step);
    strcat(fileBaseName, ".txt");

    ptr = fopen(fileBaseName, "wt");

    for (i = 0; i < Imax - 1; i++) {
        for (j = 0; j < Jmax - 1; j++) { // |E|
            fprintf(ptr, "%f\t", sqrt(pow(Ex.p[i][j][ksource], 2) +
                    pow(Ey.p[i][j][ksource], 2) + pow(Ez.p[i][j][ksource], 2)));
            //	fprintf(ptr, "%f\t", Ex.p[i][j][ksource]);//Ex
            //	fprintf(ptr, "%f\t", Ey.p[i][j][ksource]);//Ey
            //	fprintf(ptr, "%f\t", Ez.p[i][j][ksource]);//Ez
            //	fprintf(ptr, "%f\t", Hx.p[i][j][ksource]);//Hx
            //	fprintf(ptr, "%f\t", Hy.p[i][j][ksource]);//Hy
            //	fprintf(ptr, "%f\t", Hz.p[i][j][ksource]);//Hz

            // |H|
            //	fprintf(ptr, "%f\t", sqrt(pow(Hx.p[i][j][ksource], 2) +
            //		pow(Hy.p[i][j][ksource], 2) + pow( Hz.p[i][j][ksource], 2)));

        }
        fprintf(ptr, "\n");
    }
    fclose(ptr);
}

//start up

void fdtd::StartUp() {
    cout << "initializing(in Statup)..." << endl;
    initialize();
    cout << "initial pml (in Statup)" << endl;
    pml.Initial(Imax, Jmax, Kmax, 11);
    cout << "setUp (in Statup)" << endl;
    setUp();
    cout << "buildObject (in Statup)" << endl;
    buildObject();
    cout << "computing (in Statup)" << endl;
    compute();
    cout << "exit Statup" << endl;
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
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// END OF PROGRAM CPMLFDTD3D
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
