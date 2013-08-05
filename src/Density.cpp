#include "Density.h"
#include "fdtd.h"
#include "InonizationFormula.h"

Density::Density(void)
{
}


Density::~Density(void)
{
}


void Density::SetPlasmaVar(MyDataF _rei, MyDataF _vm, MyDataF _p, int _ftype) {
    rei = _rei;
    vm = _vm;
    p = _p;
    niutype = _ftype;
}

void Density::IntegerEeff() {
    unsigned i, j, k;
    unsigned io, jo, ko;
    //MyDataF vxIJK,vyIJK;
    for (i = istart, io = istart * neGrid; i <= iend; i++, io += neGrid) {
        for (j = jstart, jo = jstart * neGrid; j <= jend; j++, jo += neGrid) {
            for (k = kstart, ko = kstart * neGrid; k <= kend; k++, ko += neGrid) {
                Erms.p[io][jo][ko] = m / e * sqrt(Erms.p[io][jo][ko] / dtf * Nu_c.p[i][j][k] / 2);
            }
        }
    }
}

void Density::UpdateErms(int srcType,data3d<MyDataF> Ex,data3d<MyDataF> Ey,data3d<MyDataF> Ez) {
    unsigned i, j, k;
    unsigned io, jo, ko;
    switch (srcType) {
	case fdtd::SOURCE_GAUSSIAN:
            MyDataF vxIJK, vyIJK;
            for (i = istart, io = istart * neGrid; i <= iend; i++, io += neGrid) {
                for (j = jstart, jo = jstart * neGrid; j <= jend; j++, jo += neGrid) {
                    for (k = kstart, ko = kstart * neGrid; k <= kend; k++, ko += neGrid) {
                        vxIJK = (Vx.p[i - 1][j][k - 1] + Vx.p[i + 1][j][k - 1] + Vx.p[i - 1][j][k + 1] + Vx.p[i + 1][j][k + 1]) / 4;
                        vyIJK = (Vy.p[i][j - 1][k - 1] + Vy.p[i][j + 1][k - 1] + Vy.p[i][j - 1][k + 1] + Vy.p[i][j + 1][k + 1]) / 4;
                        Erms.p[io][jo][ko] += (Vz.p[i][j][k] * Vz.p[i][j][k] + vxIJK * vxIJK + vyIJK * vyIJK) * dt;
                    }
                }
            }
            break;
	case fdtd::SOURCE_SINE:
        default:
            MyDataF exIJK, eyIJK;
            for (i = istart, io = istart * neGrid; i <= iend; i++, io += neGrid) {
                for (j = jstart, jo = jstart * neGrid; j <= jend; j++, jo += neGrid) {
                    for (k = kstart, ko = kstart * neGrid; k <= kend; k++, ko += neGrid) {
                        exIJK = (Ex.p[i - 1][j][k - 1] + Ex.p[i + 1][j][k - 1] + Ex.p[i - 1][j][k + 1] + Ex.p[i + 1][j][k + 1]) / 4;
                        eyIJK = (Ey.p[i][j - 1][k - 1] + Ey.p[i][j + 1][k - 1] + Ey.p[i][j - 1][k + 1] + Ey.p[i][j + 1][k + 1]) / 4;
                        Erms.p[io][jo][ko] = sqrt(Ez.p[i][j][k] * Ez.p[i][j][k] + exIJK * exIJK + eyIJK * eyIJK);
                    }
                }
            }
    }
}

void Density::updateCollisionFrequency() {
    int i, j, k;
    //unsigned io, jo, ko;
    MyDataF EeffDivP;
    MyDataF DivParam = 100 * p * 133.3;
    MyDataF C1 = 5.20e8 * p;
    MyDataF C2 = 2.93e8 * p;
    MyDataF C3 = 3.24e8 * p;
    //    for (i = istart, io = istart * neGrid; i <= iend; i++, io += neGrid) {
    //        for (j = jstart, jo = jstart * neGrid; j <= jend; j++, jo += neGrid) {
    //            for (k = kstart, ko = kstart * neGrid; k <= kend; k++, ko += neGrid) {
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

void Density::InterpErms() {
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
}

void Density::UpdateDensity(int srcType) {

    int i, j, k, mt = 1;

    MyDataF Eeff, alpha_t, tau_m, kasi;
    MyDataF Ne_ijk, Neip1, Neim1, Nejm1, Nejp1, Nekp1, Nekm1;
    MyDataF Deff;
    MyDataF maxvi = 0, minvi = 0;
    MyDataF vi, va;

    unsigned ci = 0, cj = 0, ck = 0;
    Ne_pre = Ne;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) \
        schedule(dynamic) private(i,j,k,Eeff,Ne_ijk, Neip1, Neim1, Nejm1, Nejp1, Nekp1, Nekm1,vi,va,alpha_t,Deff,tau_m,kasi)
#endif
    for (i = mt; i < Ne.nx - mt; i++) {
        for (j = mt; j < Ne.ny - mt; j++) {
            for (k = mt; k < Ne.nz - mt; k++) {
                switch (srcType) {
				case fdtd::SOURCE_GAUSSIAN:
                        Eeff = Erms.p[i][j][k] / 100; //convert to V/cm
                        break;
                    default:
                        Eeff = Erms.p[i][j][k] / 100 * pow(1 / (1 + omega * omega / vm / vm), 0.5);
                }

                Ne_ijk = Ne_pre.p[i][j][k];
                Neip1 = Ne_pre.p[i + 1][j][k];
                Neim1 = Ne_pre.p[i - 1][j][k];
                Nejp1 = Ne_pre.p[i][j + 1][k];
                Nejm1 = Ne_pre.p[i][j - 1][k];
                Nekp1 = Ne_pre.p[i][j][k + 1];
                Nekm1 = Ne_pre.p[i][j][k - 1];

                switch (niutype) {
                    case MORROW_AND_LOWKE:
                        Niu_MorrowAndLowke(&vi, &va, Eeff, Ne_ijk * 1e6);
                        break;
                    case NIKONOV:
                        Niu_Nikonov(&vi, &va, Eeff, p);
                        break;
                    case KANG:
                        Niu_Kang(&vi, &va, Eeff);
                        break;
                    case ALI:
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
}

void Density::UpdateVeloity(void) {
    //return 0;
}

void Density::WallCircleBound(data3d<MyDataF> &stru) {
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

void Density::createCoeff(int srcType,data3d<MyDataF> Ex,data3d<MyDataF> Ey,data3d<MyDataF> Ez)  {
    // velocity coefficients
    //Cvxex.CreateStruct(Vx,0.0);
    //Cvyey.CreateStruct(Vy,0.0);
    //Cvzez.CreateStruct(Vz,0.0);
    // electricity coefficients
    //Cexe.CreateStruct(Ex, 0.0);
    //Ceye.CreateStruct(Ey, 0.0);
    //Ceze.CreateStruct(Ez, 0.0);
    Cexvx.CreateStruct(Ex, 0.0);
    Ceyvy.CreateStruct(Ey, 0.0);
    Cezvz.CreateStruct(Ez, 0.0);
    // beta
    beta.CreateStruct(Ne, 0.0);
    // velocity coefficients
    if (srcType == fdtd::SOURCE_GAUSSIAN) {
        Cvxex_guassian.CreateStruct(Vx, 0.0);
        Cvyey_guassian.CreateStruct(Vy, 0.0);
        Cvzez_guassian.CreateStruct(Vz, 0.0);
    }
}

void Density::initCoeff(const fdtd &mfdtd) {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Velocity Coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    a = vm * mfdtd.dt / 2;
    gamma = 1 + a;
    alpha = (1 - a) / gamma;
    Cvxex = Cvyey = Cvzez = e * mfdtd.dt / 2 / me / gamma;
    Coeff_velocity = half_e * mfdtd.dt / me;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update collision frequency
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (mfdtd.srcType == fdtd::SOURCE_GAUSSIAN) {
        Nu_c.ResetStructData(vm);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //electricity coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    updateCoeff(mfdtd);

}

void Density::updateCoeff(const fdtd &mfdtd) {

	updateBeta(mfdtd.srcType);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //electricity coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int i, j, k;
    unsigned im, jm, km;
    MyDataF tmp = eMDtDiv2DivEps0 * (1 + alpha);
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k,im,jm,km)//shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif    
    for (j = 0; j < mfdtd.Ex.ny; j++) {
        jm = j*neGrid;
        for (i = 0, im = halfNeGrid; i < mfdtd.Ex.nx; i++, im += neGrid) {
            for (k = 0, km = halfNeGrid; k < mfdtd.Ex.nz; k++, km += neGrid) {
                MyDataF kappa = (1 + beta.p[im][jm][km]);
                mfdtd.Cexe.p[i][j][k] = (1 - beta.p[im][jm][km]) / kappa;
                mfdtd.Cexhy.p[i][j][k] = -dtDivEps0DivDz / kappa;
                mfdtd.Cexhz.p[i][j][k] = dtDivEps0DivDy / kappa;
                if (mfdtd.srcType == fdtd::SOURCE_GAUSSIAN) {
                    MyDataF a = half_dt * Nu_c.p[im][jm][km];
                    MyDataF gamma_t = 1 + a;
                    MyDataF tmpc = eMDtDiv2DivEps0 * (1 + (1 - a) / gamma_t);
                    Cvxex_guassian.p[i][j][k] = Coeff_velocity / gamma_t;
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
    for (i = 0; i < mfdtd.Ey.nx; i++) {
        im = i*neGrid;
        for (j = 0, jm = halfNeGrid; j < mfdtd.Ey.ny; j++, jm += neGrid) {
            for (k = 0, km = halfNeGrid; k < mfdtd.Ey.nz; k++, km += neGrid) {
                MyDataF kappa = (1 + beta.p[im][jm][km]);
                mfdtd.Ceye.p[i][j][k] = (1 - beta.p[im][jm][km]) / kappa;
                mfdtd.Ceyhx.p[i][j][k] = dtDivEps0DivDz / kappa;
                mfdtd.Ceyhz.p[i][j][k] = -dtDivEps0DivDx / kappa;
                if (mfdtd.srcType == fdtd::SOURCE_GAUSSIAN) {
                    MyDataF a = half_dt * Nu_c.p[im][jm][km];
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
    for (i = 0; i < mfdtd.Ez.nx; i++) {
        im = i*neGrid;
        for (j = 0, jm = 0; j < mfdtd.Ez.ny; j++, jm += neGrid) {
            for (k = 0, km = 0; k < mfdtd.Ez.nz; k++, km += neGrid) {
                MyDataF kappa = (1 + beta.p[im][jm][km]);
                mfdtd.Ceze.p[i][j][k] = (1 - beta.p[im][jm][km]) / kappa;
                mfdtd.Cezhy.p[i][j][k] = dtDivEps0DivDx / kappa;
                mfdtd.Cezhx.p[i][j][k] = -dtDivEps0DivDy / kappa;
                if (mfdtd.srcType == fdtd::SOURCE_GAUSSIAN) {
                    MyDataF a = half_dt * Nu_c.p[im][jm][km];
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

void Density::updateBeta(int srcType) {
    unsigned is = istart*neGrid;
    unsigned ie = iend*neGrid;
    unsigned js = jstart*neGrid;
    unsigned je = jend*neGrid;
    unsigned ks = kstart*neGrid;
    unsigned ke = kend*neGrid;
    MyDataF temp = e2Dt2Div4DivEps0DivMe;
    if (srcType != fdtd::SOURCE_GAUSSIAN) {
        temp = temp / gamma;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)  //shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif
        for (int i = is; i < ie; i++) {
            for (unsigned j = js; j < je; j++) {
                for (unsigned k = ks; k < ke; k++) {
                    beta.p[i][j][k] = temp * Ne.p[i][j][k];
                }
            }
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif
        for (int i = is; i < ie; i++) {
            for (unsigned j = js; j < je; j++) {
                for (unsigned k = ks; k < ke; k++) {
                    beta.p[i][j][k] = temp / (1 + half_dt * Nu_c.p[i][j][k]) * Ne.p[i][j][k];
                }
            }
        }
    }
}

void Density::initDensity() {
    MyDataF tmp=pow(50e-6,3);
    int i0=isp;
    int j0=jsp+30;
    int k0=ksp;
    for(unsigned i=0;i<Ne.nx;i++){
        for(unsigned j=0;j<Ne.ny;j++){
            for(unsigned k=0;k<Ne.nz;k++){
                Ne.p[i][j][k]=Ne0*exp((pow((i-i0)*dx,2)+pow((j-j0)*dx,2)+pow((k-k0)*dx,2))/tmp);
            }
        }
    }
}
