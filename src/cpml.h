/* 
 * File:   cpml.h
 * Author: skiloop
 *
 * Created on April 3, 2013, 8:41 AM
 */

#ifndef CPML_H
#define	CPML_H

#include <math.h>
#include "datastruct.h"

template<class type1, class type2>
class cpml {
public:

    //  Fundamental Constants (MKS units)
    const static double pi;
    const static double C;
    static double mu_0;
    static double eps_0;

    // Specify Material Relative Permittivity and Conductivity
    static double epsR; //free space


    //  CPML components (Taflove 3rd Edition, Chapter 7)
    data3d<type1> psi_Ezx_1;
    data3d<type1> psi_Ezx_2;
    data3d<type1> psi_Hyx_1;
    data3d<type1> psi_Hyx_2;
    data3d<type1> psi_Ezy_1;
    data3d<type1> psi_Ezy_2;
    data3d<type1> psi_Hxy_1;
    data3d<type1> psi_Hxy_2;
    data3d<type1> psi_Hxz_1;
    data3d<type1> psi_Hxz_2;
    data3d<type1> psi_Hyz_1;
    data3d<type1> psi_Hyz_2;
    data3d<type1> psi_Exz_1;
    data3d<type1> psi_Exz_2;
    data3d<type1> psi_Eyz_1;
    data3d<type1> psi_Eyz_2;
    data3d<type1> psi_Hzx_1;
    data3d<type1> psi_Eyx_1;
    data3d<type1> psi_Hzx_2;
    data3d<type1> psi_Eyx_2;
    data3d<type1> psi_Hzy_1;
    data3d<type1> psi_Exy_1;
    data3d<type1> psi_Hzy_2;
    data3d<type1> psi_Exy_2;

    data1d<type1> be_x_1, ce_x_1, alphae_x_PML_1, sige_x_PML_1, kappae_x_PML_1;
    data1d<type1> bh_x_1, ch_x_1, alphah_x_PML_1, sigh_x_PML_1, kappah_x_PML_1;
    data1d<type1> be_x_2, ce_x_2, alphae_x_PML_2, sige_x_PML_2, kappae_x_PML_2;
    data1d<type1> bh_x_2, ch_x_2, alphah_x_PML_2, sigh_x_PML_2, kappah_x_PML_2;
    data1d<type1> be_y_1, ce_y_1, alphae_y_PML_1, sige_y_PML_1, kappae_y_PML_1;
    data1d<type1> bh_y_1, ch_y_1, alphah_y_PML_1, sigh_y_PML_1, kappah_y_PML_1;
    data1d<type1> be_y_2, ce_y_2, alphae_y_PML_2, sige_y_PML_2, kappae_y_PML_2;
    data1d<type1> bh_y_2, ch_y_2, alphah_y_PML_2, sigh_y_PML_2, kappah_y_PML_2;
    data1d<type1> be_z_1, ce_z_1, alphae_z_PML_1, sige_z_PML_1, kappae_z_PML_1;
    data1d<type1> bh_z_1, ch_z_1, alphah_z_PML_1, sigh_z_PML_1, kappah_z_PML_1;
    data1d<type1> be_z_2, ce_z_2, alphae_z_PML_2, sige_z_PML_2, kappae_z_PML_2;
    data1d<type1> bh_z_2, ch_z_2, alphah_z_PML_2, sigh_z_PML_2, kappah_z_PML_2;

    // denominators for the update equations
    data1d<type1> den_ex;
    data1d<type1> den_hx;
    data1d<type1> den_ey;
    data1d<type1> den_hy;
    data1d<type1> den_ez;
    data1d<type1> den_hz;

    cpml() {
    };

    cpml(const cpml<type1, type2>& orig) {
    };

    virtual ~cpml() {
    };
   
    void updateEz(const int k, data3d<type1> &Ez, const data3d<type1>&Hx, const data3d<type1>&Hy, const data3d<type2>&ID3, const type1* CB, const double dx, const double dy) {
        int i, j, ii, jj;
        short id;
        for (j = 1; j < Jmax - 1; ++j) {
            //............................................................
            //  PML for bottom Ez, x-direction
            //.............................................................
            for (i = 1; i < nxPML_1; ++i) {
                id = ID3.p[i][j][k];
                psi_Ezx_1.p[i][j][k] = be_x_1.p[i] * psi_Ezx_1.p[i][j][k]
                        + ce_x_1.p[i] * (Hy.p[i][j][k] - Hy.p[i - 1][j][k]) / dx;
                Ez.p[i][j][k] = Ez.p[i][j][k] + CB[id] * psi_Ezx_1.p[i][j][k];
            }
            //............................................................
            //  PML for top Ez, x-direction
            //............................................................
            ii = nxPML_2 - 1;
            for (i = Imax - nxPML_2; i < Imax - 1; ++i) {
                id = ID3.p[i][j][k];
                psi_Ezx_2.p[ii][j][k] = be_x_2.p[ii] * psi_Ezx_2.p[ii][j][k]
                        + ce_x_2.p[ii] * (Hy.p[i][j][k] - Hy.p[i - 1][j][k]) / dx;
                Ez.p[i][j][k] = Ez.p[i][j][k] + CB[id] * psi_Ezx_2.p[ii][j][k];
                ii = ii - 1;
            }
        }

        for (i = 1; i < Imax - 1; ++i) {
            //..........................................................
            //  PML for bottom Ez, y-direction
            //..........................................................
            for (j = 1; j < nyPML_1; ++j) {
                id = ID3.p[i][j][k];
                psi_Ezy_1.p[i][j][k] = be_y_1.p[j] * psi_Ezy_1.p[i][j][k]
                        + ce_y_1.p[j] * (Hx.p[i][j - 1][k] - Hx.p[i][j][k]) / dy;
                Ez.p[i][j][k] = Ez.p[i][j][k] + CB[id] * psi_Ezy_1.p[i][j][k];
            }
            //............................................................
            //  PML for top Ez, y-direction
            //............................................................
            jj = nyPML_2 - 1;
            for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {
                id = ID3.p[i][j][k];
                psi_Ezy_2.p[i][jj][k] = be_y_2.p[jj] * psi_Ezy_2.p[i][jj][k]
                        + ce_y_2.p[jj] * (Hx.p[i][j - 1][k] - Hx.p[i][j][k]) / dy;
                Ez.p[i][j][k] = Ez.p[i][j][k] + CB[id] * psi_Ezy_2.p[i][jj][k];
                jj = jj - 1;
            }
        }
    };

    void updateEyOut(data3d<type1> &Ey, const data3d<type1> &Hx, const data3d<type2>&ID2, const type1* CB, const double dz) {
        int i, j, k, kk;
        short id;
        for (i = 1; i < Imax - 1; ++i) {
            for (j = 0; j < Jmax - 1; ++j) {
                //...........................................................
                //  PML for bottom Ey, k-direction
                //...........................................................
                for (k = 0; k < nzPML_1; ++k) {
                    id = ID2.p[i][j][k];
                    psi_Eyz_1.p[i][j][k] = be_z_1.p[k] * psi_Eyz_1.p[i][j][k]
                            + ce_z_1.p[k] * (Hx.p[i][j][k + 1] - Hx.p[i][j][k]) / dz;
                    Ey.p[i][j][k] = Ey.p[i][j][k] + CB[id] * psi_Eyz_1.p[i][j][k];
                }
                //...........................................................
                //  PML for top Ey, k-direction
                //............................................................
                kk = nzPML_2 - 1;
                for (k = Kmax - nzPML_2 - 1; k < Kmax - 1; ++k) {

                    id = ID2.p[i][j][k];
                    psi_Eyz_2.p[i][j][kk] = be_z_2.p[kk] * psi_Eyz_2.p[i][j][kk]
                            + ce_z_2.p[kk] * (Hx.p[i][j][k + 1] - Hx.p[i][j][k]) / dz;
                    Ey.p[i][j][k] = Ey.p[i][j][k] + CB[id] * psi_Eyz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    void updateEyIn(const int k, data3d<type1> &Ey, const data3d<type1> &Hz, const data3d<type2>&ID2, const type1* CB, const double dx) {
        int i, j, ii;
        short id;
        for (j = 0; j < Jmax - 1; ++j) {
            //...........................................................
            //  PML for bottom Ey, i-direction
            //...........................................................
            for (i = 1; i < nxPML_1; ++i) {
                id = ID2.p[i][j][k];
                psi_Eyx_1.p[i][j][k] = be_x_1.p[i] * psi_Eyx_1.p[i][j][k]
                        + ce_x_1.p[i] * (Hz.p[i - 1][j][k] - Hz.p[i][j][k]) / dx;
                Ey.p[i][j][k] = Ey.p[i][j][k] + CB[id] * psi_Eyx_1.p[i][j][k];
            }
            //............................................................
            //  PML for top Ey, i-direction
            //............................................................
            ii = nxPML_2 - 1;
            for (i = Imax - nxPML_2; i < Imax - 1; ++i) {
                id = ID2.p[i][j][k];
                psi_Eyx_2.p[ii][j][k] = be_x_2.p[ii] * psi_Eyx_2.p[ii][j][k]
                        + ce_x_2.p[ii] * (Hz.p[i - 1][j][k] - Hz.p[i][j][k]) / dx;
                Ey.p[i][j][k] = Ey.p[i][j][k] + CB[id] * psi_Eyx_2.p[ii][j][k];
                ii = ii - 1;
            }
        }
    };

    void updateExIn(const int k, data3d<type1> &Ex, const data3d<type1> &Hz, const data3d<type2>&ID1, const type1* CB, const double dy) {
        int i, j, jj;
        short id;
        for (i = 0; i < Imax - 1; ++i) {
            //..............................................................
            //  PML for bottom Ex, j-direction
            //..............................................................
            for (j = 1; j < nyPML_1; ++j) {

                id = ID1.p[i][j][k];
                psi_Exy_1.p[i][j][k] = be_y_1.p[j] * psi_Exy_1.p[i][j][k]
                        + ce_y_1.p[j] * (Hz.p[i][j][k] - Hz.p[i][j - 1][k]) / dy;
                Ex.p[i][j][k] = Ex.p[i][j][k] + CB[id] * psi_Exy_1.p[i][j][k];
            }
            //.............................................................
            //  PML for top Ex, j-direction
            //.............................................................
            jj = nyPML_2 - 1;
            for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {

                id = ID1.p[i][j][k];
                psi_Exy_2.p[i][jj][k] = be_y_2.p[jj] * psi_Exy_2.p[i][jj][k]
                        + ce_y_2.p[jj] * (Hz.p[i][j][k] - Hz.p[i][j - 1][k]) / dy;
                Ex.p[i][j][k] = Ex.p[i][j][k] + CB[id] * psi_Exy_2.p[i][jj][k];
                jj = jj - 1;
            }
        }
    };

    void updateExOut(data3d<type1> &Ex, const data3d<type1> &Hy, const data3d<type2>&ID1, const type1* CB, const double dz) {
        int i, j, k, kk;
        short id;
        for (i = 0; i < Imax - 1; ++i) {

            for (j = 1; j < Jmax - 1; ++j) {
                //.............................................................
                //  PML for bottom Ex, k-direction
                //.............................................................
                for (k = 0; k < nzPML_1; ++k) {

                    id = ID1.p[i][j][k];
                    psi_Exz_1.p[i][j][k] = be_z_1.p[k] * psi_Exz_1.p[i][j][k]
                            + ce_z_1.p[k] * (Hy.p[i][j][k] - Hy.p[i][j][k + 1]) / dz;
                    Ex.p[i][j][k] = Ex.p[i][j][k] + CB[id] * psi_Exz_1.p[i][j][k];
                }
                //..............................................................
                //  PML for top Ex, k-direction
                //..............................................................
                kk = nzPML_2 - 1;
                for (k = Kmax - nzPML_2 - 1; k < Kmax - 1; ++k) {

                    id = ID1.p[i][j][k];
                    psi_Exz_2.p[i][j][kk] = be_z_2.p[kk] * psi_Exz_2.p[i][j][kk]
                            + ce_z_2.p[kk] * (Hy.p[i][j][k] - Hy.p[i][j][k + 1]) / dz;
                    Ex.p[i][j][k] = Ex.p[i][j][k] + CB[id] * psi_Exz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    void updateHz(const int k, data3d<type1> &Hz, const data3d<type1>&Ex, const data3d<type1>&Ey, const float DB, const double dx, const double dy) {
        int i, j, ii, jj;
        for (j = 0; j < Jmax - 1; ++j) {
            //..........................................................
            //  PML for bottom Hz, x-direction
            //..........................................................
            for (i = 0; i < nxPML_1 - 1; ++i) {

                psi_Hzx_1.p[i][j][k] = bh_x_1.p[i] * psi_Hzx_1.p[i][j][k]
                        + ch_x_1.p[i] * (Ey.p[i][j][k] - Ey.p[i + 1][j][k]) / dx;
                Hz.p[i][j][k] = Hz.p[i][j][k] + DB * psi_Hzx_1.p[i][j][k];
            }
            //..........................................................
            //  PML for top Hz, x-direction
            //..........................................................
            ii = nxPML_2 - 2;
            for (i = Imax - nxPML_2; i < Imax - 1; ++i) {

                psi_Hzx_2.p[ii][j][k] = bh_x_2.p[ii] * psi_Hzx_2.p[ii][j][k]
                        + ch_x_2.p[ii] * (Ey.p[i][j][k] - Ey.p[i + 1][j][k]) / dx;
                Hz.p[i][j][k] = Hz.p[i][j][k] + DB * psi_Hzx_2.p[ii][j][k];
                ii = ii - 1;
            }
        }

        for (i = 0; i < Imax - 1; ++i) {
            //........................................................
            //  PML for bottom Hz, y-direction
            //.........................................................
            for (j = 0; j < nyPML_1 - 1; ++j) {
                psi_Hzy_1.p[i][j][k] = bh_y_1.p[j] * psi_Hzy_1.p[i][j][k]
                        + ch_y_1.p[j] * (Ex.p[i][j + 1][k] - Ex.p[i][j][k]) / dy;
                Hz.p[i][j][k] = Hz.p[i][j][k] + DB * psi_Hzy_1.p[i][j][k];

            }
            //.........................................................
            //  PML for top Hz, y-direction
            //..........................................................
            jj = nyPML_2 - 2;
            for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {

                psi_Hzy_2.p[i][jj][k] = bh_y_2.p[jj] * psi_Hzy_2.p[i][jj][k]
                        + ch_y_2.p[jj] * (Ex.p[i][j + 1][k] - Ex.p[i][j][k]) / dy;
                Hz.p[i][j][k] = Hz.p[i][j][k] + DB * psi_Hzy_2.p[i][jj][k];
                jj = jj - 1;
            }
        }
    };

    void updateHyOut(data3d<type1> &Hy, const data3d<type1> &Ex, const float DB, const double dz) {
        int i, j, k, kk;
        for (i = 0; i < Imax - 1; ++i) {
            for (j = 0; j < Jmax - 1; ++j) {
                //.......................................................
                //  PML for bottom Hy, k-direction
                //......................................................
                for (k = 1; k < nzPML_1; ++k) {

                    psi_Hyz_1.p[i][j][k - 1] = bh_z_1.p[k - 1] * psi_Hyz_1.p[i][j][k - 1]
                            + ch_z_1.p[k - 1] * (Ex.p[i][j][k - 1] - Ex.p[i][j][k]) / dz;
                    Hy.p[i][j][k] = Hy.p[i][j][k] + DB * psi_Hyz_1.p[i][j][k - 1];
                }
                //.......................................................
                //  PML for top Hy, k-direction
                //.........................................................
                kk = nzPML_2 - 2;
                for (k = Kmax - nzPML_2; k < Kmax - 1; ++k) {
                    psi_Hyz_2.p[i][j][kk] = bh_z_2.p[kk] * psi_Hyz_2.p[i][j][kk]
                            + ch_z_2.p[kk] * (Ex.p[i][j][k - 1] - Ex.p[i][j][k]) / dz;
                    Hy.p[i][j][k] = Hy.p[i][j][k] + DB * psi_Hyz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    void updateHyIn(const int k, data3d<type1> &Hy, const data3d<type1> &Ez, const float DB, const double dx) {
        int i, j, ii;
        for (j = 0; j < Jmax - 1; ++j) {
            //.......................................................
            //  PML for bottom Hy, i-direction
            //.......................................................
            for (i = 0; i < nxPML_1 - 1; ++i) {

                psi_Hyx_1.p[i][j][k] = bh_x_1.p[i] * psi_Hyx_1.p[i][j][k]
                        + ch_x_1.p[i] * (Ez.p[i + 1][j][k] - Ez.p[i][j][k]) / dx;
                Hy.p[i][j][k] = Hy.p[i][j][k] + DB * psi_Hyx_1.p[i][j][k];
            }
            //.........................................................
            //  PML for top Hy, i-direction
            //.........................................................
            ii = nxPML_2 - 2;
            for (i = Imax - nxPML_2; i < Imax - 1; ++i) {

                psi_Hyx_2.p[ii][j][k] = bh_x_2.p[ii] * psi_Hyx_2.p[ii][j][k]
                        + ch_x_2.p[ii] * (Ez.p[i + 1][j][k] - Ez.p[i][j][k]) / dx;
                Hy.p[i][j][k] = Hy.p[i][j][k] + DB * psi_Hyx_2.p[ii][j][k];
                ii = ii - 1;
            }
        }
    };

    void updateHxOut(data3d<type1> &Hx, const data3d<type1> &Ey, const float DB, const double dz) {
        int i, j, k, kk;
        for (i = 0; i < Imax - 1; ++i) {

            for (j = 0; j < Jmax - 1; ++j) {
                //....................................................
                //  PML for bottom Hx, k-direction
                //................................................
                for (k = 1; k < nzPML_1; ++k) {

                    psi_Hxz_1.p[i][j][k - 1] = bh_z_1.p[k - 1] * psi_Hxz_1.p[i][j][k - 1]
                            + ch_z_1.p[k - 1] * (Ey.p[i][j][k] - Ey.p[i][j][k - 1]) / dz;
                    Hx.p[i][j][k] = Hx.p[i][j][k] + DB * psi_Hxz_1.p[i][j][k - 1];
                }
                //....................................................
                //  PML for top Hx, k-direction
                //...............................................
                kk = nzPML_2 - 2;
                for (k = Kmax - nzPML_2; k < Kmax - 1; ++k) {

                    psi_Hxz_2.p[i][j][kk] = bh_z_2.p[kk] * psi_Hxz_2.p[i][j][kk]
                            + ch_z_2.p[kk] * (Ey.p[i][j][k] - Ey.p[i][j][k - 1]) / dz;
                    Hx.p[i][j][k] = Hx.p[i][j][k] + DB * psi_Hxz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    void updateHxIn(const int k, data3d<type1> &Hx, const data3d<type1> &Ez, const float DB, const double dy) {
        int i, j, jj;
        for (i = 0; i < Imax - 1; ++i) {
            //...............................................
            //  PML for bottom Hx, j-direction
            //...............................................
            for (j = 0; j < nyPML_1 - 1; ++j) {
                psi_Hxy_1.p[i][j][k] = bh_y_1.p[j] * psi_Hxy_1.p[i][j][k]
                        + ch_y_1.p[j] * (Ez.p[i][j][k] - Ez.p[i][j + 1][k]) / dy;
                Hx.p[i][j][k] = Hx.p[i][j][k] + DB * psi_Hxy_1.p[i][j][k];
            }
            //....................................................
            //  PML for top Hx, j-direction
            //.....................................................
            jj = nyPML_2 - 2;
            for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {
                psi_Hxy_2.p[i][jj][k] = bh_y_2.p[jj] * psi_Hxy_2.p[i][jj][k]
                        + ch_y_2.p[jj] * (Ez.p[i][j][k] - Ez.p[i][j + 1][k]) / dy;
                Hx.p[i][j][k] = Hx.p[i][j][k] + DB * psi_Hxy_2.p[i][jj][k];
                jj = jj - 1;
            }
        }
    };

    void initParmeters(const double dx, const double dy, const double dz, int m_, int ma_) {
        m = m_;
        ma = ma_;
        sig_x_max = 0.75 * (0.8 * (m + 1) / (dx * sqrt(mu_0 / (eps_0 * epsR))));
        sig_y_max = 0.75 * (0.8 * (m + 1) / (dy * sqrt(mu_0 / (eps_0 * epsR))));
        sig_z_max = 0.75 * (0.8 * (m + 1) / (dz * sqrt(mu_0 / (eps_0 * epsR))));
        alpha_x_max = 0.03;
        alpha_y_max = alpha_x_max;
        alpha_z_max = alpha_x_max;
        kappa_x_max = 8.0;
        kappa_y_max = kappa_x_max;
        kappa_z_max = kappa_x_max;
    };

    void Initial(unsigned nx, unsigned ny, unsigned nz, unsigned ncpml) {

        //PML Layers (10 layers)
        nxPML_1 = ncpml;
        nxPML_2 = ncpml;
        nyPML_1 = ncpml;
        nyPML_2 = ncpml;
        nzPML_1 = ncpml;
        nzPML_2 = ncpml;

        //domain size
        Imax = nx;
        Jmax = ny;
        Kmax = nz;
    };

    void InitialMuEps() {
        mu_0 = 4.0 * pi * 1.0E-7;
        eps_0 = 1.0 / (C * C * mu_0);
    };

    void createCPMLArray() {
        createPsi();
        createDen();
        createCBKAP();
    };

    void initCPML(const double dt, const double dx, const double dy, const double dz) {


        initCBKAP(dt, dx, dy, dz);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  DENOMINATORS FOR FIELD UPDATES
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        initDen(dt, dx, dy, dz);

    };

private:

    // Specify the CPML Thickness in Each Direction (Value of Zero
    // Corresponds to No PML, and the Grid is Terminated with a PEC)
    // PML thickness in each direction
    int nxPML_1, nxPML_2, nyPML_1;
    int nyPML_2, nzPML_1, nzPML_2;

    // grid size corresponding to the number of Ez field components
    int Imax;
    int Jmax;
    int Kmax;

    // Specify the CPML Order and Other Parameters:
    int m, ma;

    double sig_x_max;
    double sig_y_max;
    double sig_z_max;
    double alpha_x_max;
    double alpha_y_max;
    double alpha_z_max;
    double kappa_x_max;
    double kappa_y_max;
    double kappa_z_max;

    void initCBKAP(const double dt, const double dx, const double dy, const double dz) {
        int i, j, k;
        for (i = 0; i < nxPML_1; ++i) {

            sige_x_PML_1.p[i] = sig_x_max * pow(((nxPML_1 - 1 - i)
                    / (nxPML_1 - 1.0)), m);
            alphae_x_PML_1.p[i] = alpha_x_max * pow(((float) i
                    / (nxPML_1 - 1.0)), ma);
            kappae_x_PML_1.p[i] = 1.0 + (kappa_x_max - 1.0) *
                    pow((nxPML_1 - 1 - i) / (nxPML_1 - 1.0), m);
            be_x_1.p[i] = exp(-(sige_x_PML_1.p[i] / kappae_x_PML_1.p[i] +
                    alphae_x_PML_1.p[i]) * dt / eps_0);

            if ((sige_x_PML_1.p[i] == 0.0) &&
                    (alphae_x_PML_1.p[i] == 0.0) && (i == nxPML_1 - 1)) {
                ce_x_1.p[i] = 0.0;

            } else {
                ce_x_1.p[i] = sige_x_PML_1.p[i] * (be_x_1.p[i] - 1.0) /
                        (sige_x_PML_1.p[i] + kappae_x_PML_1.p[i] * alphae_x_PML_1.p[i])
                        / kappae_x_PML_1.p[i];
            }
        }

        for (i = 0; i < nxPML_1 - 1; ++i) {

            sigh_x_PML_1.p[i] = sig_x_max * pow(((nxPML_1 - 1 - i - 0.5)
                    / (nxPML_1 - 1.0)), m);
            alphah_x_PML_1.p[i] = alpha_x_max * pow(((i + 1 - 0.5)
                    / (nxPML_1 - 1.0)), ma);
            kappah_x_PML_1.p[i] = 1.0 + (kappa_x_max - 1.0) *
                    pow(((nxPML_1 - 1 - i - 0.5) / (nxPML_1 - 1.0)), m);
            bh_x_1.p[i] = exp(-(sigh_x_PML_1.p[i] / kappah_x_PML_1.p[i] +
                    alphah_x_PML_1.p[i]) * dt / eps_0);
            ch_x_1.p[i] = sigh_x_PML_1.p[i] * (bh_x_1.p[i] - 1.0) /
                    (sigh_x_PML_1.p[i] + kappah_x_PML_1.p[i] * alphah_x_PML_1.p[i])
                    / kappah_x_PML_1.p[i];
        }

        for (i = 0; i < nxPML_2; ++i) {

            sige_x_PML_2.p[i] = sig_x_max * pow(((nxPML_2 - 1 - i)
                    / (nxPML_2 - 1.0)), m);
            alphae_x_PML_2.p[i] = alpha_x_max * pow(((float) i
                    / (nxPML_2 - 1.0)), ma);
            kappae_x_PML_2.p[i] = 1.0 + (kappa_x_max - 1.0) *
                    pow((nxPML_2 - 1 - i) / (nxPML_2 - 1.0), m);
            be_x_2.p[i] = exp(-(sige_x_PML_2.p[i] / kappae_x_PML_2.p[i] +
                    alphae_x_PML_2.p[i]) * dt / eps_0);

            if ((sige_x_PML_2.p[i] == 0.0) &&
                    (alphae_x_PML_2.p[i] == 0.0) && (i == nxPML_2 - 1)) {

                ce_x_2.p[i] = 0.0;

            } else {

                ce_x_2.p[i] = sige_x_PML_2.p[i] * (be_x_2.p[i] - 1.0) /
                        (sige_x_PML_2.p[i] + kappae_x_PML_2.p[i] * alphae_x_PML_2.p[i])
                        / kappae_x_PML_2.p[i];
            }

        }

        for (i = 0; i < nxPML_2 - 1; ++i) {

            sigh_x_PML_2.p[i] = sig_x_max * pow(((nxPML_2 - 1 - i - 0.5)
                    / (nxPML_2 - 1.0)), m);
            alphah_x_PML_2.p[i] = alpha_x_max * pow(((i + 1 - 0.5)
                    / (nxPML_2 - 1.0)), ma);
            kappah_x_PML_2.p[i] = 1.0 + (kappa_x_max - 1.0) *
                    pow(((nxPML_2 - 1 - i - 0.5) / (nxPML_2 - 1.0)), m);
            bh_x_2.p[i] = exp(-(sigh_x_PML_2.p[i] / kappah_x_PML_2.p[i] +
                    alphah_x_PML_2.p[i]) * dt / eps_0);
            ch_x_2.p[i] = sigh_x_PML_2.p[i] * (bh_x_2.p[i] - 1.0) /
                    (sigh_x_PML_2.p[i] + kappah_x_PML_2.p[i] * alphah_x_PML_2.p[i])
                    / kappah_x_PML_2.p[i];
        }

        for (j = 0; j < nyPML_1; ++j) {

            sige_y_PML_1.p[j] = sig_y_max * pow(((nyPML_1 - 1 - j)
                    / (nyPML_1 - 1.0)), m);
            alphae_y_PML_1.p[j] = alpha_y_max * pow(((float) j
                    / (nyPML_1 - 1.0)), ma);
            kappae_y_PML_1.p[j] = 1.0 + (kappa_y_max - 1.0) *
                    pow((nyPML_1 - 1 - j) / (nyPML_1 - 1.0), m);
            be_y_1.p[j] = exp(-(sige_y_PML_1.p[j] / kappae_y_PML_1.p[j] +
                    alphae_y_PML_1.p[j]) * dt / eps_0);

            if ((sige_y_PML_1.p[j] == 0.0) &&
                    (alphae_y_PML_1.p[j] == 0.0) && (j == nyPML_1 - 1)) {
                ce_y_1.p[j] = 0.0;

            } else {
                ce_y_1.p[j] = sige_y_PML_1.p[j] * (be_y_1.p[j] - 1.0) /
                        (sige_y_PML_1.p[j] + kappae_y_PML_1.p[j] * alphae_y_PML_1.p[j])
                        / kappae_y_PML_1.p[j];
            }
        }

        for (j = 0; j < nyPML_1 - 1; ++j) {

            sigh_y_PML_1.p[j] = sig_y_max * pow(((nyPML_1 - 1 - j - 0.5)
                    / (nyPML_1 - 1.0)), m);
            alphah_y_PML_1.p[j] = alpha_y_max * pow(((j + 1 - 0.5)
                    / (nyPML_1 - 1.0)), ma);
            kappah_y_PML_1.p[j] = 1.0 + (kappa_y_max - 1.0) *
                    pow(((nyPML_1 - 1 - j - 0.5) / (nyPML_1 - 1.0)), m);
            bh_y_1.p[j] = exp(-(sigh_y_PML_1.p[j] / kappah_y_PML_1.p[j] +
                    alphah_y_PML_1.p[j]) * dt / eps_0);
            ch_y_1.p[j] = sigh_y_PML_1.p[j] * (bh_y_1.p[j] - 1.0) /
                    (sigh_y_PML_1.p[j] + kappah_y_PML_1.p[j] * alphah_y_PML_1.p[j])
                    / kappah_y_PML_1.p[j];
        }

        for (j = 0; j < nyPML_2; ++j) {

            sige_y_PML_2.p[j] = sig_y_max * pow(((nyPML_2 - 1 - j)
                    / (nyPML_2 - 1.0)), m);
            alphae_y_PML_2.p[j] = alpha_y_max * pow(((float) j
                    / (nyPML_2 - 1.0)), ma);
            kappae_y_PML_2.p[j] = 1.0 + (kappa_y_max - 1.0) *
                    pow((nyPML_2 - 1 - j) / (nyPML_2 - 1.0), m);
            be_y_2.p[j] = exp(-(sige_y_PML_2.p[j] / kappae_y_PML_2.p[j] +
                    alphae_y_PML_2.p[j]) * dt / eps_0);

            if ((sige_y_PML_2.p[j] == 0.0) &&
                    (alphae_y_PML_2.p[j] == 0.0) && (j == nyPML_2 - 1)) {

                ce_y_2.p[j] = 0.0;

            } else {

                ce_y_2.p[j] = sige_y_PML_2.p[j] * (be_y_2.p[j] - 1.0) /
                        (sige_y_PML_2.p[j] + kappae_y_PML_2.p[j] * alphae_y_PML_2.p[j])
                        / kappae_y_PML_2.p[j];
            }
        }

        for (j = 0; j < nyPML_2 - 1; ++j) {

            sigh_y_PML_2.p[j] = sig_y_max * pow(((nyPML_2 - 1 - j - 0.5)
                    / (nyPML_2 - 1.0)), m);
            alphah_y_PML_2.p[j] = alpha_y_max * pow(((j + 1 - 0.5)
                    / (nyPML_2 - 1.0)), ma);
            kappah_y_PML_2.p[j] = 1.0 + (kappa_y_max - 1.0) *
                    pow(((nyPML_2 - 1 - j - 0.5) / (nyPML_2 - 1.0)), m);
            bh_y_2.p[j] = exp(-(sigh_y_PML_2.p[j] / kappah_y_PML_2.p[j] +
                    alphah_y_PML_2.p[j]) * dt / eps_0);
            ch_y_2.p[j] = sigh_y_PML_2.p[j] * (bh_y_2.p[j] - 1.0) /
                    (sigh_y_PML_2.p[j] + kappah_y_PML_2.p[j] * alphah_y_PML_2.p[j])
                    / kappah_y_PML_2.p[j];
        }

        for (k = 0; k < nzPML_1; ++k) {

            sige_z_PML_1.p[k] = sig_z_max * pow(((nzPML_1 - 1 - k)
                    / (nzPML_1 - 1.0)), m);
            alphae_z_PML_1.p[k] = alpha_z_max * pow(((float) k
                    / (nzPML_1 - 1.0)), ma);
            kappae_z_PML_1.p[k] = 1.0 + (kappa_z_max - 1.0) *
                    pow((nzPML_1 - 1 - k) / (nzPML_1 - 1.0), m);
            be_z_1.p[k] = exp(-(sige_z_PML_1.p[k] / kappae_z_PML_1.p[k] +
                    alphae_z_PML_1.p[k]) * dt / eps_0);

            if ((sige_z_PML_1.p[k] == 0.0) &&
                    (alphae_z_PML_1.p[k] == 0.0) && (k == nzPML_1 - 1)) {

                ce_z_1.p[k] = 0.0;

            } else {

                ce_z_1.p[k] = sige_z_PML_1.p[k] * (be_z_1.p[k] - 1.0) /
                        (sige_z_PML_1.p[k] + kappae_z_PML_1.p[k] * alphae_z_PML_1.p[k])
                        / kappae_z_PML_1.p[k];
            }
        }

        for (k = 0; k < nzPML_1 - 1; ++k) {

            sigh_z_PML_1.p[k] = sig_z_max * pow(((nzPML_1 - 1 - k - 0.5)
                    / (nzPML_1 - 1.0)), m);
            alphah_z_PML_1.p[k] = alpha_z_max * pow(((k + 1 - 0.5)
                    / (nzPML_1 - 1.0)), ma);
            kappah_z_PML_1.p[k] = 1.0 + (kappa_z_max - 1.0) *
                    pow(((nzPML_1 - 1 - k - 0.5) / (nzPML_1 - 1.0)), m);
            bh_z_1.p[k] = exp(-(sigh_z_PML_1.p[k] / kappah_z_PML_1.p[k] +
                    alphah_z_PML_1.p[k]) * dt / eps_0);
            ch_z_1.p[k] = sigh_z_PML_1.p[k] * (bh_z_1.p[k] - 1.0) /
                    (sigh_z_PML_1.p[k] + kappah_z_PML_1.p[k] * alphah_z_PML_1.p[k])
                    / kappah_z_PML_1.p[k];
        }

        for (k = 0; k < nzPML_2; ++k) {

            sige_z_PML_2.p[k] = sig_z_max * pow(((nzPML_2 - 1 - k)
                    / (nzPML_2 - 1.0)), m);
            alphae_z_PML_2.p[k] = alpha_z_max * pow(((float) k
                    / (nzPML_2 - 1.0)), ma);
            kappae_z_PML_2.p[k] = 1.0 + (kappa_z_max - 1.0) *
                    pow((nzPML_2 - 1 - k) / (nzPML_2 - 1.0), m);
            be_z_2.p[k] = exp(-(sige_z_PML_2.p[k] / kappae_z_PML_2.p[k] +
                    alphae_z_PML_2.p[k]) * dt / eps_0);

            if ((sige_z_PML_2.p[k] == 0.0) &&
                    (alphae_z_PML_2.p[k] == 0.0) && (k == nzPML_2 - 1)) {

                ce_z_2.p[k] = 0.0;

            } else {

                ce_z_2.p[k] = sige_z_PML_2.p[k] * (be_z_2.p[k] - 1.0) /
                        (sige_z_PML_2.p[k] + kappae_z_PML_2.p[k] * alphae_z_PML_2.p[k])
                        / kappae_z_PML_2.p[k];
            }
        }

        for (k = 0; k < nzPML_2 - 1; ++k) {

            sigh_z_PML_2.p[k] = sig_z_max * pow(((nzPML_2 - 1 - k - 0.5)
                    / (nzPML_2 - 1.0)), m);
            alphah_z_PML_2.p[k] = alpha_z_max * pow(((k + 1 - 0.5)
                    / (nzPML_2 - 1.0)), ma);
            kappah_z_PML_2.p[k] = 1.0 + (kappa_z_max - 1.0) *
                    pow(((nzPML_2 - 1 - k - 0.5) / (nzPML_2 - 1.0)), m);
            bh_z_2.p[k] = exp(-(sigh_z_PML_2.p[k] / kappah_z_PML_2.p[k] +
                    alphah_z_PML_2.p[k]) * dt / eps_0);
            ch_z_2.p[k] = sigh_z_PML_2.p[k] * (bh_z_2.p[k] - 1.0) /
                    (sigh_z_PML_2.p[k] + kappah_z_PML_2.p[k] * alphah_z_PML_2.p[k])
                    / kappah_z_PML_2.p[k];
        }
    };

    void initDen(const double dt, const double dx, const double dy, const double dz) {
        int i, j, k, ii, jj, kk;
        ii = nxPML_2 - 2;

        for (i = 0; i < Imax - 1; ++i) {

            if (i < nxPML_1 - 1) {

                den_hx.p[i] = 1.0 / (kappah_x_PML_1.p[i] * dx);

            } else if (i >= Imax - nxPML_2) {

                den_hx.p[i] = 1.0 / (kappah_x_PML_2.p[ii] * dx);
                ii = ii - 1;

            } else {

                den_hx.p[i] = 1.0 / dx;
            }
        }

        jj = nyPML_2 - 2;

        for (j = 0; j < Jmax - 1; ++j) {

            if (j < nyPML_1 - 1) {

                den_hy.p[j] = 1.0 / (kappah_y_PML_1.p[j] * dy);

            } else if (j >= Jmax - nyPML_2) {

                den_hy.p[j] = 1.0 / (kappah_y_PML_2.p[jj] * dy);
                jj = jj - 1;

            } else {

                den_hy.p[j] = 1.0 / dy;
            }
        }

        kk = nzPML_2 - 2;

        for (k = 1; k < Kmax - 1; ++k) {

            if (k < nzPML_1) {

                den_hz.p[k] = 1.0 / (kappah_z_PML_1.p[k - 1] * dz);

            } else if (k >= Kmax - nzPML_2) {

                den_hz.p[k] = 1.0 / (kappah_z_PML_2.p[kk] * dz);
                kk = kk - 1;

            } else {

                den_hz.p[k] = 1.0 / dz;
            }
        }

        ii = nxPML_2 - 1;

        for (i = 0; i < Imax - 1; ++i) {

            if (i < nxPML_1) {

                den_ex.p[i] = 1.0 / (kappae_x_PML_1.p[i] * dx);

            } else if (i >= Imax - nxPML_2) {

                den_ex.p[i] = 1.0 / (kappae_x_PML_2.p[ii] * dx);
                ii = ii - 1;

            } else {

                den_ex.p[i] = 1.0 / dx;
            }
        }

        jj = nyPML_2 - 1;

        for (j = 0; j < Jmax - 1; ++j) {

            if (j < nyPML_1) {

                den_ey.p[j] = 1.0 / (kappae_y_PML_1.p[j] * dy);

            } else if (j >= Jmax - nyPML_2) {

                den_ey.p[j] = 1.0 / (kappae_y_PML_2.p[jj] * dy);
                jj = jj - 1;

            } else {

                den_ey.p[j] = 1.0 / dy;
            }
        }

        kk = nzPML_2 - 1;

        for (k = 0; k < Kmax - 1; ++k) {

            if (k < nzPML_1) {

                den_ez.p[k] = 1.0 / (kappae_z_PML_1.p[k] * dz);

            } else if (k >= Kmax - 1 - nzPML_2) {

                den_ez.p[k] = 1.0 / (kappae_z_PML_2.p[kk] * dz);
                kk = kk - 1;

            } else {

                den_ez.p[k] = 1.0 / dz;
            }
        }
    };

    void createDen() {
        den_ex.createArray(Imax - 1, 0.0);
        den_hx.createArray(Imax - 1, 0.0);
        den_ey.createArray(Jmax - 1, 0.0);
        den_hy.createArray(Jmax - 1, 0.0);
        den_ez.createArray(Kmax - 1, 0.0);
        den_hz.createArray(Kmax - 1, 0.0);
    };

    void createCBKAP() {

        be_x_1.createArray(nxPML_1, 0.0);
        ce_x_1.createArray(nxPML_1, 0.0);
        alphae_x_PML_1.createArray(nxPML_1, 0.0);
        sige_x_PML_1.createArray(nxPML_1, 0.0);
        kappae_x_PML_1.createArray(nxPML_1, 0.0);
        bh_x_1.createArray(nxPML_1 - 1, 0.0);
        ch_x_1.createArray(nxPML_1 - 1, 0.0);
        alphah_x_PML_1.createArray(nxPML_1 - 1, 0.0);
        sigh_x_PML_1.createArray(nxPML_1 - 1, 0.0);
        kappah_x_PML_1.createArray(nxPML_1 - 1, 0.0);

        be_x_2.createArray(nxPML_2, 0.0);
        ce_x_2.createArray(nxPML_2, 0.0);
        alphae_x_PML_2.createArray(nxPML_2, 0.0);
        sige_x_PML_2.createArray(nxPML_2, 0.0);
        kappae_x_PML_2.createArray(nxPML_2, 0.0);

        bh_x_2.createArray(nxPML_2 - 1, 0.0);
        ch_x_2.createArray(nxPML_2 - 1, 0.0);
        alphah_x_PML_2.createArray(nxPML_2 - 1, 0.0);
        sigh_x_PML_2.createArray(nxPML_2 - 1, 0.0);
        kappah_x_PML_2.createArray(nxPML_2 - 1, 0.0);

        be_y_1.createArray(nxPML_1, 0.0);
        ce_y_1.createArray(nxPML_1, 0.0);
        kappae_y_PML_1.createArray(nxPML_1, 0.0);
        sige_y_PML_1.createArray(nxPML_1, 0.0);
        alphae_y_PML_1.createArray(nxPML_1, 0.0);

        bh_y_1.createArray(nyPML_1 - 1, 0.0);
        ch_y_1.createArray(nyPML_1 - 1, 0.0);
        alphah_y_PML_1.createArray(nyPML_1 - 1, 0.0);
        sigh_y_PML_1.createArray(nyPML_1 - 1, 0.0);
        kappah_y_PML_1.createArray(nyPML_1 - 1, 0.0);

        be_y_2.createArray(nyPML_2, 0.0);
        ce_y_2.createArray(nyPML_2, 0.0);
        alphae_y_PML_2.createArray(nyPML_2, 0.0);
        sige_y_PML_2.createArray(nyPML_2, 0.0);
        kappae_y_PML_2.createArray(nyPML_2, 0.0);

        bh_y_2.createArray(nyPML_2 - 1, 0.0);
        ch_y_2.createArray(nyPML_2 - 1, 0.0);
        alphah_y_PML_2.createArray(nyPML_2 - 1, 0.0);
        sigh_y_PML_2.createArray(nyPML_2 - 1, 0.0);
        kappah_y_PML_2.createArray(nyPML_2 - 1, 0.0);

        be_z_1.createArray(nzPML_1, 0.0);
        ce_z_1.createArray(nzPML_1, 0.0);
        alphae_z_PML_1.createArray(nzPML_1, 0.0);
        sige_z_PML_1.createArray(nzPML_1, 0.0);
        kappae_z_PML_1.createArray(nzPML_1, 0.0);

        bh_z_1.createArray(nzPML_1 - 1, 0.0);
        ch_z_1.createArray(nzPML_1 - 1, 0.0);
        alphah_z_PML_1.createArray(nzPML_1 - 1, 0.0);
        sigh_z_PML_1.createArray(nzPML_1 - 1, 0.0);
        kappah_z_PML_1.createArray(nzPML_1 - 1, 0.0);

        be_z_2.createArray(nzPML_2, 0.0);
        ce_z_2.createArray(nzPML_2, 0.0);
        alphae_z_PML_2.createArray(nzPML_2, 0.0);
        sige_z_PML_2.createArray(nzPML_2, 0.0);
        kappae_z_PML_2.createArray(nzPML_2, 0.0);

        bh_z_2.createArray(nzPML_2 - 1, 0.0);
        ch_z_2.createArray(nzPML_2 - 1, 0.0);
        alphah_z_PML_2.createArray(nzPML_2 - 1, 0.0);
        sigh_z_PML_2.createArray(nzPML_2 - 1, 0.0);
        kappah_z_PML_2.createArray(nzPML_2 - 1, 0.0);

    };

    void createPsi() {
        psi_Ezx_1.CreateStruct(nxPML_1, Jmax, Kmax, 0);
        psi_Ezx_2.CreateStruct(nxPML_2, Jmax, Kmax, 0);
        psi_Hyx_1.CreateStruct(nxPML_1 - 1, Jmax, Kmax, 0);
        psi_Hyx_2.CreateStruct(nxPML_2 - 1, Jmax, Kmax, 0);

        psi_Ezy_1.CreateStruct(Imax, nyPML_1, Kmax, 0);
        psi_Ezy_2.CreateStruct(Imax, nyPML_2, Kmax, 0);
        psi_Hxy_1.CreateStruct(Imax, nyPML_1 - 1, Kmax, 0);
        psi_Hxy_2.CreateStruct(Imax, nyPML_2 - 1, Kmax, 0);

        psi_Hxz_1.CreateStruct(Imax, Jmax - 1, nzPML_1 - 1, 0);
        psi_Hxz_2.CreateStruct(Imax, Jmax - 1, nzPML_2 - 1, 0);
        psi_Eyz_1.CreateStruct(Imax, Jmax - 1, nzPML_1, 0);
        psi_Eyz_2.CreateStruct(Imax, Jmax - 1, nzPML_2, 0);

        psi_Hyz_1.CreateStruct(Imax - 1, Jmax, nzPML_1 - 1, 0);
        psi_Hyz_2.CreateStruct(Imax - 1, Jmax, nzPML_2 - 1, 0);
        psi_Exz_1.CreateStruct(Imax - 1, Jmax, nzPML_1, 0);
        psi_Exz_2.CreateStruct(Imax - 1, Jmax, nzPML_2, 0);

        psi_Hzx_1.CreateStruct(nxPML_1 - 1, Jmax - 1, Kmax - 1, 0);
        psi_Hzx_2.CreateStruct(nxPML_2 - 1, Jmax - 1, Kmax - 1, 0);
        psi_Eyx_1.CreateStruct(nxPML_1, Jmax - 1, Kmax - 1, 0);
        psi_Eyx_2.CreateStruct(nxPML_2, Jmax - 1, Kmax - 1, 0);

        psi_Hzy_1.CreateStruct(Imax - 1, nyPML_1 - 1, Kmax - 1, 0);
        psi_Hzy_2.CreateStruct(Imax - 1, nyPML_2 - 1, Kmax - 1, 0);
        psi_Exy_1.CreateStruct(Imax - 1, nyPML_1, Kmax - 1, 0);
        psi_Exy_2.CreateStruct(Imax - 1, nyPML_2, Kmax - 1, 0);
    };
};
template<class type1, class type2>
double cpml<type1, type2>::eps_0 = 0;
template<class type1, class type2>
double cpml<type1, type2>::mu_0 = 0;
template<class type1, class type2>
double cpml<type1, type2>::epsR = 1.0;
template<class type1, class type2>
const double cpml<type1, type2>::C = 2.99792458E8;
template<class type1, class type2>
const double cpml<type1, type2>::pi = 3.14159265358979;
#endif	/* CPML_H */

