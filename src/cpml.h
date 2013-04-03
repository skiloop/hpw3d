#pragma once
#include "datastruct.h"
//#include "data1d.h"

extern MyDataF eps_0, epsR;
extern MyDataF mu_0;
extern MyDataF dt, dx, dy, dz;

class cpml {
public:

    cpml(int _m, int _ma) : m(_m), ma(_ma) {
    };
    ~cpml(void);
private:
    // grid size corresponding to the number of Ez field components
    unsigned Imax;
    unsigned Jmax;
    unsigned Kmax;
    //  Specify the CPML Thickness in Each Direction (Value of Zero
    //  Corresponds to No PML, and the Grid is Terminated with a PEC)
    // PML thickness in each direction
    unsigned nxPML_1, nxPML_2, nyPML_1;
    unsigned nyPML_2, nzPML_1, nzPML_2;

    ////  Specify the dipole Boundaries(A cuboidal rode- NOT as a cylinder)
    //unsigned istart, iend, jstart;
    //unsigned jend, kstart, kend;

    //  Specify the CPML Order and Other Parameters:
    int m, ma;

    MyDataF sig_x_max;
    MyDataF sig_y_max;
    MyDataF sig_z_max;
    MyDataF alpha_x_max;
    MyDataF alpha_y_max;
    MyDataF alpha_z_max;
    MyDataF kappa_x_max;
    MyDataF kappa_y_max;
    MyDataF kappa_z_max;

    //  CPML components (Taflove 3rd Edition, Chapter 7)
    data3d<MyDataF> psi_Ezx_1;
    data3d<MyDataF> psi_Ezx_2;
    data3d<MyDataF> psi_Hyx_1;
    data3d<MyDataF> psi_Hyx_2;
    data3d<MyDataF> psi_Ezy_1;
    data3d<MyDataF> psi_Ezy_2;
    data3d<MyDataF> psi_Hxy_1;
    data3d<MyDataF> psi_Hxy_2;
    data3d<MyDataF> psi_Hxz_1;
    data3d<MyDataF> psi_Hxz_2;
    data3d<MyDataF> psi_Hyz_1;
    data3d<MyDataF> psi_Hyz_2;
    data3d<MyDataF> psi_Exz_1;
    data3d<MyDataF> psi_Exz_2;
    data3d<MyDataF> psi_Eyz_1;
    data3d<MyDataF> psi_Eyz_2;
    data3d<MyDataF> psi_Hzx_1;
    data3d<MyDataF> psi_Eyx_1;
    data3d<MyDataF> psi_Hzx_2;
    data3d<MyDataF> psi_Eyx_2;
    data3d<MyDataF> psi_Hzy_1;
    data3d<MyDataF> psi_Exy_1;
    data3d<MyDataF> psi_Hzy_2;
    data3d<MyDataF> psi_Exy_2;

    data1d<MyDataF> be_x_1, ce_x_1, alphae_x_PML_1, sige_x_PML_1, kappae_x_PML_1;
    data1d<MyDataF> bh_x_1, ch_x_1, alphah_x_PML_1, sigh_x_PML_1, kappah_x_PML_1;
    data1d<MyDataF> be_x_2, ce_x_2, alphae_x_PML_2, sige_x_PML_2, kappae_x_PML_2;
    data1d<MyDataF> bh_x_2, ch_x_2, alphah_x_PML_2, sigh_x_PML_2, kappah_x_PML_2;
    data1d<MyDataF> be_y_1, ce_y_1, alphae_y_PML_1, sige_y_PML_1, kappae_y_PML_1;
    data1d<MyDataF> bh_y_1, ch_y_1, alphah_y_PML_1, sigh_y_PML_1, kappah_y_PML_1;
    data1d<MyDataF> be_y_2, ce_y_2, alphae_y_PML_2, sige_y_PML_2, kappae_y_PML_2;
    data1d<MyDataF> bh_y_2, ch_y_2, alphah_y_PML_2, sigh_y_PML_2, kappah_y_PML_2;
    data1d<MyDataF> be_z_1, ce_z_1, alphae_z_PML_1, sige_z_PML_1, kappae_z_PML_1;
    data1d<MyDataF> bh_z_1, ch_z_1, alphah_z_PML_1, sigh_z_PML_1, kappah_z_PML_1;
    data1d<MyDataF> be_z_2, ce_z_2, alphae_z_PML_2, sige_z_PML_2, kappae_z_PML_2;
    data1d<MyDataF> bh_z_2, ch_z_2, alphah_z_PML_2, sigh_z_PML_2, kappah_z_PML_2;
public:
    // denominators for the update equations
    data1d<MyDataF> den_ex;
    data1d<MyDataF> den_hx;
    data1d<MyDataF> den_ey;
    data1d<MyDataF> den_hy;
    data1d<MyDataF> den_ez;
    data1d<MyDataF> den_hz;

public:
    int Initial(unsigned nx, unsigned ny, unsigned nz, unsigned ncpml = 11);

    template<typename Type1, typename Type2>
    void UpdatePMLForHx(data3d<Type1> &hx, const data3d<Type1> &ey, const data3d<Type1> &ez, Type2 DB) {
        unsigned i, j, k;
        unsigned jj, kk;
        for (k = 1; k < Kmax - 1; ++k) {
            for (i = 0; i < Imax - 1; ++i) {
                //...............................................
                //  PML for bottom hx.p, j-direction
                //...............................................
                for (j = 0; j < nyPML_1 - 1; ++j) {
                    psi_Hxy_1.p[i][j][k] = bh_y_1.p[j] * psi_Hxy_1.p[i][j][k]
                            + ch_y_1.p[j] * (ez.p[i][j][k] - ez.p[i][j + 1][k]) / dy;
                    hx.p[i][j][k] = hx.p[i][j][k] + DB * psi_Hxy_1.p[i][j][k];
                }
                //....................................................
                //  PML for top hx.p, j-direction
                //.....................................................
                jj = nyPML_2 - 2;
                for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {
                    psi_Hxy_2.p[i][jj][k] = bh_y_2.p[jj] * psi_Hxy_2.p[i][jj][k]
                            + ch_y_2.p[jj] * (ez.p[i][j][k] - ez.p[i][j + 1][k]) / dy;
                    hx.p[i][j][k] = hx.p[i][j][k] + DB * psi_Hxy_2.p[i][jj][k];
                    jj = jj - 1;
                }
            }
        }

        for (i = 0; i < Imax - 1; ++i) {
            for (j = 0; j < Jmax - 1; ++j) {
                //....................................................
                //  PML for bottom hx.p, k-direction
                //................................................
                for (k = 1; k < nzPML_1; ++k) {
                    psi_Hxz_1.p[i][j][k - 1] = bh_z_1.p[k - 1] * psi_Hxz_1.p[i][j][k - 1]
                            + ch_z_1.p[k - 1] * (ey.p[i][j][k] - ey.p[i][j][k - 1]) / dz;
                    hx.p[i][j][k] = hx.p[i][j][k] + DB * psi_Hxz_1.p[i][j][k - 1];
                }
                //....................................................
                //  PML for top hx.p, k-direction
                //...............................................
                kk = nzPML_2 - 2;
                for (k = Kmax - nzPML_2; k < Kmax - 1; ++k) {
                    psi_Hxz_2.p[i][j][kk] = bh_z_2.p[kk] * psi_Hxz_2.p[i][j][kk]
                            + ch_z_2.p[kk] * (ey.p[i][j][k] - ey.p[i][j][k - 1]) / dz;
                    hx.p[i][j][k] = hx.p[i][j][k] + DB * psi_Hxz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForHy(data3d<Type1> &hy, const data3d<Type1> &ez, const data3d<Type1> &ex, Type2 DB) {
        unsigned i, j, k;
        unsigned ii, kk;
        for (k = 1; k < Kmax - 1; ++k) {
            for (j = 0; j < Jmax - 1; ++j) {
                //.......................................................
                //  PML for bottom hy.p, i-direction
                //.......................................................
                for (i = 0; i < nxPML_1 - 1; ++i) {
                    psi_Hyx_1.p[i][j][k] = bh_x_1.p[i] * psi_Hyx_1.p[i][j][k]
                            + ch_x_1.p[i] * (ez.p[i + 1][j][k] - ez.p[i][j][k]) / dx;
                    hy.p[i][j][k] = hy.p[i][j][k] + DB * psi_Hyx_1.p[i][j][k];
                }
                //.........................................................
                //  PML for top hy.p, i-direction
                //.........................................................
                ii = nxPML_2 - 2;
                for (i = Imax - nxPML_2; i < Imax - 1; ++i) {
                    psi_Hyx_2.p[ii][j][k] = bh_x_2.p[ii] * psi_Hyx_2.p[ii][j][k]
                            + ch_x_2.p[ii] * (ez.p[i + 1][j][k] - ez.p[i][j][k]) / dx;
                    hy.p[i][j][k] = hy.p[i][j][k] + DB * psi_Hyx_2.p[ii][j][k];
                    ii = ii - 1;
                }
            }
        }

        for (i = 0; i < Imax - 1; ++i) {

            for (j = 0; j < Jmax - 1; ++j) {
                //.......................................................
                //  PML for bottom hy.p, k-direction
                //......................................................
                for (k = 1; k < nzPML_1; ++k) {
                    psi_Hyz_1.p[i][j][k - 1] = bh_z_1.p[k - 1] * psi_Hyz_1.p[i][j][k - 1]
                            + ch_z_1.p[k - 1] * (ex.p[i][j][k - 1] - ex.p[i][j][k]) / dz;
                    hy.p[i][j][k] = hy.p[i][j][k] + DB * psi_Hyz_1.p[i][j][k - 1];
                }
                //.......................................................
                //  PML for top hy.p, k-direction
                //.........................................................
                kk = nzPML_2 - 2;
                for (k = Kmax - nzPML_2; k < Kmax - 1; ++k) {
                    psi_Hyz_2.p[i][j][kk] = bh_z_2.p[kk] * psi_Hyz_2.p[i][j][kk]
                            + ch_z_2.p[kk] * (ex.p[i][j][k - 1] - ex.p[i][j][k]) / dz;
                    hy.p[i][j][k] = hy.p[i][j][k] + DB * psi_Hyz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForHz(data3d<Type1> &hz, const data3d<Type1> &ex, const data3d<Type1> &ey, Type2 DB) {
        unsigned i, j, k;
        unsigned ii, jj;
        for (k = 0; k < Kmax - 1; ++k) {
            for (j = 0; j < Jmax - 1; ++j) {
                //..........................................................
                //  PML for bottom hz.p, x-direction
                //..........................................................
                for (i = 0; i < nxPML_1 - 1; ++i) {
                    psi_Hzx_1.p[i][j][k] = bh_x_1.p[i] * psi_Hzx_1.p[i][j][k]
                            + ch_x_1.p[i] * (ey.p[i][j][k] - ey.p[i + 1][j][k]) / dx;
                    hz.p[i][j][k] = hz.p[i][j][k] + DB * psi_Hzx_1.p[i][j][k];
                }
                //..........................................................
                //  PML for top hz.p, x-direction
                //..........................................................
                ii = nxPML_2 - 2;
                for (i = Imax - nxPML_2; i < Imax - 1; ++i) {
                    psi_Hzx_2.p[ii][j][k] = bh_x_2.p[ii] * psi_Hzx_2.p[ii][j][k]
                            + ch_x_2.p[ii] * (ey.p[i][j][k] - ey.p[i + 1][j][k]) / dx;
                    hz.p[i][j][k] = hz.p[i][j][k] + DB * psi_Hzx_2.p[ii][j][k];
                    ii = ii - 1;
                }
            }

            for (i = 0; i < Imax - 1; ++i) {
                //........................................................
                //  PML for bottom hz.p, y-direction
                //.........................................................
                for (j = 0; j < nyPML_1 - 1; ++j) {
                    psi_Hzy_1.p[i][j][k] = bh_y_1.p[j] * psi_Hzy_1.p[i][j][k]
                            + ch_y_1.p[j] * (ex.p[i][j + 1][k] - ex.p[i][j][k]) / dy;
                    hz.p[i][j][k] = hz.p[i][j][k] + DB * psi_Hzy_1.p[i][j][k];

                }
                //.........................................................
                //  PML for top hz.p, y-direction
                //..........................................................
                jj = nyPML_2 - 2;
                for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {

                    psi_Hzy_2.p[i][jj][k] = bh_y_2.p[jj] * psi_Hzy_2.p[i][jj][k]
                            + ch_y_2.p[jj] * (ex.p[i][j + 1][k] - ex.p[i][j][k]) / dy;
                    hz.p[i][j][k] = hz.p[i][j][k] + DB * psi_Hzy_2.p[i][jj][k];
                    jj = jj - 1;
                }
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEx(data3d<Type1> &ex, const data3d<Type1> &hy, const data3d<Type1> &hz, Type1* CB, Type2*** ID1) {
        unsigned i, j, k;
        unsigned jj, kk, id;
        for (k = 0; k < Kmax - 1; ++k) {
            for (i = 0; i < Imax - 1; ++i) {
                //..............................................................
                //  PML for bottom ex.p, j-direction
                //..............................................................
                for (j = 1; j < nyPML_1; ++j) {
                    id = ID1[i][j][k];
                    psi_Exy_1.p[i][j][k] = be_y_1.p[j] * psi_Exy_1.p[i][j][k]
                            + ce_y_1.p[j] * (hz.p[i][j][k] - hz.p[i][j - 1][k]) / dy;
                    ex.p[i][j][k] = ex.p[i][j][k] + CB[id] * psi_Exy_1.p[i][j][k];
                }
                //.............................................................
                //  PML for top ex.p, j-direction
                //.............................................................
                jj = nyPML_2 - 1;
                for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {
                    id = ID1[i][j][k];
                    psi_Exy_2.p[i][jj][k] = be_y_2.p[jj] * psi_Exy_2.p[i][jj][k]
                            + ce_y_2.p[jj] * (hz.p[i][j][k] - hz.p[i][j - 1][k]) / dy;
                    ex.p[i][j][k] = ex.p[i][j][k] + CB[id] * psi_Exy_2.p[i][jj][k];
                    jj = jj - 1;
                }
            }
        }

        for (i = 0; i < Imax - 1; ++i) {
            for (j = 1; j < Jmax - 1; ++j) {
                //.............................................................
                //  PML for bottom ex.p, k-direction
                //.............................................................
                for (k = 0; k < nzPML_1; ++k) {
                    id = ID1[i][j][k];
                    psi_Exz_1.p[i][j][k] = be_z_1.p[k] * psi_Exz_1.p[i][j][k]
                            + ce_z_1.p[k] * (hy.p[i][j][k] - hy.p[i][j][k + 1]) / dz;
                    ex.p[i][j][k] = ex.p[i][j][k] + CB[id] * psi_Exz_1.p[i][j][k];
                }
                //..............................................................
                //  PML for top ex.p, k-direction
                //..............................................................
                kk = nzPML_2 - 1;
                for (k = Kmax - nzPML_2 - 1; k < Kmax - 1; ++k) {
                    id = ID1[i][j][k];
                    psi_Exz_2.p[i][j][kk] = be_z_2.p[kk] * psi_Exz_2.p[i][j][kk]
                            + ce_z_2.p[kk] * (hy.p[i][j][k] - hy.p[i][j][k + 1]) / dz;
                    ex.p[i][j][k] = ex.p[i][j][k] + CB[id] * psi_Exz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEy(data3d<Type1> &ey, const data3d<Type1> &hz, const data3d<Type1> &hx, Type1* CB, Type2*** ID2) {
        unsigned i, j, k;
        unsigned ii, kk;
        unsigned id;
        for (k = 0; k < Kmax - 1; ++k) {
            for (j = 0; j < Jmax - 1; ++j) {
                //...........................................................
                //  PML for bottom ey.p, i-direction
                //...........................................................
                for (i = 1; i < nxPML_1; ++i) {
                    id = ID2.p[i][j][k];
                    psi_Eyx_1.p[i][j][k] = be_x_1.p[i] * psi_Eyx_1.p[i][j][k]
                            + ce_x_1.p[i] * (hz.p[i - 1][j][k] - hz.p[i][j][k]) / dx;
                    ey.p[i][j][k] = ey.p[i][j][k] + CB[id] * psi_Eyx_1.p[i][j][k];
                }
                //............................................................
                //  PML for top ey.p, i-direction
                //............................................................
                ii = nxPML_2 - 1;
                for (i = Imax - nxPML_2; i < Imax - 1; ++i) {
                    id = ID2.p[i][j][k];
                    psi_Eyx_2.p[ii][j][k] = be_x_2.p[ii] * psi_Eyx_2.p[ii][j][k]
                            + ce_x_2.p[ii] * (hz.p[i - 1][j][k] - hz.p[i][j][k]) / dx;
                    ey.p[i][j][k] = ey.p[i][j][k] + CB[id] * psi_Eyx_2.p[ii][j][k];
                    ii = ii - 1;
                }
            }
        }

        for (i = 1; i < Imax - 1; ++i) {

            for (j = 0; j < Jmax - 1; ++j) {
                //...........................................................
                //  PML for bottom ey.p, k-direction
                //...........................................................
                for (k = 0; k < nzPML_1; ++k) {
                    id = ID2.p[i][j][k];
                    psi_Eyz_1.p[i][j][k] = be_z_1.p[k] * psi_Eyz_1.p[i][j][k]
                            + ce_z_1.p[k] * (hx.p[i][j][k + 1] - hx.p[i][j][k]) / dz;
                    ey.p[i][j][k] = ey.p[i][j][k] + CB[id] * psi_Eyz_1.p[i][j][k];
                }
                //...........................................................
                //  PML for top ey.p, k-direction
                //............................................................
                kk = nzPML_2 - 1;
                for (k = Kmax - nzPML_2 - 1; k < Kmax - 1; ++k) {

                    id = ID2.p[i][j][k];
                    psi_Eyz_2.p[i][j][kk] = be_z_2.p[kk] * psi_Eyz_2.p[i][j][kk]
                            + ce_z_2.p[kk] * (hx.p[i][j][k + 1] - hx.p[i][j][k]) / dz;
                    ey.p[i][j][k] = ey.p[i][j][k] + CB[id] * psi_Eyz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEz(data3d<Type1> &ez, const data3d<Type1> &hx, const data3d<Type1> &hy, Type1* CB, Type2*** ID3) {
        unsigned i, j, k;
        unsigned jj, ii, id;
        for (k = 1; k < Kmax - 1; ++k) {
            for (j = 1; j < Jmax - 1; ++j) {
                //............................................................
                //  PML for bottom ez.p, x-direction
                //.............................................................
                for (i = 1; i < nxPML_1; ++i) {

                    id = ID3[i][j][k];
                    psi_Ezx_1.p[i][j][k] = be_x_1.p[i] * psi_Ezx_1.p[i][j][k]
                            + ce_x_1.p[i] * (hy.p[i][j][k] - hy.p[i - 1][j][k]) / dx;
                    ez.p[i][j][k] = ez.p[i][j][k] + CB[id] * psi_Ezx_1.p[i][j][k];
                }
                //............................................................
                //  PML for top ez.p, x-direction
                //............................................................
                ii = nxPML_2 - 1;
                for (i = Imax - nxPML_2; i < Imax - 1; ++i) {

                    id = ID3[i][j][k];
                    psi_Ezx_2.p[ii][j][k] = be_x_2.p[ii] * psi_Ezx_2.p[ii][j][k]
                            + ce_x_2.p[ii] * (hy.p[i][j][k] - hy.p[i - 1][j][k]) / dx;
                    ez.p[i][j][k] = ez.p[i][j][k] + CB[id] * psi_Ezx_2.p[ii][j][k];
                    ii = ii - 1;
                }
            }

            for (i = 1; i < Imax - 1; ++i) {
                //..........................................................
                //  PML for bottom ez.p, y-direction
                //..........................................................
                for (j = 1; j < nyPML_1; ++j) {

                    id = ID3[i][j][k];
                    psi_Ezy_1.p[i][j][k] = be_y_1.p[j] * psi_Ezy_1.p[i][j][k]
                            + ce_y_1.p[j] * (hx.p[i][j - 1][k] - hx.p[i][j][k]) / dy;
                    ez.p[i][j][k] = ez.p[i][j][k] + CB[id] * psi_Ezy_1.p[i][j][k];
                }
                //............................................................
                //  PML for top ez.p, y-direction
                //............................................................
                jj = nyPML_2 - 1;
                for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {
                    id = ID3[i][j][k];
                    psi_Ezy_2.p[i][jj][k] = be_y_2.p[jj] * psi_Ezy_2.p[i][jj][k]
                            + ce_y_2.p[jj] * (hx.p[i][j - 1][k] - hx.p[i][j][k]) / dy;
                    ez.p[i][j][k] = ez.p[i][j][k] + CB[id] * psi_Ezy_2.p[i][jj][k];
                    jj = jj - 1;
                }
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEx(data3d<Type1> &ex, const data3d<Type1> &hy, const data3d<Type1> &hz, Type1* CB, data3d<Type2> &ID1) {
        UpdatePMLForEx(ex, hy, hz, CB, ID1.p);
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEy(data3d<Type1> &ey, const data3d<Type1> &hz, const data3d<Type1> &hx, Type1* CB, data3d<Type2> &ID2) {
        UpdatePMLForEz(ey, hz, hx, CB, ID2.p);
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEz(data3d<Type1> &ez, const data3d<Type1> &hx, const data3d<Type1> &hy, Type1* CB, data3d<Type2> &ID3) {
        UpdatePMLForEz(ez, hx, hy, CB, ID3.p);
    };

    //=============================================
    // With index 
    //=============================================

    template<typename Type1, typename Type2>
    void UpdatePMLForHxIn(unsigned k, data3d<Type1> &hx, const data3d<Type1> &ey, const data3d<Type1> &ez, Type2 DB) {
        unsigned i, j;
        unsigned jj;

        for (i = 0; i < Imax - 1; ++i) {
            //...............................................
            //  PML for bottom hx.p, j-direction
            //...............................................
            for (j = 0; j < nyPML_1 - 1; ++j) {
                psi_Hxy_1.p[i][j][k] = bh_y_1.p[j] * psi_Hxy_1.p[i][j][k]
                        + ch_y_1.p[j] * (ez.p[i][j][k] - ez.p[i][j + 1][k]) / dy;
                hx.p[i][j][k] = hx.p[i][j][k] + DB * psi_Hxy_1.p[i][j][k];
            }
            //....................................................
            //  PML for top hx.p, j-direction
            //.....................................................
            jj = nyPML_2 - 2;
            for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {
                psi_Hxy_2.p[i][jj][k] = bh_y_2.p[jj] * psi_Hxy_2.p[i][jj][k]
                        + ch_y_2.p[jj] * (ez.p[i][j][k] - ez.p[i][j + 1][k]) / dy;
                hx.p[i][j][k] = hx.p[i][j][k] + DB * psi_Hxy_2.p[i][jj][k];
                jj = jj - 1;
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForHxOut(data3d<Type1> &hx, const data3d<Type1> &ey, const data3d<Type1> &ez, Type2 DB) {
        unsigned i, j, k;
        unsigned kk;
        for (i = 0; i < Imax - 1; ++i) {
            for (j = 0; j < Jmax - 1; ++j) {
                //....................................................
                //  PML for bottom hx.p, k-direction
                //................................................
                for (k = 1; k < nzPML_1; ++k) {
                    psi_Hxz_1.p[i][j][k - 1] = bh_z_1.p[k - 1] * psi_Hxz_1.p[i][j][k - 1]
                            + ch_z_1.p[k - 1] * (ey.p[i][j][k] - ey.p[i][j][k - 1]) / dz;
                    hx.p[i][j][k] = hx.p[i][j][k] + DB * psi_Hxz_1.p[i][j][k - 1];
                }
                //....................................................
                //  PML for top hx.p, k-direction
                //...............................................
                kk = nzPML_2 - 2;
                for (k = Kmax - nzPML_2; k < Kmax - 1; ++k) {
                    psi_Hxz_2.p[i][j][kk] = bh_z_2.p[kk] * psi_Hxz_2.p[i][j][kk]
                            + ch_z_2.p[kk] * (ey.p[i][j][k] - ey.p[i][j][k - 1]) / dz;
                    hx.p[i][j][k] = hx.p[i][j][k] + DB * psi_Hxz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForHyIn(unsigned k, data3d<Type1> &hy, const data3d<Type1> &ez, const data3d<Type1> &ex, Type2 DB) {
        unsigned j, i;
        unsigned ii;

        for (j = 0; j < Jmax - 1; ++j) {
            //.......................................................
            //  PML for bottom hy.p, i-direction
            //.......................................................
            for (i = 0; i < nxPML_1 - 1; ++i) {
                psi_Hyx_1.p[i][j][k] = bh_x_1.p[i] * psi_Hyx_1.p[i][j][k]
                        + ch_x_1.p[i] * (ez.p[i + 1][j][k] - ez.p[i][j][k]) / dx;
                hy.p[i][j][k] = hy.p[i][j][k] + DB * psi_Hyx_1.p[i][j][k];
            }
            //.........................................................
            //  PML for top hy.p, i-direction
            //.........................................................
            ii = nxPML_2 - 2;
            for (i = Imax - nxPML_2; i < Imax - 1; ++i) {
                psi_Hyx_2.p[ii][j][k] = bh_x_2.p[ii] * psi_Hyx_2.p[ii][j][k]
                        + ch_x_2.p[ii] * (ez.p[i + 1][j][k] - ez.p[i][j][k]) / dx;
                hy.p[i][j][k] = hy.p[i][j][k] + DB * psi_Hyx_2.p[ii][j][k];
                ii = ii - 1;
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForHyOut(data3d<Type1> &hy, const data3d<Type1> &ez, const data3d<Type1> &ex, Type2 DB) {
        unsigned i, j, k;
        unsigned kk;

        for (i = 0; i < Imax - 1; ++i) {

            for (j = 0; j < Jmax - 1; ++j) {
                //.......................................................
                //  PML for bottom hy.p, k-direction
                //......................................................
                for (k = 1; k < nzPML_1; ++k) {
                    psi_Hyz_1.p[i][j][k - 1] = bh_z_1.p[k - 1] * psi_Hyz_1.p[i][j][k - 1]
                            + ch_z_1.p[k - 1] * (ex.p[i][j][k - 1] - ex.p[i][j][k]) / dz;
                    hy.p[i][j][k] = hy.p[i][j][k] + DB * psi_Hyz_1.p[i][j][k - 1];
                }
                //.......................................................
                //  PML for top hy.p, k-direction
                //.........................................................
                kk = nzPML_2 - 2;
                for (k = Kmax - nzPML_2; k < Kmax - 1; ++k) {
                    psi_Hyz_2.p[i][j][kk] = bh_z_2.p[kk] * psi_Hyz_2.p[i][j][kk]
                            + ch_z_2.p[kk] * (ex.p[i][j][k - 1] - ex.p[i][j][k]) / dz;
                    hy.p[i][j][k] = hy.p[i][j][k] + DB * psi_Hyz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForHz(unsigned k, data3d<Type1> &hz, const data3d<Type1> &ex, const data3d<Type1> &ey, Type2 DB) {
        unsigned i, j;
        unsigned ii, jj;

        for (j = 0; j < Jmax - 1; ++j) {
            //..........................................................
            //  PML for bottom hz.p, x-direction
            //..........................................................
            for (i = 0; i < nxPML_1 - 1; ++i) {
                psi_Hzx_1.p[i][j][k] = bh_x_1.p[i] * psi_Hzx_1.p[i][j][k]
                        + ch_x_1.p[i] * (ey.p[i][j][k] - ey.p[i + 1][j][k]) / dx;
                hz.p[i][j][k] = hz.p[i][j][k] + DB * psi_Hzx_1.p[i][j][k];
            }
            //..........................................................
            //  PML for top hz.p, x-direction
            //..........................................................
            ii = nxPML_2 - 2;
            for (i = Imax - nxPML_2; i < Imax - 1; ++i) {
                psi_Hzx_2.p[ii][j][k] = bh_x_2.p[ii] * psi_Hzx_2.p[ii][j][k]
                        + ch_x_2.p[ii] * (ey.p[i][j][k] - ey.p[i + 1][j][k]) / dx;
                hz.p[i][j][k] = hz.p[i][j][k] + DB * psi_Hzx_2.p[ii][j][k];
                ii = ii - 1;
            }
        }

        for (i = 0; i < Imax - 1; ++i) {
            //........................................................
            //  PML for bottom hz.p, y-direction
            //.........................................................
            for (j = 0; j < nyPML_1 - 1; ++j) {
                psi_Hzy_1.p[i][j][k] = bh_y_1.p[j] * psi_Hzy_1.p[i][j][k]
                        + ch_y_1.p[j] * (ex.p[i][j + 1][k] - ex.p[i][j][k]) / dy;
                hz.p[i][j][k] = hz.p[i][j][k] + DB * psi_Hzy_1.p[i][j][k];

            }
            //.........................................................
            //  PML for top hz.p, y-direction
            //..........................................................
            jj = nyPML_2 - 2;
            for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {

                psi_Hzy_2.p[i][jj][k] = bh_y_2.p[jj] * psi_Hzy_2.p[i][jj][k]
                        + ch_y_2.p[jj] * (ex.p[i][j + 1][k] - ex.p[i][j][k]) / dy;
                hz.p[i][j][k] = hz.p[i][j][k] + DB * psi_Hzy_2.p[i][jj][k];
                jj = jj - 1;
            }
        }

    };

    template<typename Type1, typename Type2>
    void UpdatePMLForExIn(unsigned k, data3d<Type1> &ex, const data3d<Type1> &hy, const data3d<Type1> &hz, Type1* CB, Type2*** ID1) {
        unsigned i, j;
        unsigned jj, id;
        for (i = 0; i < Imax - 1; ++i) {
            //..............................................................
            //  PML for bottom ex.p, j-direction
            //..............................................................
            for (j = 1; j < nyPML_1; ++j) {
                id = ID1[i][j][k];
                psi_Exy_1.p[i][j][k] = be_y_1.p[j] * psi_Exy_1.p[i][j][k]
                        + ce_y_1.p[j] * (hz.p[i][j][k] - hz.p[i][j - 1][k]) / dy;
                ex.p[i][j][k] = ex.p[i][j][k] + CB[id] * psi_Exy_1.p[i][j][k];
            }
            //.............................................................
            //  PML for top ex.p, j-direction
            //.............................................................
            jj = nyPML_2 - 1;
            for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {
                id = ID1[i][j][k];
                psi_Exy_2.p[i][jj][k] = be_y_2.p[jj] * psi_Exy_2.p[i][jj][k]
                        + ce_y_2.p[jj] * (hz.p[i][j][k] - hz.p[i][j - 1][k]) / dy;
                ex.p[i][j][k] = ex.p[i][j][k] + CB[id] * psi_Exy_2.p[i][jj][k];
                jj = jj - 1;
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForExOut(data3d<Type1> &ex, const data3d<Type1> &hy, const data3d<Type1> &hz, Type1* CB, Type2*** ID1) {
        unsigned i, j, k;
        unsigned kk, id;
        for (i = 0; i < Imax - 1; ++i) {
            for (j = 1; j < Jmax - 1; ++j) {
                //.............................................................
                //  PML for bottom ex.p, k-direction
                //.............................................................
                for (k = 0; k < nzPML_1; ++k) {
                    id = ID1[i][j][k];
                    psi_Exz_1.p[i][j][k] = be_z_1.p[k] * psi_Exz_1.p[i][j][k]
                            + ce_z_1.p[k] * (hy.p[i][j][k] - hy.p[i][j][k + 1]) / dz;
                    ex.p[i][j][k] = ex.p[i][j][k] + CB[id] * psi_Exz_1.p[i][j][k];
                }
                //..............................................................
                //  PML for top ex.p, k-direction
                //..............................................................
                kk = nzPML_2 - 1;
                for (k = Kmax - nzPML_2 - 1; k < Kmax - 1; ++k) {
                    id = ID1[i][j][k];
                    psi_Exz_2.p[i][j][kk] = be_z_2.p[kk] * psi_Exz_2.p[i][j][kk]
                            + ce_z_2.p[kk] * (hy.p[i][j][k] - hy.p[i][j][k + 1]) / dz;
                    ex.p[i][j][k] = ex.p[i][j][k] + CB[id] * psi_Exz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEyIn(unsigned k, data3d<Type1> &ey, const data3d<Type1> &hz, const data3d<Type1> &hx, Type1* CB, Type2*** ID2) {
        unsigned i, j;
        unsigned ii;
        unsigned id;

        for (j = 0; j < Jmax - 1; ++j) {
            //...........................................................
            //  PML for bottom ey.p, i-direction
            //...........................................................
            for (i = 1; i < nxPML_1; ++i) {
                id = ID2[i][j][k];
                psi_Eyx_1.p[i][j][k] = be_x_1.p[i] * psi_Eyx_1.p[i][j][k]
                        + ce_x_1.p[i] * (hz.p[i - 1][j][k] - hz.p[i][j][k]) / dx;
                ey.p[i][j][k] = ey.p[i][j][k] + CB[id] * psi_Eyx_1.p[i][j][k];
            }
            //............................................................
            //  PML for top ey.p, i-direction
            //............................................................
            ii = nxPML_2 - 1;
            for (i = Imax - nxPML_2; i < Imax - 1; ++i) {
                id = ID2[i][j][k];
                psi_Eyx_2.p[ii][j][k] = be_x_2.p[ii] * psi_Eyx_2.p[ii][j][k]
                        + ce_x_2.p[ii] * (hz.p[i - 1][j][k] - hz.p[i][j][k]) / dx;
                ey.p[i][j][k] = ey.p[i][j][k] + CB[id] * psi_Eyx_2.p[ii][j][k];
                ii = ii - 1;
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEyOut(data3d<Type1> &ey, const data3d<Type1> &hz, const data3d<Type1> &hx, Type1* CB, Type2*** ID2) {
        unsigned i, j, k;
        unsigned kk;
        Type2 id;
        for (i = 1; i < Imax - 1; ++i) {
            for (j = 0; j < Jmax - 1; ++j) {
                //...........................................................
                //  PML for bottom ey.p, k-direction
                //...........................................................
                for (k = 0; k < nzPML_1; ++k) {
                    id = ID2[i][j][k];
                    psi_Eyz_1.p[i][j][k] = be_z_1.p[k] * psi_Eyz_1.p[i][j][k]
                            + ce_z_1.p[k] * (hx.p[i][j][k + 1] - hx.p[i][j][k]) / dz;
                    ey.p[i][j][k] = ey.p[i][j][k] + CB[id] * psi_Eyz_1.p[i][j][k];
                }
                //...........................................................
                //  PML for top ey.p, k-direction
                //............................................................
                kk = nzPML_2 - 1;
                for (k = Kmax - nzPML_2 - 1; k < Kmax - 1; ++k) {

                    id = ID2[i][j][k];
                    psi_Eyz_2.p[i][j][kk] = be_z_2.p[kk] * psi_Eyz_2.p[i][j][kk]
                            + ce_z_2.p[kk] * (hx.p[i][j][k + 1] - hx.p[i][j][k]) / dz;
                    ey.p[i][j][k] = ey.p[i][j][k] + CB[id] * psi_Eyz_2.p[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEz(unsigned k, data3d<Type1> &ez, const data3d<Type1> &hx, const data3d<Type1> &hy, Type1* CB, Type2*** ID3) {
        unsigned i, j;
        unsigned jj, ii;
        Type2 id;

        for (j = 1; j < Jmax - 1; ++j) {
            //............................................................
            //  PML for bottom ez.p, x-direction
            //.............................................................
            for (i = 1; i < nxPML_1; ++i) {

                id = ID3[i][j][k];
                psi_Ezx_1.p[i][j][k] = be_x_1.p[i] * psi_Ezx_1.p[i][j][k]
                        + ce_x_1.p[i] * (hy.p[i][j][k] - hy.p[i - 1][j][k]) / dx;
                ez.p[i][j][k] = ez.p[i][j][k] + CB[id] * psi_Ezx_1.p[i][j][k];
            }
            //............................................................
            //  PML for top ez.p, x-direction
            //............................................................
            ii = nxPML_2 - 1;
            for (i = Imax - nxPML_2; i < Imax - 1; ++i) {

                id = ID3[i][j][k];
                psi_Ezx_2.p[ii][j][k] = be_x_2.p[ii] * psi_Ezx_2.p[ii][j][k]
                        + ce_x_2.p[ii] * (hy.p[i][j][k] - hy.p[i - 1][j][k]) / dx;
                ez.p[i][j][k] = ez.p[i][j][k] + CB[id] * psi_Ezx_2.p[ii][j][k];
                ii = ii - 1;
            }
        }

        for (i = 1; i < Imax - 1; ++i) {
            //..........................................................
            //  PML for bottom ez.p, y-direction
            //..........................................................
            for (j = 1; j < nyPML_1; ++j) {
                id = ID3[i][j][k];
                psi_Ezy_1.p[i][j][k] = be_y_1.p[j] * psi_Ezy_1.p[i][j][k]
                        + ce_y_1.p[j] * (hx.p[i][j - 1][k] - hx.p[i][j][k]) / dy;
                ez.p[i][j][k] = ez.p[i][j][k] + CB[id] * psi_Ezy_1.p[i][j][k];
            }
            //............................................................
            //  PML for top ez.p, y-direction
            //............................................................
            jj = nyPML_2 - 1;
            for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {
                id = ID3[i][j][k];
                psi_Ezy_2.p[i][jj][k] = be_y_2.p[jj] * psi_Ezy_2.p[i][jj][k]
                        + ce_y_2.p[jj] * (hx.p[i][j - 1][k] - hx.p[i][j][k]) / dy;
                ez.p[i][j][k] = ez.p[i][j][k] + CB[id] * psi_Ezy_2.p[i][jj][k];
                jj = jj - 1;
            }
        }

    };

    template<typename Type1, typename Type2>
    void UpdatePMLForExIn(unsigned k, data3d<Type1> &ex, const data3d<Type1> &hy, const data3d<Type1> &hz, Type1* CB, data3d<Type2> &ID1) {
        UpdatePMLForExIn(k, ex, hy, hz, CB, ID1.p);
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEyIn(unsigned k, data3d<Type1> &ey, const data3d<Type1> &hz, const data3d<Type1> &hx, Type1* CB, data3d<Type2> &ID2) {
        UpdatePMLForEyIn(k, ey, hz, hx, CB, ID2.p);
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForExOut(data3d<Type1> &ex, const data3d<Type1> &hy, const data3d<Type1> &hz, Type1* CB, data3d<Type2> &ID1) {
        UpdatePMLForExOut(ex, hy, hz, CB, ID1.p);
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEyOut(data3d<Type1> &ey, const data3d<Type1> &hz, const data3d<Type1> &hx, Type1* CB, data3d<Type2> &ID2) {
        UpdatePMLForEyOut(ey, hz, hx, CB, ID2.p);
    };

    template<typename Type1, typename Type2>
    void UpdatePMLForEz(unsigned k, data3d<Type1> &ez, const data3d<Type1> &hx, const data3d<Type1> &hy, Type1* CB, data3d<Type2> &ID3) {
        UpdatePMLForEz(k, ez, hx, hy, CB, ID3.p);
    };

private:
    void InitializeCPML();
    int CreatePMLArrays();
};

