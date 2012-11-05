
#include <exception>
#include <iostream>
#include <math.h>

#include "cpml2.h"

using namespace std;

cpml::~cpml(void) {

}

int cpml::Initial(unsigned nx, unsigned ny, unsigned nz, unsigned ncpml) {

    //PML Layers (10 layers)
    nxPML_1 = ncpml;
    nxPML_2 = ncpml;
    nyPML_1 = ncpml;
    nyPML_2 = ncpml;
    nzPML_1 = ncpml;
    nzPML_2 = ncpml;

    Imax = nx;
    Jmax = ny;
    Kmax = nz;

    CreatePMLArrays();
    InitializeCPML();

    return 0;
}

int cpml::CreatePMLArrays() {
    try
    {
        unsigned i;

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
        //psi_Eyz_1.CreateStruct(Imax, Jmax - 1, nzPML_1, 0);
        //psi_Eyz_2.CreateStruct(Imax, Jmax - 1, nzPML_2, 0);
        
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

        be_x_1.createArray(nxPML_1,0.0);
        ce_x_1.createArray(nxPML_1,0.0);
        alphae_x_PML_1.createArray(nxPML_1,0.0);
        sige_x_PML_1.createArray(nxPML_1,0.0);
        kappae_x_PML_1.createArray(nxPML_1,0.0);
        bh_x_1.createArray(nxPML_1-1,0.0);
        ch_x_1.createArray(nxPML_1-1,0.0);
        alphah_x_PML_1.createArray(nxPML_1-1,0.0);
        sigh_x_PML_1.createArray(nxPML_1-1,0.0);
        kappah_x_PML_1.createArray(nxPML_1-1,0.0);


        be_x_2.createArray(nxPML_2,0.0);
        ce_x_2.createArray(nxPML_2,0.0);
        alphae_x_PML_2.createArray(nxPML_2,0.0);
        sige_x_PML_2.createArray(nxPML_2,0.0);
        kappae_x_PML_2.createArray(nxPML_2,0.0);

        bh_x_2.createArray(nxPML_2-1,0.0);
        ch_x_2.createArray(nxPML_2-1,0.0);
        alphah_x_PML_2.createArray(nxPML_2-1,0.0);
        sigh_x_PML_2.createArray(nxPML_2-1,0.0);
        kappah_x_PML_2.createArray(nxPML_2-1,0.0);


        be_y_1.createArray(nxPML_1,0.0);
        ce_y_1.createArray(nxPML_1,0.0);
        kappae_y_PML_1.createArray(nxPML_1,0.0);
        sige_y_PML_1.createArray(nxPML_1,0.0);
        alphae_y_PML_1.createArray(nxPML_1,0.0);


        bh_y_1.createArray(nyPML_1 - 1,0.0);
        ch_y_1.createArray(nyPML_1 - 1,0.0);
        alphah_y_PML_1.createArray(nyPML_1 - 1,0.0);
        sigh_y_PML_1.createArray(nyPML_1 - 1,0.0);
        kappah_y_PML_1.createArray(nyPML_1 - 1,0.0);

        be_y_2.createArray(nyPML_2,0.0);
        ce_y_2.createArray(nyPML_2,0.0);
        alphae_y_PML_2.createArray(nyPML_2,0.0);
        sige_y_PML_2.createArray(nyPML_2,0.0);
        kappae_y_PML_2.createArray(nyPML_2,0.0);

        bh_y_2.createArray(nyPML_2 - 1,0.0);
        ch_y_2.createArray(nyPML_2 - 1,0.0);
        alphah_y_PML_2.createArray(nyPML_2 - 1,0.0);
        sigh_y_PML_2.createArray(nyPML_2 - 1,0.0);
        kappah_y_PML_2.createArray(nyPML_2 - 1,0.0);

        be_z_1.createArray(nzPML_1,0.0);
        ce_z_1.createArray(nzPML_1,0.0);
        alphae_z_PML_1.createArray(nzPML_1,0.0);
        sige_z_PML_1.createArray(nzPML_1,0.0);
        kappae_z_PML_1.createArray(nzPML_1,0.0);

        bh_z_1.createArray(nzPML_1 - 1,0.0);
        ch_z_1.createArray(nzPML_1 - 1,0.0);
        alphah_z_PML_1.createArray(nzPML_1 - 1,0.0);
        sigh_z_PML_1.createArray(nzPML_1 - 1,0.0);
        kappah_z_PML_1.createArray(nzPML_1 - 1,0.0);

        be_z_2.createArray(nzPML_2,0.0);
        ce_z_2.createArray(nzPML_2,0.0);
        alphae_z_PML_2.createArray(nzPML_2,0.0);
        sige_z_PML_2.createArray(nzPML_2,0.0);
         kappae_z_PML_2.createArray(nzPML_2,0.0);

        bh_z_2.createArray(nzPML_2 - 1,0.0);
        ch_z_2.createArray(nzPML_2 - 1,0.0);
        alphah_z_PML_2.createArray(nzPML_2 - 1,0.0);
        sigh_z_PML_2.createArray(nzPML_2 - 1,0.0);
        kappah_z_PML_2.createArray(nzPML_2 - 1,0.0);

        den_ex.createArray(Imax - 1,0.0);

        den_hx.createArray(Imax - 1,0.0);

        den_ey.createArray(Jmax - 1,0.0);

        den_hy.createArray(Jmax - 1,0.0);

        den_ez.createArray(Kmax - 1,0.0);

        den_hz.createArray(Kmax - 1,0.0);
    }

    catch(exception & e) {
        cerr << "Error:" << e.what() << endl;
        return -1;
    }
}

void cpml::InitializeCPML() {
    unsigned i, j, k;
    unsigned ii, jj, kk;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  PML parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    sig_x_max = 0.75 * (0.8 * (m + 1) / (dx * sqrt(mu_0 / (eps_0 * epsR))));
    sig_y_max = 0.75 * (0.8 * (m + 1) / (dy * sqrt(mu_0 / (eps_0 * epsR))));
    sig_z_max = 0.75 * (0.8 * (m + 1) / (dz * sqrt(mu_0 / (eps_0 * epsR))));
    alpha_x_max = 0.24;
    alpha_y_max = alpha_x_max;
    alpha_z_max = alpha_x_max;
    kappa_x_max = 15.0;
    kappa_y_max = kappa_x_max;
    kappa_z_max = kappa_x_max;

    for (i = 0; i < nxPML_1; ++i) {

        sige_x_PML_1.p[i] = sig_x_max * pow(((nxPML_1 - 1 - i)
                / (nxPML_1 - 1.0)), m);
        alphae_x_PML_1.p[i] = alpha_x_max * pow(((MyDataF) i
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
        alphae_x_PML_2.p[i] = alpha_x_max * pow(((MyDataF) i
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
        alphae_y_PML_1.p[j] = alpha_y_max * pow(((MyDataF) j
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
        alphae_y_PML_2.p[j] = alpha_y_max * pow(((MyDataF) j
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
        alphae_z_PML_1.p[k] = alpha_z_max * pow(((MyDataF) k
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
        alphae_z_PML_2.p[k] = alpha_z_max * pow(((MyDataF) k
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

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  DENOMINATORS FOR FIELD UPDATES
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
}
/*
void cpml::UpdatePMLForHx(data3d<MyDataF> &hx, const data3d<MyDataF> &ey, const data3d<MyDataF> &ez, MyDataF DB) {
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
}

void cpml::UpdatePMLForHy(data3d<MyDataF> &hy, const data3d<MyDataF> &ez, const data3d<MyDataF> &ex, MyDataF DB) {
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
}

void cpml::UpdatePMLForHz(data3d<MyDataF> &hz, const data3d<MyDataF> &ex, const data3d<MyDataF> &ey, MyDataF DB) {
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
}

void cpml::UpdatePMLForEx(data3d<MyDataF> &ex, const data3d<MyDataF> &hy, const data3d<MyDataF> &hz, MyDataF *CB, unsigned*** ID1) {
    unsigned i, j, k;
    unsigned jj, kk, id;
    for (k = 0; k < Kmax - 1; ++k) {
        for (i = 0; i < Imax - 1; ++i) {
            //..............................................................
            //  PML for bottom ex.p, j-direction
            //..............................................................
            for (j = 1; j < nyPML_1; ++j) {
                id = ID1.p[i][j][k];
                psi_Exy_1.p[i][j][k] = be_y_1.p[j] * psi_Exy_1.p[i][j][k]
                        + ce_y_1.p[j] * (hz.p[i][j][k] - hz.p[i][j - 1][k]) / dy;
                ex.p[i][j][k] = ex.p[i][j][k] + CB.p[id] * psi_Exy_1.p[i][j][k];
            }
            //.............................................................
            //  PML for top ex.p, j-direction
            //.............................................................
            jj = nyPML_2 - 1;
            for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {
                id = ID1.p[i][j][k];
                psi_Exy_2.p[i][jj][k] = be_y_2.p[jj] * psi_Exy_2.p[i][jj][k]
                        + ce_y_2.p[jj] * (hz.p[i][j][k] - hz.p[i][j - 1][k]) / dy;
                ex.p[i][j][k] = ex.p[i][j][k] + CB.p[id] * psi_Exy_2.p[i][jj][k];
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
                id = ID1.p[i][j][k];
                psi_Exz_1.p[i][j][k] = be_z_1.p[k] * psi_Exz_1.p[i][j][k]
                        + ce_z_1.p[k] * (hy.p[i][j][k] - hy.p[i][j][k + 1]) / dz;
                ex.p[i][j][k] = ex.p[i][j][k] + CB.p[id] * psi_Exz_1.p[i][j][k];
            }
            //..............................................................
            //  PML for top ex.p, k-direction
            //..............................................................
            kk = nzPML_2 - 1;
            for (k = Kmax - nzPML_2 - 1; k < Kmax - 1; ++k) {
                id = ID1.p[i][j][k];
                psi_Exz_2.p[i][j][kk] = be_z_2.p[kk] * psi_Exz_2.p[i][j][kk]
                        + ce_z_2.p[kk] * (hy.p[i][j][k] - hy.p[i][j][k + 1]) / dz;
                ex.p[i][j][k] = ex.p[i][j][k] + CB.p[id] * psi_Exz_2.p[i][j][kk];
                kk = kk - 1;
            }
        }
    }
}

void cpml::UpdatePMLForEy(data3d<MyDataF> &ey, const data3d<MyDataF> &hz, const data3d<MyDataF> &hx, MyDataF* CB, unsigned*** ID2) {
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
                ey.p[i][j][k] = ey.p[i][j][k] + CB.p[id] * psi_Eyx_1.p[i][j][k];
            }
            //............................................................
            //  PML for top ey.p, i-direction
            //............................................................
            ii = nxPML_2 - 1;
            for (i = Imax - nxPML_2; i < Imax - 1; ++i) {
                id = ID2.p[i][j][k];
                psi_Eyx_2.p[ii][j][k] = be_x_2.p[ii] * psi_Eyx_2.p[ii][j][k]
                        + ce_x_2.p[ii] * (hz.p[i - 1][j][k] - hz.p[i][j][k]) / dx;
                ey.p[i][j][k] = ey.p[i][j][k] + CB.p[id] * psi_Eyx_2.p[ii][j][k];
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
                ey.p[i][j][k] = ey.p[i][j][k] + CB.p[id] * psi_Eyz_1.p[i][j][k];
            }
            //...........................................................
            //  PML for top ey.p, k-direction
            //............................................................
            kk = nzPML_2 - 1;
            for (k = Kmax - nzPML_2 - 1; k < Kmax - 1; ++k) {

                id = ID2.p[i][j][k];
                psi_Eyz_2.p[i][j][kk] = be_z_2.p[kk] * psi_Eyz_2.p[i][j][kk]
                        + ce_z_2.p[kk] * (hx.p[i][j][k + 1] - hx.p[i][j][k]) / dz;
                ey.p[i][j][k] = ey.p[i][j][k] + CB.p[id] * psi_Eyz_2.p[i][j][kk];
                kk = kk - 1;
            }
        }
    }
}

void cpml::UpdatePMLForEz(data3d<MyDataF> &ez, const data3d<MyDataF> &hx, const data3d<MyDataF> &hy, MyDataF* CB, unsigned*** ID3) {
    unsigned i, j, k;
    unsigned jj, ii, id;
    for (k = 1; k < Kmax - 1; ++k) {
        for (j = 1; j < Jmax - 1; ++j) {
            //............................................................
            //  PML for bottom ez.p, x-direction
            //.............................................................
            for (i = 1; i < nxPML_1; ++i) {

                id = ID3.p[i][j][k];
                psi_Ezx_1.p[i][j][k] = be_x_1.p[i] * psi_Ezx_1.p[i][j][k]
                        + ce_x_1.p[i] * (hy.p[i][j][k] - hy.p[i - 1][j][k]) / dx;
                ez.p[i][j][k] = ez.p[i][j][k] + CB.p[id] * psi_Ezx_1.p[i][j][k];
            }
            //............................................................
            //  PML for top ez.p, x-direction
            //............................................................
            ii = nxPML_2 - 1;
            for (i = Imax - nxPML_2; i < Imax - 1; ++i) {

                id = ID3.p[i][j][k];
                psi_Ezx_2.p[ii][j][k] = be_x_2.p[ii] * psi_Ezx_2.p[ii][j][k]
                        + ce_x_2.p[ii] * (hy.p[i][j][k] - hy.p[i - 1][j][k]) / dx;
                ez.p[i][j][k] = ez.p[i][j][k] + CB.p[id] * psi_Ezx_2.p[ii][j][k];
                ii = ii - 1;
            }
        }

        for (i = 1; i < Imax - 1; ++i) {
            //..........................................................
            //  PML for bottom ez.p, y-direction
            //..........................................................
            for (j = 1; j < nyPML_1; ++j) {

                id = ID3.p[i][j][k];
                psi_Ezy_1.p[i][j][k] = be_y_1.p[j] * psi_Ezy_1.p[i][j][k]
                        + ce_y_1.p[j] * (hx.p[i][j - 1][k] - hx.p[i][j][k]) / dy;
                ez.p[i][j][k] = ez.p[i][j][k] + CB.p[id] * psi_Ezy_1.p[i][j][k];
            }
            //............................................................
            //  PML for top ez.p, y-direction
            //............................................................
            jj = nyPML_2 - 1;
            for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {
                id = ID3.p[i][j][k];
                psi_Ezy_2.p[i][jj][k] = be_y_2.p[jj] * psi_Ezy_2.p[i][jj][k]
                        + ce_y_2.p[jj] * (hx.p[i][j - 1][k] - hx.p[i][j][k]) / dy;
                ez.p[i][j][k] = ez.p[i][j][k] + CB.p[id] * psi_Ezy_2.p[i][jj][k];
                jj = jj - 1;
            }
        }
    }
}

 */
