/* 
 * File:   cpmld.h
 * Author: skiloop
 *
 * Created on April 3, 2013, 8:41 AM
 */

#ifndef CPMLD_H
#define	CPMLD_H

#include "datastruct.h"

class cpmld {
public:
    cpmld();
    cpmld(const cpmld& orig);
    virtual ~cpmld();
    //  CPML components (Taflove 3rd Edition, Chapter 7)
    data3d<float> psi_Ezx_1;
    data3d<float> psi_Ezx_2;
    data3d<float> psi_Hyx_1;
    data3d<float> psi_Hyx_2;
    data3d<float> psi_Ezy_1;
    data3d<float> psi_Ezy_2;
    data3d<float> psi_Hxy_1;
    data3d<float> psi_Hxy_2;
    data3d<float> psi_Hxz_1;
    data3d<float> psi_Hxz_2;
    data3d<float> psi_Hyz_1;
    data3d<float> psi_Hyz_2;
    data3d<float> psi_Exz_1;
    data3d<float> psi_Exz_2;
    data3d<float> psi_Eyz_1;
    data3d<float> psi_Eyz_2;
    data3d<float> psi_Hzx_1;
    data3d<float> psi_Eyx_1;
    data3d<float> psi_Hzx_2;
    data3d<float> psi_Eyx_2;
    data3d<float> psi_Hzy_1;
    data3d<float> psi_Exy_1;
    data3d<float> psi_Hzy_2;
    data3d<float> psi_Exy_2;

    data1d<float> be_x_1, ce_x_1, alphae_x_PML_1, sige_x_PML_1, kappae_x_PML_1;
    data1d<float> bh_x_1, ch_x_1, alphah_x_PML_1, sigh_x_PML_1, kappah_x_PML_1;
    data1d<float> be_x_2, ce_x_2, alphae_x_PML_2, sige_x_PML_2, kappae_x_PML_2;
    data1d<float> bh_x_2, ch_x_2, alphah_x_PML_2, sigh_x_PML_2, kappah_x_PML_2;
    data1d<float> be_y_1, ce_y_1, alphae_y_PML_1, sige_y_PML_1, kappae_y_PML_1;
    data1d<float> bh_y_1, ch_y_1, alphah_y_PML_1, sigh_y_PML_1, kappah_y_PML_1;
    data1d<float> be_y_2, ce_y_2, alphae_y_PML_2, sige_y_PML_2, kappae_y_PML_2;
    data1d<float> bh_y_2, ch_y_2, alphah_y_PML_2, sigh_y_PML_2, kappah_y_PML_2;
    data1d<float> be_z_1, ce_z_1, alphae_z_PML_1, sige_z_PML_1, kappae_z_PML_1;
    data1d<float> bh_z_1, ch_z_1, alphah_z_PML_1, sigh_z_PML_1, kappah_z_PML_1;
    data1d<float> be_z_2, ce_z_2, alphae_z_PML_2, sige_z_PML_2, kappae_z_PML_2;
    data1d<float> bh_z_2, ch_z_2, alphah_z_PML_2, sigh_z_PML_2, kappah_z_PML_2;

    // denominators for the update equations
    data1d<float> den_ex;
    data1d<float> den_hx;
    data1d<float> den_ey;
    data1d<float> den_hy;
    data1d<float> den_ez;
    data1d<float> den_hz;

    // Specify the CPML Thickness in Each Direction (Value of Zero
    // Corresponds to No PML, and the Grid is Terminated with a PEC)
    // PML thickness in each direction
    int nxPML_1, nxPML_2, nyPML_1;
    int nyPML_2, nzPML_1, nzPML_2;

    // grid size corresponding to the number of Ez field components
    int Imax;
    int Jmax;
    int Kmax;

    //  Fundamental Constants (MKS units)
    static const double pi = 3.14159265358979;
    static const double C = 2.99792458E8;
    static double mu_0;
    static double eps_0;
    // Specify Material Relative Permittivity and Conductivity
    static double epsR; //free space

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

    void Initial(unsigned nx, unsigned ny, unsigned nz, unsigned ncpml);
    static void InitialMuEps();
    void // <editor-fold defaultstate="collapsed" desc="comment">
    createCPMLArray()// </editor-fold>
    ;
    void createPsi();
    void createCBKAP();
    void createDen();
    void initCPML(float dt, float dx, float dy, float dz);
    void initPsi();
    void initCBKAP(float dt, float dx, float dy, float dz);
    void initDen(float dt, float dx, float dy, float dz);
    void initParmeters(float dx, float dy, float dz, int m_, int ma_);
private:

};

#endif	/* CPMLD_H */

