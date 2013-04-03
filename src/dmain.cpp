
/*******************************************************************
//     3-D FDTD code with CPML absorbing boundary conditions - Dipole
// ******************************************************************
//
//     This C code implements the finite-difference time-domain
//     solution of Maxwell's curl equations over a three-dimensional
//     Cartesian space lattice.  The grid is terminated by CPML ABCs.
//     The thickness of the PML in each Cartesian direction can be 
//     varied independently.
//
// ******************************************************************/
//Header files (Libraries to be included)
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

#include"datastruct.h"
#include"cpmld.h"

//  Fundamental Constants (MKS units)
double pi = 3.14159265358979;
double C = 2.99792458E8;
double mu_0;
double eps_0;

//  Specify Material Relative Permittivity and Conductivity
double epsR = 1.0; //free space

//  Specify Number of Time Steps and Grid Size Parameters
int nMax = 500; // total number of time steps

// grid size corresponding to the number of Ez field components
int Imax = 51;
int Jmax = 126;
int Kmax = 26;

//  Specify Grid Cell Size in Each Direction and Calculate the
//  Resulting Courant-Stable Time Step
double dx = 1.0E-3;
double dy = 1.0E-3;
double dz = 1.0E-3; // cell size in each direction
// time step increment
double dt;

//  Specify the Impulsive Source (Differentiated Gaussian) parameters
double tw = 53.0E-12; //pulse width
double tO; //delay
double source; //Differentiated Gaussian source
double amp = 1000; // Amplitude

//Specify the Time Step at which the data has to be saved for Visualization
int save_modulus = 10;

//  Specify the dipole Boundaries(A cuboidal rode- NOT as a cylinder)
int istart, iend, jstart;
int jend, kstart, kend;

//Output recording point
int ksource = 12;

//  Specify the CPML Order and Other Parameters:
int m = 3, ma = 1;

//double pml.sig_x_max;
//double pml.sig_y_max;
//double pml.sig_z_max;
//double pml.alpha_x_max;
//double pml.alpha_y_max;
//double pml.alpha_z_max;
//double pml.kappa_x_max;
//double pml.kappa_y_max;
//double pml.kappa_z_max;

//Loop indices
int i, j, ii, jj, k, kk, n;


// H & E Field components
data3d<float> Hx;
data3d<float> Hy;
data3d<float> Hz;
data3d<float> Ex;
data3d<float> Ey;
data3d<float> Ez;

data3d<short> ID1; //medium definition array for Ex
data3d<short> ID2; //medium definition array for Ey
data3d<short> ID3; //medium definition array for Ez

// cpml
cpmld<float, short> pml;
//Max number of materials allowed
int numMaterials = 50;

//permittivity, permeability and conductivity of diffrent materials
double *epsilon;
double *mu;
double *sigma;

//E field update coefficients
float *CA;
float *CB;

//H field update coefficients
float DA;
float DB;

//Function prototype definitions
void initialize(); //Memeory initialization
void setUp(); //Coefficients, parameters etc will get computed
void initializeCPML(); //CPML coefficient computation
void compute(); //E & H Field update equation
void buildObject(); //Creates the object geometry
void yeeCube(int, int, int, short); //Sets material properties to a cell
void writeField(int); //Writes output
void buildSphere(); //Builds a spherical object
void buildDipole(); //Builds a dipole

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int main() {

    initialize();
    setUp();
    buildObject();
    initializeCPML();
    compute();
    return 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void initialize() {

    mu_0 = 4.0 * pi * 1.0E-7;
    eps_0 = 1.0 / (C * C * mu_0);

    pml.InitialMuEps();

    pml.Initial(Imax, Jmax, Kmax, 11);
    //PML Layers (10 layers)
    //    pml.nxPML_1 = 11;
    //    pml.nxPML_2 = 11;
    //    pml.nyPML_1 = 11;
    //    pml.nyPML_2 = 11;
    //    pml.nzPML_1 = 11;
    //    pml.nzPML_2 = 11;

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

    CA = (float *) malloc((numMaterials) * sizeof (float));

    for (i = 0; i < numMaterials; i++) {

        CA[i] = 0.0;
    }

    CB = (float *) malloc((numMaterials) * sizeof (float));

    for (i = 0; i < numMaterials; i++) {

        CB[i] = 0.0;
    }

    //    pml.psi_Ezx_1.CreateStruct(pml.nxPML_1, Jmax, Kmax, 0);
    //    pml.psi_Ezx_2.CreateStruct(pml.nxPML_2, Jmax, Kmax, 0);
    //    pml.psi_Hyx_1.CreateStruct(pml.nxPML_1 - 1, Jmax, Kmax, 0);
    //    pml.psi_Hyx_2.CreateStruct(pml.nxPML_2 - 1, Jmax, Kmax, 0);
    //
    //    pml.psi_Ezy_1.CreateStruct(Imax, pml.nyPML_1, Kmax, 0);
    //    pml.psi_Ezy_2.CreateStruct(Imax, pml.nyPML_2, Kmax, 0);
    //    pml.psi_Hxy_1.CreateStruct(Imax, pml.nyPML_1 - 1, Kmax, 0);
    //    pml.psi_Hxy_2.CreateStruct(Imax, pml.nyPML_2 - 1, Kmax, 0);
    //
    //    pml.psi_Hxz_1.CreateStruct(Imax, Jmax - 1, pml.nzPML_1 - 1, 0);
    //    pml.psi_Hxz_2.CreateStruct(Imax, Jmax - 1, pml.nzPML_2 - 1, 0);
    //    pml.psi_Eyz_1.CreateStruct(Imax, Jmax - 1, pml.nzPML_1, 0);
    //    pml.psi_Eyz_2.CreateStruct(Imax, Jmax - 1, pml.nzPML_2, 0);
    //
    //    pml.psi_Hyz_1.CreateStruct(Imax - 1, Jmax, pml.nzPML_1 - 1, 0);
    //    pml.psi_Hyz_2.CreateStruct(Imax - 1, Jmax, pml.nzPML_2 - 1, 0);
    //    pml.psi_Exz_1.CreateStruct(Imax - 1, Jmax, pml.nzPML_1, 0);
    //    pml.psi_Exz_2.CreateStruct(Imax - 1, Jmax, pml.nzPML_2, 0);
    //
    //    pml.psi_Hzx_1.CreateStruct(pml.nxPML_1 - 1, Jmax - 1, Kmax - 1, 0);
    //    pml.psi_Hzx_2.CreateStruct(pml.nxPML_2 - 1, Jmax - 1, Kmax - 1, 0);
    //    pml.psi_Eyx_1.CreateStruct(pml.nxPML_1, Jmax - 1, Kmax - 1, 0);
    //    pml.psi_Eyx_2.CreateStruct(pml.nxPML_2, Jmax - 1, Kmax - 1, 0);
    //
    //    pml.psi_Hzy_1.CreateStruct(Imax - 1, pml.nyPML_1 - 1, Kmax - 1, 0);
    //    pml.psi_Hzy_2.CreateStruct(Imax - 1, pml.nyPML_2 - 1, Kmax - 1, 0);
    //    pml.psi_Exy_1.CreateStruct(Imax - 1, pml.nyPML_1, Kmax - 1, 0);
    //    pml.psi_Exy_2.CreateStruct(Imax - 1, pml.nyPML_2, Kmax - 1, 0);

    Ez.CreateStruct(Imax, Jmax, Kmax, 0);
    Ey.CreateStruct(Imax, Jmax - 1, Kmax - 1, 0);
    Ex.CreateStruct(Imax - 1, Jmax, Kmax - 1, 0);

    Hx.CreateStruct(Imax, Jmax - 1, Kmax, 0);
    Hy.CreateStruct(Imax - 1, Jmax, Kmax, 0);
    Hz.CreateStruct((Imax - 1), (Jmax - 1), (Kmax - 1), 0);

    ID1.CreateStruct(Imax, Jmax, Kmax, 0);
    ID2.CreateStruct(Imax, Jmax, Kmax, 0);
    ID3.CreateStruct(Imax, Jmax, Kmax, 0);

    //    pml.be_x_1.createArray(pml.nxPML_1, 0.0);
    //    pml.ce_x_1.createArray(pml.nxPML_1, 0.0);
    //    pml.alphae_x_PML_1.createArray(pml.nxPML_1, 0.0);
    //    pml.sige_x_PML_1.createArray(pml.nxPML_1, 0.0);
    //    pml.kappae_x_PML_1.createArray(pml.nxPML_1, 0.0);
    //    pml.bh_x_1.createArray(pml.nxPML_1 - 1, 0.0);
    //    pml.ch_x_1.createArray(pml.nxPML_1 - 1, 0.0);
    //    pml.alphah_x_PML_1.createArray(pml.nxPML_1 - 1, 0.0);
    //    pml.sigh_x_PML_1.createArray(pml.nxPML_1 - 1, 0.0);
    //    pml.kappah_x_PML_1.createArray(pml.nxPML_1 - 1, 0.0);
    //
    //    pml.be_x_2.createArray(pml.nxPML_2, 0.0);
    //    pml.ce_x_2.createArray(pml.nxPML_2, 0.0);
    //    pml.alphae_x_PML_2.createArray(pml.nxPML_2, 0.0);
    //    pml.sige_x_PML_2.createArray(pml.nxPML_2, 0.0);
    //    pml.kappae_x_PML_2.createArray(pml.nxPML_2, 0.0);
    //
    //    pml.bh_x_2.createArray(pml.nxPML_2 - 1, 0.0);
    //    pml.ch_x_2.createArray(pml.nxPML_2 - 1, 0.0);
    //    pml.alphah_x_PML_2.createArray(pml.nxPML_2 - 1, 0.0);
    //    pml.sigh_x_PML_2.createArray(pml.nxPML_2 - 1, 0.0);
    //    pml.kappah_x_PML_2.createArray(pml.nxPML_2 - 1, 0.0);
    //
    //    pml.be_y_1.createArray(pml.nxPML_1, 0.0);
    //    pml.ce_y_1.createArray(pml.nxPML_1, 0.0);
    //    pml.kappae_y_PML_1.createArray(pml.nxPML_1, 0.0);
    //    pml.sige_y_PML_1.createArray(pml.nxPML_1, 0.0);
    //    pml.alphae_y_PML_1.createArray(pml.nxPML_1, 0.0);
    //
    //
    //    pml.bh_y_1.createArray(pml.nyPML_1 - 1, 0.0);
    //    pml.ch_y_1.createArray(pml.nyPML_1 - 1, 0.0);
    //    pml.alphah_y_PML_1.createArray(pml.nyPML_1 - 1, 0.0);
    //    pml.sigh_y_PML_1.createArray(pml.nyPML_1 - 1, 0.0);
    //    pml.kappah_y_PML_1.createArray(pml.nyPML_1 - 1, 0.0);
    //
    //    pml.be_y_2.createArray(pml.nyPML_2, 0.0);
    //    pml.ce_y_2.createArray(pml.nyPML_2, 0.0);
    //    pml.alphae_y_PML_2.createArray(pml.nyPML_2, 0.0);
    //    pml.sige_y_PML_2.createArray(pml.nyPML_2, 0.0);
    //    pml.kappae_y_PML_2.createArray(pml.nyPML_2, 0.0);
    //
    //    pml.bh_y_2.createArray(pml.nyPML_2 - 1, 0.0);
    //    pml.ch_y_2.createArray(pml.nyPML_2 - 1, 0.0);
    //    pml.alphah_y_PML_2.createArray(pml.nyPML_2 - 1, 0.0);
    //    pml.sigh_y_PML_2.createArray(pml.nyPML_2 - 1, 0.0);
    //    pml.kappah_y_PML_2.createArray(pml.nyPML_2 - 1, 0.0);
    //
    //    pml.be_z_1.createArray(pml.nzPML_1, 0.0);
    //    pml.ce_z_1.createArray(pml.nzPML_1, 0.0);
    //    pml.alphae_z_PML_1.createArray(pml.nzPML_1, 0.0);
    //    pml.sige_z_PML_1.createArray(pml.nzPML_1, 0.0);
    //    pml.kappae_z_PML_1.createArray(pml.nzPML_1, 0.0);
    //
    //    pml.bh_z_1.createArray(pml.nzPML_1 - 1, 0.0);
    //    pml.ch_z_1.createArray(pml.nzPML_1 - 1, 0.0);
    //    pml.alphah_z_PML_1.createArray(pml.nzPML_1 - 1, 0.0);
    //    pml.sigh_z_PML_1.createArray(pml.nzPML_1 - 1, 0.0);
    //    pml.kappah_z_PML_1.createArray(pml.nzPML_1 - 1, 0.0);
    //
    //    pml.be_z_2.createArray(pml.nzPML_2, 0.0);
    //    pml.ce_z_2.createArray(pml.nzPML_2, 0.0);
    //    pml.alphae_z_PML_2.createArray(pml.nzPML_2, 0.0);
    //    pml.sige_z_PML_2.createArray(pml.nzPML_2, 0.0);
    //    pml.kappae_z_PML_2.createArray(pml.nzPML_2, 0.0);
    //
    //    pml.bh_z_2.createArray(pml.nzPML_2 - 1, 0.0);
    //    pml.ch_z_2.createArray(pml.nzPML_2 - 1, 0.0);
    //    pml.alphah_z_PML_2.createArray(pml.nzPML_2 - 1, 0.0);
    //    pml.sigh_z_PML_2.createArray(pml.nzPML_2 - 1, 0.0);
    //    pml.kappah_z_PML_2.createArray(pml.nzPML_2 - 1, 0.0);

    pml.createCPMLArray();
    //    pml.den_ex.createArray(Imax - 1, 0.0);
    //    pml.den_hx.createArray(Imax - 1, 0.0);
    //    pml.den_ey.createArray(Jmax - 1, 0.0);
    //    pml.den_hy.createArray(Jmax - 1, 0.0);
    //    pml.den_ez.createArray(Kmax - 1, 0.0);
    //    pml.den_hz.createArray(Kmax - 1, 0.0);

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void setUp() {

    //Time step
    dt = 0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) +
            1.0 / (dz * dz)));
    //delay
    tO = 4.0 * tw;

    //  Specify the dipole size 
    istart = 24;
    iend = 26;
    jstart = 55;
    jend = 71;
    kstart = 11;
    kend = 13;

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
    //        pml.sig_x_max = 0.75 * (0.8 * (m + 1) / (dx * sqrt(mu_0 / (eps_0 * epsR))));
    //        pml.sig_y_max = 0.75 * (0.8 * (m + 1) / (dy * sqrt(mu_0 / (eps_0 * epsR))));
    //        pml.sig_z_max = 0.75 * (0.8 * (m + 1) / (dz * sqrt(mu_0 / (eps_0 * epsR))));
    //        pml.alpha_x_max = 0.24;
    //        pml.alpha_y_max = pml.alpha_x_max;
    //        pml.alpha_z_max = pml.alpha_x_max;
    //        pml.kappa_x_max = 15.0;
    //        pml.kappa_y_max = pml.kappa_x_max;
    //        pml.kappa_z_max = pml.kappa_x_max;
    printf("\nTIme step = %e", dt);
    printf("\n Number of steps = %d", nMax);
    printf("\n Total Simulation time = %e Seconds", nMax * dt);

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SET CPML PARAMETERS IN EACH DIRECTION
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void initializeCPML() {
    pml.initCPML(dt, dx, dy, dz);
    //    for (i = 0; i < pml.nxPML_1; ++i) {
    //
    //        pml.sige_x_PML_1.p[i] = pml.sig_x_max * pow(((pml.nxPML_1 - 1 - i)
    //                / (pml.nxPML_1 - 1.0)), m);
    //        pml.alphae_x_PML_1.p[i] = pml.alpha_x_max * pow(((float) i
    //                / (pml.nxPML_1 - 1.0)), ma);
    //        pml.kappae_x_PML_1.p[i] = 1.0 + (pml.kappa_x_max - 1.0) *
    //                pow((pml.nxPML_1 - 1 - i) / (pml.nxPML_1 - 1.0), m);
    //        pml.be_x_1.p[i] = exp(-(pml.sige_x_PML_1.p[i] / pml.kappae_x_PML_1.p[i] +
    //                pml.alphae_x_PML_1.p[i]) * dt / eps_0);
    //
    //        if ((pml.sige_x_PML_1.p[i] == 0.0) &&
    //                (pml.alphae_x_PML_1.p[i] == 0.0) && (i == pml.nxPML_1 - 1)) {
    //            pml.ce_x_1.p[i] = 0.0;
    //
    //        } else {
    //            pml.ce_x_1.p[i] = pml.sige_x_PML_1.p[i] * (pml.be_x_1.p[i] - 1.0) /
    //                    (pml.sige_x_PML_1.p[i] + pml.kappae_x_PML_1.p[i] * pml.alphae_x_PML_1.p[i])
    //                    / pml.kappae_x_PML_1.p[i];
    //        }
    //    }
    //
    //    for (i = 0; i < pml.nxPML_1 - 1; ++i) {
    //
    //        pml.sigh_x_PML_1.p[i] = pml.sig_x_max * pow(((pml.nxPML_1 - 1 - i - 0.5)
    //                / (pml.nxPML_1 - 1.0)), m);
    //        pml.alphah_x_PML_1.p[i] = pml.alpha_x_max * pow(((i + 1 - 0.5)
    //                / (pml.nxPML_1 - 1.0)), ma);
    //        pml.kappah_x_PML_1.p[i] = 1.0 + (pml.kappa_x_max - 1.0) *
    //                pow(((pml.nxPML_1 - 1 - i - 0.5) / (pml.nxPML_1 - 1.0)), m);
    //        pml.bh_x_1.p[i] = exp(-(pml.sigh_x_PML_1.p[i] / pml.kappah_x_PML_1.p[i] +
    //                pml.alphah_x_PML_1.p[i]) * dt / eps_0);
    //        pml.ch_x_1.p[i] = pml.sigh_x_PML_1.p[i] * (pml.bh_x_1.p[i] - 1.0) /
    //                (pml.sigh_x_PML_1.p[i] + pml.kappah_x_PML_1.p[i] * pml.alphah_x_PML_1.p[i])
    //                / pml.kappah_x_PML_1.p[i];
    //    }
    //
    //    for (i = 0; i < pml.nxPML_2; ++i) {
    //
    //        pml.sige_x_PML_2.p[i] = pml.sig_x_max * pow(((pml.nxPML_2 - 1 - i)
    //                / (pml.nxPML_2 - 1.0)), m);
    //        pml.alphae_x_PML_2.p[i] = pml.alpha_x_max * pow(((float) i
    //                / (pml.nxPML_2 - 1.0)), ma);
    //        pml.kappae_x_PML_2.p[i] = 1.0 + (pml.kappa_x_max - 1.0) *
    //                pow((pml.nxPML_2 - 1 - i) / (pml.nxPML_2 - 1.0), m);
    //        pml.be_x_2.p[i] = exp(-(pml.sige_x_PML_2.p[i] / pml.kappae_x_PML_2.p[i] +
    //                pml.alphae_x_PML_2.p[i]) * dt / eps_0);
    //
    //        if ((pml.sige_x_PML_2.p[i] == 0.0) &&
    //                (pml.alphae_x_PML_2.p[i] == 0.0) && (i == pml.nxPML_2 - 1)) {
    //
    //            pml.ce_x_2.p[i] = 0.0;
    //
    //        } else {
    //
    //            pml.ce_x_2.p[i] = pml.sige_x_PML_2.p[i] * (pml.be_x_2.p[i] - 1.0) /
    //                    (pml.sige_x_PML_2.p[i] + pml.kappae_x_PML_2.p[i] * pml.alphae_x_PML_2.p[i])
    //                    / pml.kappae_x_PML_2.p[i];
    //        }
    //
    //    }
    //
    //    for (i = 0; i < pml.nxPML_2 - 1; ++i) {
    //
    //        pml.sigh_x_PML_2.p[i] = pml.sig_x_max * pow(((pml.nxPML_2 - 1 - i - 0.5)
    //                / (pml.nxPML_2 - 1.0)), m);
    //        pml.alphah_x_PML_2.p[i] = pml.alpha_x_max * pow(((i + 1 - 0.5)
    //                / (pml.nxPML_2 - 1.0)), ma);
    //        pml.kappah_x_PML_2.p[i] = 1.0 + (pml.kappa_x_max - 1.0) *
    //                pow(((pml.nxPML_2 - 1 - i - 0.5) / (pml.nxPML_2 - 1.0)), m);
    //        pml.bh_x_2.p[i] = exp(-(pml.sigh_x_PML_2.p[i] / pml.kappah_x_PML_2.p[i] +
    //                pml.alphah_x_PML_2.p[i]) * dt / eps_0);
    //        pml.ch_x_2.p[i] = pml.sigh_x_PML_2.p[i] * (pml.bh_x_2.p[i] - 1.0) /
    //                (pml.sigh_x_PML_2.p[i] + pml.kappah_x_PML_2.p[i] * pml.alphah_x_PML_2.p[i])
    //                / pml.kappah_x_PML_2.p[i];
    //    }
    //
    //    for (j = 0; j < pml.nyPML_1; ++j) {
    //
    //        pml.sige_y_PML_1.p[j] = pml.sig_y_max * pow(((pml.nyPML_1 - 1 - j)
    //                / (pml.nyPML_1 - 1.0)), m);
    //        pml.alphae_y_PML_1.p[j] = pml.alpha_y_max * pow(((float) j
    //                / (pml.nyPML_1 - 1.0)), ma);
    //        pml.kappae_y_PML_1.p[j] = 1.0 + (pml.kappa_y_max - 1.0) *
    //                pow((pml.nyPML_1 - 1 - j) / (pml.nyPML_1 - 1.0), m);
    //        pml.be_y_1.p[j] = exp(-(pml.sige_y_PML_1.p[j] / pml.kappae_y_PML_1.p[j] +
    //                pml.alphae_y_PML_1.p[j]) * dt / eps_0);
    //
    //        if ((pml.sige_y_PML_1.p[j] == 0.0) &&
    //                (pml.alphae_y_PML_1.p[j] == 0.0) && (j == pml.nyPML_1 - 1)) {
    //            pml.ce_y_1.p[j] = 0.0;
    //
    //        } else {
    //            pml.ce_y_1.p[j] = pml.sige_y_PML_1.p[j] * (pml.be_y_1.p[j] - 1.0) /
    //                    (pml.sige_y_PML_1.p[j] + pml.kappae_y_PML_1.p[j] * pml.alphae_y_PML_1.p[j])
    //                    / pml.kappae_y_PML_1.p[j];
    //        }
    //    }
    //
    //    for (j = 0; j < pml.nyPML_1 - 1; ++j) {
    //
    //        pml.sigh_y_PML_1.p[j] = pml.sig_y_max * pow(((pml.nyPML_1 - 1 - j - 0.5)
    //                / (pml.nyPML_1 - 1.0)), m);
    //        pml.alphah_y_PML_1.p[j] = pml.alpha_y_max * pow(((j + 1 - 0.5)
    //                / (pml.nyPML_1 - 1.0)), ma);
    //        pml.kappah_y_PML_1.p[j] = 1.0 + (pml.kappa_y_max - 1.0) *
    //                pow(((pml.nyPML_1 - 1 - j - 0.5) / (pml.nyPML_1 - 1.0)), m);
    //        pml.bh_y_1.p[j] = exp(-(pml.sigh_y_PML_1.p[j] / pml.kappah_y_PML_1.p[j] +
    //                pml.alphah_y_PML_1.p[j]) * dt / eps_0);
    //        pml.ch_y_1.p[j] = pml.sigh_y_PML_1.p[j] * (pml.bh_y_1.p[j] - 1.0) /
    //                (pml.sigh_y_PML_1.p[j] + pml.kappah_y_PML_1.p[j] * pml.alphah_y_PML_1.p[j])
    //                / pml.kappah_y_PML_1.p[j];
    //    }
    //
    //    for (j = 0; j < pml.nyPML_2; ++j) {
    //
    //        pml.sige_y_PML_2.p[j] = pml.sig_y_max * pow(((pml.nyPML_2 - 1 - j)
    //                / (pml.nyPML_2 - 1.0)), m);
    //        pml.alphae_y_PML_2.p[j] = pml.alpha_y_max * pow(((float) j
    //                / (pml.nyPML_2 - 1.0)), ma);
    //        pml.kappae_y_PML_2.p[j] = 1.0 + (pml.kappa_y_max - 1.0) *
    //                pow((pml.nyPML_2 - 1 - j) / (pml.nyPML_2 - 1.0), m);
    //        pml.be_y_2.p[j] = exp(-(pml.sige_y_PML_2.p[j] / pml.kappae_y_PML_2.p[j] +
    //                pml.alphae_y_PML_2.p[j]) * dt / eps_0);
    //
    //        if ((pml.sige_y_PML_2.p[j] == 0.0) &&
    //                (pml.alphae_y_PML_2.p[j] == 0.0) && (j == pml.nyPML_2 - 1)) {
    //
    //            pml.ce_y_2.p[j] = 0.0;
    //
    //        } else {
    //
    //            pml.ce_y_2.p[j] = pml.sige_y_PML_2.p[j] * (pml.be_y_2.p[j] - 1.0) /
    //                    (pml.sige_y_PML_2.p[j] + pml.kappae_y_PML_2.p[j] * pml.alphae_y_PML_2.p[j])
    //                    / pml.kappae_y_PML_2.p[j];
    //        }
    //    }
    //
    //    for (j = 0; j < pml.nyPML_2 - 1; ++j) {
    //
    //        pml.sigh_y_PML_2.p[j] = pml.sig_y_max * pow(((pml.nyPML_2 - 1 - j - 0.5)
    //                / (pml.nyPML_2 - 1.0)), m);
    //        pml.alphah_y_PML_2.p[j] = pml.alpha_y_max * pow(((j + 1 - 0.5)
    //                / (pml.nyPML_2 - 1.0)), ma);
    //        pml.kappah_y_PML_2.p[j] = 1.0 + (pml.kappa_y_max - 1.0) *
    //                pow(((pml.nyPML_2 - 1 - j - 0.5) / (pml.nyPML_2 - 1.0)), m);
    //        pml.bh_y_2.p[j] = exp(-(pml.sigh_y_PML_2.p[j] / pml.kappah_y_PML_2.p[j] +
    //                pml.alphah_y_PML_2.p[j]) * dt / eps_0);
    //        pml.ch_y_2.p[j] = pml.sigh_y_PML_2.p[j] * (pml.bh_y_2.p[j] - 1.0) /
    //                (pml.sigh_y_PML_2.p[j] + pml.kappah_y_PML_2.p[j] * pml.alphah_y_PML_2.p[j])
    //                / pml.kappah_y_PML_2.p[j];
    //    }
    //
    //    for (k = 0; k < pml.nzPML_1; ++k) {
    //
    //        pml.sige_z_PML_1.p[k] = pml.sig_z_max * pow(((pml.nzPML_1 - 1 - k)
    //                / (pml.nzPML_1 - 1.0)), m);
    //        pml.alphae_z_PML_1.p[k] = pml.alpha_z_max * pow(((float) k
    //                / (pml.nzPML_1 - 1.0)), ma);
    //        pml.kappae_z_PML_1.p[k] = 1.0 + (pml.kappa_z_max - 1.0) *
    //                pow((pml.nzPML_1 - 1 - k) / (pml.nzPML_1 - 1.0), m);
    //        pml.be_z_1.p[k] = exp(-(pml.sige_z_PML_1.p[k] / pml.kappae_z_PML_1.p[k] +
    //                pml.alphae_z_PML_1.p[k]) * dt / eps_0);
    //
    //        if ((pml.sige_z_PML_1.p[k] == 0.0) &&
    //                (pml.alphae_z_PML_1.p[k] == 0.0) && (k == pml.nzPML_1 - 1)) {
    //
    //            pml.ce_z_1.p[k] = 0.0;
    //
    //        } else {
    //
    //            pml.ce_z_1.p[k] = pml.sige_z_PML_1.p[k] * (pml.be_z_1.p[k] - 1.0) /
    //                    (pml.sige_z_PML_1.p[k] + pml.kappae_z_PML_1.p[k] * pml.alphae_z_PML_1.p[k])
    //                    / pml.kappae_z_PML_1.p[k];
    //        }
    //    }
    //
    //    for (k = 0; k < pml.nzPML_1 - 1; ++k) {
    //
    //        pml.sigh_z_PML_1.p[k] = pml.sig_z_max * pow(((pml.nzPML_1 - 1 - k - 0.5)
    //                / (pml.nzPML_1 - 1.0)), m);
    //        pml.alphah_z_PML_1.p[k] = pml.alpha_z_max * pow(((k + 1 - 0.5)
    //                / (pml.nzPML_1 - 1.0)), ma);
    //        pml.kappah_z_PML_1.p[k] = 1.0 + (pml.kappa_z_max - 1.0) *
    //                pow(((pml.nzPML_1 - 1 - k - 0.5) / (pml.nzPML_1 - 1.0)), m);
    //        pml.bh_z_1.p[k] = exp(-(pml.sigh_z_PML_1.p[k] / pml.kappah_z_PML_1.p[k] +
    //                pml.alphah_z_PML_1.p[k]) * dt / eps_0);
    //        pml.ch_z_1.p[k] = pml.sigh_z_PML_1.p[k] * (pml.bh_z_1.p[k] - 1.0) /
    //                (pml.sigh_z_PML_1.p[k] + pml.kappah_z_PML_1.p[k] * pml.alphah_z_PML_1.p[k])
    //                / pml.kappah_z_PML_1.p[k];
    //    }
    //
    //    for (k = 0; k < pml.nzPML_2; ++k) {
    //
    //        pml.sige_z_PML_2.p[k] = pml.sig_z_max * pow(((pml.nzPML_2 - 1 - k)
    //                / (pml.nzPML_2 - 1.0)), m);
    //        pml.alphae_z_PML_2.p[k] = pml.alpha_z_max * pow(((float) k
    //                / (pml.nzPML_2 - 1.0)), ma);
    //        pml.kappae_z_PML_2.p[k] = 1.0 + (pml.kappa_z_max - 1.0) *
    //                pow((pml.nzPML_2 - 1 - k) / (pml.nzPML_2 - 1.0), m);
    //        pml.be_z_2.p[k] = exp(-(pml.sige_z_PML_2.p[k] / pml.kappae_z_PML_2.p[k] +
    //                pml.alphae_z_PML_2.p[k]) * dt / eps_0);
    //
    //        if ((pml.sige_z_PML_2.p[k] == 0.0) &&
    //                (pml.alphae_z_PML_2.p[k] == 0.0) && (k == pml.nzPML_2 - 1)) {
    //
    //            pml.ce_z_2.p[k] = 0.0;
    //
    //        } else {
    //
    //            pml.ce_z_2.p[k] = pml.sige_z_PML_2.p[k] * (pml.be_z_2.p[k] - 1.0) /
    //                    (pml.sige_z_PML_2.p[k] + pml.kappae_z_PML_2.p[k] * pml.alphae_z_PML_2.p[k])
    //                    / pml.kappae_z_PML_2.p[k];
    //        }
    //    }
    //
    //    for (k = 0; k < pml.nzPML_2 - 1; ++k) {
    //
    //        pml.sigh_z_PML_2.p[k] = pml.sig_z_max * pow(((pml.nzPML_2 - 1 - k - 0.5)
    //                / (pml.nzPML_2 - 1.0)), m);
    //        pml.alphah_z_PML_2.p[k] = pml.alpha_z_max * pow(((k + 1 - 0.5)
    //                / (pml.nzPML_2 - 1.0)), ma);
    //        pml.kappah_z_PML_2.p[k] = 1.0 + (pml.kappa_z_max - 1.0) *
    //                pow(((pml.nzPML_2 - 1 - k - 0.5) / (pml.nzPML_2 - 1.0)), m);
    //        pml.bh_z_2.p[k] = exp(-(pml.sigh_z_PML_2.p[k] / pml.kappah_z_PML_2.p[k] +
    //                pml.alphah_z_PML_2.p[k]) * dt / eps_0);
    //        pml.ch_z_2.p[k] = pml.sigh_z_PML_2.p[k] * (pml.bh_z_2.p[k] - 1.0) /
    //                (pml.sigh_z_PML_2.p[k] + pml.kappah_z_PML_2.p[k] * pml.alphah_z_PML_2.p[k])
    //                / pml.kappah_z_PML_2.p[k];
    //    }
    //
    //    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //    //  DENOMINATORS FOR FIELD UPDATES
    //    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //    ii = pml.nxPML_2 - 2;
    //
    //    for (i = 0; i < Imax - 1; ++i) {
    //
    //        if (i < pml.nxPML_1 - 1) {
    //
    //            pml.den_hx.p[i] = 1.0 / (pml.kappah_x_PML_1.p[i] * dx);
    //
    //        } else if (i >= Imax - pml.nxPML_2) {
    //
    //            pml.den_hx.p[i] = 1.0 / (pml.kappah_x_PML_2.p[ii] * dx);
    //            ii = ii - 1;
    //
    //        } else {
    //
    //            pml.den_hx.p[i] = 1.0 / dx;
    //        }
    //    }
    //
    //    jj = pml.nyPML_2 - 2;
    //
    //    for (j = 0; j < Jmax - 1; ++j) {
    //
    //        if (j < pml.nyPML_1 - 1) {
    //
    //            pml.den_hy.p[j] = 1.0 / (pml.kappah_y_PML_1.p[j] * dy);
    //
    //        } else if (j >= Jmax - pml.nyPML_2) {
    //
    //            pml.den_hy.p[j] = 1.0 / (pml.kappah_y_PML_2.p[jj] * dy);
    //            jj = jj - 1;
    //
    //        } else {
    //
    //            pml.den_hy.p[j] = 1.0 / dy;
    //        }
    //    }
    //
    //    kk = pml.nzPML_2 - 2;
    //
    //    for (k = 1; k < Kmax - 1; ++k) {
    //
    //        if (k < pml.nzPML_1) {
    //
    //            pml.den_hz.p[k] = 1.0 / (pml.kappah_z_PML_1.p[k - 1] * dz);
    //
    //        } else if (k >= Kmax - pml.nzPML_2) {
    //
    //            pml.den_hz.p[k] = 1.0 / (pml.kappah_z_PML_2.p[kk] * dz);
    //            kk = kk - 1;
    //
    //        } else {
    //
    //            pml.den_hz.p[k] = 1.0 / dz;
    //        }
    //    }
    //
    //    ii = pml.nxPML_2 - 1;
    //
    //    for (i = 0; i < Imax - 1; ++i) {
    //
    //        if (i < pml.nxPML_1) {
    //
    //            pml.den_ex.p[i] = 1.0 / (pml.kappae_x_PML_1.p[i] * dx);
    //
    //        } else if (i >= Imax - pml.nxPML_2) {
    //
    //            pml.den_ex.p[i] = 1.0 / (pml.kappae_x_PML_2.p[ii] * dx);
    //            ii = ii - 1;
    //
    //        } else {
    //
    //            pml.den_ex.p[i] = 1.0 / dx;
    //        }
    //    }
    //
    //    jj = pml.nyPML_2 - 1;
    //
    //    for (j = 0; j < Jmax - 1; ++j) {
    //
    //        if (j < pml.nyPML_1) {
    //
    //            pml.den_ey.p[j] = 1.0 / (pml.kappae_y_PML_1.p[j] * dy);
    //
    //        } else if (j >= Jmax - pml.nyPML_2) {
    //
    //            pml.den_ey.p[j] = 1.0 / (pml.kappae_y_PML_2.p[jj] * dy);
    //            jj = jj - 1;
    //
    //        } else {
    //
    //            pml.den_ey.p[j] = 1.0 / dy;
    //        }
    //    }
    //
    //    kk = pml.nzPML_2 - 1;
    //
    //    for (k = 0; k < Kmax - 1; ++k) {
    //
    //        if (k < pml.nzPML_1) {
    //
    //            pml.den_ez.p[k] = 1.0 / (pml.kappae_z_PML_1.p[k] * dz);
    //
    //        } else if (k >= Kmax - 1 - pml.nzPML_2) {
    //
    //            pml.den_ez.p[k] = 1.0 / (pml.kappae_z_PML_2.p[kk] * dz);
    //            kk = kk - 1;
    //
    //        } else {
    //
    //            pml.den_ez.p[k] = 1.0 / dz;
    //        }
    //    }
}

void compute() {

    short id;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  BEGIN TIME STEP
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    printf("\nBegin time-stepping...\n");

    for (n = 1; n <= nMax; ++n) {

        printf("Ez at time step %d at (25, 40, 12) :  %f\n", n, Ez.p[25][40][12]);

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
            //            for (i = 0; i < Imax - 1; ++i) {
            //                //...............................................
            //                //  PML for bottom Hx, j-direction
            //                //...............................................
            //                for (j = 0; j < pml.nyPML_1 - 1; ++j) {
            //                    pml.psi_Hxy_1.p[i][j][k] = pml.bh_y_1.p[j] * pml.psi_Hxy_1.p[i][j][k]
            //                            + pml.ch_y_1.p[j] * (Ez.p[i][j][k] - Ez.p[i][j + 1][k]) / dy;
            //                    Hx.p[i][j][k] = Hx.p[i][j][k] + DB * pml.psi_Hxy_1.p[i][j][k];
            //                }
            //                //....................................................
            //                //  PML for top Hx, j-direction
            //                //.....................................................
            //                jj = pml.nyPML_2 - 2;
            //                for (j = Jmax - pml.nyPML_2; j < Jmax - 1; ++j) {
            //                    pml.psi_Hxy_2.p[i][jj][k] = pml.bh_y_2.p[jj] * pml.psi_Hxy_2.p[i][jj][k]
            //                            + pml.ch_y_2.p[jj] * (Ez.p[i][j][k] - Ez.p[i][j + 1][k]) / dy;
            //                    Hx.p[i][j][k] = Hx.p[i][j][k] + DB * pml.psi_Hxy_2.p[i][jj][k];
            //                    jj = jj - 1;
            //                }
            //            }
        }
        pml.updateHxOut(Hx, Ey, DB, dz);

        //        for (i = 0; i < Imax - 1; ++i) {
        //
        //            for (j = 0; j < Jmax - 1; ++j) {
        //                //....................................................
        //                //  PML for bottom Hx, k-direction
        //                //................................................
        //                for (k = 1; k < pml.nzPML_1; ++k) {
        //
        //                    pml.psi_Hxz_1.p[i][j][k - 1] = pml.bh_z_1.p[k - 1] * pml.psi_Hxz_1.p[i][j][k - 1]
        //                            + pml.ch_z_1.p[k - 1] * (Ey.p[i][j][k] - Ey.p[i][j][k - 1]) / dz;
        //                    Hx.p[i][j][k] = Hx.p[i][j][k] + DB * pml.psi_Hxz_1.p[i][j][k - 1];
        //                }
        //                //....................................................
        //                //  PML for top Hx, k-direction
        //                //...............................................
        //                kk = pml.nzPML_2 - 2;
        //                for (k = Kmax - pml.nzPML_2; k < Kmax - 1; ++k) {
        //
        //                    pml.psi_Hxz_2.p[i][j][kk] = pml.bh_z_2.p[kk] * pml.psi_Hxz_2.p[i][j][kk]
        //                            + pml.ch_z_2.p[kk] * (Ey.p[i][j][k] - Ey.p[i][j][k - 1]) / dz;
        //                    Hx.p[i][j][k] = Hx.p[i][j][k] + DB * pml.psi_Hxz_2.p[i][j][kk];
        //                    kk = kk - 1;
        //                }
        //            }
        //        }

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

            //            for (j = 0; j < Jmax - 1; ++j) {
            //                //.......................................................
            //                //  PML for bottom Hy, i-direction
            //                //.......................................................
            //                for (i = 0; i < pml.nxPML_1 - 1; ++i) {
            //
            //                    pml.psi_Hyx_1.p[i][j][k] = pml.bh_x_1.p[i] * pml.psi_Hyx_1.p[i][j][k]
            //                            + pml.ch_x_1.p[i] * (Ez.p[i + 1][j][k] - Ez.p[i][j][k]) / dx;
            //                    Hy.p[i][j][k] = Hy.p[i][j][k] + DB * pml.psi_Hyx_1.p[i][j][k];
            //                }
            //                //.........................................................
            //                //  PML for top Hy, i-direction
            //                //.........................................................
            //                ii = pml.nxPML_2 - 2;
            //                for (i = Imax - pml.nxPML_2; i < Imax - 1; ++i) {
            //
            //                    pml.psi_Hyx_2.p[ii][j][k] = pml.bh_x_2.p[ii] * pml.psi_Hyx_2.p[ii][j][k]
            //                            + pml.ch_x_2.p[ii] * (Ez.p[i + 1][j][k] - Ez.p[i][j][k]) / dx;
            //                    Hy.p[i][j][k] = Hy.p[i][j][k] + DB * pml.psi_Hyx_2.p[ii][j][k];
            //                    ii = ii - 1;
            //                }
            //            }
        }
        pml.updateHyOut(Hy, Ex, DB, dz);
        //        for (i = 0; i < Imax - 1; ++i) {
        //            for (j = 0; j < Jmax - 1; ++j) {
        //                //.......................................................
        //                //  PML for bottom Hy, k-direction
        //                //......................................................
        //                for (k = 1; k < pml.nzPML_1; ++k) {
        //
        //                    pml.psi_Hyz_1.p[i][j][k - 1] = pml.bh_z_1.p[k - 1] * pml.psi_Hyz_1.p[i][j][k - 1]
        //                            + pml.ch_z_1.p[k - 1] * (Ex.p[i][j][k - 1] - Ex.p[i][j][k]) / dz;
        //                    Hy.p[i][j][k] = Hy.p[i][j][k] + DB * pml.psi_Hyz_1.p[i][j][k - 1];
        //                }
        //                //.......................................................
        //                //  PML for top Hy, k-direction
        //                //.........................................................
        //                kk = pml.nzPML_2 - 2;
        //                for (k = Kmax - pml.nzPML_2; k < Kmax - 1; ++k) {
        //                    pml.psi_Hyz_2.p[i][j][kk] = pml.bh_z_2.p[kk] * pml.psi_Hyz_2.p[i][j][kk]
        //                            + pml.ch_z_2.p[kk] * (Ex.p[i][j][k - 1] - Ex.p[i][j][k]) / dz;
        //                    Hy.p[i][j][k] = Hy.p[i][j][k] + DB * pml.psi_Hyz_2.p[i][j][kk];
        //                    kk = kk - 1;
        //                }
        //            }
        //        }

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

            //            for (j = 0; j < Jmax - 1; ++j) {
            //                //..........................................................
            //                //  PML for bottom Hz, x-direction
            //                //..........................................................
            //                for (i = 0; i < pml.nxPML_1 - 1; ++i) {
            //
            //                    pml.psi_Hzx_1.p[i][j][k] = pml.bh_x_1.p[i] * pml.psi_Hzx_1.p[i][j][k]
            //                            + pml.ch_x_1.p[i] * (Ey.p[i][j][k] - Ey.p[i + 1][j][k]) / dx;
            //                    Hz.p[i][j][k] = Hz.p[i][j][k] + DB * pml.psi_Hzx_1.p[i][j][k];
            //                }
            //                //..........................................................
            //                //  PML for top Hz, x-direction
            //                //..........................................................
            //                ii = pml.nxPML_2 - 2;
            //                for (i = Imax - pml.nxPML_2; i < Imax - 1; ++i) {
            //
            //                    pml.psi_Hzx_2.p[ii][j][k] = pml.bh_x_2.p[ii] * pml.psi_Hzx_2.p[ii][j][k]
            //                            + pml.ch_x_2.p[ii] * (Ey.p[i][j][k] - Ey.p[i + 1][j][k]) / dx;
            //                    Hz.p[i][j][k] = Hz.p[i][j][k] + DB * pml.psi_Hzx_2.p[ii][j][k];
            //                    ii = ii - 1;
            //                }
            //            }
            //
            //            for (i = 0; i < Imax - 1; ++i) {
            //                //........................................................
            //                //  PML for bottom Hz, y-direction
            //                //.........................................................
            //                for (j = 0; j < pml.nyPML_1 - 1; ++j) {
            //                    pml.psi_Hzy_1.p[i][j][k] = pml.bh_y_1.p[j] * pml.psi_Hzy_1.p[i][j][k]
            //                            + pml.ch_y_1.p[j] * (Ex.p[i][j + 1][k] - Ex.p[i][j][k]) / dy;
            //                    Hz.p[i][j][k] = Hz.p[i][j][k] + DB * pml.psi_Hzy_1.p[i][j][k];
            //
            //                }
            //                //.........................................................
            //                //  PML for top Hz, y-direction
            //                //..........................................................
            //                jj = pml.nyPML_2 - 2;
            //                for (j = Jmax - pml.nyPML_2; j < Jmax - 1; ++j) {
            //
            //                    pml.psi_Hzy_2.p[i][jj][k] = pml.bh_y_2.p[jj] * pml.psi_Hzy_2.p[i][jj][k]
            //                            + pml.ch_y_2.p[jj] * (Ex.p[i][j + 1][k] - Ex.p[i][j][k]) / dy;
            //                    Hz.p[i][j][k] = Hz.p[i][j][k] + DB * pml.psi_Hzy_2.p[i][jj][k];
            //                    jj = jj - 1;
            //                }
            //            }
        }

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
            pml.updateExIn(k, Ex, Hz, ID1, CB, dy);
            //            for (i = 0; i < Imax - 1; ++i) {
            //                //..............................................................
            //                //  PML for bottom Ex, j-direction
            //                //..............................................................
            //                for (j = 1; j < pml.nyPML_1; ++j) {
            //
            //                    id = ID1.p[i][j][k];
            //                    pml.psi_Exy_1.p[i][j][k] = pml.be_y_1.p[j] * pml.psi_Exy_1.p[i][j][k]
            //                            + pml.ce_y_1.p[j] * (Hz.p[i][j][k] - Hz.p[i][j - 1][k]) / dy;
            //                    Ex.p[i][j][k] = Ex.p[i][j][k] + CB[id] * pml.psi_Exy_1.p[i][j][k];
            //                }
            //                //.............................................................
            //                //  PML for top Ex, j-direction
            //                //.............................................................
            //                jj = pml.nyPML_2 - 1;
            //                for (j = Jmax - pml.nyPML_2; j < Jmax - 1; ++j) {
            //
            //                    id = ID1.p[i][j][k];
            //                    pml.psi_Exy_2.p[i][jj][k] = pml.be_y_2.p[jj] * pml.psi_Exy_2.p[i][jj][k]
            //                            + pml.ce_y_2.p[jj] * (Hz.p[i][j][k] - Hz.p[i][j - 1][k]) / dy;
            //                    Ex.p[i][j][k] = Ex.p[i][j][k] + CB[id] * pml.psi_Exy_2.p[i][jj][k];
            //                    jj = jj - 1;
            //                }
            //            }
        }

        pml.updateExOut(Ex, Hy, ID1, CB, dz);
        //        for (i = 0; i < Imax - 1; ++i) {
        //
        //            for (j = 1; j < Jmax - 1; ++j) {
        //                //.............................................................
        //                //  PML for bottom Ex, k-direction
        //                //.............................................................
        //                for (k = 0; k < pml.nzPML_1; ++k) {
        //
        //                    id = ID1.p[i][j][k];
        //                    pml.psi_Exz_1.p[i][j][k] = pml.be_z_1.p[k] * pml.psi_Exz_1.p[i][j][k]
        //                            + pml.ce_z_1.p[k] * (Hy.p[i][j][k] - Hy.p[i][j][k + 1]) / dz;
        //                    Ex.p[i][j][k] = Ex.p[i][j][k] + CB[id] * pml.psi_Exz_1.p[i][j][k];
        //                }
        //                //..............................................................
        //                //  PML for top Ex, k-direction
        //                //..............................................................
        //                kk = pml.nzPML_2 - 1;
        //                for (k = Kmax - pml.nzPML_2 - 1; k < Kmax - 1; ++k) {
        //
        //                    id = ID1.p[i][j][k];
        //                    pml.psi_Exz_2.p[i][j][kk] = pml.be_z_2.p[kk] * pml.psi_Exz_2.p[i][j][kk]
        //                            + pml.ce_z_2.p[kk] * (Hy.p[i][j][k] - Hy.p[i][j][k + 1]) / dz;
        //                    Ex.p[i][j][k] = Ex.p[i][j][k] + CB[id] * pml.psi_Exz_2.p[i][j][kk];
        //                    kk = kk - 1;
        //                }
        //            }
        //        }

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

            pml.updateEyIn(k, Ey, Hz, ID2, CB, dx);
            //            for (j = 0; j < Jmax - 1; ++j) {
            //                //...........................................................
            //                //  PML for bottom Ey, i-direction
            //                //...........................................................
            //                for (i = 1; i < pml.nxPML_1; ++i) {
            //                    id = ID2.p[i][j][k];
            //                    pml.psi_Eyx_1.p[i][j][k] = pml.be_x_1.p[i] * pml.psi_Eyx_1.p[i][j][k]
            //                            + pml.ce_x_1.p[i] * (Hz.p[i - 1][j][k] - Hz.p[i][j][k]) / dx;
            //                    Ey.p[i][j][k] = Ey.p[i][j][k] + CB[id] * pml.psi_Eyx_1.p[i][j][k];
            //                }
            //                //............................................................
            //                //  PML for top Ey, i-direction
            //                //............................................................
            //                ii = pml.nxPML_2 - 1;
            //                for (i = Imax - pml.nxPML_2; i < Imax - 1; ++i) {
            //                    id = ID2.p[i][j][k];
            //                    pml.psi_Eyx_2.p[ii][j][k] = pml.be_x_2.p[ii] * pml.psi_Eyx_2.p[ii][j][k]
            //                            + pml.ce_x_2.p[ii] * (Hz.p[i - 1][j][k] - Hz.p[i][j][k]) / dx;
            //                    Ey.p[i][j][k] = Ey.p[i][j][k] + CB[id] * pml.psi_Eyx_2.p[ii][j][k];
            //                    ii = ii - 1;
            //                }
            //            }
        }
        pml.updateEyOut(Ey, Hx, ID2, CB, dz);
        //
        //        for (i = 1; i < Imax - 1; ++i) {
        //            for (j = 0; j < Jmax - 1; ++j) {
        //                //...........................................................
        //                //  PML for bottom Ey, k-direction
        //                //...........................................................
        //                for (k = 0; k < pml.nzPML_1; ++k) {
        //                    id = ID2.p[i][j][k];
        //                    pml.psi_Eyz_1.p[i][j][k] = pml.be_z_1.p[k] * pml.psi_Eyz_1.p[i][j][k]
        //                            + pml.ce_z_1.p[k] * (Hx.p[i][j][k + 1] - Hx.p[i][j][k]) / dz;
        //                    Ey.p[i][j][k] = Ey.p[i][j][k] + CB[id] * pml.psi_Eyz_1.p[i][j][k];
        //                }
        //                //...........................................................
        //                //  PML for top Ey, k-direction
        //                //............................................................
        //                kk = pml.nzPML_2 - 1;
        //                for (k = Kmax - pml.nzPML_2 - 1; k < Kmax - 1; ++k) {
        //
        //                    id = ID2.p[i][j][k];
        //                    pml.psi_Eyz_2.p[i][j][kk] = pml.be_z_2.p[kk] * pml.psi_Eyz_2.p[i][j][kk]
        //                            + pml.ce_z_2.p[kk] * (Hx.p[i][j][k + 1] - Hx.p[i][j][k]) / dz;
        //                    Ey.p[i][j][k] = Ey.p[i][j][k] + CB[id] * pml.psi_Eyz_2.p[i][j][kk];
        //                    kk = kk - 1;
        //                }
        //            }
        //        }

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
            pml.updateEz(k, Ez, Hx, Hy, ID3, CB, dx, dy);

            //            for (j = 1; j < Jmax - 1; ++j) {
            //                //............................................................
            //                //  PML for bottom Ez, x-direction
            //                //.............................................................
            //                for (i = 1; i < pml.nxPML_1; ++i) {
            //                    id = ID3.p[i][j][k];
            //                    pml.psi_Ezx_1.p[i][j][k] = pml.be_x_1.p[i] * pml.psi_Ezx_1.p[i][j][k]
            //                            + pml.ce_x_1.p[i] * (Hy.p[i][j][k] - Hy.p[i - 1][j][k]) / dx;
            //                    Ez.p[i][j][k] = Ez.p[i][j][k] + CB[id] * pml.psi_Ezx_1.p[i][j][k];
            //                }
            //                //............................................................
            //                //  PML for top Ez, x-direction
            //                //............................................................
            //                ii = pml.nxPML_2 - 1;
            //                for (i = Imax - pml.nxPML_2; i < Imax - 1; ++i) {
            //                    id = ID3.p[i][j][k];
            //                    pml.psi_Ezx_2.p[ii][j][k] = pml.be_x_2.p[ii] * pml.psi_Ezx_2.p[ii][j][k]
            //                            + pml.ce_x_2.p[ii] * (Hy.p[i][j][k] - Hy.p[i - 1][j][k]) / dx;
            //                    Ez.p[i][j][k] = Ez.p[i][j][k] + CB[id] * pml.psi_Ezx_2.p[ii][j][k];
            //                    ii = ii - 1;
            //                }
            //            }
            //
            //            for (i = 1; i < Imax - 1; ++i) {
            //                //..........................................................
            //                //  PML for bottom Ez, y-direction
            //                //..........................................................
            //                for (j = 1; j < pml.nyPML_1; ++j) {
            //                    id = ID3.p[i][j][k];
            //                    pml.psi_Ezy_1.p[i][j][k] = pml.be_y_1.p[j] * pml.psi_Ezy_1.p[i][j][k]
            //                            + pml.ce_y_1.p[j] * (Hx.p[i][j - 1][k] - Hx.p[i][j][k]) / dy;
            //                    Ez.p[i][j][k] = Ez.p[i][j][k] + CB[id] * pml.psi_Ezy_1.p[i][j][k];
            //                }
            //                //............................................................
            //                //  PML for top Ez, y-direction
            //                //............................................................
            //                jj = pml.nyPML_2 - 1;
            //                for (j = Jmax - pml.nyPML_2; j < Jmax - 1; ++j) {
            //                    id = ID3.p[i][j][k];
            //                    pml.psi_Ezy_2.p[i][jj][k] = pml.be_y_2.p[jj] * pml.psi_Ezy_2.p[i][jj][k]
            //                            + pml.ce_y_2.p[jj] * (Hx.p[i][j - 1][k] - Hx.p[i][j][k]) / dy;
            //                    Ez.p[i][j][k] = Ez.p[i][j][k] + CB[id] * pml.psi_Ezy_2.p[i][jj][k];
            //                    jj = jj - 1;
            //                }
            //            }
        }


        //-----------------------------------------------------------
        //   Apply a point source (Soft)
        //-----------------------------------------------------------
        i = 25;
        j = 63;
        k = 12;
        source = amp * -2.0 * ((n * dt - tO) / tw)
                * exp(-pow(((n * dt - tO) / tw), 2)); //Differentiated Gaussian pulse

        Ez.p[i][j][k] = Ez.p[i][j][k] - CB[ID3.p[i][j][k]] * source;


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

void buildObject() {

    //buildSphere();
    buildDipole();
}

//Builds a sphere (Sample code - NOt used in this program)

void buildSphere() {

    double dist; //distance
    double rad = 8; //(double)Imax / 5.0; // sphere radius
    double sc = (double) Imax / 2.0; //sphere centre
    //double rad2 = 0.3; //(double)Imax / 5.0 - 3.0; // sphere radius

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

void buildDipole() {

    int centre = (jstart + jend) / 2;

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

void yeeCube(int I, int J, int K, short mType) {

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

void writeField(int iteration) {

    FILE *ptr;
    char step[10];
    char fileBaseName[] = "E_Field_";
    sprintf(step, "%d", iteration);
    strcat(fileBaseName, step);
    strcat(fileBaseName, ".txt");

    ptr = fopen(fileBaseName, "wt");

    for (i = 0; i < Imax - 1; i++) {

        for (j = 0; j < Jmax - 1; j++) {
            // |E|
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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// END OF PROGRAM CPMLFDTD3D
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
