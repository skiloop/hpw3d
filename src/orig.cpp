
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

//  Specify the CPML Thickness in Each Direction (Value of Zero
//  Corresponds to No PML, and the Grid is Terminated with a PEC)
// PML thickness in each direction
int nxPML_1, nxPML_2, nyPML_1;
int nyPML_2, nzPML_1, nzPML_2;

//  Specify the CPML Order and Other Parameters:
int m = 3, ma = 1;

double sig_x_max;
double sig_y_max;
double sig_z_max;
double alpha_x_max;
double alpha_y_max;
double alpha_z_max;
double kappa_x_max;
double kappa_y_max;
double kappa_z_max;

//Loop indices
int i, j, ii, jj, k, kk, n;


// H & E Field components
float ***Hx;
float ***Hy;
float ***Hz;
float ***Ex;
float ***Ey;
float ***Ez;

short ***ID1; //medium definition array for Ex
short ***ID2; //medium definition array for Ey
short ***ID3; //medium definition array for Ez

//  CPML components (Taflove 3rd Edition, Chapter 7)
float ***psi_Ezx_1;
float ***psi_Ezx_2;
float ***psi_Hyx_1;
float ***psi_Hyx_2;
float ***psi_Ezy_1;
float ***psi_Ezy_2;
float ***psi_Hxy_1;
float ***psi_Hxy_2;
float ***psi_Hxz_1;
float ***psi_Hxz_2;
float ***psi_Hyz_1;
float ***psi_Hyz_2;
float ***psi_Exz_1;
float ***psi_Exz_2;
float ***psi_Eyz_1;
float ***psi_Eyz_2;
float ***psi_Hzx_1;
float ***psi_Eyx_1;
float ***psi_Hzx_2;
float ***psi_Eyx_2;
float ***psi_Hzy_1;
float ***psi_Exy_1;
float ***psi_Hzy_2;
float ***psi_Exy_2;

float *be_x_1, *ce_x_1, *alphae_x_PML_1, *sige_x_PML_1, *kappae_x_PML_1;
float *bh_x_1, *ch_x_1, *alphah_x_PML_1, *sigh_x_PML_1, *kappah_x_PML_1;
float *be_x_2, *ce_x_2, *alphae_x_PML_2, *sige_x_PML_2, *kappae_x_PML_2;
float *bh_x_2, *ch_x_2, *alphah_x_PML_2, *sigh_x_PML_2, *kappah_x_PML_2;
float *be_y_1, *ce_y_1, *alphae_y_PML_1, *sige_y_PML_1, *kappae_y_PML_1;
float *bh_y_1, *ch_y_1, *alphah_y_PML_1, *sigh_y_PML_1, *kappah_y_PML_1;
float *be_y_2, *ce_y_2, *alphae_y_PML_2, *sige_y_PML_2, *kappae_y_PML_2;
float *bh_y_2, *ch_y_2, *alphah_y_PML_2, *sigh_y_PML_2, *kappah_y_PML_2;
float *be_z_1, *ce_z_1, *alphae_z_PML_1, *sige_z_PML_1, *kappae_z_PML_1;
float *bh_z_1, *ch_z_1, *alphah_z_PML_1, *sigh_z_PML_1, *kappah_z_PML_1;
float *be_z_2, *ce_z_2, *alphae_z_PML_2, *sige_z_PML_2, *kappae_z_PML_2;
float *bh_z_2, *ch_z_2, *alphah_z_PML_2, *sigh_z_PML_2, *kappah_z_PML_2;

// denominators for the update equations
float *den_ex;
float *den_hx;
float *den_ey;
float *den_hy;
float *den_ez;
float *den_hz;

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
void UpdatePMLForHxIn(unsigned k, float***Hx, float*** Ey, float*** Ez, float DB) {
    for (i = 0; i < Imax - 1; ++i) {
        //...............................................
        //  PML for bottom Hx, j-direction
        //...............................................
        for (j = 0; j < nyPML_1 - 1; ++j) {
            psi_Hxy_1[i][j][k] = bh_y_1[j] * psi_Hxy_1[i][j][k]
                    + ch_y_1[j] * (Ez[i][j][k] - Ez[i][j + 1][k]) / dy;
            Hx[i][j][k] = Hx[i][j][k] + DB * psi_Hxy_1[i][j][k];
        }
        //....................................................
        //  PML for top Hx, j-direction
        //.....................................................
        jj = nyPML_2 - 2;
        for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {
            psi_Hxy_2[i][jj][k] = bh_y_2[jj] * psi_Hxy_2[i][jj][k]
                    + ch_y_2[jj] * (Ez[i][j][k] - Ez[i][j + 1][k]) / dy;
            Hx[i][j][k] = Hx[i][j][k] + DB * psi_Hxy_2[i][jj][k];
            jj = jj - 1;
        }
    }
}
void UpdatePMLForHxOut(float***Hx, float*** Ey, float*** Ez, float DB){
	unsigned i,j,k;
	unsigned kk;
        for (i = 0; i < Imax - 1; ++i) {

            for (j = 0; j < Jmax - 1; ++j) {
                //....................................................
                //  PML for bottom Hx, k-direction
                //................................................
                for (k = 1; k < nzPML_1; ++k) {

                    psi_Hxz_1[i][j][k - 1] = bh_z_1[k - 1] * psi_Hxz_1[i][j][k - 1]
                            + ch_z_1[k - 1] * (Ey[i][j][k] - Ey[i][j][k - 1]) / dz;
                    Hx[i][j][k] = Hx[i][j][k] + DB * psi_Hxz_1[i][j][k - 1];
                }
                //....................................................
                //  PML for top Hx, k-direction
                //...............................................
                kk = nzPML_2 - 2;
                for (k = Kmax - nzPML_2; k < Kmax - 1; ++k) {

                    psi_Hxz_2[i][j][kk] = bh_z_2[kk] * psi_Hxz_2[i][j][kk]
                            + ch_z_2[kk] * (Ey[i][j][k] - Ey[i][j][k - 1]) / dz;
                    Hx[i][j][k] = Hx[i][j][k] + DB * psi_Hxz_2[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
}
void UpdatePMLForHyIn(unsigned k,float*** Hy, float*** Ez, float*** Ex, float DB){
	unsigned i,j;
	unsigned ii;

            for (j = 0; j < Jmax - 1; ++j) {
                //.......................................................
                //  PML for bottom Hy, i-direction
                //.......................................................
                for (i = 0; i < nxPML_1 - 1; ++i) {

                    psi_Hyx_1[i][j][k] = bh_x_1[i] * psi_Hyx_1[i][j][k]
                            + ch_x_1[i] * (Ez[i + 1][j][k] - Ez[i][j][k]) / dx;
                    Hy[i][j][k] = Hy[i][j][k] + DB * psi_Hyx_1[i][j][k];
                }
                //.........................................................
                //  PML for top Hy, i-direction
                //.........................................................
                ii = nxPML_2 - 2;
                for (i = Imax - nxPML_2; i < Imax - 1; ++i) {

                    psi_Hyx_2[ii][j][k] = bh_x_2[ii] * psi_Hyx_2[ii][j][k]
                            + ch_x_2[ii] * (Ez[i + 1][j][k] - Ez[i][j][k]) / dx;
                    Hy[i][j][k] = Hy[i][j][k] + DB * psi_Hyx_2[ii][j][k];
                    ii = ii - 1;
                }
            }
}
void UpdatePMLForHyOut(float***Hy, float*** Ez, float*** Ex, float DB){
	unsigned i,j,k;
	unsigned kk;

        for (i = 0; i < Imax - 1; ++i) {
            for (j = 0; j < Jmax - 1; ++j) {
                //.......................................................
                //  PML for bottom Hy, k-direction
                //......................................................
                for (k = 1; k < nzPML_1; ++k) {

                    psi_Hyz_1[i][j][k - 1] = bh_z_1[k - 1] * psi_Hyz_1[i][j][k - 1]
                            + ch_z_1[k - 1] * (Ex[i][j][k - 1] - Ex[i][j][k]) / dz;
                    Hy[i][j][k] = Hy[i][j][k] + DB * psi_Hyz_1[i][j][k - 1];
                }
                //.......................................................
                //  PML for top Hy, k-direction
                //.........................................................
                kk = nzPML_2 - 2;
                for (k = Kmax - nzPML_2; k < Kmax - 1; ++k) {

                    psi_Hyz_2[i][j][kk] = bh_z_2[kk] * psi_Hyz_2[i][j][kk]
                            + ch_z_2[kk] * (Ex[i][j][k - 1] - Ex[i][j][k]) / dz;
                    Hy[i][j][k] = Hy[i][j][k] + DB * psi_Hyz_2[i][j][kk];
                    kk = kk - 1;
                }
            }
        }

}
void UpdatePMLForExIn(unsigned k, float***Ex, float*** Hy, float*** Hz, float *CB,short ***ID1){
	unsigned i,j;
	unsigned jj;
	short id;
            for (i = 0; i < Imax - 1; ++i) {
                //..............................................................
                //  PML for bottom Ex, j-direction
                //..............................................................
                for (j = 1; j < nyPML_1; ++j) {

                    id = ID1[i][j][k];
                    psi_Exy_1[i][j][k] = be_y_1[j] * psi_Exy_1[i][j][k]
                            + ce_y_1[j] * (Hz[i][j][k] - Hz[i][j - 1][k]) / dy;
                    Ex[i][j][k] = Ex[i][j][k] + CB[id] * psi_Exy_1[i][j][k];
                }
                //.............................................................
                //  PML for top Ex, j-direction
                //.............................................................
                jj = nyPML_2 - 1;
                for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {

                    id = ID1[i][j][k];
                    psi_Exy_2[i][jj][k] = be_y_2[jj] * psi_Exy_2[i][jj][k]
                            + ce_y_2[jj] * (Hz[i][j][k] - Hz[i][j - 1][k]) / dy;
                    Ex[i][j][k] = Ex[i][j][k] + CB[id] * psi_Exy_2[i][jj][k];
                    jj = jj - 1;
                }
            }
}
void UpdatePMLForExOut(float***Ex, float*** Hy, float*** Hz, float *CB,short ***ID1){
unsigned i,j,k;
unsigned kk;
	short id;
        for (i = 0; i < Imax - 1; ++i) {

            for (j = 1; j < Jmax - 1; ++j) {
                //.............................................................
                //  PML for bottom Ex, k-direction
                //.............................................................
                for (k = 0; k < nzPML_1; ++k) {

                    id = ID1[i][j][k];
                    psi_Exz_1[i][j][k] = be_z_1[k] * psi_Exz_1[i][j][k]
                            + ce_z_1[k] * (Hy[i][j][k] - Hy[i][j][k + 1]) / dz;
                    Ex[i][j][k] = Ex[i][j][k] + CB[id] * psi_Exz_1[i][j][k];
                }
                //..............................................................
                //  PML for top Ex, k-direction
                //..............................................................
                kk = nzPML_2 - 1;
                for (k = Kmax - nzPML_2 - 1; k < Kmax - 1; ++k) {

                    id = ID1[i][j][k];
                    psi_Exz_2[i][j][kk] = be_z_2[kk] * psi_Exz_2[i][j][kk]
                            + ce_z_2[kk] * (Hy[i][j][k] - Hy[i][j][k + 1]) / dz;
                    Ex[i][j][k] = Ex[i][j][k] + CB[id] * psi_Exz_2[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
}
void UpdatePMLForEyOut(float***Ey, float*** Hz, float*** Hx, float *CB,short ***ID2){
	unsigned i,j,k;
	unsigned kk;
	short id;
        for (i = 1; i < Imax - 1; ++i) {
            for (j = 0; j < Jmax - 1; ++j) {
                //...........................................................
                //  PML for bottom Ey, k-direction
                //...........................................................
                for (k = 0; k < nzPML_1; ++k) {
                    id = ID2[i][j][k];
                    psi_Eyz_1[i][j][k] = be_z_1[k] * psi_Eyz_1[i][j][k]
                            + ce_z_1[k] * (Hx[i][j][k + 1] - Hx[i][j][k]) / dz;
                    Ey[i][j][k] = Ey[i][j][k] + CB[id] * psi_Eyz_1[i][j][k];
                }
                //...........................................................
                //  PML for top Ey, k-direction
                //............................................................
                kk = nzPML_2 - 1;
                for (k = Kmax - nzPML_2 - 1; k < Kmax - 1; ++k) {

                    id = ID2[i][j][k];
                    psi_Eyz_2[i][j][kk] = be_z_2[kk] * psi_Eyz_2[i][j][kk]
                            + ce_z_2[kk] * (Hx[i][j][k + 1] - Hx[i][j][k]) / dz;
                    Ey[i][j][k] = Ey[i][j][k] + CB[id] * psi_Eyz_2[i][j][kk];
                    kk = kk - 1;
                }
            }
        }
}
void UpdatePMLForEyIn(unsigned k,float***Ey, float*** Hz, float*** Hx, float *CB,short ***ID2){
	unsigned i,j;
	unsigned ii;
	short id;
            for (j = 0; j < Jmax - 1; ++j) {
                //...........................................................
                //  PML for bottom Ey, i-direction
                //...........................................................
                for (i = 1; i < nxPML_1; ++i) {

                    id = ID2[i][j][k];
                    psi_Eyx_1[i][j][k] = be_x_1[i] * psi_Eyx_1[i][j][k]
                            + ce_x_1[i] * (Hz[i - 1][j][k] - Hz[i][j][k]) / dx;
                    Ey[i][j][k] = Ey[i][j][k] + CB[id] * psi_Eyx_1[i][j][k];
                }
                //............................................................
                //  PML for top Ey, i-direction
                //............................................................
                ii = nxPML_2 - 1;
                for (i = Imax - nxPML_2; i < Imax - 1; ++i) {

                    id = ID2[i][j][k];
                    psi_Eyx_2[ii][j][k] = be_x_2[ii] * psi_Eyx_2[ii][j][k]
                            + ce_x_2[ii] * (Hz[i - 1][j][k] - Hz[i][j][k]) / dx;
                    Ey[i][j][k] = Ey[i][j][k] + CB[id] * psi_Eyx_2[ii][j][k];
                    ii = ii - 1;
                }
            }
}
void UpdatePMLForHz(unsigned k,float***Hz, float*** Ex, float*** Ey, float DB){
	unsigned i,j;
	unsigned ii,jj;
            for (j = 0; j < Jmax - 1; ++j) {
                //..........................................................
                //  PML for bottom Hz, x-direction
                //..........................................................
                for (i = 0; i < nxPML_1 - 1; ++i) {

                    psi_Hzx_1[i][j][k] = bh_x_1[i] * psi_Hzx_1[i][j][k]
                            + ch_x_1[i] * (Ey[i][j][k] - Ey[i + 1][j][k]) / dx;
                    Hz[i][j][k] = Hz[i][j][k] + DB * psi_Hzx_1[i][j][k];
                }
                //..........................................................
                //  PML for top Hz, x-direction
                //..........................................................
                ii = nxPML_2 - 2;
                for (i = Imax - nxPML_2; i < Imax - 1; ++i) {

                    psi_Hzx_2[ii][j][k] = bh_x_2[ii] * psi_Hzx_2[ii][j][k]
                            + ch_x_2[ii] * (Ey[i][j][k] - Ey[i + 1][j][k]) / dx;
                    Hz[i][j][k] = Hz[i][j][k] + DB * psi_Hzx_2[ii][j][k];
                    ii = ii - 1;
                }
            }

            for (i = 0; i < Imax - 1; ++i) {
                //........................................................
                //  PML for bottom Hz, y-direction
                //.........................................................
                for (j = 0; j < nyPML_1 - 1; ++j) {
                    psi_Hzy_1[i][j][k] = bh_y_1[j] * psi_Hzy_1[i][j][k]
                            + ch_y_1[j] * (Ex[i][j + 1][k] - Ex[i][j][k]) / dy;
                    Hz[i][j][k] = Hz[i][j][k] + DB * psi_Hzy_1[i][j][k];

                }
                //.........................................................
                //  PML for top Hz, y-direction
                //..........................................................
                jj = nyPML_2 - 2;
                for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {

                    psi_Hzy_2[i][jj][k] = bh_y_2[jj] * psi_Hzy_2[i][jj][k]
                            + ch_y_2[jj] * (Ex[i][j + 1][k] - Ex[i][j][k]) / dy;
                    Hz[i][j][k] = Hz[i][j][k] + DB * psi_Hzy_2[i][jj][k];
                    jj = jj - 1;
                }
            }
}
void UpdatePMLForEz(unsigned k,float***Ez, float*** Hx, float*** Hy, float *CB,short ***ID3){

	unsigned i,j;
	unsigned ii,jj;
	short id;
            for (j = 1; j < Jmax - 1; ++j) {
                //............................................................
                //  PML for bottom Ez, x-direction
                //.............................................................
                for (i = 1; i < nxPML_1; ++i) {
                    id = ID3[i][j][k];
                    psi_Ezx_1[i][j][k] = be_x_1[i] * psi_Ezx_1[i][j][k]
                            + ce_x_1[i] * (Hy[i][j][k] - Hy[i - 1][j][k]) / dx;
                    Ez[i][j][k] = Ez[i][j][k] + CB[id] * psi_Ezx_1[i][j][k];
                }
                //............................................................
                //  PML for top Ez, x-direction
                //............................................................
                ii = nxPML_2 - 1;
                for (i = Imax - nxPML_2; i < Imax - 1; ++i) {
                    id = ID3[i][j][k];
                    psi_Ezx_2[ii][j][k] = be_x_2[ii] * psi_Ezx_2[ii][j][k]
                            + ce_x_2[ii] * (Hy[i][j][k] - Hy[i - 1][j][k]) / dx;
                    Ez[i][j][k] = Ez[i][j][k] + CB[id] * psi_Ezx_2[ii][j][k];
                    ii = ii - 1;
                }
            }

            for (i = 1; i < Imax - 1; ++i) {
                //..........................................................
                //  PML for bottom Ez, y-direction
                //..........................................................
                for (j = 1; j < nyPML_1; ++j) {

                    id = ID3[i][j][k];
                    psi_Ezy_1[i][j][k] = be_y_1[j] * psi_Ezy_1[i][j][k]
                            + ce_y_1[j] * (Hx[i][j - 1][k] - Hx[i][j][k]) / dy;
                    Ez[i][j][k] = Ez[i][j][k] + CB[id] * psi_Ezy_1[i][j][k];
                }
                //............................................................
                //  PML for top Ez, y-direction
                //............................................................
                jj = nyPML_2 - 1;
                for (j = Jmax - nyPML_2; j < Jmax - 1; ++j) {

                    id = ID3[i][j][k];
                    psi_Ezy_2[i][jj][k] = be_y_2[jj] * psi_Ezy_2[i][jj][k]
                            + ce_y_2[jj] * (Hx[i][j - 1][k] - Hx[i][j][k]) / dy;
                    Ez[i][j][k] = Ez[i][j][k] + CB[id] * psi_Ezy_2[i][jj][k];
                    jj = jj - 1;
                }
            }
}

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

    //PML Layers (10 layers)
    nxPML_1 = 11;
    nxPML_2 = 11;
    nyPML_1 = 11;
    nyPML_2 = 11;
    nzPML_1 = 11;
    nzPML_2 = 11;

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

    Ez = (float ***) malloc(Imax * sizeof (float **));

    for (i = 0; i < Imax; i++) {

        Ez[i] = (float **) malloc(Jmax * sizeof (float *));

        for (j = 0; j < Jmax; j++) {

            Ez[i][j] = (float *) malloc(Kmax * sizeof (float));

            for (k = 0; k < Kmax; k++) {

                Ez[i][j][k] = 0.0;
            }
        }
    }

    Ey = (float ***) malloc((Imax) * sizeof (float **));

    for (i = 0; i < Imax; i++) {

        Ey[i] = (float **) malloc((Jmax - 1) * sizeof (float *));

        for (j = 0; j < Jmax - 1; j++) {

            Ey[i][j] = (float *) malloc((Kmax - 1) * sizeof (float));

            for (k = 0; k < Kmax - 1; k++) {

                Ey[i][j][k] = 0.0;
            }
        }
    }

    Ex = (float ***) malloc((Imax - 1) * sizeof (float **));

    for (i = 0; i < Imax - 1; i++) {

        Ex[i] = (float **) malloc(Jmax * sizeof (float *));

        for (j = 0; j < Jmax; j++) {

            Ex[i][j] = (float *) malloc((Kmax - 1) * sizeof (float));

            for (k = 0; k < Kmax - 1; k++) {

                Ex[i][j][k] = 0.0;
            }
        }
    }

    Hx = (float ***) malloc(Imax * sizeof (float **));

    for (i = 0; i < Imax; i++) {

        Hx[i] = (float **) malloc((Jmax - 1) * sizeof (float *));

        for (j = 0; j < Jmax - 1; j++) {

            Hx[i][j] = (float *) malloc(Kmax * sizeof (float));

            for (k = 0; k < Kmax; k++) {

                Hx[i][j][k] = 0.0;
            }
        }
    }

    Hy = (float ***) malloc((Imax - 1) * sizeof (float **));

    for (i = 0; i < Imax - 1; i++) {

        Hy[i] = (float **) malloc(Jmax * sizeof (float *));

        for (j = 0; j < Jmax; j++) {

            Hy[i][j] = (float *) malloc(Kmax * sizeof (float));

            for (k = 0; k < Kmax; k++) {

                Hy[i][j][k] = 0.0;
            }
        }
    }

    Hz = (float ***) malloc((Imax - 1) * sizeof (float **));

    for (i = 0; i < Imax - 1; i++) {

        Hz[i] = (float **) malloc((Jmax - 1) * sizeof (float *));

        for (j = 0; j < Jmax - 1; j++) {

            Hz[i][j] = (float *) malloc((Kmax - 1) * sizeof (float));

            for (k = 0; k < Kmax - 1; k++) {

                Hz[i][j][k] = 0.0;
            }
        }
    }

    ID1 = (short ***) malloc(Imax * sizeof (short **));

    for (i = 0; i < Imax; i++) {

        ID1[i] = (short **) malloc(Jmax * sizeof (short *));

        for (j = 0; j < Jmax; j++) {

            ID1[i][j] = (short *) malloc(Kmax * sizeof (short));

            for (k = 0; k < Kmax; k++) {

                ID1[i][j][k] = 0;
            }
        }
    }

    ID2 = (short ***) malloc(Imax * sizeof (short **));

    for (i = 0; i < Imax; i++) {

        ID2[i] = (short **) malloc(Jmax * sizeof (short *));

        for (j = 0; j < Jmax; j++) {

            ID2[i][j] = (short *) malloc(Kmax * sizeof (short));

            for (k = 0; k < Kmax; k++) {

                ID2[i][j][k] = 0;
            }
        }
    }

    ID3 = (short ***) malloc(Imax * sizeof (short **));

    for (i = 0; i < Imax; i++) {

        ID3[i] = (short **) malloc(Jmax * sizeof (short *));

        for (j = 0; j < Jmax; j++) {

            ID3[i][j] = (short *) malloc(Kmax * sizeof (short));

            for (k = 0; k < Kmax; k++) {

                ID3[i][j][k] = 0;
            }
        }
    }

    psi_Ezx_1 = (float ***) malloc(nxPML_1 * sizeof (float **));

    for (i = 0; i < nxPML_1; i++) {

        psi_Ezx_1[i] = (float **) malloc(Jmax * sizeof (float *));

        for (j = 0; j < Jmax; j++) {

            psi_Ezx_1[i][j] = (float *) malloc(Kmax * sizeof (float));

            for (k = 0; k < Kmax; k++) {

                psi_Ezx_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Ezx_2 = (float ***) malloc(nxPML_2 * sizeof (float **));

    for (i = 0; i < nxPML_2; i++) {

        psi_Ezx_2[i] = (float **) malloc(Jmax * sizeof (float *));

        for (j = 0; j < Jmax; j++) {

            psi_Ezx_2[i][j] = (float *) malloc(Kmax * sizeof (float));

            for (k = 0; k < Kmax; k++) {

                psi_Ezx_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hyx_1 = (float ***) malloc((nxPML_1 - 1) * sizeof (float **));

    for (i = 0; i < nxPML_1 - 1; i++) {

        psi_Hyx_1[i] = (float **) malloc(Jmax * sizeof (float *));

        for (j = 0; j < Jmax; j++) {

            psi_Hyx_1[i][j] = (float *) malloc(Kmax * sizeof (float));

            for (k = 0; k < Kmax; k++) {

                psi_Hyx_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hyx_2 = (float ***) malloc((nxPML_2 - 1) * sizeof (float **));

    for (i = 0; i < nxPML_2 - 1; i++) {

        psi_Hyx_2[i] = (float **) malloc(Jmax * sizeof (float *));

        for (j = 0; j < Jmax; j++) {

            psi_Hyx_2[i][j] = (float *) malloc(Kmax * sizeof (float));

            for (k = 0; k < Kmax; k++) {

                psi_Hyx_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Ezy_1 = (float ***) malloc(Imax * sizeof (float **));

    for (i = 0; i < Imax; i++) {

        psi_Ezy_1[i] = (float **) malloc(nyPML_1 * sizeof (float *));

        for (j = 0; j < nyPML_1; j++) {

            psi_Ezy_1[i][j] = (float *) malloc(Kmax * sizeof (float));

            for (k = 0; k < Kmax; k++) {

                psi_Ezy_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Ezy_2 = (float ***) malloc(Imax * sizeof (float **));

    for (i = 0; i < Imax; i++) {

        psi_Ezy_2[i] = (float **) malloc(nyPML_2 * sizeof (float *));

        for (j = 0; j < nyPML_2; j++) {

            psi_Ezy_2[i][j] = (float *) malloc(Kmax * sizeof (float));

            for (k = 0; k < Kmax; k++) {

                psi_Ezy_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hxy_1 = (float ***) malloc(Imax * sizeof (float **));

    for (i = 0; i < Imax; i++) {

        psi_Hxy_1[i] = (float **) malloc((nyPML_1 - 1) * sizeof (float *));

        for (j = 0; j < nyPML_1 - 1; j++) {

            psi_Hxy_1[i][j] = (float *) malloc(Kmax * sizeof (float));

            for (k = 0; k < Kmax; k++) {

                psi_Hxy_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hxy_2 = (float ***) malloc(Imax * sizeof (float **));

    for (i = 0; i < Imax; i++) {

        psi_Hxy_2[i] = (float **) malloc((nyPML_2 - 1) * sizeof (float *));

        for (j = 0; j < nyPML_2 - 1; j++) {

            psi_Hxy_2[i][j] = (float *) malloc(Kmax * sizeof (float));

            for (k = 0; k < Kmax; k++) {

                psi_Hxy_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hxz_1 = (float ***) malloc(Imax * sizeof (float **));

    for (i = 0; i < Imax; i++) {

        psi_Hxz_1[i] = (float **) malloc((Jmax - 1) * sizeof (float *));

        for (j = 0; j < Jmax - 1; j++) {

            psi_Hxz_1[i][j] = (float *) malloc((nzPML_1 - 1) * sizeof (float));

            for (k = 0; k < nzPML_1 - 1; k++) {

                psi_Hxz_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hxz_2 = (float ***) malloc(Imax * sizeof (float **));

    for (i = 0; i < Imax; i++) {

        psi_Hxz_2[i] = (float **) malloc((Jmax - 1) * sizeof (float *));

        for (j = 0; j < Jmax - 1; j++) {

            psi_Hxz_2[i][j] = (float *) malloc((nzPML_2 - 1) * sizeof (float));

            for (k = 0; k < nzPML_2 - 1; k++) {

                psi_Hxz_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hyz_1 = (float ***) malloc((Imax - 1) * sizeof (float **));

    for (i = 0; i < Imax - 1; i++) {

        psi_Hyz_1[i] = (float **) malloc(Jmax * sizeof (float *));

        for (j = 0; j < Jmax; j++) {

            psi_Hyz_1[i][j] = (float *) malloc((nzPML_1 - 1) * sizeof (float));

            for (k = 0; k < nzPML_1 - 1; k++) {

                psi_Hyz_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hyz_2 = (float ***) malloc((Imax - 1) * sizeof (float **));

    for (i = 0; i < Imax - 1; i++) {

        psi_Hyz_2[i] = (float **) malloc(Jmax * sizeof (float *));

        for (j = 0; j < Jmax; j++) {

            psi_Hyz_2[i][j] = (float *) malloc((nzPML_2 - 1) * sizeof (float));

            for (k = 0; k < nzPML_2 - 1; k++) {

                psi_Hyz_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Exz_1 = (float ***) malloc((Imax - 1) * sizeof (float **));

    for (i = 0; i < Imax - 1; i++) {

        psi_Exz_1[i] = (float **) malloc(Jmax * sizeof (float *));

        for (j = 0; j < Jmax; j++) {

            psi_Exz_1[i][j] = (float *) malloc(nzPML_1 * sizeof (float));

            for (k = 0; k < nzPML_1; k++) {

                psi_Exz_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Exz_2 = (float ***) malloc((Imax - 1) * sizeof (float **));

    for (i = 0; i < Imax - 1; i++) {

        psi_Exz_2[i] = (float **) malloc(Jmax * sizeof (float *));

        for (j = 0; j < Jmax; j++) {

            psi_Exz_2[i][j] = (float *) malloc(nzPML_2 * sizeof (float));

            for (k = 0; k < nzPML_2; k++) {

                psi_Exz_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Eyz_1 = (float ***) malloc((Imax) * sizeof (float **));

    for (i = 0; i < Imax; i++) {

        psi_Eyz_1[i] = (float **) malloc((Jmax - 1) * sizeof (float *));

        for (j = 0; j < Jmax - 1; j++) {

            psi_Eyz_1[i][j] = (float *) malloc(nzPML_1 * sizeof (float));

            for (k = 0; k < nzPML_1; k++) {

                psi_Eyz_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Eyz_2 = (float ***) malloc((Imax) * sizeof (float **));

    for (i = 0; i < Imax; i++) {

        psi_Eyz_2[i] = (float **) malloc((Jmax - 1) * sizeof (float *));

        for (j = 0; j < Jmax - 1; j++) {

            psi_Eyz_2[i][j] = (float *) malloc(nzPML_2 * sizeof (float));

            for (k = 0; k < nzPML_2; k++) {

                psi_Eyz_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hzx_1 = (float ***) malloc((nxPML_1 - 1) * sizeof (float **));

    for (i = 0; i < nxPML_1 - 1; i++) {

        psi_Hzx_1[i] = (float **) malloc((Jmax - 1) * sizeof (float *));

        for (j = 0; j < Jmax - 1; j++) {

            psi_Hzx_1[i][j] = (float *) malloc((Kmax - 1) * sizeof (float));

            for (k = 0; k < Kmax - 1; k++) {

                psi_Hzx_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hzx_2 = (float ***) malloc((nxPML_2 - 1) * sizeof (float **));

    for (i = 0; i < nxPML_2 - 1; i++) {

        psi_Hzx_2[i] = (float **) malloc((Jmax - 1) * sizeof (float *));

        for (j = 0; j < Jmax - 1; j++) {

            psi_Hzx_2[i][j] = (float *) malloc((Kmax - 1) * sizeof (float));

            for (k = 0; k < Kmax - 1; k++) {

                psi_Hzx_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Eyx_1 = (float ***) malloc(nxPML_1 * sizeof (float **));

    for (i = 0; i < nxPML_1; i++) {

        psi_Eyx_1[i] = (float **) malloc((Jmax - 1) * sizeof (float *));

        for (j = 0; j < Jmax - 1; j++) {

            psi_Eyx_1[i][j] = (float *) malloc((Kmax - 1) * sizeof (float));

            for (k = 0; k < Kmax - 1; k++) {

                psi_Eyx_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Eyx_2 = (float ***) malloc(nxPML_2 * sizeof (float **));

    for (i = 0; i < nxPML_2; i++) {

        psi_Eyx_2[i] = (float **) malloc((Jmax - 1) * sizeof (float *));

        for (j = 0; j < Jmax - 1; j++) {

            psi_Eyx_2[i][j] = (float *) malloc((Kmax - 1) * sizeof (float));

            for (k = 0; k < Kmax - 1; k++) {

                psi_Eyx_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Hzy_1 = (float ***) malloc((Imax - 1) * sizeof (float **));

    for (i = 0; i < Imax - 1; i++) {

        psi_Hzy_1[i] = (float **) malloc((nyPML_1 - 1) * sizeof (float *));

        for (j = 0; j < nyPML_1 - 1; j++) {

            psi_Hzy_1[i][j] = (float *) malloc((Kmax - 1) * sizeof (float));

            for (k = 0; k < Kmax - 1; k++) {

                psi_Hzy_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Hzy_2 = (float ***) malloc((Imax - 1) * sizeof (float **));

    for (i = 0; i < Imax - 1; i++) {

        psi_Hzy_2[i] = (float **) malloc((nyPML_2 - 1) * sizeof (float *));

        for (j = 0; j < nyPML_2 - 1; j++) {

            psi_Hzy_2[i][j] = (float *) malloc((Kmax - 1) * sizeof (float));

            for (k = 0; k < Kmax - 1; k++) {

                psi_Hzy_2[i][j][k] = 0.0;
            }
        }
    }

    psi_Exy_1 = (float ***) malloc((Imax - 1) * sizeof (float **));

    for (i = 0; i < Imax - 1; i++) {

        psi_Exy_1[i] = (float **) malloc(nyPML_1 * sizeof (float *));

        for (j = 0; j < nyPML_1; j++) {

            psi_Exy_1[i][j] = (float *) malloc((Kmax - 1) * sizeof (float));

            for (k = 0; k < Kmax - 1; k++) {

                psi_Exy_1[i][j][k] = 0.0;
            }
        }
    }

    psi_Exy_2 = (float ***) malloc((Imax - 1) * sizeof (float **));

    for (i = 0; i < Imax - 1; i++) {

        psi_Exy_2[i] = (float **) malloc(nyPML_2 * sizeof (float *));

        for (j = 0; j < nyPML_2; j++) {

            psi_Exy_2[i][j] = (float *) malloc((Kmax - 1) * sizeof (float));

            for (k = 0; k < Kmax - 1; k++) {

                psi_Exy_2[i][j][k] = 0.0;
            }
        }
    }

    den_ex = (float *) malloc((Imax - 1) * sizeof (float));
    for (i = 0; i < Imax - 1; i++) {

        den_ex[i] = 0.0;
    }

    den_hx = (float *) malloc((Imax - 1) * sizeof (float));
    for (i = 0; i < Imax - 1; i++) {

        den_hx[i] = 0.0;
    }

    den_ey = (float *) malloc((Jmax - 1) * sizeof (float));
    for (i = 0; i < Jmax - 1; i++) {

        den_ey[i] = 0.0;
    }

    den_hy = (float *) malloc((Jmax - 1) * sizeof (float));
    for (i = 0; i < Jmax - 1; i++) {

        den_hy[i] = 0.0;
    }

    den_ez = (float *) malloc((Kmax - 1) * sizeof (float));
    for (i = 0; i < Kmax - 1; i++) {

        den_ez[i] = 0.0;
    }

    den_hz = (float *) malloc((Kmax - 1) * sizeof (float));
    for (i = 0; i < Kmax - 1; i++) {

        den_hz[i] = 0.0;
    }

    be_x_1 = (float *) malloc((nxPML_1) * sizeof (float));
    for (i = 0; i < nxPML_1; i++) {

        be_x_1[i] = 0.0;
    }

    ce_x_1 = (float *) malloc((nxPML_1) * sizeof (float));
    for (i = 0; i < nxPML_1; i++) {

        ce_x_1[i] = 0.0;
    }

    alphae_x_PML_1 = (float *) malloc((nxPML_1) * sizeof (float));
    for (i = 0; i < nxPML_1; i++) {

        alphae_x_PML_1[i] = 0.0;
    }

    sige_x_PML_1 = (float *) malloc((nxPML_1) * sizeof (float));
    for (i = 0; i < nxPML_1; i++) {

        sige_x_PML_1[i] = 0.0;
    }

    kappae_x_PML_1 = (float *) malloc((nxPML_1) * sizeof (float));
    for (i = 0; i < nxPML_1; i++) {

        kappae_x_PML_1[i] = 0.0;
    }

    bh_x_1 = (float *) malloc((nxPML_1 - 1) * sizeof (float));
    for (i = 0; i < nxPML_1 - 1; i++) {

        bh_x_1[i] = 0.0;
    }

    ch_x_1 = (float *) malloc((nxPML_1 - 1) * sizeof (float));
    for (i = 0; i < nxPML_1 - 1; i++) {

        ch_x_1[i] = 0.0;
    }

    alphah_x_PML_1 = (float *) malloc((nxPML_1 - 1) * sizeof (float));
    for (i = 0; i < nxPML_1 - 1; i++) {

        alphah_x_PML_1[i] = 0.0;
    }

    sigh_x_PML_1 = (float *) malloc((nxPML_1 - 1) * sizeof (float));
    for (i = 0; i < nxPML_1 - 1; i++) {

        sigh_x_PML_1[i] = 0.0;
    }

    kappah_x_PML_1 = (float *) malloc((nxPML_1 - 1) * sizeof (float));
    for (i = 0; i < nxPML_1 - 1; i++) {

        kappah_x_PML_1[i] = 0.0;
    }

    be_x_2 = (float *) malloc((nxPML_2) * sizeof (float));
    for (i = 0; i < nxPML_2; i++) {

        be_x_2[i] = 0.0;
    }

    ce_x_2 = (float *) malloc((nxPML_2) * sizeof (float));
    for (i = 0; i < nxPML_2; i++) {

        ce_x_2[i] = 0.0;
    }

    alphae_x_PML_2 = (float *) malloc((nxPML_2) * sizeof (float));
    for (i = 0; i < nxPML_2; i++) {

        alphae_x_PML_2[i] = 0.0;
    }


    sige_x_PML_2 = (float *) malloc((nxPML_2) * sizeof (float));
    for (i = 0; i < nxPML_2; i++) {

        sige_x_PML_2[i] = 0.0;
    }


    kappae_x_PML_2 = (float *) malloc((nxPML_2) * sizeof (float));
    for (i = 0; i < nxPML_2; i++) {

        kappae_x_PML_2[i] = 0.0;
    }


    bh_x_2 = (float *) malloc((nxPML_2 - 1) * sizeof (float));
    for (i = 0; i < nxPML_2 - 1; i++) {

        bh_x_2[i] = 0.0;
    }


    ch_x_2 = (float *) malloc((nxPML_2 - 1) * sizeof (float));
    for (i = 0; i < nxPML_2 - 1; i++) {

        ch_x_2[i] = 0.0;
    }

    alphah_x_PML_2 = (float *) malloc((nxPML_2 - 1) * sizeof (float));
    for (i = 0; i < nxPML_2 - 1; i++) {

        alphah_x_PML_2[i] = 0.0;
    }

    sigh_x_PML_2 = (float *) malloc((nxPML_2 - 1) * sizeof (float));
    for (i = 0; i < nxPML_2 - 1; i++) {
        sigh_x_PML_2[i] = 0.0;
    }

    kappah_x_PML_2 = (float *) malloc((nxPML_2 - 1) * sizeof (float));
    for (i = 0; i < nxPML_1 - 1; i++) {

        kappah_x_PML_2[i] = 0.0;
    }

    be_y_1 = (float *) malloc((nyPML_1) * sizeof (float));
    for (i = 0; i < nyPML_1; i++) {

        be_y_1[i] = 0.0;
    }

    ce_y_1 = (float *) malloc((nyPML_1) * sizeof (float));
    for (i = 0; i < nyPML_1; i++) {

        ce_y_1[i] = 0.0;
    }

    alphae_y_PML_1 = (float *) malloc((nyPML_1) * sizeof (float));
    for (i = 0; i < nyPML_1; i++) {

        alphae_y_PML_1[i] = 0.0;
    }

    sige_y_PML_1 = (float *) malloc((nyPML_1) * sizeof (float));
    for (i = 0; i < nyPML_1; i++) {

        sige_y_PML_1[i] = 0.0;
    }

    kappae_y_PML_1 = (float *) malloc((nyPML_1) * sizeof (float));
    for (i = 0; i < nyPML_1; i++) {

        kappae_y_PML_1[i] = 0.0;
    }

    bh_y_1 = (float *) malloc((nyPML_1 - 1) * sizeof (float));
    for (i = 0; i < nyPML_1 - 1; i++) {

        bh_y_1[i] = 0.0;
    }

    ch_y_1 = (float *) malloc((nyPML_1 - 1) * sizeof (float));
    for (i = 0; i < nyPML_1 - 1; i++) {

        ch_y_1[i] = 0.0;
    }

    alphah_y_PML_1 = (float *) malloc((nyPML_1 - 1) * sizeof (float));
    for (i = 0; i < nyPML_1 - 1; i++) {

        alphah_y_PML_1[i] = 0.0;
    }

    sigh_y_PML_1 = (float *) malloc((nyPML_1 - 1) * sizeof (float));
    for (i = 0; i < nyPML_1 - 1; i++) {

        sigh_y_PML_1[i] = 0.0;
    }

    kappah_y_PML_1 = (float *) malloc((nyPML_1 - 1) * sizeof (float));
    for (i = 0; i < nyPML_1 - 1; i++) {

        kappah_y_PML_1[i] = 0.0;
    }

    be_y_2 = (float *) malloc((nyPML_2) * sizeof (float));
    for (i = 0; i < nyPML_2; i++) {

        be_y_2[i] = 0.0;
    }

    ce_y_2 = (float *) malloc((nyPML_2) * sizeof (float));
    for (i = 0; i < nyPML_2; i++) {

        ce_y_2[i] = 0.0;
    }

    alphae_y_PML_2 = (float *) malloc((nyPML_2) * sizeof (float));
    for (i = 0; i < nyPML_2; i++) {

        alphae_y_PML_2[i] = 0.0;
    }

    sige_y_PML_2 = (float *) malloc((nyPML_2) * sizeof (float));
    for (i = 0; i < nyPML_2; i++) {

        sige_y_PML_2[i] = 0.0;
    }

    kappae_y_PML_2 = (float *) malloc((nyPML_2) * sizeof (float));
    for (i = 0; i < nyPML_2; i++) {

        kappae_y_PML_2[i] = 0.0;
    }

    bh_y_2 = (float *) malloc((nyPML_2 - 1) * sizeof (float));
    for (i = 0; i < nyPML_2 - 1; i++) {

        bh_y_2[i] = 0.0;
    }

    ch_y_2 = (float *) malloc((nyPML_2 - 1) * sizeof (float));
    for (i = 0; i < nyPML_2 - 1; i++) {

        ch_y_2[i] = 0.0;
    }

    alphah_y_PML_2 = (float *) malloc((nyPML_2 - 1) * sizeof (float));
    for (i = 0; i < nyPML_2 - 1; i++) {

        alphah_y_PML_2[i] = 0.0;
    }

    sigh_y_PML_2 = (float *) malloc((nyPML_2 - 1) * sizeof (float));
    for (i = 0; i < nyPML_2 - 1; i++) {

        sigh_y_PML_2[i] = 0.0;
    }

    kappah_y_PML_2 = (float *) malloc((nyPML_2 - 1) * sizeof (float));
    for (i = 0; i < nyPML_1 - 1; i++) {

        kappah_y_PML_2[i] = 0.0;
    }

    be_z_1 = (float *) malloc((nzPML_1) * sizeof (float));
    for (i = 0; i < nzPML_1; i++) {

        be_z_1[i] = 0.0;
    }

    ce_z_1 = (float *) malloc((nzPML_1) * sizeof (float));
    for (i = 0; i < nzPML_1; i++) {

        ce_z_1[i] = 0.0;
    }

    alphae_z_PML_1 = (float *) malloc((nzPML_1) * sizeof (float));
    for (i = 0; i < nzPML_1; i++) {

        alphae_z_PML_1[i] = 0.0;
    }

    sige_z_PML_1 = (float *) malloc((nzPML_1) * sizeof (float));
    for (i = 0; i < nzPML_1; i++) {

        sige_z_PML_1[i] = 0.0;
    }

    kappae_z_PML_1 = (float *) malloc((nzPML_1) * sizeof (float));
    for (i = 0; i < nzPML_1; i++) {

        kappae_z_PML_1[i] = 0.0;
    }

    bh_z_1 = (float *) malloc((nzPML_1 - 1) * sizeof (float));
    for (i = 0; i < nzPML_1 - 1; i++) {

        bh_z_1[i] = 0.0;
    }

    ch_z_1 = (float *) malloc((nzPML_1 - 1) * sizeof (float));
    for (i = 0; i < nzPML_1 - 1; i++) {

        ch_z_1[i] = 0.0;
    }

    alphah_z_PML_1 = (float *) malloc((nzPML_1 - 1) * sizeof (float));
    for (i = 0; i < nzPML_1 - 1; i++) {

        alphah_z_PML_1[i] = 0.0;
    }

    sigh_z_PML_1 = (float *) malloc((nzPML_1 - 1) * sizeof (float));
    for (i = 0; i < nzPML_1 - 1; i++) {

        sigh_z_PML_1[i] = 0.0;
    }

    kappah_z_PML_1 = (float *) malloc((nzPML_1 - 1) * sizeof (float));
    for (i = 0; i < nzPML_1 - 1; i++) {

        kappah_z_PML_1[i] = 0.0;
    }

    be_z_2 = (float *) malloc((nzPML_2) * sizeof (float));
    for (i = 0; i < nzPML_2; i++) {

        be_z_2[i] = 0.0;
    }

    ce_z_2 = (float *) malloc((nzPML_2) * sizeof (float));
    for (i = 0; i < nzPML_2; i++) {

        ce_z_2[i] = 0.0;
    }

    alphae_z_PML_2 = (float *) malloc((nzPML_2) * sizeof (float));
    for (i = 0; i < nzPML_2; i++) {

        alphae_z_PML_2[i] = 0.0;
    }

    sige_z_PML_2 = (float *) malloc((nzPML_2) * sizeof (float));
    for (i = 0; i < nzPML_2; i++) {

        sige_z_PML_2[i] = 0.0;
    }

    kappae_z_PML_2 = (float *) malloc((nzPML_2) * sizeof (float));
    for (i = 0; i < nzPML_2; i++) {

        kappae_z_PML_2[i] = 0.0;
    }

    bh_z_2 = (float *) malloc((nzPML_2 - 1) * sizeof (float));
    for (i = 0; i < nzPML_2 - 1; i++) {

        bh_z_2[i] = 0.0;
    }

    ch_z_2 = (float *) malloc((nzPML_2 - 1) * sizeof (float));
    for (i = 0; i < nzPML_2 - 1; i++) {

        ch_z_2[i] = 0.0;
    }

    alphah_z_PML_2 = (float *) malloc((nzPML_2 - 1) * sizeof (float));
    for (i = 0; i < nzPML_2 - 1; i++) {

        alphah_z_PML_2[i] = 0.0;
    }

    sigh_z_PML_2 = (float *) malloc((nzPML_2 - 1) * sizeof (float));
    for (i = 0; i < nzPML_2 - 1; i++) {

        sigh_z_PML_2[i] = 0.0;
    }


    kappah_z_PML_2 = (float *) malloc((nzPML_2 - 1) * sizeof (float));
    for (i = 0; i < nzPML_1 - 1; i++) {

        kappah_z_PML_2[i] = 0.0;
    }

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

    sig_x_max = 0.75 * (0.8 * (m + 1) / (dx * sqrt(mu_0 / (eps_0 * epsR))));
    sig_y_max = 0.75 * (0.8 * (m + 1) / (dy * sqrt(mu_0 / (eps_0 * epsR))));
    sig_z_max = 0.75 * (0.8 * (m + 1) / (dz * sqrt(mu_0 / (eps_0 * epsR))));
    alpha_x_max = 0.24;
    alpha_y_max = alpha_x_max;
    alpha_z_max = alpha_x_max;
    kappa_x_max = 15.0;
    kappa_y_max = kappa_x_max;
    kappa_z_max = kappa_x_max;
    printf("\nTIme step = %e", dt);
    printf("\n Number of steps = %d", nMax);
    printf("\n Total Simulation time = %e Seconds", nMax * dt);

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SET CPML PARAMETERS IN EACH DIRECTION
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void initializeCPML() {

    for (i = 0; i < nxPML_1; ++i) {

        sige_x_PML_1[i] = sig_x_max * pow(((nxPML_1 - 1 - i)
                / (nxPML_1 - 1.0)), m);
        alphae_x_PML_1[i] = alpha_x_max * pow(((float) i
                / (nxPML_1 - 1.0)), ma);
        kappae_x_PML_1[i] = 1.0 + (kappa_x_max - 1.0) *
                pow((nxPML_1 - 1 - i) / (nxPML_1 - 1.0), m);
        be_x_1[i] = exp(-(sige_x_PML_1[i] / kappae_x_PML_1[i] +
                alphae_x_PML_1[i]) * dt / eps_0);

        if ((sige_x_PML_1[i] == 0.0) &&
                (alphae_x_PML_1[i] == 0.0) && (i == nxPML_1 - 1)) {

            ce_x_1[i] = 0.0;

        } else {

            ce_x_1[i] = sige_x_PML_1[i] * (be_x_1[i] - 1.0) /
                    (sige_x_PML_1[i] + kappae_x_PML_1[i] * alphae_x_PML_1[i])
                    / kappae_x_PML_1[i];
        }
    }

    for (i = 0; i < nxPML_1 - 1; ++i) {

        sigh_x_PML_1[i] = sig_x_max * pow(((nxPML_1 - 1 - i - 0.5)
                / (nxPML_1 - 1.0)), m);
        alphah_x_PML_1[i] = alpha_x_max * pow(((i + 1 - 0.5)
                / (nxPML_1 - 1.0)), ma);
        kappah_x_PML_1[i] = 1.0 + (kappa_x_max - 1.0) *
                pow(((nxPML_1 - 1 - i - 0.5) / (nxPML_1 - 1.0)), m);
        bh_x_1[i] = exp(-(sigh_x_PML_1[i] / kappah_x_PML_1[i] +
                alphah_x_PML_1[i]) * dt / eps_0);
        ch_x_1[i] = sigh_x_PML_1[i] * (bh_x_1[i] - 1.0) /
                (sigh_x_PML_1[i] + kappah_x_PML_1[i] * alphah_x_PML_1[i])
                / kappah_x_PML_1[i];
    }

    for (i = 0; i < nxPML_2; ++i) {

        sige_x_PML_2[i] = sig_x_max * pow(((nxPML_2 - 1 - i)
                / (nxPML_2 - 1.0)), m);
        alphae_x_PML_2[i] = alpha_x_max * pow(((float) i
                / (nxPML_2 - 1.0)), ma);
        kappae_x_PML_2[i] = 1.0 + (kappa_x_max - 1.0) *
                pow((nxPML_2 - 1 - i) / (nxPML_2 - 1.0), m);
        be_x_2[i] = exp(-(sige_x_PML_2[i] / kappae_x_PML_2[i] +
                alphae_x_PML_2[i]) * dt / eps_0);

        if ((sige_x_PML_2[i] == 0.0) &&
                (alphae_x_PML_2[i] == 0.0) && (i == nxPML_2 - 1)) {

            ce_x_2[i] = 0.0;

        } else {

            ce_x_2[i] = sige_x_PML_2[i] * (be_x_2[i] - 1.0) /
                    (sige_x_PML_2[i] + kappae_x_PML_2[i] * alphae_x_PML_2[i])
                    / kappae_x_PML_2[i];
        }

    }

    for (i = 0; i < nxPML_2 - 1; ++i) {

        sigh_x_PML_2[i] = sig_x_max * pow(((nxPML_2 - 1 - i - 0.5)
                / (nxPML_2 - 1.0)), m);
        alphah_x_PML_2[i] = alpha_x_max * pow(((i + 1 - 0.5)
                / (nxPML_2 - 1.0)), ma);
        kappah_x_PML_2[i] = 1.0 + (kappa_x_max - 1.0) *
                pow(((nxPML_2 - 1 - i - 0.5) / (nxPML_2 - 1.0)), m);
        bh_x_2[i] = exp(-(sigh_x_PML_2[i] / kappah_x_PML_2[i] +
                alphah_x_PML_2[i]) * dt / eps_0);
        ch_x_2[i] = sigh_x_PML_2[i] * (bh_x_2[i] - 1.0) /
                (sigh_x_PML_2[i] + kappah_x_PML_2[i] * alphah_x_PML_2[i])
                / kappah_x_PML_2[i];
    }

    for (j = 0; j < nyPML_1; ++j) {

        sige_y_PML_1[j] = sig_y_max * pow(((nyPML_1 - 1 - j)
                / (nyPML_1 - 1.0)), m);
        alphae_y_PML_1[j] = alpha_y_max * pow(((float) j
                / (nyPML_1 - 1.0)), ma);
        kappae_y_PML_1[j] = 1.0 + (kappa_y_max - 1.0) *
                pow((nyPML_1 - 1 - j) / (nyPML_1 - 1.0), m);
        be_y_1[j] = exp(-(sige_y_PML_1[j] / kappae_y_PML_1[j] +
                alphae_y_PML_1[j]) * dt / eps_0);

        if ((sige_y_PML_1[j] == 0.0) &&
                (alphae_y_PML_1[j] == 0.0) && (j == nyPML_1 - 1)) {

            ce_y_1[j] = 0.0;

        } else {

            ce_y_1[j] = sige_y_PML_1[j] * (be_y_1[j] - 1.0) /
                    (sige_y_PML_1[j] + kappae_y_PML_1[j] * alphae_y_PML_1[j])
                    / kappae_y_PML_1[j];
        }
    }

    for (j = 0; j < nyPML_1 - 1; ++j) {

        sigh_y_PML_1[j] = sig_y_max * pow(((nyPML_1 - 1 - j - 0.5)
                / (nyPML_1 - 1.0)), m);
        alphah_y_PML_1[j] = alpha_y_max * pow(((j + 1 - 0.5)
                / (nyPML_1 - 1.0)), ma);
        kappah_y_PML_1[j] = 1.0 + (kappa_y_max - 1.0) *
                pow(((nyPML_1 - 1 - j - 0.5) / (nyPML_1 - 1.0)), m);
        bh_y_1[j] = exp(-(sigh_y_PML_1[j] / kappah_y_PML_1[j] +
                alphah_y_PML_1[j]) * dt / eps_0);
        ch_y_1[j] = sigh_y_PML_1[j] * (bh_y_1[j] - 1.0) /
                (sigh_y_PML_1[j] + kappah_y_PML_1[j] * alphah_y_PML_1[j])
                / kappah_y_PML_1[j];
    }

    for (j = 0; j < nyPML_2; ++j) {

        sige_y_PML_2[j] = sig_y_max * pow(((nyPML_2 - 1 - j)
                / (nyPML_2 - 1.0)), m);
        alphae_y_PML_2[j] = alpha_y_max * pow(((float) j
                / (nyPML_2 - 1.0)), ma);
        kappae_y_PML_2[j] = 1.0 + (kappa_y_max - 1.0) *
                pow((nyPML_2 - 1 - j) / (nyPML_2 - 1.0), m);
        be_y_2[j] = exp(-(sige_y_PML_2[j] / kappae_y_PML_2[j] +
                alphae_y_PML_2[j]) * dt / eps_0);

        if ((sige_y_PML_2[j] == 0.0) &&
                (alphae_y_PML_2[j] == 0.0) && (j == nyPML_2 - 1)) {

            ce_y_2[j] = 0.0;

        } else {

            ce_y_2[j] = sige_y_PML_2[j] * (be_y_2[j] - 1.0) /
                    (sige_y_PML_2[j] + kappae_y_PML_2[j] * alphae_y_PML_2[j])
                    / kappae_y_PML_2[j];
        }
    }

    for (j = 0; j < nyPML_2 - 1; ++j) {

        sigh_y_PML_2[j] = sig_y_max * pow(((nyPML_2 - 1 - j - 0.5)
                / (nyPML_2 - 1.0)), m);
        alphah_y_PML_2[j] = alpha_y_max * pow(((j + 1 - 0.5)
                / (nyPML_2 - 1.0)), ma);
        kappah_y_PML_2[j] = 1.0 + (kappa_y_max - 1.0) *
                pow(((nyPML_2 - 1 - j - 0.5) / (nyPML_2 - 1.0)), m);
        bh_y_2[j] = exp(-(sigh_y_PML_2[j] / kappah_y_PML_2[j] +
                alphah_y_PML_2[j]) * dt / eps_0);
        ch_y_2[j] = sigh_y_PML_2[j] * (bh_y_2[j] - 1.0) /
                (sigh_y_PML_2[j] + kappah_y_PML_2[j] * alphah_y_PML_2[j])
                / kappah_y_PML_2[j];
    }

    for (k = 0; k < nzPML_1; ++k) {

        sige_z_PML_1[k] = sig_z_max * pow(((nzPML_1 - 1 - k)
                / (nzPML_1 - 1.0)), m);
        alphae_z_PML_1[k] = alpha_z_max * pow(((float) k
                / (nzPML_1 - 1.0)), ma);
        kappae_z_PML_1[k] = 1.0 + (kappa_z_max - 1.0) *
                pow((nzPML_1 - 1 - k) / (nzPML_1 - 1.0), m);
        be_z_1[k] = exp(-(sige_z_PML_1[k] / kappae_z_PML_1[k] +
                alphae_z_PML_1[k]) * dt / eps_0);

        if ((sige_z_PML_1[k] == 0.0) &&
                (alphae_z_PML_1[k] == 0.0) && (k == nzPML_1 - 1)) {

            ce_z_1[k] = 0.0;

        } else {

            ce_z_1[k] = sige_z_PML_1[k] * (be_z_1[k] - 1.0) /
                    (sige_z_PML_1[k] + kappae_z_PML_1[k] * alphae_z_PML_1[k])
                    / kappae_z_PML_1[k];
        }
    }

    for (k = 0; k < nzPML_1 - 1; ++k) {

        sigh_z_PML_1[k] = sig_z_max * pow(((nzPML_1 - 1 - k - 0.5)
                / (nzPML_1 - 1.0)), m);
        alphah_z_PML_1[k] = alpha_z_max * pow(((k + 1 - 0.5)
                / (nzPML_1 - 1.0)), ma);
        kappah_z_PML_1[k] = 1.0 + (kappa_z_max - 1.0) *
                pow(((nzPML_1 - 1 - k - 0.5) / (nzPML_1 - 1.0)), m);
        bh_z_1[k] = exp(-(sigh_z_PML_1[k] / kappah_z_PML_1[k] +
                alphah_z_PML_1[k]) * dt / eps_0);
        ch_z_1[k] = sigh_z_PML_1[k] * (bh_z_1[k] - 1.0) /
                (sigh_z_PML_1[k] + kappah_z_PML_1[k] * alphah_z_PML_1[k])
                / kappah_z_PML_1[k];
    }

    for (k = 0; k < nzPML_2; ++k) {

        sige_z_PML_2[k] = sig_z_max * pow(((nzPML_2 - 1 - k)
                / (nzPML_2 - 1.0)), m);
        alphae_z_PML_2[k] = alpha_z_max * pow(((float) k
                / (nzPML_2 - 1.0)), ma);
        kappae_z_PML_2[k] = 1.0 + (kappa_z_max - 1.0) *
                pow((nzPML_2 - 1 - k) / (nzPML_2 - 1.0), m);
        be_z_2[k] = exp(-(sige_z_PML_2[k] / kappae_z_PML_2[k] +
                alphae_z_PML_2[k]) * dt / eps_0);

        if ((sige_z_PML_2[k] == 0.0) &&
                (alphae_z_PML_2[k] == 0.0) && (k == nzPML_2 - 1)) {

            ce_z_2[k] = 0.0;

        } else {

            ce_z_2[k] = sige_z_PML_2[k] * (be_z_2[k] - 1.0) /
                    (sige_z_PML_2[k] + kappae_z_PML_2[k] * alphae_z_PML_2[k])
                    / kappae_z_PML_2[k];
        }
    }

    for (k = 0; k < nzPML_2 - 1; ++k) {

        sigh_z_PML_2[k] = sig_z_max * pow(((nzPML_2 - 1 - k - 0.5)
                / (nzPML_2 - 1.0)), m);
        alphah_z_PML_2[k] = alpha_z_max * pow(((k + 1 - 0.5)
                / (nzPML_2 - 1.0)), ma);
        kappah_z_PML_2[k] = 1.0 + (kappa_z_max - 1.0) *
                pow(((nzPML_2 - 1 - k - 0.5) / (nzPML_2 - 1.0)), m);
        bh_z_2[k] = exp(-(sigh_z_PML_2[k] / kappah_z_PML_2[k] +
                alphah_z_PML_2[k]) * dt / eps_0);
        ch_z_2[k] = sigh_z_PML_2[k] * (bh_z_2[k] - 1.0) /
                (sigh_z_PML_2[k] + kappah_z_PML_2[k] * alphah_z_PML_2[k])
                / kappah_z_PML_2[k];
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  DENOMINATORS FOR FIELD UPDATES
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ii = nxPML_2 - 2;

    for (i = 0; i < Imax - 1; ++i) {

        if (i < nxPML_1 - 1) {

            den_hx[i] = 1.0 / (kappah_x_PML_1[i] * dx);

        } else if (i >= Imax - nxPML_2) {

            den_hx[i] = 1.0 / (kappah_x_PML_2[ii] * dx);
            ii = ii - 1;

        } else {

            den_hx[i] = 1.0 / dx;
        }
    }

    jj = nyPML_2 - 2;

    for (j = 0; j < Jmax - 1; ++j) {

        if (j < nyPML_1 - 1) {

            den_hy[j] = 1.0 / (kappah_y_PML_1[j] * dy);

        } else if (j >= Jmax - nyPML_2) {

            den_hy[j] = 1.0 / (kappah_y_PML_2[jj] * dy);
            jj = jj - 1;

        } else {

            den_hy[j] = 1.0 / dy;
        }
    }

    kk = nzPML_2 - 2;

    for (k = 1; k < Kmax - 1; ++k) {

        if (k < nzPML_1) {

            den_hz[k] = 1.0 / (kappah_z_PML_1[k - 1] * dz);

        } else if (k >= Kmax - nzPML_2) {

            den_hz[k] = 1.0 / (kappah_z_PML_2[kk] * dz);
            kk = kk - 1;

        } else {

            den_hz[k] = 1.0 / dz;
        }
    }

    ii = nxPML_2 - 1;

    for (i = 0; i < Imax - 1; ++i) {

        if (i < nxPML_1) {

            den_ex[i] = 1.0 / (kappae_x_PML_1[i] * dx);

        } else if (i >= Imax - nxPML_2) {

            den_ex[i] = 1.0 / (kappae_x_PML_2[ii] * dx);
            ii = ii - 1;

        } else {

            den_ex[i] = 1.0 / dx;
        }
    }

    jj = nyPML_2 - 1;

    for (j = 0; j < Jmax - 1; ++j) {

        if (j < nyPML_1) {

            den_ey[j] = 1.0 / (kappae_y_PML_1[j] * dy);

        } else if (j >= Jmax - nyPML_2) {

            den_ey[j] = 1.0 / (kappae_y_PML_2[jj] * dy);
            jj = jj - 1;

        } else {

            den_ey[j] = 1.0 / dy;
        }
    }

    kk = nzPML_2 - 1;

    for (k = 0; k < Kmax - 1; ++k) {

        if (k < nzPML_1) {

            den_ez[k] = 1.0 / (kappae_z_PML_1[k] * dz);

        } else if (k >= Kmax - 1 - nzPML_2) {

            den_ez[k] = 1.0 / (kappae_z_PML_2[kk] * dz);
            kk = kk - 1;

        } else {

            den_ez[k] = 1.0 / dz;
        }
    }
}

void compute() {

    short id;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  BEGIN TIME STEP
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    printf("\nBegin time-stepping...\n");

    for (n = 1; n <= nMax; ++n) {

        printf("Ez at time step %d at (25, 40, 12) :  %f\n", n, Ez[25][40][12]);

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Hx
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 1; k < Kmax - 1; ++k) {
            for (i = 0; i < Imax - 1; ++i) {
                for (j = 0; j < Jmax - 1; ++j) {
                    Hx[i][j][k] = DA * Hx[i][j][k] + DB *
                            ((Ez[i][j][k] - Ez[i][j + 1][k]) * den_hy[j] +
                            (Ey[i][j][k] - Ey[i][j][k - 1]) * den_hz[k]);
                }
            }
            // update PML for Hx inner
            UpdatePMLForHxIn(k, Hx, Ey, Ez, DB);
        }
        UpdatePMLForHxOut(Hx,Ey,Ez,DB);

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Hy
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 1; k < Kmax - 1; ++k) {

            for (i = 0; i < Imax - 1; ++i) {

                for (j = 0; j < Jmax - 1; ++j) {

                    Hy[i][j][k] = DA * Hy[i][j][k] + DB *
                            ((Ez[i + 1][j][k] - Ez[i][j][k]) * den_hx[i] +
                            (Ex[i][j][k - 1] - Ex[i][j][k]) * den_hz[k]);
                }
            }
	    UpdatePMLForHyIn(k,Hy,Ez,Ex,DB);
        }

	UpdatePMLForHyOut(Hy,Ez,Ex,DB);

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Hz
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 0; k < Kmax - 1; ++k) {
            for (i = 0; i < Imax - 1; ++i) {
                for (j = 0; j < Jmax - 1; ++j) {
                    Hz[i][j][k] = DA * Hz[i][j][k] + DB
                            * ((Ey[i][j][k] - Ey[i + 1][j][k]) * den_hx[i] +
                            (Ex[i][j + 1][k] - Ex[i][j][k]) * den_hy[j]);
                }
            }
	    UpdatePMLForHz(k,Hz,Ex,Ey,DB);
        }

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Ex
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 0; k < Kmax - 1; ++k) {

            for (i = 0; i < Imax - 1; ++i) {

                for (j = 1; j < Jmax - 1; ++j) {

                    id = ID1[i][j][k];
                    if (id == 1) { // PEC

                        Ex[i][j][k] = 0;

                    } else {

                        Ex[i][j][k] = CA[id] * Ex[i][j][k] + CB[id] *
                                ((Hz[i][j][k] - Hz[i][j - 1][k]) * den_ey[j] +
                                (Hy[i][j][k] - Hy[i][j][k + 1]) * den_ez[k]);
                    }
                }
            }

	    UpdatePMLForExIn(k,Ex,Hy,Hz,CB,ID1);
        }

	    UpdatePMLForExOut(Ex,Hy,Hz,CB,ID1);

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Ey
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 0; k < Kmax - 1; ++k) {

            for (i = 1; i < Imax - 1; ++i) {

                for (j = 0; j < Jmax - 1; ++j) {

                    id = ID2[i][j][k];
                    if (id == 1) { // PEC
                        Ey[i][j][k] = 0;
                    } else {
                        Ey[i][j][k] = CA[id] * Ey[i][j][k] + CB[id] *
                                ((Hz[i - 1][j][k] - Hz[i][j][k]) * den_ex[i] +
                                (Hx[i][j][k + 1] - Hx[i][j][k]) * den_ez[k]);
                    }
                }
            }
UpdatePMLForEyIn(k,Ey,Hz,Hx,CB,ID2);
        }

UpdatePMLForEyOut(Ey,Hz,Hx,CB,ID2);

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //  UPDATE Ez
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (k = 1; k < Kmax - 1; ++k) {
            for (i = 1; i < Imax - 1; ++i) {
                for (j = 1; j < Jmax - 1; ++j) {
                    id = ID3[i][j][k];
                    if (id == 1) { // PEC
                        Ez[i][j][k] = 0;
                    } else {
                        Ez[i][j][k] = CA[id] * Ez[i][j][k] + CB[id]
                                * ((Hy[i][j][k] - Hy[i - 1][j][k]) * den_ex[i] +
                                (Hx[i][j - 1][k] - Hx[i][j][k]) * den_ey[j]);
                    }
                }
            }
		UpdatePMLForEz(k,Ez,Hx,Hy,CB,ID3);
        }


        //-----------------------------------------------------------
        //   Apply a point source (Soft)
        //-----------------------------------------------------------
        i = 25;
        j = 63;
        k = 12;
        source = amp * -2.0 * ((n * dt - tO) / tw)
                * exp(-pow(((n * dt - tO) / tw), 2)); //Differentiated Gaussian pulse

        Ez[i][j][k] = Ez[i][j][k] - CB[ID3[i][j][k]] * source;


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
    ID1[I][J][K] = mType;
    ID1[I][J][K + 1] = mType;
    ID1[I][J + 1][K + 1] = mType;
    ID1[I][J + 1][K] = mType;

    //set face 2 (for EY)
    ID2[I][J][K] = mType;
    ID2[I + 1][J][K] = mType;
    ID2[I + 1][J][K + 1] = mType;
    ID2[I][J][K + 1] = mType;

    //set face 3 (for EZ)
    ID3[I][J][K] = mType;
    ID3[I + 1][J][K] = mType;
    ID3[I + 1][J + 1][K] = mType;
    ID3[I][J + 1][K] = mType;
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
            fprintf(ptr, "%f\t", sqrt(pow(Ex[i][j][ksource], 2) +
                    pow(Ey[i][j][ksource], 2) + pow(Ez[i][j][ksource], 2)));

            //	fprintf(ptr, "%f\t", Ex[i][j][ksource]);//Ex
            //	fprintf(ptr, "%f\t", Ey[i][j][ksource]);//Ey
            //	fprintf(ptr, "%f\t", Ez[i][j][ksource]);//Ez
            //	fprintf(ptr, "%f\t", Hx[i][j][ksource]);//Hx
            //	fprintf(ptr, "%f\t", Hy[i][j][ksource]);//Hy
            //	fprintf(ptr, "%f\t", Hz[i][j][ksource]);//Hz

            // |H|
            //	fprintf(ptr, "%f\t", sqrt(pow(Hx[i][j][ksource], 2) +
            //		pow(Hy[i][j][ksource], 2) + pow( Hz[i][j][ksource], 2)));

        }
        fprintf(ptr, "\n");
    }

    fclose(ptr);

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// END OF PROGRAM CPMLFDTD3D
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
