
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
#include<iostream>
#include<stdlib.h>

#include "cpml.h"

using namespace std;
//  Fundamental Constants (MKS units)
MyDataF pi = 3.14159265358979;
MyDataF C = 2.99792458E8;
MyDataF mu_0;
MyDataF eps_0;

//  Specify Material Relative Permittivity and Conductivity
MyDataF epsR = 1.0; //free space

//  Specify Number of Time Steps and Grid Size Parameters
int nMax = 500; // total number of time steps

// grid size corresponding to the number of Ez field components
int Imax = 51;
int Jmax = 126;
int Kmax = 26;

//  Specify Grid Cell Size in Each Direction and Calculate the
//  Resulting Courant-Stable Time Step
MyDataF dx = 1.0E-3;
MyDataF dy = 1.0E-3;
MyDataF dz = 1.0E-3; // cell size in each direction
// time step increment
MyDataF dt;

//  Specify the Impulsive Source (Differentiated Gaussian) parameters
MyDataF tw = 53.0E-12; //pulse width
MyDataF tO; //delay
MyDataF source; //Differentiated Gaussian source
MyDataF amp = 1000; // Amplitude

//Specify the Time Step at which the data has to be saved for Visualization
int save_modulus = 10;

//  Specify the dipole Boundaries(A cuboidal rode- NOT as a cylinder)
int istart, iend, jstart;
int jend, kstart, kend;

//Output recording point
int ksource = 12;

//  Specify the CPML Order and Other Parameters:
int m = 3, ma = 1;

//Loop indices
int i, j, ii, jj, k, kk, n;


// H & E Field components
data3d<MyDataF> Hx;
data3d<MyDataF> Hy;
data3d<MyDataF> Hz;
data3d<MyDataF> Ex;
data3d<MyDataF> Ey;
data3d<MyDataF> Ez;

data3d<short> ID1; //medium definition array for Ex
data3d<short> ID2; //medium definition array for Ey
data3d<short> ID3; //medium definition array for Ez

//Max number of materials allowed
int numMaterials = 50;

//permittivity, permeability and conductivity of diffrent materials
MyDataF *epsilon;
MyDataF *mu;
MyDataF *sigma;

//E field update coefficients
MyDataF *CA;
MyDataF *CB;

//H field update coefficients
MyDataF DA;
MyDataF DB;

//Function prototype definitions
void initialize(); //Memeory initialization
void setUp(); //Coefficients, parameters etc will get computed
void compute(cpml &pml); //E & H Field update equation
void buildObject(); //Creates the object geometry
void yeeCube(int, int, int, unsigned); //Sets material properties to a cell
void writeField(int); //Writes output
void buildSphere(); //Builds a spherical object
void buildDipole(); //Builds a dipole
void putvars();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int main() {
    cpml pml(m, ma);
    initialize();
    pml.Initial(Imax, Jmax, Kmax, 11);
    setUp();
    buildObject();
    putvars();
    compute(pml);
    return 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void initialize() {

    mu_0 = 4.0 * pi * 1.0E-7;
    eps_0 = 1.0 / (C * C * mu_0);

    //Dynamic memory allocation
    epsilon = (MyDataF *) malloc((numMaterials) * sizeof (MyDataF));

    for (i = 0; i < numMaterials; i++) {

        epsilon[i] = eps_0;
    }

    mu = (MyDataF *) malloc((numMaterials) * sizeof (MyDataF));

    for (i = 0; i < numMaterials; i++) {

        mu[i] = mu_0;
    }

    sigma = (MyDataF *) malloc((numMaterials) * sizeof (MyDataF));

    for (i = 0; i < numMaterials; i++) {

        sigma[i] = 0.0;
    }

    CA = (MyDataF *) malloc((numMaterials) * sizeof (MyDataF));

    for (i = 0; i < numMaterials; i++) {

        CA[i] = 0.0;
    }

    CB = (MyDataF *) malloc((numMaterials) * sizeof (MyDataF));

    for (i = 0; i < numMaterials; i++) {

        CB[i] = 0.0;
    }

    Ez.CreateStruct(Imax, Jmax, Kmax,0);
    Ey.CreateStruct(Imax, Jmax - 1, Kmax - 1,0);
    Ex.CreateStruct(Imax - 1, Jmax, Kmax - 1,0);

    Hx.CreateStruct(Imax, Jmax - 1, Kmax,0);
    Hy.CreateStruct(Imax - 1, Jmax, Kmax,0);
    Hz.CreateStruct((Imax - 1), (Jmax - 1), (Kmax - 1),0);

    ID1.CreateStruct(Imax,Jmax,Kmax,0);
    ID2.CreateStruct(Imax,Jmax,Kmax,0);
    ID3.CreateStruct(Imax,Jmax,Kmax,0);

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

    printf("\nTIme step = %e", dt);
    printf("\n Number of steps = %d", nMax);
    printf("\n Total Simulation time = %e Seconds", nMax * dt);

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SET CPML PARAMETERS IN EACH DIRECTION
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void compute(cpml &pml) {

    unsigned id;
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
            pml.UpdatePMLForHxIn(k,Hx,Ey,Ez,DB);
        }
        pml.UpdatePMLForHxOut(Hx, Ey, Ez, DB);
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
            pml.UpdatePMLForHyIn(k,Hy, Ez, Ex, DB);
        }
        pml.UpdatePMLForHyOut(Hy, Ez, Ex, DB);
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
            pml.UpdatePMLForHz(k,Hz, Ex, Ey, DB);
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
            pml.UpdatePMLForExIn(k,Ex, Hy, Hz, CB, ID1.p);
        }
        pml.UpdatePMLForExOut(Ex, Hy, Hz, CB, ID1.p);
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
            pml.UpdatePMLForEyIn(k,Ey, Hz, Hx, CB, ID2.p);
        }
        pml.UpdatePMLForEyOut(Ey, Hz, Hx, CB, ID2.p);

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
            pml.UpdatePMLForEz(k,Ez, Hx, Hy, CB, ID3.p);
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

    MyDataF dist; //distance
    MyDataF rad = 8; //(MyDataF)Imax / 5.0; // sphere radius
    MyDataF sc = (MyDataF) Imax / 2.0; //sphere centre
    //MyDataF rad2 = 0.3; //(MyDataF)Imax / 5.0 - 3.0; // sphere radius

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

void yeeCube(int I, int J, int K, unsigned mType) {

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
    char fileBaseName[100] = "E_Field_";
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

void putvars() {
    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;
    cout << "dz = " << dz << endl;
    cout << "(Imax,Jmax,Kmax) = (" << Imax << "," << Jmax << "," << Kmax << ")" << endl;
    // time step increment
    cout << "dt = " << dt << endl;

    //  Specify the Impulsive Source (Differentiated Gaussian) parameters
    cout << "tw = " << tw << endl; //pulse width
    cout << "tO = " << tO << endl; //delay
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
    cout << " m  = " << m << endl;
    cout << " ma = " << ma << endl;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// END OF PROGRAM CPMLFDTD3D
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
