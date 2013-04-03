#pragma once

#include "cpml.h"

class fdtd {
public:
#ifdef WITH_DENSITY
    fdtd(unsigned _nmax = 500, unsigned _imax = 51, unsigned _jmax = 126, unsigned _kmax = 26,
            MyDataF _tw = 53.0e-12, MyDataF _dx = 1e-3, MyDataF _dy = 1e-3, MyDataF _dz = 1e-3,
            unsigned _amp = 1000, unsigned _savemodulus = 10, unsigned _ksource = 12,
            unsigned _m = 3, unsigned _ma = 1, unsigned _nmatrial = 50, unsigned _neGrid = 10);
    
    void SetPlasmaVar(MyDataF _rei, MyDataF _vm, MyDataF _p, int _ftype);
#else
        fdtd(unsigned _nmax = 500, unsigned _imax = 51, unsigned _jmax = 126, unsigned _kmax = 26,
            MyDataF _tw = 53.0e-12, MyDataF _dx = 1e-3, MyDataF _dy = 1e-3, MyDataF _dz = 1e-3,
            unsigned _amp = 1000, unsigned _savemodulus = 10, unsigned _ksource = 12,
            unsigned _m = 3, unsigned _ma = 1, unsigned _nmatrial = 50);
#endif
    ~fdtd(void);
    //Function prototype definitions
    void initialize(); //Memeory initialization
    void setUp(); //Coefficients, parameters etc will get computed
    void compute(); //E & H Field update equation
    void StartUp();    
private:
    void buildObject(); //Creates the object geometry
    void yeeCube(unsigned, unsigned, unsigned, unsigned); //Sets material properties to a cell
    void writeField(unsigned); //Writes output
    void buildSphere(); //Builds a spherical object
    void buildDipole(); //Builds a dipole
    void putvars();
private:
    //  Specify Number of Time Steps and Grid Size Parameters
    unsigned nMax; // total number of time steps

    // grid size corresponding to the number of Ez field components
    unsigned Imax;
    unsigned Jmax;
    unsigned Kmax;
    //  Specify the Impulsive Source (Differentiated Gaussian) parameters
    MyDataF tw; //pulse width
    MyDataF dt, dx, dy, dz;
    
    MyDataF t0; //delay
    MyDataF source; //Differentiated Gaussian source
    MyDataF amp; // Amplitude
    MyDataF omega; //wave angle frequency

    // Specify the Time Step at which the data has to be saved for Visualization
    unsigned save_modulus;
    // Output recording pounsigned
    unsigned ksource;
    // Specify the CPML Order and Other Parameters:
    unsigned m, ma;

    //  Specify the dipole Boundaries(A cuboidal rode- NOT as a cylinder)
    unsigned istart, iend, jstart;
    unsigned jend, kstart, kend;
    
    // source position
    unsigned isp,jsp,ksp;

    //Max number of materials allowed
    unsigned numMaterials;

    // H & E Field components
    data3d<MyDataF> Hx;
    data3d<MyDataF> Hy;
    data3d<MyDataF> Hz;
    data3d<MyDataF> Ex;
    data3d<MyDataF> Ey;
    data3d<MyDataF> Ez;

    data3d<unsigned> ID1; //medium definition array for Ex
    data3d<unsigned> ID2; //medium definition array for Ey
    data3d<unsigned> ID3; //medium definition array for Ez

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
    
#ifdef WITH_DENSITY

    //how many fine grids per coarse grid 
    unsigned neGrid;
    int niutype;
    //Fine Grid size
    MyDataF dsf;
    //time step of plasma
    MyDataF dtf;
    
    //plasma variables
    MyDataF vm;
    MyDataF p;
    MyDataF De;
    MyDataF Da;
    MyDataF rei;
    MyDataF mu_i;
    MyDataF mu_e; 
    
    //Plasma
    data3d<MyDataF> Ne, Ne_pre;
    //
    data3d<MyDataF> Erms;
    //intial plasma value
    MyDataF Ne0;
    
    //
    data3d<MyDataF> Vx;
    data3d<MyDataF> Vy;
    data3d<MyDataF> Vz;
    
    int UpdateErms(void);
    int UpdateDensity(void);
    int UpdateVeloity(void);
    void WallCircleBound(data3d<MyDataF> &stru);
#endif
    cpml<MyDataF,unsigned int> pml;

};

