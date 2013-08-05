#pragma once

#include "datastruct.h"
#include "fdtd.h"

class Density
{
public:
	Density(void);
	~Density(void);
	void SetPlasmaVar(MyDataF _rei, MyDataF _vm, MyDataF _p, int _ftype);
private:
	//  Specify the dipole Boundaries(A cuboidal rode- NOT as a cylinder)
    unsigned istart, iend, jstart;
    unsigned jend, kstart, kend;
    //how many fine grids per coarse grid 
    unsigned neGrid;
    //initial plasma value
    MyDataF Ne0;


    int niutype;
    //Fine Grid size
    MyDataF dsf;
    //time step of plasma
    MyDataF dtf;
    int neSkipStep;

    // temporary variables that often used
    MyDataF half_dt;
    MyDataF half_e;
	MyDataF eMDtDiv2DivEps0;
	MyDataF dtDivEps0DivDx;
	MyDataF e2Dt2Div4DivEps0DivMe;
	MyDataF dtDivEps0DivDy;
	MyDataF dtDivEps0DivDz;
    MyDataF Coeff_velocity;
    unsigned halfNeGrid;

    //plasma variables
    MyDataF vm; //collision frequency
    MyDataF p; // air pressure
    MyDataF De;
    MyDataF Da;
    MyDataF rei;
    MyDataF mu_i;
    MyDataF mu_e;
    MyDataF gamma;
    MyDataF a;
    MyDataF alpha;

    // update coefficients
    //    MyDataF Chxey,Chxez,Chyez,Chyex,Chzex,Chzey;
    data3d<MyDataF> Cezvz, Ceyvy, Cexvx;
    data3d<MyDataF> Cvxex_guassian, Cvyey_guassian, Cvzez_guassian;
    MyDataF Cvxex, Cvyey, Cvzez;

    //Plasma
    data3d<MyDataF> Ne, Ne_pre;
    //
    data3d<MyDataF> Erms;
    // Beta
    data3d<MyDataF> beta;

    // collision frequency
    data3d<MyDataF> Nu_c;
    //
    data3d<MyDataF> Vx;
    data3d<MyDataF> Vy;
    data3d<MyDataF> Vz;

    //initials
    void initCoeff(const fdtd &mfdtd);
    void initDensity(const fdtd &mfdtd);
    void createCoeff(int srcType,data3d<MyDataF> Ex,data3d<MyDataF> Ey,data3d<MyDataF> Ez) ;
    void updateCoeff(const fdtd &mfdtd);
    void updateBeta(int srcType);

    // Erms or Eeff operation
    void IntegerEeff();
    void UpdateErms(int srcType,data3d<MyDataF> Ex,data3d<MyDataF> Ey,data3d<MyDataF> Ez) ;
    void updateCollisionFrequency();
    void InterpErms();
    void UpdateDensity(int srcType);
    void UpdateVeloity(void);
    void WallCircleBound(data3d<MyDataF> &stru);
};

