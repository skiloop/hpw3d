
#include<iostream>
//#define WITH_DENSITY
#include "fdtd.h"
MyDataF eps_0, epsR;
MyDataF mu_0;
MyDataF dx, dy, dz;
MyDataF pi, C;
MyDataF me, e;
MyDataF tw;
MyDataF omika;
MyDataF T;// ns
MyDataF Amp;
unsigned pmlw;
void initComData();

int main() {

    initComData();
    fdtd hpw(500,50,100,26,tw,dx,dy,dz,Amp,10,12,4,1,pmlw);
#ifdef WITH_DENSITY
    hpw.SetPlasmaVar(0,760*5.3E9,760,0);
#endif
    //hpw.initialize();
    hpw.StartUp();
    return 0;
}

void initComData() {
    pi = 3.14159265358979;
    C = 2.99792458E8;
    mu_0 = 4.0 * pi * 1.0E-7;
    eps_0 = 1.0 / (C * C * mu_0);
    epsR = 1.0;
    me=9.110e-31;
    e=1.602e-19;
    dx=dy=dz=1e-3;
    Amp=1e-10;
    tw=2e-9;
    pmlw=12;
    omika=2*pi*C/150/dx;
}
