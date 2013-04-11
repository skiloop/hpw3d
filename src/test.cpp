
#include<iostream>
//#define WITH_DENSITY
#include "fdtd.h"
MyDataF eps_0, epsR;
MyDataF mu_0;
MyDataF dx, dy, dz;
MyDataF pi, C;
MyDataF me, e;
MyDataF T;
MyDataF omika;
void initComData();

int main() {

    initComData();
    fdtd hpw(500,150,150,150,53.0e-12,dx,dy,dz,1000,10,12,3,1,6);
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
    omika=2*pi*C/100/dx;
}
