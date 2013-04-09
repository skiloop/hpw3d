
#include<iostream>
//#define WITH_DENSITY
#include "fdtd.h"
MyDataF eps_0, epsR;
MyDataF mu_0;
MyDataF dt, dx, dy, dz;
MyDataF pi, C;
MyDataF me, e;
MyDataF T;
void initComData();

int main() {

    initComData();
    fdtd hpw;
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
}
