/* 
 * File:   sourceTest.cpp
 * Author: skiloop
 *
 * Created on 2013年12月16日, 下午11:32
 */

#include <cstdlib>
#include <iostream>
#include "../src/SineWaveSource.h"
#include "../src/GaussianWaveSource.h"
#include "../src/CosineGaussianWave.h"



using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    if (argc < 2)return 0;

    int n;
    n = atoi(argv[1]);

    if (n < 0)n = 100;
    MyDataF f = 1e9;
    MyDataF T = 1.0 / f;
    MyDataF omega = 2.0 * M_PI / T;
    SineWaveSource sineSource(omega);
    GaussianWaveSource gsource(f);
    CosineGaussianWave cosGaussian(f, 0.5*f);
    double dt = T / n;
    for (int i = 0; i <= 50 * n; i++) {
        MyDataF v = cosGaussian.valueAtTime(dt * i);
        cout << v << endl;
    }
    return 0;
}

