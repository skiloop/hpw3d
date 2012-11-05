
#include "InonizationFormula.h"
#include <math.h>
//////////////////////////////////////////////////////////////////////////
//int sign(double val)
//if val>0 return 1
//if val<0 return -1
//return 0 if val==0

int sign(double val) {
    return (val < 0 ? -1 : (val > 0 ? 1 : 0));
}

//////////////////////////////////////////////////////////////////////////
// Nikonov formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Nikonov formula

double Alpha_Nikonov(double E, double P) {
    double EDivP = fabs(E) / P;
    if (EDivP < 108.0) {
        return 3.9 * P * exp(-213 / EDivP);
    }
    else {
        return 14.5 * P * exp(-316 / EDivP);
    }
}
// Calculate Eta by Nikonov formula

double Eta_Nikonov(double E, double P) {
    double EDivP;
    EDivP = fabs(E) / P;
    if (EDivP < 50.0) {
        double val1 = 4.47e-3 * (EDivP)*(EDivP);
        if (EDivP >= 10.0) {
            return val1;
        }
        else {
            double val2 = 4.47 / EDivP;
            return (val1 > val2 ? val1 : val2);
        }
    }
    else {
        double temp = sqrt(EDivP);
        return (EDivP <= 90 ? 1.58 * temp : 142 / temp);
    }

}
// Calculate We by Nikonov formula

double We_Nikonov(double E, double P) {
    return -0.0382 * E - 2.9e5 * E / P;
}
// Niu_a

double Niu_a_Nikonov(double E, double P) {
    return Eta_Nikonov(E, P) * fabs(We_Nikonov(E, P));
}
// Niu_i

double Niu_i_Nikonov(double E, double P) {
    return Alpha_Nikonov(E, P) * fabs(We_Nikonov(E, P));
}

//////////////////////////////////////////////////////////////////////////
//Morrow and Lowke formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Morrow and Lowke formula

double Alpha_MorrowAndLowke(double E, double N) {
    double edn = fabs(E) / N / 1e-15;

    if (edn > 1.5) {
        return 2e-16 * N * exp(-7.248 / edn);
    }
    else {
        return 6.619e-17 * N * exp(-5.593 / edn);
    }
}
// Calculate Eta by Morrow and Lowke formula

double Eta_MorrowAndLowke(double E, double N) {
    double edn;
    edn = fabs(E) / N / 1e-16;

    if (edn > 0.12) {
        if (edn < 10.50)
            return N * ((6.089e-20 * edn - 2.893e-19) + N * 4.47778e-59 * pow(edn * 1e-16, -1.2749));
        else
            return N * ((8.889e-21 * edn + 2.567e-19) + N * 4.47778e-59 * pow(edn * 1e-16, -1.2749));
    }
    else {
        if (edn < 0)
            return 0;
        else
            return 106.81;
    }
}
// Calculate We by Morrow and Lowke formula

double We_MorrowAndLowke(double E, double N) {
    double edn = fabs(E) / N / 1e-16;

    if (edn > 1.0) {
        if (edn <= 20.00)
            return -sign(E)*(1.03 * edn + 1.3)*1e-10;
        else
            return -sign(E)*(7.4e5 * edn + 7.1e6)*1e-16;
    }
    else {
        if (edn <= 0.26)
            return -sign(E)*(6.87e6 * edn + 3.38e4)*1e-16;
        else
            return -sign(E)*(7.2973 * edn + 16.3)*1e-11;
        ;
    }
}
// Niu_a

double Niu_a_MorrowAndLowke(double E, double N) {
    return Eta_MorrowAndLowke(E, N) * fabs(We_MorrowAndLowke(E, N));
}
// Niu_i

double Niu_i_MorrowAndLowke(double E, double N) {
    return Alpha_MorrowAndLowke(E, N) * fabs(We_MorrowAndLowke(E, N));
}

//////////////////////////////////////////////////////////////////////////
//Kang formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Kang formula

double Alpha_Kang(double E) {
    return 0.0035 * exp(-1.65e5 / E);
}
// Calculate Eta by Kang formula

double Eta_Kang(double E) {
    return 0.15 * exp(-2.5e4 / E);
}
// Calculate We by Kang formula

double We_Kang(double E) {
    return -6060 * pow(E, 0.75);
}



// Niu_a

double Niu_a_Kang(double E) {
    return Eta_Kang(E) * fabs(We_Kang(E));
}
// Niu_i

double Niu_i_Kang(double E) {
    return Alpha_Kang(E) * fabs(We_Kang(E));
}

//Calculate Niu_i and Niu_a together

void Niu_Kang(double *pNiu_i, double *pNiu_a, double E) {
    double We;
    We = fabs(We_Kang(E));
    *pNiu_a = Eta_Kang(E) * We;
    *pNiu_i = Alpha_Kang(E) * We;
}

void Niu_Nikonov(double *pNiu_i, double *pNiu_a, double E, double P) {
    double we;
    we = fabs(We_Nikonov(E, P));
    *pNiu_a = Eta_Nikonov(E, P) * we;
    *pNiu_i = Alpha_Nikonov(E, P) * we;
}

void Niu_MorrowAndLowke(double *pNiu_i, double *pNiu_a, double E, double N) {
    double we;
    we = fabs(We_MorrowAndLowke(E, N));
    *pNiu_a = Eta_MorrowAndLowke(E, N) * we;
    *pNiu_i = Alpha_MorrowAndLowke(E, N) * we;
}
