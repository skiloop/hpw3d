
#include <math.h>
#include "InonizationFormula.h"
//////////////////////////////////////////////////////////////////////////
//int sign(MyDataF val)
//if val>0 return 1
//if val<0 return -1
//return 0 if val==0

int sign(MyDataF val) {
    return (val < 0 ? -1 : (val > 0 ? 1 : 0));
}

//////////////////////////////////////////////////////////////////////////
// Nikonov formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Nikonov formula

MyDataF Alpha_Nikonov(MyDataF E, MyDataF P) {
    MyDataF EDivP = fabs(E) / P;
    if (EDivP < 108.0) {
        return 3.9 * P * exp(-213 / EDivP);
    } else {
        return 14.5 * P * exp(-316 / EDivP);
    }
}
// Calculate Eta by Nikonov formula

MyDataF Eta_Nikonov(MyDataF E, MyDataF P) {
    MyDataF EDivP;
    EDivP = fabs(E) / P;
    if (EDivP < 50.0) {
        MyDataF val1 = 4.47e-3 * (EDivP)*(EDivP);
        if (EDivP >= 10.0) {
            return val1;
        } else {
            MyDataF val2 = 4.47 / EDivP;
            return (val1 > val2 ? val1 : val2);
        }
    } else {
        MyDataF temp = sqrt(EDivP);
        return (EDivP <= 90 ? 1.58 * temp : 142 / temp);
    }

}
// Calculate We by Nikonov formula

MyDataF We_Nikonov(MyDataF E, MyDataF P) {
    return -0.0382 * E - 2.9e5 * E / P;
}
// Niu_a

MyDataF Niu_a_Nikonov(MyDataF E, MyDataF P) {
    return Eta_Nikonov(E, P) * fabs(We_Nikonov(E, P));
}
// Niu_i

MyDataF Niu_i_Nikonov(MyDataF E, MyDataF P) {
    return Alpha_Nikonov(E, P) * fabs(We_Nikonov(E, P));
}

//////////////////////////////////////////////////////////////////////////
//Morrow and Lowke formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Morrow and Lowke formula

MyDataF Alpha_MorrowAndLowke(MyDataF E, MyDataF N) {
    MyDataF edn = fabs(E) / N / 1e-15;

    if (edn > 1.5) {
        return 2e-16 * N * exp(-7.248 / edn);
    } else {
        return 6.619e-17 * N * exp(-5.593 / edn);
    }
}
// Calculate Eta by Morrow and Lowke formula

MyDataF Eta_MorrowAndLowke(MyDataF E, MyDataF N) {
    MyDataF edn;
    edn = fabs(E) / N / 1e-16;

    if (edn > 0.12) {
        if (edn < 10.50)
            return N * ((6.089e-20 * edn - 2.893e-19) + N * 4.47778e-59 * pow(edn * 1e-16, -1.2749));
        else
            return N * ((8.889e-21 * edn + 2.567e-19) + N * 4.47778e-59 * pow(edn * 1e-16, -1.2749));
    } else {
        if (edn < 0)
            return 0;
        else
            return 106.81;
    }
}
// Calculate We by Morrow and Lowke formula

MyDataF We_MorrowAndLowke(MyDataF E, MyDataF N) {
    MyDataF edn = fabs(E) / N / 1e-16;

    if (edn > 1.0) {
        if (edn <= 20.00)
            return -sign(E)*(1.03 * edn + 1.3)*1e-10;
        else
            return -sign(E)*(7.4e5 * edn + 7.1e6)*1e-16;
    } else {
        if (edn <= 0.26)
            return -sign(E)*(6.87e6 * edn + 3.38e4)*1e-16;
        else
            return -sign(E)*(7.2973 * edn + 16.3)*1e-11;
        ;
    }
}
// Niu_a

MyDataF Niu_a_MorrowAndLowke(MyDataF E, MyDataF N) {
    return Eta_MorrowAndLowke(E, N) * fabs(We_MorrowAndLowke(E, N));
}
// Niu_i

MyDataF Niu_i_MorrowAndLowke(MyDataF E, MyDataF N) {
    return Alpha_MorrowAndLowke(E, N) * fabs(We_MorrowAndLowke(E, N));
}

//////////////////////////////////////////////////////////////////////////
//Kang formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Kang formula

MyDataF Alpha_Kang(MyDataF E) {
    return 0.0035 * exp(-1.65e5 / E);
}
// Calculate Eta by Kang formula

MyDataF Eta_Kang(MyDataF E) {
    return 0.15 * exp(-2.5e4 / E);
}
// Calculate We by Kang formula

MyDataF We_Kang(MyDataF E) {
    return -6060 * pow(E, 0.75);
}



// Niu_a

MyDataF Niu_a_Kang(MyDataF E) {
    return Eta_Kang(E) * fabs(We_Kang(E));
}
// Niu_i

MyDataF Niu_i_Kang(MyDataF E) {
    return Alpha_Kang(E) * fabs(We_Kang(E));
}

//Calculate Niu_i and Niu_a together

void Niu_Kang(MyDataF *pNiu_i, MyDataF *pNiu_a, MyDataF E) {
    MyDataF We;
    We = fabs(We_Kang(E));
    *pNiu_a = Eta_Kang(E) * We;
    *pNiu_i = Alpha_Kang(E) * We;
}

void Niu_Nikonov(MyDataF *pNiu_i, MyDataF *pNiu_a, MyDataF E, MyDataF P) {
    MyDataF we;
    we = fabs(We_Nikonov(E, P));
    *pNiu_a = Eta_Nikonov(E, P) * we;
    *pNiu_i = Alpha_Nikonov(E, P) * we;
}

void Niu_MorrowAndLowke(MyDataF *pNiu_i, MyDataF *pNiu_a, MyDataF E, MyDataF N) {
    MyDataF we;
    we = fabs(We_MorrowAndLowke(E, N));
    *pNiu_a = Eta_MorrowAndLowke(E, N) * we;
    *pNiu_i = Alpha_MorrowAndLowke(E, N) * we;
}
