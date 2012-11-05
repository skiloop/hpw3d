
#ifndef INONIZATION_FORMULA_H
#define INONIZATION_FORMULA_H

#pragma once
//////////////////////////////////////////////////////////////////////////
// Nikonov formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Nikonov formula
double Alpha_Nikonov(double E, double P);
// Calculate Eta by Nikonov formula
double Eta_Nikonov(double E, double P);
// Calculate We by Nikonov formula
double We_Nikonov(double E, double P);
// Niu_a
double Niu_a_Nikonov(double E, double P);
// Niu_i
double Niu_i_Nikonov(double E, double P);

//////////////////////////////////////////////////////////////////////////
//Morrow and Lowke formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Morrow and Lowke formula
double Alpha_MorrowAndLowke(double E, double N);
// Calculate Eta by Morrow and Lowke formula
double Eta_MorrowAndLowke(double E, double N);
// Calculate We by Morrow and Lowke formula
double We_MorrowAndLowke(double E, double N);
// Niu_a
double Niu_a_MorrowAndLowke(double E, double N);
// Niu_i
double Niu_i_MorrowAndLowke(double E, double N);

//////////////////////////////////////////////////////////////////////////
//Kang formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Kang formula
double Alpha_Kang(double E);
// Calculate Eta by Kang formula
double Eta_Kang(double E);
// Calculate We by Kang formula
double We_Kang(double E);
// Niu_a
double Niu_a_Kang(double E);
// Niu_i
double Niu_i_Kang(double E);

//Calculate Niu_i and Niu_a together
void Niu_Kang(double *pNiu_i, double *pNiu_a, double E);
void Niu_Nikonov(double *pNiu_i, double *pNiu_a, double E, double P);
void Niu_MorrowAndLowke(double *pNiu_i, double *pNiu_a, double E, double N);

#endif
