
#ifndef INONIZATION_FORMULA_H
#define INONIZATION_FORMULA_H

#pragma once
#include "common.h"

//////////////////////////////////////////////////////////////////////////
// Nikonov formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Nikonov formula
MyDataF Alpha_Nikonov(MyDataF E, MyDataF P);
// Calculate Eta by Nikonov formula
MyDataF Eta_Nikonov(MyDataF E, MyDataF P);
// Calculate We by Nikonov formula
MyDataF We_Nikonov(MyDataF E, MyDataF P);
// Niu_a
MyDataF Niu_a_Nikonov(MyDataF E, MyDataF P);
// Niu_i
MyDataF Niu_i_Nikonov(MyDataF E, MyDataF P);

//////////////////////////////////////////////////////////////////////////
//Morrow and Lowke formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Morrow and Lowke formula
MyDataF Alpha_MorrowAndLowke(MyDataF E, MyDataF N);
// Calculate Eta by Morrow and Lowke formula
MyDataF Eta_MorrowAndLowke(MyDataF E, MyDataF N);
// Calculate We by Morrow and Lowke formula
MyDataF We_MorrowAndLowke(MyDataF E, MyDataF N);
// Niu_a
MyDataF Niu_a_MorrowAndLowke(MyDataF E, MyDataF N);
// Niu_i
MyDataF Niu_i_MorrowAndLowke(MyDataF E, MyDataF N);

//////////////////////////////////////////////////////////////////////////
//Kang formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Kang formula
MyDataF Alpha_Kang(MyDataF E);
// Calculate Eta by Kang formula
MyDataF Eta_Kang(MyDataF E);
// Calculate We by Kang formula
MyDataF We_Kang(MyDataF E);
// Niu_a
MyDataF Niu_a_Kang(MyDataF E);
// Niu_i
MyDataF Niu_i_Kang(MyDataF E);

//Calculate Niu_i and Niu_a together
void Niu_Kang(MyDataF *pNiu_i, MyDataF *pNiu_a, MyDataF E);
void Niu_Nikonov(MyDataF *pNiu_i, MyDataF *pNiu_a, MyDataF E, MyDataF P);
void Niu_MorrowAndLowke(MyDataF *pNiu_i, MyDataF *pNiu_a, MyDataF E, MyDataF N);

#endif
