//
// Created by andrey on 28.02.18.
//

#ifndef LAB6_SOLVEMETHOD_H
#define LAB6_SOLVEMETHOD_H

#define N 500
#define T0 0
#define TAU 0.055371
#define EPS 1e-4

#include <iostream>
#include <vector>
using namespace std;

class SolveMethod { //interface
public:
    virtual vector<vector<double>> solve(vector<double>(*F)(vector<double>, double),const vector<double> &initVariables, int n)=0;//n - кол-во шагов


};


#endif //LAB6_SOLVEMETHOD_H
