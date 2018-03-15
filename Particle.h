
// Created by andrey on 11.03.18.
//

#ifndef GRAVITATIONALPROBLEM_PARTICLE_H
#define GRAVITATIONALPROBLEM_PARTICLE_H

#include<vector>
#include <iostream>
using namespace std;
struct Point{
public:
    Point(double x, double y):x(x),y(y){}
    double x;
    double y;
};

class Particle {
public:
    double mass;
    Point* s;// coordinate
    Point* v;// speed
    Point *f = new Point(0,0); //сила с которой действуют на эту частицу
    Particle(Point* s,Point* v, double mass):s(s),v(v), mass(mass){}

    void calcAndSetSpeedAndCoord(double tau,double t){
        vector<double> variables = {s->x,v->x};
        auto nVariables = calculateNextValue(NewtonEq,tau,variables,t, mass , f->x);
        s->x = nVariables.at(0);
        v->x = nVariables.at(1);


        variables = {s->y,v->y};
        nVariables = calculateNextValue(NewtonEq,tau,variables,t, mass , f->y);
        s->y = nVariables.at(0);
        v->y = nVariables.at(1);

//        s->x+=v->x*t+f->x*t*t/(2*mass);
//        s->y+=v->y*t+f->y*t*t/(2*mass);
//
//        v->x += f->x * t/mass;
//        v->y += f->y * t/mass;
    }


    vector<double> calculateNextValue(vector<double>(*F)(vector<double>,double,double, double),double tau,vector<double> variables,
                                      double t, double mass, double f){
        vector<double> nVariables(variables.size());
        auto k1=F(variables,t,f,mass);

        vector<double> variablesK2(variables.size());
        for(int k =0;k<variables.size();k++){
            variablesK2[k]=variables[k]+tau/2 * k1[k];
        }
        auto k2=F(variablesK2,t+tau/2,f ,mass);

//        //если PK2
//        if(P==2){
//            for (int j = 0; j < variables.size(); j++){
//                nVariables[j] =variables[j]+tau*k2[j]; // сейчас кол-во элементов равно i
//            }
//            return nVariables;
//        }
//        // конец если PK2

        vector<double> variablesK3(variables.size());
        for(int k =0;k<variables.size();k++){
            variablesK3[k]=variables[k]+tau/2 * k2[k];
        }
        auto k3=F(variablesK3,t+tau/2,f,mass);

        vector<double> variablesK4(variables.size());
        for(int k =0;k<variables.size();k++){
            variablesK4[k]=variables[k]+tau * k3[k];
        }
        auto k4=F(variablesK4,t+tau, f,mass);

        vector<double> K(k4.size());
        for (int l = 0; l < K.size(); ++l) {
            K[l]=1.0/6.0 * (k1[l]+2*k2[l]+2*k3[l]+k4[l]);
        }

        for (int j = 0; j < variables.size(); j++){
            nVariables[j] =variables[j]+tau*K[j]; // сейчас кол-во элементов равно i
        }
        return nVariables;

    }

private:
    static vector<double> NewtonEq(vector<double> variables, double t,double f, double mass) {
        return vector<double>({
                                       variables.at(1),
                                       f/mass
                               });
    }

};


#endif //GRAVITATIONALPROBLEM_PARTICLE_H
