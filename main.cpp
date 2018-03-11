#include <iostream>
#include "Particle.h"
#include <vector>
#include <cmath>
using namespace std;
# define N 100
# define TAU 0.1
#define T0 0
#define COUNT_ITER 100
const double G = 6.67408e-11;
vector<Particle *> getInitRandom(){
    vector<Particle *> particles(N);  // освобождать память, выделенную под объекты?

    double maxValueCoord = 100;
    double maxValueSpeed = 10;
    double maxValueMass = 5;

    for(int i=0;i<N;i++){
        auto mass = (double)rand()/maxValueMass;

        auto sx = (double)rand()/maxValueCoord;
        auto sy = (double)rand()/maxValueCoord;

        auto vx = (double)rand()/maxValueSpeed;
        auto vy = (double)rand()/maxValueSpeed;

        particles[i] = new Particle(new Point(sx,sy),new Point(vx,vy),mass);

    }
    return particles;
}

void calcForces(vector<Particle *> particles) {
    for (int i = 0; i < N; i++) {
        Point * f = new Point(0.0,0.0);
        auto pi = particles[i];
        for (int j = 0; j < N; j++) {
            if (i == j) continue;

            auto pj = particles[j];

            double dx = pj->s->x - pi->s->x;
            double dy = pj->s->y - pi->s->y;
            double r = dx*dx + dy*dy;
            auto Gmij =G* pi->mass * pj->mass;
            f->x += Gmij*dx /(r*sqrt(r));
            f->y += Gmij*dy /(r*sqrt(r));
        }
        particles[i]->f=f;
        delete f;
    }

}

void move(vector<Particle *> particles, double t){
    for(int i=0;i<N;i++){
        particles[i]->calcAndSetSpeedAndCoord(t);
    }
}

int main() {
    srand (time(NULL));
    auto particles = getInitRandom();

    double t = T0;
    double tau=TAU;
    for (int i = 0; i < COUNT_ITER; i++) {
        calcForces(particles);
        move(particles,t);
        t+=tau;
    }


    return 0;
}