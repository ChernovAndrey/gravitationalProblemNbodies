#include <iostream>
#include "Particle.h"
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;
# define N 2
# define TAU 0.01
#define T0 0
#define COUNT_ITER 1000
const double G = 6.67408e-11;


void writingDataInFile(double x, double y , int k){ // k - номер тела
    const string pathToFile = "/home/andrey/CLionProjects/gravitationalProblem/files/"+to_string(k)+".txt";
    ofstream file(pathToFile, ios_base::app);
    file<<x<<" "<<y<<endl;
    file.close();
}

vector<Particle *> getInitFor2Particle(){
    double x1 = 1;
    double y1 = 0;
    double x2 = -1;
    double y2 = 0;
    vector<Particle *> particles(N);  // освобождать память, выделенную под объекты?
    particles[0] = new Particle(new Point(x1,y1),new Point(0,5.775e-6),2);
    particles[1] = new Particle(new Point(x2,y2),new Point(0,-5.775e-6),2);

    writingDataInFile(x1,y1,0);
    writingDataInFile(x2,y2,1);

    return particles;
}


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
            double r = sqrt(dx*dx + dy*dy);
            auto Gmij =G* pi->mass * pj->mass;
            f->x += Gmij*dx /(r*r*r);
            f->y += Gmij*dy /(r*r*r);
        }
        particles[i]->f=f;
        delete f;
    }

}

void move(vector<Particle *> particles, double t){
    for(int i=0;i<N;i++){
        particles[i]->calcAndSetSpeedAndCoord(TAU,t);
        if(N==2){
            writingDataInFile(particles[i]->s->x, particles[i]->s->y, i);
        }
    }
}

void cleanFiles(){
    for(int i=0;i<N;i++){
        const string pathToFile = "/home/andrey/CLionProjects/gravitationalProblem/files/"+to_string(i)+".txt";
        ofstream file(pathToFile, ios_base::out | ios_base::trunc);
        file.close();
    }
}
int main() {
    cleanFiles();
    srand (time(NULL));
    //auto particles = getInitRandom();
    auto particles = getInitFor2Particle();
    double t = T0;
    double tau=TAU;
    for (int i = 0; i < COUNT_ITER; i++) {
        calcForces(particles);
        move(particles,t);
        t+=tau;
    }


    return 0;
}