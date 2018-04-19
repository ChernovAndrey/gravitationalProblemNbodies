#include <iostream>
#include "Particle.h"
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;
# define N 3
# define TAU 0.0001
#define T0 0
//#define TF 0.034492738
#define TF 12.58984937
//const double G = 1.18555535802194e-04;
const double G = 498217402368e-12;

// время в днях

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
    vector<Particle *> particles(2);
    particles[0] = new Particle(new Point(x1,y1),new Point(0,498960e-6),2);
    particles[1] = new Particle(new Point(x2,y2),new Point(0,-498960e-6),2);

    writingDataInFile(x1,y1,0);
    writingDataInFile(x2,y2,1);

    return particles;
}


vector<Particle *> getInitFor3Particle(){
    double x1 = 1;
    double y1 = 0;
    double x2 = -1;
    double y2 = 0;
    double x3 = 0;
    double y3 = 0;
    vector<Particle *> particles(3);
    particles[0] = new Particle(new Point(x1,y1),new Point(0,498960e-6),2);
    particles[1] = new Particle(new Point(x2,y2),new Point(0,-498960e-6),2);
    particles[2] = new Particle(new Point(x3,y3),new Point(0,0),2);

    writingDataInFile(x1,y1,0);
    writingDataInFile(x2,y2,1);
    //writingDataInFile(0,0,2);

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

void calcForces(vector<Particle *> &particles) {
    for (int i = 0; i < N; i++) {
        Point * f = new Point(0.0,0.0);
        auto pi = particles[i];
        for (int j = 0; j < N; j++) {
            if (i == j) continue;

            auto pj = particles[j];

            double dx = pj->s->x - pi->s->x;
            double dy = pj->s->y - pi->s->y;
            double r = sqrt(dx*dx + dy*dy);
            double Gmij =G * pi->mass * pj->mass;
            f->x += Gmij*dx /(r*r*r);
            f->y += Gmij*dy /(r*r*r);
   //         cout<< f->x <<" " << f->y<<endl;
        }
        particles[i]->f=f;
        cout<<"f:"<<endl;
        cout<<particles[i]->f->x<<'\t'<<particles[i]->f->y<<endl;
    }

}

void moving(vector<Particle *> &particles, double t){
    for(int i=0;i<particles.size();i++){
        //cout<<particles[i]->f->y<<endl;
        particles[i]->calcAndSetSpeedAndCoord(TAU,t);
        writingDataInFile(particles[i]->s->x, particles[i]->s->y, i);
    }
}

void cleanFiles(int n ){
    for(int i=0;i<n;i++){
        const string pathToFile = "/home/andrey/CLionProjects/gravitationalProblem/files/"+to_string(i)+".txt";
        ofstream file(pathToFile, ios_base::out | ios_base::trunc);
        file.close();
    }
}

void deleteVectorObjects(vector<Particle *> v ){
    for (int i = 0; i < v.size(); ++i) {
        delete v[i];
    }
}

int main() {
    cleanFiles(N);
    srand (time(NULL));
    //auto particles = getInitRandom();
    auto particles = getInitFor3Particle();
    double t = T0;
    double tau=TAU;
    double tf=TF;
    int countIter = (tf-t)/tau;
    for (int i = 0; i < countIter; i++) {
        calcForces(particles);
     //   cout<<particles[i]->f->y<<endl;
        moving(particles,t);
        t+=tau;
    }
    deleteVectorObjects(particles);
    return 0;
}