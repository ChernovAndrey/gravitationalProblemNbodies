#define N 2
# define TAU 0.0001
#define T0 0
//#define TF 0.034492738
#define TF 12.58984937
#include "structs.h"
#include<iostream>
#include <cmath>
#include <fstream>
using namespace std;



//variables - это массив в котором хранится только s и v/
 __device__ double*  NewtonEq(double *variables, double t,double f, double mass) {
        double *res = (double *)malloc(2 * sizeof(double));
        res[0] = variables[1];
        res[1]=f/mass;
        return res;
    }

__device__ double* calculateNextValue(double tau, double *variables,  double t, double mass, double f){
        // vector<double> nVariables(variables.size());
        int size = 2;
        double *nVariables = (double *)malloc(size * sizeof(double));
        double *variablesK2 = (double *)malloc(size * sizeof(double));
        double *variablesK3 = (double *)malloc(2 * sizeof(double));
        double *variablesK4 = (double *)malloc(2 * sizeof(double));
        double *K = (double *)malloc(2 * sizeof(double));
                
        double* k1=NewtonEq(variables,t,f,mass);

        for(int k =0;k<size;k++){
            variablesK2[k]=variables[k]+tau/2 * k1[k];
        }
        double* k2=NewtonEq(variablesK2,t+tau/2,f ,mass);

        for(int k =0;k<size;k++){
            variablesK3[k]=variables[k]+tau/2 * k2[k];
        }
        double* k3 = k3=NewtonEq(variablesK3,t+tau/2,f,mass);

        for(int k =0;k<size;k++){
            variablesK4[k]=variables[k]+tau * k3[k];
        }
        double* k4=NewtonEq(variablesK4,t+tau, f,mass);

        for (int l = 0; l < size; ++l) {
            K[l]=1.0/6.0 * (k1[l]+2*k2[l]+2*k3[l]+k4[l]);
        }

        for (int j = 0; j < size; j++){
            nVariables[j] =variables[j]+tau*K[j]; // сейчас кол-во элементов равно i
        }
        free(variablesK2);free(variablesK3);free(variablesK4);
        free(K);
        free(k1);free(k2);free(k3);free(k4);
        return nVariables;

    }

__device__ void calcAndSetSpeedAndCoord(Particle particle,double tau,double t){
        double * variables = (double*)malloc(2 * sizeof(double));
        variables[0]=particle.s->x;
        variables[1]=particle.v->x;
        double* nVariables = calculateNextValue(tau,variables,t, particle.mass , particle.f->x);
        particle.s->x = nVariables[0];
        particle.v->x = nVariables[1];


        variables[0]=particle.s->y;
        variables[1]=particle.v->y;
        nVariables = calculateNextValue(tau,variables,t, particle.mass , particle.f->y);
        particle.s->y = nVariables[0];
        particle.v->y = nVariables[1];
    }

__global__ void calculate(Particle *particles, double* tau, double* t) {
        const double G = 498217402368e-12;
        int i = blockIdx.x;
        Point *f = (struct Point*)malloc (sizeof (struct Point));
        Particle pi = particles[i];
        for (int j = 0; j < N; j++) {
            if (i == j) continue;

            Particle pj = particles[j];
            double dx = pj.s->x - pi.s->x;
            double dy = pj.s->y - pi.s->y;
            double r = sqrt(dx*dx + dy*dy);
            double Gmij =G * pi.mass * pj.mass;
            f->x += Gmij*dx /(r*r*r);
            f->y += Gmij*dy /(r*r*r);
        }
        particles[i].f=f;
        particles[i].f->x +=100;
        calcAndSetSpeedAndCoord(particles[i],*tau,*t);
        //free(f);
}

// void writingDataInFile(double x, double y , int k){ // k - номер тела
//     const string pathToFile = "/home/andrey/CLionProjects/gravitationalProblem/cuda/files/"+to_string(k)+".txt";
//     ofstream file(pathToFile, ios_base::app);
//     file<<x<<" "<<y<<endl;
//     file.close();
// }

__global__ void add(int *a, int *b, int *c){
    printf("Hello, world from the device!\n");
    *c= *a + *b;
}

// int main(){
//     // cleanFiles(N);
//     srand (time(NULL));
//     //auto particles = getInitRandom();
//     Particle* particles = getInitFor2Particle();
//     double t = T0;
//     double tau=TAU;
//     double tf=TF;
//     double *dev_tau;
//     double *dev_T;
//     Particle* dev_particles;
//     int countIter = (tf-t)/tau;
    
//     int size= N * sizeof(Particle);
//     cudaMalloc((void**)&dev_particles, size);
//     cudaMalloc((void**)&dev_tau, sizeof(double));
//     cudaMalloc((void**)&dev_T, sizeof(double));
    
//     cudaMemcpy(dev_particles, particles, size, cudaMemcpyHostToDevice);
//     cudaMemcpy(dev_tau, &tau, sizeof(double), cudaMemcpyHostToDevice);
        
//     cout <<particles[0].f->x << '\t'<<particles[0].s->y <<endl;
//     for (int i = 0; i < countIter; i++) {
//         cudaMemcpy(dev_T, &t, sizeof(double), cudaMemcpyHostToDevice);
//         cout<<"T="<<*dev_T<<'\t'<<t<<endl;
//      //   *dev_T+=0.01;
//         calculate<<<N,1>>>(dev_particles,dev_tau,dev_T);
//         t+=tau;
//         cudaMemcpy(particles, dev_particles, size, cudaMemcpyDeviceToHost);
//         cout <<particles[0].f->x << '\t'<<particles[0].s->y <<endl;
//         // writingDataInFile(particles[0].s->x,particles[0].s->y,0);
//         // writingDataInFile(particles[1].s->x,particles[1].s->y,1);     
//     }

//     free(particles);
//     cudaFree(dev_particles);cudaFree(dev_T);cudaFree(dev_tau);
//     return 0;
// }

int main(void){
    int a, b, c;
    int *dev_a, *dev_b, *dev_c;
    int size= sizeof(int);
    cudaMalloc((void**)&dev_a, size);
    cudaMalloc((void**)&dev_b, size);
    cudaMalloc((void**)&dev_c, size);
    a=3;
    b=1;
    cudaMemcpy(dev_a,&a,size, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b,&b,size, cudaMemcpyHostToDevice);

    add<<< 1,1 >>>(dev_a,dev_b,dev_c);
    cudaDeviceSynchronize();
    cudaMemcpy(&c,dev_c,size, cudaMemcpyDeviceToHost);

    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);
    cout<<c<<endl;
    return 0;
}

