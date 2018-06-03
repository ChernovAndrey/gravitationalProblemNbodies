
#include<iostream>
using namespace std;
#include <fstream>
#include"structs.h"
#include<string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
__global__ void add( int *a, int *b, int *c ) {
printf("%d\n", *a);
*c = *a + *b;
}

 __device__ double*  NewtonEq(double *variables, double t,double f, double mass) {
        double *res = (double *)malloc(2 * sizeof(double));
        // double res[2];
        res[0] = variables[1];
        res[1]=f/mass;
        return res;
    }

__device__ void calculateNextValue(double tau, double *variables,  double t, double mass, double f){
        // vector<double> nVariables(variables.size());
        const int size = 2;
        // double *variablesK2 = (double *)malloc(size * sizeof(double));
        // double *variablesK3 = (double *)malloc(2 * sizeof(double));
        // double *variablesK4 = (double *)malloc(2 * sizeof(double));
        // double *K = (double *)malloc(2 * sizeof(double));        
        double k1[2];
        double k2[2];
        double k3[2];
        double k4[2];	
        // double* k1=NewtonEq(variables,t,f,mass);
        k1[0]=variables[1];
        k1[1]=f/mass;

        // double nVariables[2];
        double variablesK2[2];
        double variablesK3[2];
        double variablesK4[2];
        double K[2];
        for(int k =0;k<size;k++){
            variablesK2[k]=variables[k]+tau/2 * k1[k];
        }
        // double* k2=NewtonEq(variablesK2,t+tau/2,f ,mass);
        k2[0]=variablesK2[1];
        k2[1]=k1[1];

        for(int k =0;k<size;k++){
            variablesK3[k]=variables[k]+tau/2 * k2[k];
        }
        // double* k3 = k3=NewtonEq(variablesK3,t+tau/2,f,mass);
        k3[0]=variablesK3[1];
        k3[1]=k1[1];
        for(int k =0;k<size;k++){
            variablesK4[k]=variables[k]+tau * k3[k];
        }
        // double* k4=NewtonEq(variablesK4,t+tau, f,mass);
        k4[0]=variablesK4[1];
        k4[1]=k1[1];
        for (int l = 0; l < size; ++l) {
            K[l]=1.0/6.0 * (k1[l]+2*k2[l]+2*k3[l]+k4[l]);
        }

        for (int j = 0; j < size; j++){
            variables[j] =variables[j]+tau*K[j]; // сейчас кол-во элементов равно i
        }
    //новые значения в variables
    }

__device__ void calcAndSetSpeedAndCoord(Particle &particle,double tau,double t){
        //double * variables = (double*)malloc(2 * sizeof(double));
		double variables[2];        
        variables[0]=particle.s.x;
        variables[1]=particle.v.x;
        calculateNextValue(tau,variables,t, particle.mass , particle.f.x);
        particle.s.x = variables[0];
        particle.v.x = variables[1];


        variables[0]=particle.s.y;
        variables[1]=particle.v.y;
        calculateNextValue(tau,variables,t, particle.mass , particle.f.y);
        particle.s.y = variables[0];
        particle.v.y = variables[1];
        // delete nVariables;
    }

// string convertIntToString(int a){
// 	// char *intStr = itoa(a);
// 	// return  string(intStr);

// 	stringstream ss;
//  	ss << a;
// 	return ss.str();
// }

// void writingDataInFile(double x, double y , int k){ // k - номер тела
//     string pathToFile = "files/"+convertIntToString(k)+".txt";
//     char *cpath = new char[pathToFile.length() + 1];
// 	strcpy(cpath, pathToFile.c_str());
// 	// do stuff
//     FILE* f = fopen(cpath, "a");
// 	fprintf(f, "%f %f", x,y);
// 	fclose(f);
// 	delete [] cpath;
//     //ofstream file(pathToFile, ios_base::app);
//     // file<<x<<" "<<y<<endl;
//     // file.close();
// }



__device__ void  calcForces(Particle *particles, int i, double G, int *N){
	Point f;
	Particle pi = particles[i];
	for (int j = 0; j < *N; j++) {
	    // printf("%s\n", "start i iteration");
	    if (i == j) continue;

	    Particle pj = particles[j];
	    double dx = pj.s.x - pi.s.x;
	    double dy = pj.s.y - pi.s.y;
	    double r = sqrt(dx*dx + dy*dy);
	    double Gmij =G * pi.mass * pj.mass;
	    f.x += Gmij*dx /(r*r*r);
	    f.y += Gmij*dy /(r*r*r);
	}
	particles[i].f=f;

}

__global__ void calculate(Particle *particles,double *t, double *tau, int *N){
	int i =  threadIdx.x + blockIdx.x * blockDim.x;;
	// printf("blockdim=%d\n",blockDim.x);
	const double G = 498217402368e-12;
	// printf("N=%d\n", *N);
	// printf("%s\n","start calculate" );
    // printf("in gpu start %f  %f\n",particles[i].s.x,particles[i].s.y);
	// Point f;
	// Particle pi = particles[i];
	// for (int j = 0; j < *N; j++) {
	//     // printf("%s\n", "start i iteration");
	//     if (i == j) continue;

	//     Particle pj = particles[j];
	//     double dx = pj.s.x - pi.s.x;
	//     double dy = pj.s.y - pi.s.y;
	//     double r = sqrt(dx*dx + dy*dy);
	//     double Gmij =G * pi.mass * pj.mass;
	//     f.x += Gmij*dx /(r*r*r);
	//     f.y += Gmij*dy /(r*r*r);
	// }
	// particles[i].f=f;
	 calcForces(particles,i,G,N);
	__syncthreads();
	calcAndSetSpeedAndCoord(particles[i],*tau,*t);
	// printf("f=%f %f\n",f.x,f.y);
	__syncthreads();
}

__global__ void printParticles(Particle *particles){
	printf("sx = %f\n", particles[0].s.x);
}

void copyParticlesDeviceToHost(Particle *&particles, Particle *&dev_particles, int size){
	// cudaMemcpy(&particles, &dev_particles, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(&particles[0].s, &dev_particles[0].s, sizeof(struct Point), cudaMemcpyDeviceToHost);
    // cudaMemcpy(&particles[0].v, &dev_particles[0].v, sizeof(struct Point), cudaMemcpyDeviceToHost);
    // cudaMemcpy(&particles[0].f, &dev_particles[0].f, sizeof(struct Point), cudaMemcpyDeviceToHost);
    cudaMemcpy(&particles[1].s, &dev_particles[1].s, sizeof(struct Point), cudaMemcpyDeviceToHost);
    // cudaMemcpy(&particles[1].v, &dev_particles[1].v, sizeof(struct Point), cudaMemcpyDeviceToHost);
    // cudaMemcpy(&particles[1].f, &dev_particles[1].f, sizeof(struct Point), cudaMemcpyDeviceToHost);
}


int main( void ) {
// cleanFiles(N);
srand (time(NULL));

Particle *particles=NULL;
Particle *dev_particles=NULL;
int countBlocks = 13;
int countThreads = 64;
int N=countBlocks*countThreads;
getInitRandom(particles,dev_particles,N);
// getInitFor2Particle(particles,dev_particles);
double t = 0.0;
// double tau=0.0001;
double tau =1.0;
//double tf = 0.001;
// double tf=12.58984937;
double tf=10.0;
double *dev_tau;
double *dev_T;

int *dev_N;
int countIter = (tf-t)/tau;
int size= N*sizeof (struct Particle);
//int size= N * sizeof(Particle);
cudaMalloc((void**)&dev_tau, sizeof(double));
cudaMalloc((void**)&dev_T, sizeof(double));
cudaMalloc((void**)&dev_N,sizeof(double));
   
cudaMemcpy(dev_N,&N,sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(dev_tau, &tau, sizeof(double), cudaMemcpyHostToDevice);
// cudaMemcpy(dev_T, &t, sizeof(double), cudaMemcpyHostToDevice);

cudaEvent_t start, stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);

cudaEventRecord(start);

cout<<"count iter="<<countIter<<endl;
for (int i = 0; i < countIter; i++) {
  //   cout<<"start iter in host"<<endl;
     cudaMemcpy(dev_T, &t, sizeof(double), cudaMemcpyHostToDevice);
     // cout<<"T="<<*dev_T<<'\t'<<t<<endl;
  //   *dev_T+=0.01;
     calculate<<<countBlocks,countThreads>>>(dev_particles,dev_T,dev_tau,dev_N);
     t+=tau;
     copyParticlesDeviceToHost(particles,dev_particles,size);
     // cout <<particles[0].s.x << " "<<particles[0].s.y <<endl;
     // cout <<particles[1].s.x << " "<<particles[1].s.y <<endl;
     // cout<<"finish iter="<<i<<endl;
      // writingDataInFile(particles[0].s.x,particles[0].s.y,0);
      // writingDataInFile(particles[1].s.x,particles[1].s.y,1);
      // if (i%1000==0){
      // 	cout<<"finish iter="<<i<<endl;
      // }     
 }
cudaEventRecord(stop);
cudaEventSynchronize(stop);
float milliseconds = 0;
cudaEventElapsedTime(&milliseconds, start, stop);    
cout<<"time execute:"<<milliseconds<<endl;
cudaEventDestroy (start);
cudaEventDestroy (stop);
cudaFree( dev_particles );
cudaFree( dev_tau );
cudaFree( dev_T );
free(particles); // все надо освободить
return 0;
}