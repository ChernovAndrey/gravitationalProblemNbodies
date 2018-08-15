
#include<iostream>
using namespace std;
#include <fstream>
#include"structs.h"
#include<string>


__device__ Particle  **vecParticlesK;
__device__ int dM; // кол-во необходимых массивов.

__device__ void  calcForces(Particle *particles, int i, int N,double G){ //i-номер тела
	Point f={0,0};
	Particle pi = particles[i];
	// printf("i = %d ,  N= %d\n",i,N);
	for (int j = 0; j < N; j++) {
	    // printf("%s\n", "start i iteration");
	    if (i == j) continue;

	    Particle pj = particles[j];
	    double dx = pj.s.x - pi.s.x;
	    double dy = pj.s.y - pi.s.y;
	    // printf("dx = %f ,  dy=%f\n",dx,dy);    
	    double r = sqrt(dx*dx + dy*dy);
	    double Gmij =G * pi.mass * pj.mass;
	    f.x += Gmij*dx /(r*r*r);
	    f.y += Gmij*dy /(r*r*r);
	    // printf("f = %f    %f\n",f.x,f.y);
    
	}
	particles[i].f=f;

}

__device__ void copyParticlesToVecParticlesK(Particle * particles, int i){
	for(int j=0;j<dM;j++){
		vecParticlesK[j][i] = particles[i];
	}
}
__device__ void moveParticle(Particle *&particles, Particle *ParticlesK, double step, int i) { // i - номер тела
            particles[i].s.x += step * ParticlesK[i].v.x;//подставляем непосредственно уравнение Ньютона.
            particles[i].s.y += step * ParticlesK[i].v.y;//подставляем непосредственно уравнение Ньютона.

            particles[i].v.x += step * ParticlesK[i].f.x /
                     ParticlesK[i].mass;//подставляем непосредственно уравнение Ньютона.
            particles[i].v.y += step * ParticlesK[i].f.y /
                    ParticlesK[i].mass;//подставляем непосредственно уравнение Ньютона.
       
}

// __device__ void moveParticle(Particle *&particles, double step[], int i) { // i - номер тела
// 		for (int j = 0; j < vecParticlesK.size(); ++j) {
//             particles[i].s.x += step[j] * vecParticlesK[j][i]->v.x;//подставляем непосредственно уравнение Ньютона.
//             particles[i].s.y += step[j] * vecParticlesK[j][i]->v.y;//подставляем непосредственно уравнение Ньютона.

//             particles[i].v.x += step[j] * vecParticlesK[j][i]->f.x /
//                     vecParticlesK[j][i]->mass;//подставляем непосредственно уравнение Ньютона.
//             particles[i].v.y += step[j] * vecParticlesK[j][i]->f.y /
//                     vecParticlesK[j][i]->mass;//подставляем непосредственно уравнение Ньютона.
//      	   }
// }


__device__ void moveParticle(Particle *&particles, Particle* vecParticlesK[], double step[], int i, int sizeVecParticles) { // i - номер тела
		for (int j = 0; j < sizeVecParticles; ++j) {
            particles[i].s.x += step[j] * vecParticlesK[j][i].v.x;//подставляем непосредственно уравнение Ньютона.
            particles[i].s.y += step[j] * vecParticlesK[j][i].v.y;//подставляем непосредственно уравнение Ньютона.

            particles[i].v.x += step[j] * vecParticlesK[j][i].f.x /
                    vecParticlesK[j][i].mass;//подставляем непосредственно уравнение Ньютона.
            particles[i].v.y += step[j] * vecParticlesK[j][i].f.y /
                    vecParticlesK[j][i].mass;//подставляем непосредственно уравнение Ньютона.
     	   }
}



__device__ void calculateRK1(Particle *&particles, double tau,int i, int N, double G) { // i - номер тела.

    calcForces(particles,i,N,G);
    __syncthreads();
    moveParticle(particles,particles,tau,i);
}


__device__ void calculateRK2(Particle  *&particles, double tau, int i, int N, double G) {
	calcForces(particles,i,N,G);
	 __syncthreads();
	copyParticlesToVecParticlesK(particles,i);
    __syncthreads();
    moveParticle(vecParticlesK[0],vecParticlesK[0],tau/2.0,i);
    __syncthreads();
    calcForces(vecParticlesK[0],i,N,G);
    __syncthreads();
    moveParticle(particles,vecParticlesK[0],tau,i);
}



__device__ void calculateRK3(Particle *&particles, double tau, int i, int N, double G) {

    calcForces(particles,i,N,G);
	__syncthreads();
	copyParticlesToVecParticlesK(particles,i);
    __syncthreads();
    moveParticle(vecParticlesK[0],vecParticlesK[0],tau/2.0,i);
    __syncthreads();
    calcForces(vecParticlesK[0],i,N,G);
	__syncthreads();
	moveParticle(vecParticlesK[1],(Particle*[]){particles,vecParticlesK[0]},(double[]){-tau,2*tau},i,2);
    __syncthreads();
    calcForces(vecParticlesK[1],i,N,G);
    __syncthreads();
    moveParticle(particles,(Particle*[]){particles,vecParticlesK[0],vecParticlesK[1]},(double[]){tau/6.0,4*tau/6.0,tau/6.0},i,3);

}

__device__ void calculateRK4(Particle *&particles, double tau, int i, int N, double G) {

    calcForces(particles,i,N,G);
    __syncthreads();
   	copyParticlesToVecParticlesK(particles,i);
   	__syncthreads();
    moveParticle(vecParticlesK[0],vecParticlesK[0],tau/2.0,i);
	__syncthreads();
    calcForces(vecParticlesK[0],i,N,G);
	__syncthreads();
	moveParticle(vecParticlesK[1],vecParticlesK[0],tau/2.0,i);
    __syncthreads();
    calcForces(vecParticlesK[1],i,N,G);
	__syncthreads();
	moveParticle(vecParticlesK[2],vecParticlesK[1],tau,i);
    __syncthreads();
    calcForces(vecParticlesK[2],i,N,G);
    __syncthreads();
    moveParticle(particles,(Particle*[]){particles,vecParticlesK[0],vecParticlesK[1],vecParticlesK[2]},(double[]){tau/6.0,2*tau/6.0,2*tau/6.0,tau/6.0},i,4);

}

__global__ void calculate(Particle *particles,double *t, double *tau, int *N){
	int i =  threadIdx.x + blockIdx.x * blockDim.x;;
	// int i = threadIdx.x;
	// printf("blockdim=%d\n",blockDim.x);
	const double G = 498217402368e-12;
	// const double G = 4.0;
	calculateRK1(particles,*tau,i,*N,G);
	__syncthreads();
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

__global__ void createVecParticlesK(int *M, int *N){// M - кол-во необходимых массивов.
	 vecParticlesK = new Particle*[*M];
	 for(int i=0;i<*M;i++){
	 	vecParticlesK[i] = new Particle[*N];
	 }
	 dM = *M;
}


__global__ void deleteVecParticlesK(int *M, int *N){// M - кол-во необходимых массивов.
	 for(int i=0;i<*M;i++){
	 	delete[] vecParticlesK[i];
	 }
	 delete[] vecParticlesK;
}

int main( void ) {
// cleanFiles(N);
srand (time(NULL));

Particle *particles=NULL;
Particle *dev_particles=NULL;
int countBlocks = 13;
int countThreads = 128;
int N=countBlocks*countThreads;
// int N=2;
int M=3;
getInitRandom(particles,dev_particles,N);
// getInitFor2Particle(particles,dev_particles);
double t = 0.0;
// double tau=0.0001;
double tau =0.001;
//double tf = 0.001;
// double tf=12.58984937;
double tf = 0.01;
// double tf=2*M_PI;
double *dev_tau;
double *dev_T;

int *dev_N;
int *dev_M;
int countIter = (tf-t)/tau;
int size= N*sizeof (struct Particle);
//int size= N * sizeof(Particle);
cudaMalloc((void**)&dev_tau, sizeof(double));
cudaMalloc((void**)&dev_T, sizeof(double));
cudaMalloc((void**)&dev_N,sizeof(int));
cudaMalloc((void**)&dev_M,sizeof(int));

cudaMemcpy(dev_N,&N,sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(dev_M,&M,sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(dev_tau, &tau, sizeof(double), cudaMemcpyHostToDevice);
createVecParticlesK<<<1,1>>>(dev_M,dev_N);

cudaEvent_t start, stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);

cudaEventRecord(start);

cout<<"count iter="<<countIter<<endl;


cout.precision(14);
for (int i = 0; i < countIter; i++) {
	 t+=tau;	
     cudaMemcpy(dev_T, &t, sizeof(double), cudaMemcpyHostToDevice);
     calculate<<<countBlocks,countThreads>>>(dev_particles,dev_T,dev_tau,dev_N);
     // calculate<<<1,2>>>(dev_particles,dev_T,dev_tau,dev_N);
     copyParticlesDeviceToHost(particles,dev_particles,size);
     // cout << particles[0].s.x << "   "<< particles[0].s.y <<endl;
     // cout << particles[1].s.x << "   "<< particles[1].s.y <<endl;
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
deleteVecParticlesK<<<1,1>>>(dev_M,dev_N);
return 0;
}