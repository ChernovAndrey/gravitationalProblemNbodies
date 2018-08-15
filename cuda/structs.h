#include <iostream>
#include <string>
#include <typeinfo>


using namespace std;
struct Point{
// public:	
    double x;
    double y;
};

struct Particle {
// public:	
    double mass;
    Point s;// coordinate
    Point v;// speed
    Point f; //сила с которой действуют на эту частицу

};


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

void getInitFor2Particle(Particle* &particles,Particle* &dev_particles){ 
    double x1 =1;
    double y1 = 0;
    double x2 = -1;
    double y2 = 0;
    
    double vx1 = 0;
    double vy1 = 1.0; 
    double vx2 = 0;
    double vy2 = -1.0;
    double m = 1.0;
    particles = (struct Particle*)malloc (2*sizeof (struct Particle));
 //    particles[0].s = (struct Point*)malloc(sizeof(struct Point));
 //    particles[0].v = (struct Point*)malloc(sizeof(struct Point));
	// particles[0].f = (struct Point*)malloc(sizeof(struct Point));
	// particles[1].s = (struct Point*)malloc(sizeof(struct Point));
	// particles[1].v = (struct Point*)malloc(sizeof(struct Point));
	// particles[1].f = (struct Point*)malloc(sizeof(struct Point));
    
    // particles[0].mass = (double*)malloc(sizeof(double));
    // particles[1].mass = (double*)malloc(sizeof(double));

	particles[0].s.x = x1;
	particles[0].s.y = y1;
	particles[0].v.x = vx1;
	particles[0].v.y = vy1;
	particles[0].mass = m;

	particles[1].s.x = x2;
	particles[1].s.y = y2;
	particles[1].v.x = vx2;
	particles[1].v.y = vy2;
	particles[1].mass = m;
       
   // cout<<"particles s[0] "<< particles[0].s.x<<endl;
    //cuda

    
    cudaMalloc(&dev_particles, 2*sizeof (struct Particle));
    // cout<<"after cuda malloc"<<endl;
    // cudaMalloc(&particles[0].s,sizeof(struct Point));
    // cudaMalloc(&particles[0].v,sizeof(struct Point));
    // cudaMalloc(&particles[0].f,sizeof(struct Point));
    // cudaMalloc(&particles[1].s,sizeof(struct Point));
    // cudaMalloc(&particles[1].v,sizeof(struct Point));
    // cudaMalloc(&particles[1].f,sizeof(struct Point));

    // cudaMemcpy(&dev_particles,&particles,sizeof(2*sizeof (struct Particle)),cudaMemcpyHostToDevice);
    // cout<<"after copy on device"<<endl;
    cudaMemcpy(&dev_particles[0].s, &particles[0].s, sizeof(struct Point), cudaMemcpyHostToDevice);
    cudaMemcpy(&dev_particles[0].v, &particles[0].v, sizeof(struct Point), cudaMemcpyHostToDevice);
    cudaMemcpy(&dev_particles[0].f, &particles[0].f, sizeof(struct Point), cudaMemcpyHostToDevice);
    cudaMemcpy(&dev_particles[1].s, &particles[1].s, sizeof(struct Point), cudaMemcpyHostToDevice);
    cudaMemcpy(&dev_particles[1].v, &particles[1].v, sizeof(struct Point), cudaMemcpyHostToDevice);
    cudaMemcpy(&dev_particles[1].f, &particles[1].f, sizeof(struct Point), cudaMemcpyHostToDevice);
    cudaMemcpy(&dev_particles[0].mass, &particles[0].mass, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&dev_particles[1].mass, &particles[1].mass, sizeof(double), cudaMemcpyHostToDevice);

}



// Функция, генерирующая случайное действительное число от min до max
double randomDouble(double min, double max)
{
    return (double)(rand())/RAND_MAX*(max - min) + min;
}

void getInitRandom(Particle* &particles,Particle* &dev_particles,int N){
    particles = (struct Particle*)malloc (N*sizeof (struct Particle));
    cudaMalloc(&dev_particles, N*sizeof (struct Particle));
 
    double maxValueCoord = 100;
    double minValueCoord = -100;
    double maxValueSpeed = 5e-5;
    double minValueSpeed = 5e-6;
    double maxValueMass = 5.0;
    double minValueMass = 1.0;
   
    for(int i=0;i<N;i++){
        double mass = randomDouble(minValueMass,maxValueMass);

        double sx = randomDouble(minValueCoord,maxValueCoord);
        double sy = randomDouble(minValueCoord,maxValueCoord);

        double vx = randomDouble(minValueSpeed,maxValueSpeed);
        double vy = randomDouble(minValueSpeed,maxValueSpeed);


        particles[i].s.x = sx;
        particles[i].s.y = sy;
        particles[i].v.x = vx;
        particles[i].v.y = vy;
        particles[i].mass = mass;
        // cout<<"s="<<particles[i].s.x<<" "<<particles[i].s.y<<endl;
        cudaMemcpy(&dev_particles[i].s, &particles[i].s, sizeof(struct Point), cudaMemcpyHostToDevice);
        cudaMemcpy(&dev_particles[i].v, &particles[i].v, sizeof(struct Point), cudaMemcpyHostToDevice);
        cudaMemcpy(&dev_particles[i].f, &particles[i].f, sizeof(struct Point), cudaMemcpyHostToDevice);
        cudaMemcpy(&dev_particles[i].mass, &particles[i].mass, sizeof(double), cudaMemcpyHostToDevice);
    

        // particles[i] = new Particle(new Point(sx,sy),new Point(vx,vy),mass);

    }
}