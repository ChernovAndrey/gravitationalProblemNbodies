#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
#include <tuple>
// для эталона tau =0.0001;T0 0; TF 1.0
using namespace std;
# define N 13*512
//#define N 2
# define TAU 0.001
//# define TAU 10
#define T0 0
//#define TF 0.034492738
//# define TF 12.58984937
//# define TF 12.589849537
//# define TF 2*M_PI
#define TF 0.01

//const double G = 1.18555535802194e-04;
const double G = 498217402368e-12; // дни сантиметры килограммы.
//const double G = 4.0; // дни сантиметры килограммы.
//const double G = 0.004981367808;
//const double G = 498136780.8;
// время в днях
bool parallelism_enabled = false;

struct Point {
public:
//    Point(){}
//    Point(double x, double y):x(x),y(y){}
    double x;
    double y;
};
struct Particle {
public:
    double mass;
    Point s;// coordinate
    Point v;// speed
    Point f = {0, 0}; //сила с которой действуют на эту частицу
};


void deleteParticles(vector<Particle *> v) {
    for (int i = 0; i < v.size(); ++i) {
        delete v[i];
    }
}

void writingDataInFile(double x, double y, int k) { // k - номер тела
    const string pathToFile = "/home/andrey/CLionProjects/gravitationalProblem/files/" + to_string(k) + ".txt";
    ofstream file(pathToFile, ios_base::app);
    file.precision(14);
    file << x << " " << y << endl;
    file.close();
}

vector<Particle *> getInitFor2Particle() {
    double x1 = 1;
    double y1 = 0;
    double x2 = -1;
    double y2 = 0;
    vector<Particle *> particles(2);
    particles[0] = new Particle();
    particles[1] = new Particle();
    particles[0]->s = {x1, y1};
    particles[0]->v = {0, 1.0};
    particles[0]->mass = 1.0;
    particles[1]->s = {x2, y2};
    particles[1]->v = {0, -1.0};
    particles[1]->mass = 1.0;
//    particles[0] = new Particle(new Point(x1,y1),new Point(0,1.0),1.0);
//    particles[1] = new Particle(new Point(x2,y2),new Point(0,-1.0),1.0);

    writingDataInFile(x1, y1, 0);
    writingDataInFile(x2, y2, 1);

    return particles;
}

//
//vector<Particle *> getInitFor3Particle(){
//    double x1 = 1;
//    double y1 = 0;
//    double x2 = -1;
//    double y2 = 0;
//    double x3 = 0;
//    double y3 = 0;
//    vector<Particle *> particles(3);
//    particles[0] = new Particle(new Point(x1,y1),new Point(0,498960e-6),2);
//    particles[1] = new Particle(new Point(x2,y2),new Point(0,-498960e-6),2);
//    particles[2] = new Particle(new Point(x3,y3),new Point(0,0),2);
//
//    writingDataInFile(x1,y1,0);
//    writingDataInFile(x2,y2,1);
//    //writingDataInFile(0,0,2);
//
//    return particles;
//}

// Функция, генерирующая случайное действительное число от min до max
double randomDouble(double min, double max) {
    return (double) (rand()) / RAND_MAX * (max - min) + min;
}

vector<Particle *> getInitRandom(){
    vector<Particle *> particles(N);  // освобождать память, выделенную под объекты?

    double maxValueCoord = 100;
    double minValueCoord = -100;
    double maxValueSpeed = 5e-5;
    double minValueSpeed = 5e-6;
    double maxValueMass = 5.0;
    double minValueMass = 1.0;

    for(int i=0;i<N;i++){
        auto mass = randomDouble(minValueMass,maxValueMass);

        auto sx = randomDouble(minValueCoord,maxValueCoord);
        auto sy = randomDouble(minValueCoord,maxValueCoord);

        auto vx = randomDouble(minValueSpeed,maxValueSpeed);
        auto vy = randomDouble(minValueSpeed,maxValueSpeed);

        particles[i] = new Particle();

        particles[i]->s={sx,sy};
        particles[i]->s={sx,sy};
        particles[i]->v={vx,vy};
        particles[i]->mass=mass;

    }
    return particles;
}

void calcForces(vector<Particle *> &particles) {
#pragma omp parallel for if(parallelism_enabled)
    for (int i = 0; i < N; i++) {
        Point f = {0, 0};
        auto pi = particles[i];
        for (int j = 0; j < N; j++) {
            if (i == j) continue;

            auto pj = particles[j];

            double dx = pj->s.x - pi->s.x;
            double dy = pj->s.y - pi->s.y;
            double r = sqrt(dx * dx + dy * dy);
            double Gmij = G * pi->mass * pj->mass;
            f.x += Gmij * dx / (r * r * r);
            f.y += Gmij * dy / (r * r * r);
        }
        particles[i]->f = f;
    }

}
//
//vector<double> NewtonEq(double v,double f, double mass) {
//    return vector<double>({
//                                  v,
//                                  f/mass
//                          });
//}

vector<Particle *> createCopyParticles(const vector<Particle *> &particles) {
    vector<Particle *> particles2(particles.size());
    for (int j = 0; j < particles2.size(); ++j) {
        particles2[j] = new Particle();
        particles2[j]->s = particles[j]->s;
        particles2[j]->v = particles[j]->v;
        particles2[j]->mass = particles[j]->mass;
        particles2[j]->f= particles[j]->f;
    }
    return particles2;
}

//cдвигаем частицы. particles- старые значения. particlesK- нужны для уравнения Ньютона. ответ пишется в particles.
void moveParticle(vector<Particle *> &particles, const vector<vector<Particle *>> &vecParticlesK, vector<double> step) {
    for (int i = 0; i < particles.size(); i++) {
        for (int j = 0; j < vecParticlesK.size(); ++j) {
            particles[i]->s.x += step[j] * vecParticlesK[j][i]->v.x;//подставляем непосредственно уравнение Ньютона.
            particles[i]->s.y += step[j] * vecParticlesK[j][i]->v.y;//подставляем непосредственно уравнение Ньютона.

            particles[i]->v.x += step[j] * vecParticlesK[j][i]->f.x /
                    vecParticlesK[j][i]->mass;//подставляем непосредственно уравнение Ньютона.
            particles[i]->v.y += step[j] * vecParticlesK[j][i]->f.y /
                    vecParticlesK[j][i]->mass;//подставляем непосредственно уравнение Ньютона.
        }
    }
}


void calculateRK1(vector<Particle *> &particles, double tau) {

    calcForces(particles);
//    cout<<"s0 "<<particles[0]->s.x <<'\t'<< particles[0]->s.y <<endl;
//    cout<<"s1 "<<particles[1]->s.x <<'\t'<< particles[1]->s.y <<endl;
//    cout<<"f0 "<<particles[0]->f.x <<'\t'<< particles[0]->f.y <<endl;
//    cout<<"f1 "<<particles[1]->f.x <<'\t'<< particles[1]->f.y <<endl;
    moveParticle(particles,{particles},{tau});

}


void calculateRK2(vector<Particle *> &particles, double tau) {

    calcForces(particles);
    auto particlesK2 = createCopyParticles(particles);
    moveParticle(particlesK2,{particlesK2},{tau/2.0});
    calcForces(particlesK2);
    moveParticle(particles,{particlesK2},{tau});

    deleteParticles(particlesK2);
}

void calculateRK3(vector<Particle *> &particles, double tau) {

    calcForces(particles);

    auto particlesK2 = createCopyParticles(particles);
    moveParticle(particlesK2,{particlesK2},{tau/2.0});
    calcForces(particlesK2);

    auto particlesK3 = createCopyParticles(particles);
    moveParticle(particlesK3,{particles,particlesK2},{-tau,2*tau});
    calcForces(particlesK3);

    moveParticle(particles,{particles,particlesK2,particlesK3},{tau/6.0,4*tau/6.0,tau/6.0});

    deleteParticles(particlesK2);
    deleteParticles(particlesK3);
}


void calculateRK4(vector<Particle *> &particles, double tau) {

    calcForces(particles);

    auto particlesK2 = createCopyParticles(particles);
    moveParticle(particlesK2,{particlesK2},{tau/2.0});
    calcForces(particlesK2);

    auto particlesK3 = createCopyParticles(particles);
    moveParticle(particlesK3,{particlesK2},{tau/2.0});
    calcForces(particlesK3);

    auto particlesK4 = createCopyParticles(particles);
    moveParticle(particlesK4,{particlesK3},{tau});
    calcForces(particlesK4);

    moveParticle(particles,{particles,particlesK2,particlesK3,particlesK4},{tau/6.0,2*tau/6.0,2*tau/6.0,tau/6.0});

    deleteParticles(particlesK2);
    deleteParticles(particlesK3);
    deleteParticles(particlesK4);
}




void cleanFiles(int n) {
    for (int i = 0; i < n; i++) {
        const string pathToFile = "/home/andrey/CLionProjects/gravitationalProblem/files/" + to_string(i) + ".txt";
        ofstream file(pathToFile, ios_base::out | ios_base::trunc);
        file.close();
    }
}



//pair<double,double> getAccuracy2Particles(double t){
////    double period = 12.589849537;//0.0117833
//    double period = 12.5807;
//    double x = cos(2*M_PI*t/period);
//    double y = sin(2*M_PI*t/period);
//    writingDataInFile(x, y, 2);
//    return make_pair(x,y);
//}


pair<double, double> getAccuracy2Particles2(double t) {
//    double period = 12.589849537;//0.0117833
    double period = 2 * M_PI;
    double x = cos(t);
    double y = sin(t);
    writingDataInFile(x, y, 2);
    return make_pair(x, y);
}

double getError2Particles(double x, double y, double t) {
    double xAccur;
    double yAccur;
    tie(xAccur, yAccur) = getAccuracy2Particles2(t);
//    cout << "compare" << endl;
//    cout << "x " << xAccur << " " << x << endl;
//    cout << "y " << yAccur << " " << y << endl;
    cout.precision(12);

    return max(abs(xAccur - x), abs(yAccur - y));
}
//читаем из файла два решения и сравниваем их. tau1=2*tau2
void compareSolution(){
    ifstream file1("/home/andrey/CLionProjects/gravitationalProblem/files/0.txt");
    ifstream file2("/home/andrey/CLionProjects/gravitationalProblem/files/1.txt");
    int n = 100;
    double error=0.0;
    double x1,y1,x2,y2;
    for (int i = 0; i < n; i++) {
        file2 >> x2;
        file2 >> y2;
        file1>>x1;
        file1>>y1;
        error = max(abs(x1-x2),abs(y1-y2));
    }

    cout<<"error files "<< error<<endl;
    file1.close();
    file2.close();
}
int main() {
    cout<<"N="<<N<<endl;
//    compareSolution();
//
//    cout<<"coef: "<<0.000240217/8.06703e-05;
//    return 0;
//    cout.precision(12);
//    srand(time(NULL));
//    cout<<scientific;
    cleanFiles(2);
    auto particles = getInitRandom();
//    auto particles = getInitFor2Particle();
    double t = T0;
    double tau = TAU;
    double tf = TF;
    int countIter = (tf - t) / tau;
    cout << "count Iter=" << countIter << endl;
    auto begin = chrono::high_resolution_clock::now();
    double error = 0;

    for (int i = 0; i < countIter; i++) {
        t += tau;
        calculateRK4(particles, tau);
//
//        writingDataInFile(particles[0]->s.x, particles[0]->s.y, 0);
//        writingDataInFile(particles[1]->s.x, particles[1]->s.y, 1);
//        auto curError = getError2Particles(particles[0]->s.x, particles[0]->s.y, t);
//        if (curError > error) {
//            error = curError;
//        }
    }
//    for (int j = 0; j < particles.size(); ++j) {
//        writingDataInFile(particles[j]->s.x, particles[j]->s.y, 1);
//    }
    auto end = chrono::high_resolution_clock::now();
    cout << "ERROR: " << error << endl;
    auto timens = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();
    cout << timens/(1e6) << "ms" << endl;
    deleteParticles(particles);
    return 0;
}