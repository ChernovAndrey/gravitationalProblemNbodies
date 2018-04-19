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
    Point *s;// coordinate
    Point* v;// speed
    Point *f; //сила с которой действуют на эту частицу

};

Particle* getInitFor2Particle(){ 
    double x1 =1;
    double y1 = 0;
    double x2 = -1;
    double y2 = 0;
    
    double vx1 = 0;
    double vy1 = 498960e-6; 
    double vx2 = 0;
    double vy2 = -498960e-6;
    Particle * particles = (struct Particle*)malloc (2*sizeof (struct Particle));
    particles[0].s = (struct Point*)malloc(sizeof(struct Point));
    particles[0].v = (struct Point*)malloc(sizeof(struct Point));
	particles[0].f = (struct Point*)malloc(sizeof(struct Point));
	particles[1].s = (struct Point*)malloc(sizeof(struct Point));
	particles[1].v = (struct Point*)malloc(sizeof(struct Point));
	particles[1].f = (struct Point*)malloc(sizeof(struct Point));


	particles[0].s->x = x1;
	particles[0].s->y = y1;
	particles[0].v->x = vx1;
	particles[0].v->y = vy1;
	particles[0].mass = 2;

	particles[1].s->x = x2;
	particles[1].s->y = y2;
	particles[1].v->x = vx2;
	particles[1].v->y = vy2;
	particles[1].mass = 2;
	return particles;
}
