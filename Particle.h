
// Created by andrey on 11.03.18.
//

#ifndef GRAVITATIONALPROBLEM_PARTICLE_H
#define GRAVITATIONALPROBLEM_PARTICLE_H

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

    void calcAndSetSpeedAndCoord(double t){
        s->x+=v->x*t+f->x*t*t/(2*mass);
        s->y+=v->y*t+f->y*t*t/(2*mass);

        v->x += f->x * t/mass;
        v->y += f->y * t/mass;
    }

};


#endif //GRAVITATIONALPROBLEM_PARTICLE_H
