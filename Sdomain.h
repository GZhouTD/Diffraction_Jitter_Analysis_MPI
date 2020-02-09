//
// Created by GQZhou on 4/5/2019.
//

#ifndef YURI_SDOMAIN_H
#define YURI_SDOMAIN_H
#include "Readin.h"
#include "scicont.h"
#include "math.h"
#include "mpi.h"
#include <time.h>
#include <random>

typedef struct{
    vector<complex<double> > R00;
    vector<complex<double> > R0H;
    vector<complex<double> > R01;
}RT;

typedef struct{
    double pos;
    double bragg;
    double asym;
}PTB;


typedef struct{
    vector<double> n_gx;
    vector<double> n_gh;
    vector<double> n_asym;
    vector<double> n_strain;
    vector<vector<double> > efieldT_real;
    vector<vector<double> > efieldT_imag;
}CON;

typedef struct{
    vector<complex<double> >  rfield;
    vector<complex<double> >  tfield;
    vector<double> absfield;
    vector<complex<double> >  R00;
    vector<complex<double> >  R0H;
    vector<double> abs;
    int rayid;
}OUT;

class Sdomain {
public:
    Sdomain();
    virtual ~Sdomain();
    PTB g_sampling(Crystal, Jitter, int);
    Shape genshape(Shape,Jitter, int);
    OUT simulate(Crystal, CON, PTB, vector<double>,vector<double>, int);
    CON prepare(Shape, PTB, vector<double>, vector<vector<complex<double> > >);
 //   RT interaction(Crystal,vector<double>);
};


#endif //YURI_SDOMAIN_H
