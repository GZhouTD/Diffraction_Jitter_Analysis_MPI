//
// Created by GQZhou on 4/4/2019.
//

#ifndef YURI_READIN_H
#define YURI_READIN_H
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "scicont.h"
#include "mpi.h"

using namespace std;
typedef struct {
    double thickness;
    double bragg;
    double dbragg;
    double asymmetry;
    double dasym;
    double d;
    double strain;
    double photon_en;
    double xr0;
    double xi0;
    double xrh;
    double xih;
    complex<double> ele_sus0;
    complex<double> ele_susH;
    complex<double> ele_susHbar;
}Crystal;

typedef struct {
    vector<double> g_x;
    vector<double> g_h;
    vector<double> g_asym;
    vector<double> g_strain;
}Shape;

typedef struct {
    int sample;
    int seed_init;
    double pos_init;
    double pos_rms;
    double angle_init;
    double angle_rms;
    double asymangle_rms;
    double shape_rms;
}Jitter;

class Readin {
public:
    Readin();
    virtual ~Readin();
    Crystal crystalreader(string);
    vector<double> freqreader(string);
    vector<double> xreader(string);
    vector<double> strainreader(string);
    vector<vector<complex<double> > > fieldreader(string);
    Shape shapereader(string);
    Jitter jitterreader(string);
};

#endif //YURI_READIN_H
