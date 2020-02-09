//
// Created by GQZhou on 4/5/2019.
//

#ifndef YURI_WRITEOUT_H
#define YURI_WRITEOUT_H

#include <fstream>
#include <iostream>
#include "Sdomain.h"
#include "mpi.h"

using namespace std;

class Writeout {
public:
    Writeout();
    virtual ~Writeout();
    bool writer(RT, string);
    bool out_init(int, string);
    bool input_disp(Crystal, Jitter, Shape, vector<double>,vector<double>,vector<double>, vector<vector<complex<double> > >);
};


#endif //YURI_WRITEOUT_H
