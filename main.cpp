#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "Readin.h"
#include <math.h>
#include "scicont.h"
#include "Sdomain.h"
#include "Writeout.h"
#include "mpi.h"
#include "spline.h"
#include <iomanip>
#include <unistd.h>


using namespace std;

int main(int argv,char* argc[]) {
    MPI_Init(&argv,&argc);  // MPI Init
    int rank_i, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_i);  // get rank
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // get MPI size
    //define vars
    vector<double> freq;
    Crystal crystal;
    Jitter jitter;
    Shape shape;
    PTB ptb;
    vector<double> x;
    vector<double> strain;
    bool jgment, disp,pa, wrt;
    vector<vector<complex<double> > > e_field;
    Shape s_shape;
    CON con;
    OUT out;
    //Instantiation
    Sdomain *sdomain = new Sdomain;
    Readin *readin = new Readin;
    Writeout *writeout = new Writeout;

    // read in information
    crystal = readin->crystalreader("crystal.in");
    freq = readin->freqreader("freq.in");
    x = readin->xreader("x.in");
    shape = readin->shapereader("shape.in");
    jitter = readin->jitterreader("jitter.in");
    e_field = readin->fieldreader("laser_real.in");
  //  cout<<e_field.size()<<e_field[0].size()<<endl;
    //define output name;
    string outname = "diffraction";

    if (rank_i==0){
        cout << "Current working path:  " << getcwd(NULL, 0) << endl;
        cout<<"Number of processors: "<<mpi_size<<endl;
        disp = writeout->input_disp(crystal,jitter,shape,freq,x,e_field); //display the inputs.
        jgment = writeout->out_init(x.size(),outname); //generate the outputs file to be inject.
    }
    MPI_Barrier(MPI_COMM_WORLD);  //processors synced

    // sample loop
    for (int smid=0; smid < jitter.sample; smid++){

        ptb = sdomain->g_sampling(crystal,jitter,smid); // generate the samples{}
        s_shape = sdomain->genshape(shape,jitter,smid); //generate the shape changering
        con = sdomain->prepare(s_shape,ptb,x);  // generate a condition for simulation

        if (rank_i==0){
            pa = writeout->pos_angle(ptb.pos,ptb.bragg);  //record the position and input angle
        }
        for (int rayid=0; rayid < x.size(); rayid++) {
            if (rayid % mpi_size == rank_i) {
                if (rayid%5==0 && smid%10==0){
                  cout<<"Sample ID: "<<smid<<" Ray ID: "<<rayid<<endl;  // screen display to show progress
                }
                out = sdomain->simulate(crystal,con,ptb,x,freq,rayid, e_field);  //real simulation !!!!
                wrt = writeout->writer(out,outname);  // generate the output files
                // cout<<"processor "<<rank_i<<"is working: "<<pow(smid,2)<<endl;
                //    ptb = sdomain->sampling(crystal,jitter,shape_i,smid);
            }
        }
    }
  //  rt = sdomain->interaction(crystal,freq);
            //->writer(rt,"out.dat");
    MPI_Finalize();  // MPI stop
    return 0;
}
