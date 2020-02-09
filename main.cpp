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
    bool jgment, disp, wrt;
    vector<vector<complex<double> > > e_field;
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
        disp = writeout->input_disp(crystal,jitter,shape,freq,x,e_field);
        jgment = writeout->out_init(x.size(),outname);
    }
    MPI_Barrier(MPI_COMM_WORLD);  //processors synced
    Shape s_shape;
    CON con;
    OUT out;
    for (int smid=0; smid < jitter.sample; smid++){
        ptb = sdomain->g_sampling(crystal,jitter,smid);
        s_shape = sdomain->genshape(shape,jitter,smid);

        con = sdomain->prepare(s_shape,ptb,x,e_field);

        for (int rayid=0; rayid < x.size(); rayid++) {
            if (smid % mpi_size == rank_i) {
                cout<<rayid<<endl;

                out = sdomain->simulate(crystal,con,ptb,x,freq,rayid);
                wrt = writeout->writer(out,outname);
                // cout<<"processor "<<rank_i<<"is working: "<<pow(smid,2)<<endl;
                //    ptb = sdomain->sampling(crystal,jitter,shape_i,smid);

            }
        }
    }
  //  rt = sdomain->interaction(crystal,freq);
            //->writer(rt,"out.dat");
    MPI_Finalize();
    return 0;
}
