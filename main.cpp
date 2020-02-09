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
    int rank_i, mpi_size;
    MPI_Init(&argv,&argc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_i);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    vector<double> freq;
    Crystal crystal;
    Jitter jitter;
    Shape shape;
    vector<double> x;
    vector<double> strain;
    bool jgment, disp;
    RT rt;
    PTB ptb;
    Sdomain *sdomain = new Sdomain;
    Readin *readin = new Readin;
    Writeout *writeout = new Writeout;
    crystal = readin->crystalreader("crystal.in");
    vector<vector<complex<double> > > e_field;

    double wavelength = h_Plank * c_speed / crystal.photon_en / e_charge;
    double w0 = 2 * pi * c_speed / wavelength;
    double cf = w0 / 2 / pi;
    freq = readin->freqreader("freq.in");
    x = readin->xreader("x.in");
    shape = readin->shapereader("shape.in");
    jitter = readin->jitterreader("jitter.in");
    e_field = readin->fieldreader("laser_real.in");
    if (rank_i==0){
        //      cout << "Current working path:  " << getcwd(NULL, 0) << endl;
        //       cout<<"Number of processors: "<<mpi_size<<endl;
   //     disp = writeout->input_disp(crystal,jitter,shape,freq,x,strain,e_field);
        jgment = writeout->out_init(x.size()/10,"diffraction");
       // cout<<jgment<<endl;
    }
    //   cout<<"Processor "<<rank_i<<" is running well!"<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
    Shape s_shape;
    CON con;
    for (int smid=0; smid < jitter.sample; smid++){
        ptb = sdomain->g_sampling(crystal,jitter,smid);
        s_shape = sdomain->genshape(shape,jitter,smid);
        con = sdomain->prepare(s_shape,ptb,x,e_field);
   //     cout<< con.efieldT_real.size()<<"  "<<con.efieldT_real[0].size()<<endl;
        for (int rayid=0; rayid < x.size(); rayid++) {
            if (smid % mpi_size == rank_i) {

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
