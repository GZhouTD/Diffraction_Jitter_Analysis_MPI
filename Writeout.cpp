//
// Created by GQZhou on 4/5/2019.
//

#include "Writeout.h"
#include <iomanip>
#include <math.h>
Writeout::Writeout() = default;
Writeout::~Writeout() = default;


/*bool Writeout::writer(RT rt ,const string fname) {
    ofstream f2w(fname);
    if (!f2w) {
        cout << "Unable to write to disk";
        return false;
    }
    for (int it = 0; it < rt.R00.size(); it++) {
        f2w<<fixed<<setprecision(9)<<
        pow(abs(*(rt.R00.begin()+it)),2)<<"    "<<
        pow(abs(*(rt.R0H.begin()+it)),2)<<"    "<<
        pow(abs(*(rt.R01.begin()+it)),2)<<"    "<<
        real(*(rt.R00.begin()+it))<<"    "<<imag(*(rt.R00.begin()+it))<<"    "<<
        real(*(rt.R0H.begin()+it))<<"    "<<imag(*(rt.R0H.begin()+it))<<"    "<<
        real(*(rt.R01.begin()+it))<<"    "<<imag(*(rt.R01.begin()+it))<<endl;
    }
    f2w.close();
    return true;
}*/


bool Writeout::out_init(int num, string fname) {
    string fname_fit;
    string form1 = "_";
    string form3 = ".out";
    for (int fit = 0; fit < num; fit++){
        fname_fit = fname+form1+to_string(fit)+form3;
        ofstream f2w(fname_fit);
        if (!f2w) {
            cout << "Unable to write to disk";
            return false;
        }
        f2w<<"";
        f2w.close();
    };
    ofstream f2w_n("pos_angle.out");
    if (!f2w_n) {
        cout << "Unable to write to disk";
        return false;
    }
    f2w_n<<"";
    f2w_n.close();
    return true;
}

bool Writeout::pos_angle(double pos, double angle) {
    ofstream f2w("pos_angle.out", ios::app);
    if (!f2w) {
        cout << "Unable to write to disk";
        return false;
    }
    f2w<<fixed<<setprecision(9)<<pos<<"    "<<angle<<endl;
    f2w.close();
    return true;
}

bool Writeout::writer(OUT out, const string fname) {
    string fname_fit;
    string form1 = "_";
    string form3 = ".out";
    fname_fit =  fname+form1+to_string(out.rayid)+form3;
   // cout<<fname_fit<<endl;
    ofstream f2w(fname_fit, ios::app);
    if (!f2w) {
        cout << "Unable to write to disk";
        return false;
    }
/*    for (int it = 0; it < out.R0H.size(); it++) {
        f2w<<fixed<<setprecision(9)<<abs(out.R0H[it])<<"    ";
        if (it==out.R0H.size()-1){
            f2w<<endl;
        }
    }*/
    for (int it = 0; it < out.rfield.size(); it++) {
        f2w<<fixed<<setprecision(9)<<out.R0H[it].real()<<"    ";
     /*   if (it==out.rfield.size()-1){
            f2w<<endl;
        }*/
    }
    for (int it = 0; it < out.rfield.size(); it++) {
        f2w << fixed << setprecision(9) << out.R0H[it].imag() << "    ";
        if (it == out.rfield.size() - 1) {
            f2w << endl;
        }
    }
    f2w.close();
return true;
}


bool Writeout::input_disp(Crystal crystal, Jitter jitter, Shape shape, vector<double> freq,\
        vector<double> x, vector<vector<complex<double> > > efield){

    cout<<" ########Crytal Info. ########"<<endl;
    cout<<"cry_thickness: "<<crystal.thickness<<endl;
    cout<<"cry_bragg: "<<crystal.bragg<<endl;
    cout<<"cry_asymmetry: "<<crystal.asymmetry<<endl;
    cout<<"cry_d: "<<crystal.d<<endl;
    cout<<"pho_energy: "<<crystal.photon_en<<endl;
    cout<<"xr0: "<<crystal.xr0<<endl;
    cout<<"xi0: "<<crystal.xi0<<endl;
    cout<<"xrh: "<<crystal.xrh<<endl;
    cout<<"xih: "<<crystal.xih<<endl;
    cout<<"                  "<<endl;
    cout<<"##### Jitter Info. Info. #####"<<endl;
    cout<<"sample: "<<jitter.sample<<endl;
    cout<<"pos_init: "<<jitter.pos_init<<endl;
    cout<<"pos_rms: "<<jitter.pos_rms<<endl;
    cout<<"angle_init: "<<jitter.angle_init<<endl;
    cout<<"angle_rms: "<<jitter.angle_rms<<endl;
    cout<<"asymangle_rms: "<<jitter.asymangle_rms<<endl;
    cout<<"shape_rms: "<<jitter.shape_rms<<endl;
    cout<<"seed_init: "<<jitter.seed_init<<endl;
    cout<<"                  "<<endl;
    cout<<"##### Shape Input Info. ######"<<endl;
    cout<<"shape.x size: "<<shape.g_x.size()<<endl;
    cout<<"shape.h size: "<<shape.g_h.size()<<endl;
    cout<<"shape.asym size: "<<shape.g_asym.size()<<endl;
    cout<<"shape.x[500]: "<<shape.g_x[500]<<endl;
    cout<<"shape.h[900]: "<<shape.g_h[900]<<endl;
    cout<<"shape.asym[700]: "<<shape.g_asym[700]<<endl;
    cout<<"                  "<<endl;
    cout<<"##### Freq. Input Info. ######"<<endl;
    cout<<"freq.size: "<<freq.size()<<endl;
    cout<<"freq[100]: "<<freq[100]<<endl;
    cout<<"                  "<<endl;
    cout<<"####### X Input Info. ########"<<endl;
    cout<<"x.size: "<<x.size()<<endl;
    cout<<"x[100]: "<<x[100]<<endl;
    cout<<"                  "<<endl;
    cout<<"##### Strain Input Info. #####"<<endl;
    cout<<"strain.size: "<<shape.g_strain.size()<<endl;
    cout<<"strain[100]: "<<shape.g_strain[100]<<endl;
    cout<<"                  "<<endl;
    cout<<"##### Efield Input Info. #####"<<endl;
    cout<<"efield.xsize: "<<efield.size()<<"  efield.fsize: "<<efield[0].size()<<endl;
    cout<<"efield[50][777]: "<<efield[0][777]<<endl;
    cout<<"                  "<<endl;
    return  true;
}


