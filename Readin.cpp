//
// Created by GQZhou on 4/4/2019.
//
#include "Readin.h"
#include <unistd.h>
#include <math.h>
#include <iomanip>
#include "scicont.h"



Readin::Readin() = default;

Readin::~Readin() = default;


vector<string> &split(const string &str, char delim, vector<string> &elems, bool skip_empty = true) {
    istringstream iss(str);
    for (string item; getline(iss, item, delim); )
        if (skip_empty && item.empty()) continue;
        else elems.push_back(item);
    return elems;
}

Crystal Readin::crystalreader(const string fname) {
    Crystal crystal;
    string buffer;
    string num;
    string key;
    int rank_i;
    int n;
    const double pi = M_PI;
    ifstream f2r(fname);
    double tmp;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_i);
    if (!f2r) {
        cout << "Unable to open the crystal file";
        exit(1); // terminate with error
    }
    while (getline(f2r, buffer)) {
        if ((!buffer.compare(0, 1, "#") || buffer.length() < 1)) { continue; }

        if (buffer.compare("=") != 0) {
            n = buffer.find("=");

            key = buffer.substr(0, n - 1);
            num = buffer.substr(n + 1);
            tmp = atof(num.c_str());

            if (!key.find("cry_thickness")) { crystal.thickness = tmp; }
            if (!key.find("cry_bragg")) { crystal.bragg = tmp; }
            if (!key.find("cry_asymmetry")) { crystal.asymmetry = tmp; }
            if (!key.find("cry_d")) { crystal.d = tmp; }
            if (!key.find("pho_energy")) { crystal.photon_en = tmp; }
            if (!key.find("xr0")) { crystal.xr0 = tmp; }
            if (!key.find("xi0")) { crystal.xi0 = tmp; }
            if (!key.find("xrh")) { crystal.xrh = tmp; }
            if (!key.find("xih")) { crystal.xih = tmp; }
        }
    }
    crystal.dbragg = 0.0;
    crystal.dasym = 0.0;
    crystal.strain = 0.0;
    crystal.ele_sus0 = crystal.xr0+I*crystal.xi0;
    crystal.ele_susH = crystal.xrh+ I*crystal.xih;
    crystal.ele_susHbar = (pow(crystal.xrh, 2) - pow(crystal.xih, 2)+I*2.0 * abs(crystal.xrh * crystal.xih) * cos(pi))/crystal.ele_susH;

    f2r.close();
    return crystal;
}


Jitter Readin::jitterreader(const string fname) {
    Jitter jitter;
    string buffer;
    double num;
    string key;
    vector<string> values;
    ifstream f2r(fname);
    if (!f2r) {
        cout << "Unable to open the jitter file";
        exit(1); // terminate with error
    }
    while (getline(f2r, buffer)) {
        split(buffer, '=', values);

        //    cout<<n<<endl;
        key = values[0];
        num = atof(values[1].c_str());
        if (!key.find("sample")) { jitter.sample = int(num); }
        if (!key.find("seed_init")) { jitter.seed_init = int(num); }
        if (!key.find("pos_init")) { jitter.pos_init = num; }
        if (!key.find("pos_rms")) { jitter.pos_rms = num; }
        if (!key.find("angle_init")) { jitter.angle_init = num; }
        if (!key.find("angle_rms")) { jitter.angle_rms = num; }
        if (!key.find("asymangle_rms")) { jitter.asymangle_rms = num; }
        if (!key.find("shape_rms")) { jitter.shape_rms = num; }
        if (!key.find("strain_rms")) { jitter.strain_rms = num; }
        vector<string>().swap(values);
    }

    f2r.close();
    return jitter;
}

vector<double> Readin::freqreader(const string fname) {
    fstream f2r(fname);
    string buffer;
    vector<double> freq;
    double tmp;
    if (!f2r) {
        cout << "Unable to open the frequency file";
        exit(1); // terminate with error
    }
    while (getline(f2r, buffer)) {
        tmp = atof(buffer.c_str());
        freq.push_back(tmp);
    }
    f2r.close();
    return freq;
}


vector<double> Readin::xreader(const string fname) {
    fstream f2r(fname);
    string buffer;
    vector<double> x;
    double tmp;
    if (!f2r) {
        cout << "Unable to open the frequency file";
        exit(1); // terminate with error
    }
    while (getline(f2r, buffer)) {
        tmp = atof(buffer.c_str());
        x.push_back(tmp);
    }
    f2r.close();
    return x;
}

/*vector<double> Readin::strainreader(const string fname) {
    fstream f2r(fname);
    string buffer;
    vector<double> strain;
    double tmp;
    if (!f2r) {
        cout << "Unable to open the strain file";
        exit(1); // terminate with error
    }
    while (getline(f2r, buffer)) {
        tmp = atof(buffer.c_str());
        strain.push_back(tmp);
    }
    f2r.close();
    return strain;
}*/

Shape Readin::shapereader(const string fname) {
    fstream f2r(fname);
    string buffer;
    Shape shape;
    vector<string> values;
    double tmp1, tmp2,tmp3, tmp4;
    if (!f2r) {
        cout << "Unable to open the shape file";
        exit(1); // terminate with error
    }
    while (getline(f2r, buffer)) {
        split(buffer, ' ', values);
        tmp1 = atof(values[0].c_str());
        tmp2 = atof(values[1].c_str());
        tmp3 = atof(values[2].c_str());
        tmp4 = atof(values[3].c_str());
        shape.g_x.push_back(tmp1);
        shape.g_h.push_back(tmp2);
        shape.g_asym.push_back(tmp3);
        shape.g_strain.push_back(tmp4);
        vector<string>().swap(values);
    }
    f2r.close();
    return shape;
}

vector<vector<complex<double> > > Readin::fieldreader(const string fname) {
    fstream f2r(fname);
    string buffer;
    vector<string> values;
    vector<vector<complex<double> > > e_field;
    complex<double> tmp;
    if (!f2r) {
        cout << "Unable to open the e-field file";
        exit(1); // terminate with error
    }
    while(getline(f2r, buffer)) {

        split(buffer, ' ', values);
        int len = values.size();
        vector<complex<double> > e_ray(len);
        for (auto it = 0; it < len; it++) {
            tmp = atof(values[it].c_str())+0.0*I;
            e_ray[it]=tmp;
        }
        e_field.push_back(e_ray);
        vector<string>().swap(values);
    }
    f2r.close();
    return e_field;
}