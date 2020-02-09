//
// Created by GQZhou on 4/5/2019.
//

#include "Sdomain.h"
#include <iomanip>
#include "spline.h"



Sdomain::Sdomain() = default;

Sdomain::~Sdomain() = default;


PTB Sdomain::g_sampling(Crystal crystal, Jitter jitter, int smid){
    PTB ptb;
    int rank_i;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_i);
 //   int seed = time(NULL)+rank_i*10000+smid*100000+smid*rank_i;
    int seed = smid+jitter.seed_init;
    std::default_random_engine generator(seed);
    std::normal_distribution<double> pos_dist(jitter.pos_init, jitter.pos_rms);
    std::normal_distribution<double> angle_dist(crystal.bragg, jitter.angle_rms);
   // std::normal_distribution<double> asym_dist(crystal.asymmetry, jitter.asymangle_rms);
//    std::normal_distribution<double> shape_dist(shape_init, jitter.shape_rms);
    ptb.pos =pos_dist(generator);
    ptb.bragg=angle_dist(generator);
    ptb.asym=0;//asym_dist(generator);
 //   ptb.shape=shape_dist(generator);
    return ptb;
}


Shape Sdomain::genshape(Shape shape, Jitter jitter, int smid){
    Shape s_shape;
    s_shape.g_x = shape.g_x;
    int seed = smid+jitter.seed_init*10;
    std::default_random_engine generator(seed);
    double c_shape;
    double c_asym;
    double c_strain;
    for (int it=0; it<shape.g_x.size();it++){
        std::normal_distribution<double> shape_dist(shape.g_h[it], jitter.shape_rms);
        c_shape=shape_dist(generator);
        s_shape.g_h.push_back(c_shape);
        std::normal_distribution<double> asym_dist(shape.g_asym[it], jitter.asymangle_rms);
        c_asym=asym_dist(generator);
        s_shape.g_asym.push_back(c_asym);
        std::normal_distribution<double> strain_dist(shape.g_strain[it], jitter.strain_rms);
        c_strain=strain_dist(generator);
        s_shape.g_strain.push_back(c_strain);
       // cout<<c_strain<<endl;
    }

    return s_shape;
}


CON Sdomain::prepare(Shape s_shape, PTB ptb,vector<double> x, vector<vector<complex<double> > >  efield) {
    CON con;

    double trans = x[x.size()-1]-x[0];
    double trans_project, mesh;
    if (con.n_gx.size()>0){
        vector<double>().swap(con.n_gx);
        vector<double>().swap(con.n_gh);
        vector<double>().swap(con.n_asym);
    }
    trans_project = trans * sin(pi/180*ptb.bragg);
    tk::spline s;
    tk::spline asym;
    tk::spline strain;
    asym.set_points(s_shape.g_x, s_shape.g_asym);
    s.set_points(s_shape.g_x,s_shape.g_h);    // currently it is required that X is already sorted
 //   cout<< s_shape.g_x.size()<<"  "<<s_shape.g_strain.size()<<endl;
    strain.set_points(s_shape.g_x,s_shape.g_strain);
    mesh = trans_project/x.size();
/*    for (int i=0; i<s_shape.g_strain.size();i++){
        cout<<s_shape.g_strain[i]<<endl;
    }*/
    for (int mesh_i=0; mesh_i<x.size()+2; mesh_i++){
        con.n_gx.push_back(ptb.pos-trans_project/2-mesh+mesh_i*mesh);
        con.n_gh.push_back(s(con.n_gx[mesh_i]));
        con.n_asym.push_back(asym(con.n_gx[mesh_i]));
        con.n_strain.push_back(strain(con.n_gx[mesh_i]));
      //  cout<<"sss   "<<strain(con.n_gx[mesh_i])<<endl;
        //    cout<< n_gx[mesh_i]<<"    "<<n_gh[mesh_i]<<"   "<<n_asym[mesh_i]<<endl;
    }
  //  cout<< con.n_gx.size()<<"  "<<con.n_gh.size()<<"   "<<con.n_asym.size()<<endl;
    int dim1 = efield.size();
    int dim2 = efield[0].size();
    vector<vector<double> > tmp_real(dim2, vector<double> (dim1));
    vector<vector<double> > tmp_imag(dim2, vector<double> (dim1));
 //   cout<< tmp.size()<< "  "<< tmp[0].size()<< endl;
    for (int  i =0; i<dim1; i++){
        for (int j=0; j< dim2; j++){
            tmp_real[j][i] = efield[i][j].real();
            tmp_imag[j][i] = efield[i][j].imag();
        }
    }
    con.efieldT_real.swap(tmp_real);
    con.efieldT_imag.swap(tmp_imag);
/*    for (int j=0; j<con.efieldT_real[0].size();j++){
        cout<<setprecision(9)<<con.efieldT_real[0][j]<<endl;
    }*/
    return con;
}


OUT Sdomain::simulate(Crystal crystal, CON con, PTB ptb, vector<double> x, vector<double> freq,int rayid){
    OUT out;
    out.rayid = rayid;
    if ( (out.absfield.size()>0)||(out.abs.size()>0)||(out.R0H.size()>0)||(out.rfield.size()>0)\
    ||(out.tfield.size()>0)||(out.R00.size()>0) ){
        vector<complex<double> >().swap(out.rfield);
        vector<complex<double> >().swap(out.tfield);
        vector<double>().swap(out.absfield);;
        vector<complex<double> >().swap(out.R00);
        vector<complex<double> >().swap(out.R0H);
        vector<double>().swap(out.abs);
    }
    int len = con.efieldT_real.size();
    double x0 = con.n_gx[rayid+1];
    double mesh = con.n_gx[rayid+2]-con.n_gx[rayid+1];
    double shape_prime = (con.n_gh[rayid+2]-con.n_gh[rayid])/2/mesh;
  //  cout<<shape_prime<<endl;
    double eray_prime = tan((180-ptb.bragg)*pi/180);
 //   cout<<eray_prime<<endl;
    double dot = 1*1+shape_prime*eray_prime;
  //  cout<<dot<<endl;
    double costheta = dot/sqrt(1+pow(shape_prime,2))/sqrt(1+pow(eray_prime,2));
   // cout<<costheta<<endl;
    double r_theta = acos(abs(costheta));
    double r_asym = (crystal.asymmetry+con.n_asym[rayid+1])*pi/180;
    complex<double> in_laser;
    double lr;
    double lr2d;
    double alpha;
    complex<double> y;
    complex<double> Y1;
    complex<double> Y2;
    complex<double> R1;
    complex<double> R2;
    complex<double> tR00;
    complex<double> tR0H;
  //  double r_theta = (crystal.bragg+crystal.dbragg)*pi/180; //rad
    double gamma0 = cos(r_theta+r_asym-pi/2);
    double gammaH = cos(r_theta+r_asym+pi/2);
    double act_d = crystal.d*(1.0+con.n_strain[rayid+1]);
   // cout<<crystal.d<<endl;
    double asy_fac = gamma0/gammaH;
    double wavelength = h_Plank * c_speed / crystal.photon_en / e_charge;
    double ang_freq = 2*pi*c_speed/wavelength;
    double wave_num = ang_freq/c_speed;
    complex<double> extin_len = sqrt(gamma0*abs(gammaH))/(wave_num*sqrt(crystal.ele_susH*crystal.ele_susHbar));
    complex<double> A = crystal.thickness/extin_len;
    complex<double> C = exp(I*crystal.ele_sus0*wave_num*crystal.thickness/(2*gamma0));
    complex<double> G = sqrt(abs(asy_fac)*(crystal.ele_susH*crystal.ele_susHbar))/crystal.ele_susHbar;
  //  cout<<"running well"<<endl;
    for (int i =0; i<len; i++){
        tk::spline s_real;
        tk::spline s_imag;
/*        for (int j=0; j<con.efieldT_real[i].size();j++){
            cout<<setprecision(9)<<x[j]<<"  "<<con.efieldT_imag[i][j]<<endl;
        }
        exit(1);*/
        s_real.set_points(x, con.efieldT_real[i]);
     //   cout<<"running"<<endl;
        s_imag.set_points(x, con.efieldT_imag[i]);
        in_laser = s_real(x0)+s_imag(x0)*I;

        lr = c_speed*1e10/freq[i];
     //   cout<<in_laser<<endl;
        lr2d = lr/act_d;
     //   cout<<act_d<<endl;
        alpha = pow(lr2d,2)-2.0*lr2d*sin(r_theta);
       // cout<< r_theta<<endl;
        y = wave_num*extin_len/(2*gamma0)*(asy_fac*alpha+crystal.ele_sus0*(1-asy_fac));
        Y1 = -y-sqrt(pow(y,2)+asy_fac/abs(asy_fac));
        Y2 = -y+sqrt(pow(y,2)+asy_fac/abs(asy_fac));
        R1 = G*Y1;
        R2 = G*Y2;
       /// cout<<R1<<"  "<<R2<<endl;
        tR00 = exp(I*(crystal.ele_sus0*wave_num*crystal.thickness/2.0/gamma0+A/2.0*Y1))*(R2-R1)/(R2-R1*exp(I*A/2.0*(Y1-Y2)));
        tR0H = R1*R2*(1.0-exp(I*A/2.0*(Y1-Y2)))/(R2-R1*exp(I*A/2.0*(Y1-Y2)));
      //  cout<<tR0H<<endl;
        out.R0H.push_back(tR0H);
        out.R00.push_back(tR00);
        out.abs.push_back((1-pow(abs(tR0H),2)-pow(abs(tR00),2)));
        out.rfield.push_back(tR0H*in_laser);
        out.tfield.push_back(tR00*in_laser);
        out.absfield.push_back((1-pow(abs(tR0H),2)-pow(abs(tR00),2))*pow(abs(in_laser),2));
    }
    return out;
}


/*RT Sdomain::interaction(Crystal crystal, vector<double> freq) {
    RT rt;
    string buffer;
    string num;
    string key;
    double lr;
    double lr2d;
    double alpha;
    complex<double> y;
    complex<double> Y1;
    complex<double> Y2;
    complex<double> R1;
    complex<double> R2;
    complex<double> tR00;
    complex<double> tR0H;
    double r_theta = (crystal.bragg+crystal.dbragg)*pi/180; //rad
    double r_asym = (crystal.asymmetry+crystal.dasym)*pi/180;
    double gamma0 = cos(r_theta+r_asym-pi/2);
    double gammaH = cos(r_theta+r_asym+pi/2);
    double act_d = crystal.d*(1.0+crystal.strain);
    double asy_fac = gamma0/gammaH;
    double wavelength = h_Plank * c_speed / crystal.photon_en / e_charge;
    double ang_freq = 2*pi*c_speed/wavelength;
    double wave_num = ang_freq/c_speed;
    complex<double> extin_len = sqrt(gamma0*abs(gammaH))/(wave_num*sqrt(crystal.ele_susH*crystal.ele_susHbar));
    complex<double> A = crystal.thickness/extin_len;
    complex<double> C = exp(I*crystal.ele_sus0*wave_num*crystal.thickness/(2*gamma0));
    complex<double> G = sqrt(abs(asy_fac)*(crystal.ele_susH*crystal.ele_susHbar))/crystal.ele_susHbar;
    for (auto it = 0; it < freq.size(); it++) {
        lr = c_speed*1e10/freq[it];
        lr2d = lr/act_d;
        alpha = pow(lr2d,2)-2.0*lr2d*sin(r_theta);
        y = wave_num*extin_len/(2*gamma0)*(asy_fac*alpha+crystal.ele_sus0*(1-asy_fac));
        Y1 = -y-sqrt(pow(y,2)+asy_fac/abs(asy_fac));
        Y2 = -y+sqrt(pow(y,2)+asy_fac/abs(asy_fac));
        R1 = G*Y1;
        R2 = G*Y2;
        tR00 = exp(I*(crystal.ele_sus0*wave_num*crystal.thickness/2.0/gamma0+A/2.0*Y1))*(R2-R1)/(R2-R1*exp(I*A/2.0*(Y1-Y2)));
        tR0H = R1*R2*(1.0-exp(I*A/2.0*(Y1-Y2)))/(R2-R1*exp(I*A/2.0*(Y1-Y2)));
        rt.R00.push_back(tR00);
        rt.R0H.push_back(tR0H);
        rt.R01.push_back(tR00-C);
    }
    return rt;
}*/
