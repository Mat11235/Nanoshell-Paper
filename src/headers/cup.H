#include <complex>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_complex.h>
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#include <vector>


#define cc  299792458.
/** c [m/s] Speed of light **/
#define h  6.626068e-34
/** h [m²kg/s] **/

#define he 6.62606957e-27 	
/** h [erg s] serve h nel sistema di Gauss**/

#define j2eV 6.24150636309e18
#define eV2j 1.60217733000103e-19
#define eV2erg 1.60217657e-12
#define erg2eV 6.24150932e11

//Para la impresion de tau1, tau2

#define hbar 1.0545718e-34       // [J·s]
#define ps 1e12                  // 1 ps = 1e-12 s

double eV2nm(double ome){
  double lam;
  lam=cc/(ome*eV2j/h)*1.e9;
  return lam;
  }

  struct OmegaPoint {
    double t_ps;    // tiempo donde paraste (ps)
    double dip_ss;  // |dip| (o el valor que quieras)
};

class nanosphere{
    private:
	gsl_interp_accel *acc;
	gsl_interp_accel *ach;
	gsl_spline *reeps;
	gsl_spline *imeps;
	gsl_spline *reeph;
	gsl_spline *imeph;
	double *omem;
	int rows, spln;
	char *wrk;
    public:
	std::complex<double> img;
	double act;
	double pi, vF, enha;
	double A, r1, r2, rap;
	double G, nOme;
	double ome_p, Dome, lam_0, eps_0;
	double ome_p2, ome_0, gamd, kap;
	double ome_o2, ome, f, n, enne, omeB;
	double Ome_p, Gam_d, Gam, eps_inf, Ome_L, Gam_L, Deps, DELTAOme;
	std::complex<double> eps3, ceps_inf;
	void init(){
	  /** all values are in MKS **/
	  img = std::complex<double> (0., 1.);
	  pi = acos(-1.);
	  }

//Introduccion de funciones numerical
	void dflt(double* tildeN, double t_ps, double tpump, int* pip){
		if(t_ps >= tpump && *pip ==0) {
			*tildeN = 1.;
			//G = 0.1414;
		 	//GG = set_GG(G, tau2);
			*pip = 1;
		}
	}
	
	void squarewave(double* tildeN, double Tperiod, double duty, double t_ps, double tpump, std::complex<double>* GG, double G, double tau2){
		if (t_ps >= tpump) {
			    // onda cuadrada: tildeN = 1 durante el pulso, luego 0
			    double local_t = fmod(t_ps - tpump, Tperiod);
			    if (local_t < duty * Tperiod) {
			        *tildeN = 1.;  // bomba encendida
				//*GG = set_GG(0.1414, tau2);
			    } else {
			        *tildeN = 0;  // bomba apagada
			    *GG = set_GG(G, tau2);
			    }
			} else {
			    *tildeN = 0;  // antes de tpump, no hay bomba
			}
	}

	void sinusoidalwave(double* tildeN, double Tperiod, double t_ps, std::complex<double>* GG, double G, double tau2) {
	     //Rango [0,1] y N(0)=0
	    *tildeN = 0.5 * (1.0 - std::cos(2.0 * M_PI * t_ps / Tperiod));
	     //(equivalente a 0.5*(1 + sin(2π t/T - π/2)))
	}

		    // if (*tildeN > 0.6) {
	    //     *GG = set_GG(0.15, tau2);
	    // } else {
	    //     *GG = set_GG(G, tau2);
	    // }

//Terminan las funciones numerical

//Comienzan funcion de ISS

std::complex<double> gammaD (std::complex<double>** A, double Ome){
    std::complex<double> gamdD = (A[1][1]-img*Ome)*(A[2][2]-img*Ome)-A[1][2]*A[2][1];
    return gamdD;
    }

std::complex<double> gamma1 (std::complex<double>** A, double Ome){
    std::complex<double> gam1 = (A[1][2]*A[2][0]-A[1][0]*(A[2][2]-img*Ome))/gammaD(A, Ome);
    return gam1;
    }

std::complex<double> gamma2 (std::complex<double>** A, double Ome){
    std::complex<double> gam2 = (A[1][0]*A[2][1]-A[2][0]*(A[1][1]-img*Ome))/gammaD(A, Ome);
    return gam2;
    }

    
std::complex<double> gimmeXi (std::complex<double>* p0, std::complex<double> gam1, std::complex<double> gam2, double tau1){
    std::complex<double> xi = tau1*(p0[0]+p0[1]*gam1+p0[2]*gam2);
    return xi;
    }

double find_Omega(double omeeV, char* hst, char* sol=0, double rho=0){
    std::complex<double> GG, OmeH, GamP, OmeP, **A, *p0, *p1, *p2, gam1, gam2, xi, *vOme, *kap, *odv;
    double tau1, tau2, gamd, ome, ome0, eps_b, eps_s, *D, **R, **I, *c;
    double Ome, equ1, equ2, kamp, dmem, *nv;
    
    eps_b=this->set_host(hst);
    eps_s=this->set_host(sol);
    
    A = new std::complex<double>*[3];
    for(int i = 0; i <= 2; i++)
        A[i] = new std::complex<double>[3];

    c = new double [5];

    D = new double [5];

    R = new double*[3];
    for(int i = 0; i < 3; i++)
        R[i] = new double[4];

    I = new double*[3];
    for(int i = 0; i < 3; i++)
        I[i] = new double[4];

    vOme = new std::complex<double>[4];
    kap  = new std::complex<double>[3];
    nv   = new double[4];
    odv  = new std::complex<double> [3];

    // normalized variables
    nv = this->normalized_variables ();
    tau1 = nv[0];
    tau2 = nv[1];
    gamd = nv[2];
    ome0 = nv[3];
    GG   = this->set_GG(this->G, tau2);

    ome  = omeeV/this->Ome_p;

    odv  = this->set_ome_dep_vrbls(ome, ome0, tau2, gamd);
    OmeH = odv[0];
    OmeP = odv[1];
    GamP = odv[2];
    
    p0 = pcfc0(this->eps_inf, eps_b, eps_s, rho);
    p1 = pcfc1(this->eps_inf, eps_b, eps_s, rho);
    p2 = pcfc2(this->eps_inf, eps_b, eps_s, rho);

    A   = coefficients(OmeH, OmeP, GG, GamP, p0, p1, p2);
        
    D    = Dsqr(p1[1], p1[2], p2[1], p2[2], ome, gamd);
    R    = Rnum(p1, p2, ome, gamd);
    I    = Inum(p1, p2, ome, gamd);
    c    = cffcOme (D, R, I, p0, tau2, ome, ome0);
    

    
    vOme = solve4(c);
    kap  = eigenvalues(A);

    Ome=0.;

    kamp = 66.;
    dmem = 0.1;
    for (int i = 0; i < 3; i++) if (kap[i].real() > 0) kamp = kap[i].imag();
    if (kamp != 66){
        for (int i = 0; i < 5; i++){
            if (vOme[i].imag()==0){
                if (fabs(vOme[i].real()-kamp)<dmem) {
                    Ome = vOme[i].real();
                    dmem = fabs(vOme[i].real()-kamp);         
                    }
                }
            }
        }   
    gam1 = gamma1(A, Ome);
    gam2 = gamma2(A, Ome);
    xi   = gimmeXi(p0, gam1, gam2, tau1);
    equ1 = norm(gammaD(A, Ome))*(real(xi)+tau2*(Ome-ome+ome0)*imag(xi))/tau1;
    equ2 = c[0]+c[1]*Ome+c[2]*pow(Ome,2)+c[3]*pow(Ome,3)+c[4]*pow(Ome,4);
    if (abs(equ1-equ2)>=1.e-18 && fabs(Ome)> 1.e-30){
        std::cout.precision(20);
        std::cout<<"> Fatal! The error on the coefficient equation is larger than "<<abs(equ1-equ2)<<std::endl;
        std::cout<<"> ome ="<<omeeV<<std::endl;
        exit(-66);
        }
    
    
    
    return Ome;
    }

double* Dsqr(std::complex<double> p11, std::complex<double> p12, std::complex<double> p21, std::complex<double> p22, std::complex<double> ome, std::complex<double> gamd) {
    double* D;
    D = new double[5];

    D[0] = real(pow(ome, 8.) +
           (8. * pow(gamd, 2.) + 2. * p22 + 2. * p11) * pow(ome, 6.) +
           (16. * pow(gamd, 4.) + (8. * p22 + 8. * p11) * pow(gamd, 2.) + pow(p22, 2.) + 4. * p11 * p22 - 2. * p12 * p21 + pow(p11, 2.)) * pow(ome, 4.) +
           ((4. * pow(p22, 2.) + 8. * p12 * p21 + 4. * pow(p11, 2.)) * pow(gamd, 2.) + 2. * p11 * pow(p22, 2.) +
            (2. * pow(p11, 2.) - 2. * p12 * p21) * p22 - 2. * p11 * p12 * p21) * pow(ome, 2.) +
           pow(p11, 2.) * pow(p22, 2.) - 2. * p11 * p12 * p21 * p22 + pow(p12, 2.) * pow(p21, 2.));

    D[1] = real(-8. * pow(ome, 7.) +
           (-48. * pow(gamd, 2.) - 12. * p22 - 12. * p11) * pow(ome, 5.) +
           (-64. * pow(gamd, 4.) + (-32. * p22 - 32. * p11) * pow(gamd, 2.) - 4. * pow(p22, 2.) - 16. * p11 * p22 + 8. * p12 * p21 - 4. * pow(p11, 2.)) * pow(ome, 3.) +
           ((-8. * pow(p22, 2.) - 16. * p12 * p21 - 8. * pow(p11, 2.)) * pow(gamd, 2.) - 4. * p11 * pow(p22, 2.) +
            (4. * p12 * p21 - 4. * pow(p11, 2.)) * p22 + 4. * p11 * p12 * p21) * ome);

    D[2] = real(24. * pow(ome, 6.) +
           (104. * pow(gamd, 2.) + 24. * p22 + 24. * p11) * pow(ome, 4.) +
           (96. * pow(gamd, 4.) + (40. * p22 + 40. * p11) * pow(gamd, 2.) +
            4. * pow(p22, 2.) + 16. * p11 * p22 - 8. * p12 * p21 + 4. * pow(p11, 2.)) * pow(ome, 2.) +
           (4. * pow(p22, 2.) + 8. * p12 * p21 + 4. * pow(p11, 2.)) * pow(gamd, 2.));

    D[3] = real(-32. * pow(ome, 5.) +
           (-96. * pow(gamd, 2.) - 16. * p22 - 16. * p11) * pow(ome, 3.) +
           ((-16. * p22 - 16. * p11) * pow(gamd, 2.) - 64. * pow(gamd, 4.)) * ome);

    D[4] = real(16. * pow(ome, 4.) + 32. * pow(gamd, 2.) * pow(ome, 2.) + 16. * pow(gamd, 4.));

    return D;
}

double** Rnum(std::complex<double>* p1, std::complex<double>* p2, std::complex<double> ome, std::complex<double> gamd) {
    double** R;
    R = new double*[3];
    for (int i = 0; i < 3; i++)
        R[i] = new double[4];

    R[1][0] = real(-p1[0] * pow(ome, 6.) +
              (-4. * p1[0] * pow(gamd, 2.) - 2. * p1[0] * p2[2] + p1[2] * p2[0] - p1[0] * p1[1]) * pow(ome, 4.) +
              ((-4. * p1[2] * p2[0] - 4. * p1[0] * p1[1]) * pow(gamd, 2.) - p1[0] * pow(p2[2], 2.) +
               (p1[2] * p2[0] - 2. * p1[0] * p1[1]) * p2[2] + p1[0] * p1[2] * p2[1] + p1[1] * p1[2] * p2[0]) * pow(ome, 2.) -
              p1[0] * p1[1] * pow(p2[2], 2.) + (p1[0] * p1[2] * p2[1] + p1[1] * p1[2] * p2[0]) * p2[2] - pow(p1[2], 2.) * p2[0] * p2[1]);

    R[1][1] = real(6. * p1[0] * pow(ome, 5.) +
              (16. * p1[0] * pow(gamd, 2.) + 8. * p1[0] * p2[2] - 4. * p1[2] * p2[0] + 4. * p1[0] * p1[1]) * pow(ome, 3.) +
              ((8. * p1[2] * p2[0] + 8. * p1[0] * p1[1]) * pow(gamd, 2.) + 2. * p1[0] * pow(p2[2], 2.) +
               (4. * p1[0] * p1[1] - 2. * p1[2] * p2[0]) * p2[2] - 2. * p1[0] * p1[2] * p2[1] - 2. * p1[1] * p1[2] * p2[0]) * ome);

    R[1][2] = real(-12. * p1[0] * pow(ome, 4.) +
              (-20. * p1[0] * pow(gamd, 2.) - 8. * p1[0] * p2[2] + 4. * p1[2] * p2[0] - 4. * p1[0] * p1[1]) * pow(ome, 2.) +
              (-4. * p1[2] * p2[0] - 4. * p1[0] * p1[1]) * pow(gamd, 2.));

    R[1][3] = real(8. * p1[0] * pow(ome, 3.) + 8. * p1[0] * pow(gamd, 2.) * ome);

    R[2][0] = real(-p2[0] * pow(ome, 6.) +
              (-4. * p2[0] * pow(gamd, 2.) - p2[0] * p2[2] + p1[0] * p2[1] - 2. * p1[1] * p2[0]) * pow(ome, 4.) +
              ((-4. * p2[0] * p2[2] - 4. * p1[0] * p2[1]) * pow(gamd, 2.) + (p1[0] * p2[1] - 2. * p1[1] * p2[0]) * p2[2] +
               (p1[2] * p2[0] + p1[0] * p1[1]) * p2[1] - pow(p1[1], 2.) * p2[0]) * pow(ome, 2.) +
              (p1[0] * p1[1] * p2[1] - pow(p1[1], 2.) * p2[0]) * p2[2] - p1[0] * p1[2] * pow(p2[1], 2.) + p1[1] * p1[2] * p2[0] * p2[1]);

    R[2][1] = real(6. * p2[0] * pow(ome, 5.) +
              (16. * p2[0] * pow(gamd, 2.) + 4. * p2[0] * p2[2] - 4. * p1[0] * p2[1] + 8. * p1[1] * p2[0]) * pow(ome, 3.) +
              ((8. * p2[0] * p2[2] + 8. * p1[0] * p2[1]) * pow(gamd, 2.) + (4. * p1[1] * p2[0] - 2. * p1[0] * p2[1]) * p2[2] +
               (-2. * p1[2] * p2[0] - 2. * p1[0] * p1[1]) * p2[1] + 2. * pow(p1[1], 2.) * p2[0]) * ome);

    R[2][2] = real(-12. * p2[0] * pow(ome, 4.) +
              (-20. * p2[0] * pow(gamd, 2.) - 4. * p2[0] * p2[2] + 4. * p1[0] * p2[1] - 8. * p1[1] * p2[0]) * pow(ome, 2.) +
              (-4. * p2[0] * p2[2] - 4. * p1[0] * p2[1]) * pow(gamd, 2.));

    R[2][3] = real(8. * p2[0] * pow(ome, 3.) + 8. * p2[0] * pow(gamd, 2.) * ome);

    return R;
}

double** Inum(std::complex<double>* p1, std::complex<double>* p2, std::complex<double> ome, std::complex<double> gamd) {
    double** I;
    I = new double*[3];
    for (int i = 0; i < 3; i++)
        I[i] = new double[4];

    I[1][0] = real(2. * p1[0] * gamd * pow(ome, 5.) +
              (8. * p1[0] * pow(gamd, 3.) + (4. * p1[0] * p2[2] - 4. * p1[2] * p2[0]) * gamd) * pow(ome, 3.) +
              (2. * p1[0] * pow(p2[2], 2.) - 2. * p1[2] * p2[0] * p2[2] + 2. * p1[0] * p1[2] * p2[1] - 2. * p1[1] * p1[2] * p2[0]) * gamd * ome);

    I[1][1] = real(-10. * p1[0] * gamd * pow(ome, 4.) +
              ((12. * p1[2] * p2[0] - 12. * p1[0] * p2[2]) * gamd - 24. * p1[0] * pow(gamd, 3.)) * pow(ome, 2.) +
              (-2. * p1[0] * pow(p2[2], 2.) + 2. * p1[2] * p2[0] * p2[2] - 2. * p1[0] * p1[2] * p2[1] + 2. * p1[1] * p1[2] * p2[0]) * gamd);

    I[1][2] = real(16. * p1[0] * gamd * pow(ome, 3.) +
              (24. * p1[0] * pow(gamd, 3.) + (8. * p1[0] * p2[2] - 8. * p1[2] * p2[0]) * gamd) * ome);

    I[1][3] = real(-8. * p1[0] * gamd * pow(ome, 2.) - 8. * p1[0] * pow(gamd, 3.));

    I[2][0] = real(2. * p2[0] * gamd * pow(ome, 5.) +
              (8. * p2[0] * pow(gamd, 3.) + (4. * p1[1] * p2[0] - 4. * p1[0] * p2[1]) * gamd) * pow(ome, 3.) +
              (-2. * p1[0] * p2[1] * p2[2] + (2. * p1[2] * p2[0] - 2. * p1[0] * p1[1]) * p2[1] + 2. * pow(p1[1], 2.) * p2[0]) * gamd * ome);

    I[2][1] = real(-10. * p2[0] * gamd * pow(ome, 4.) +
              ((12. * p1[0] * p2[1] - 12. * p1[1] * p2[0]) * gamd - 24. * p2[0] * pow(gamd, 3.)) * pow(ome, 2.) +
              (2. * p1[0] * p2[1] * p2[2] + (2. * p1[0] * p1[1] - 2. * p1[2] * p2[0]) * p2[1] - 2. * pow(p1[1], 2.) * p2[0]) * gamd);

    I[2][2] = real(16. * p2[0] * gamd * pow(ome, 3.) +
              (24. * p2[0] * pow(gamd, 3.) + (8. * p1[1] * p2[0] - 8. * p1[0] * p2[1]) * gamd) * ome);

    I[2][3] = real(-8. * p2[0] * gamd * pow(ome, 2.) - 8. * p2[0] * pow(gamd, 3.));

    return I;
}

double *cffcOme (double *D, double **R, double **I, std::complex<double>* p0, double tau2, double ome, double ome0){
    double *c;
    c = new double[5];

    c[0] = real(-(p0[2]*I[2][0]+p0[1]*I[1][0])*(ome-ome0)*tau2+p0[2]*R[2][0]+p0[1]*R[1][0]+D[0]*p0[0]);
    c[1] = real(-(p0[2]*I[2][1]+p0[1]*I[1][1])*(ome-ome0)*tau2+p0[2]*(I[2][0]*tau2+R[2][1])+
                                                          p0[1]*(I[1][0]*tau2+R[1][1])+p0[0]*D[1]);
    c[2] = real(-(p0[2]*I[2][2]+p0[1]*I[1][2])*(ome-ome0)*tau2+p0[2]*(I[2][1]*tau2+R[2][2])+
                                                          p0[1]*(I[1][1]*tau2+R[1][2])+p0[0]*D[2]);
    c[3] = real(-(p0[2]*I[2][3]+p0[1]*I[1][3])*(ome-ome0)*tau2+p0[2]*(I[2][2]*tau2+R[2][3])+
                                                          p0[1]*(I[1][2]*tau2+R[1][3])+p0[0]*D[3]);
    c[4] = real((p0[2]*I[2][3]+p0[1]*I[1][3])*tau2+p0[0]*D[4]);

    return c;
    }

std::complex<double>* simple_intensity_steady_state(char* mdl, char* mtl, char* hst, double omeeV, char* sol,double rho, double E0){
    std::complex<double> GG, OmeH, GamP, OmeP, gam1, gam2, xi, **A, cq0s, *vec, *odv;
    double q0s, p3s, tau1, tau2, gamd, ome, ome0, ome1, ome2, Ome, *nv;
    double eps_b, eps_s, *fro;
	std::complex<double> *p0, *p1, *p2, *p3;

	Ome = find_Omega(omeeV, hst, sol, rho);
	ome1 = -1e9;
	ome2 = +1e9;

    eps_b=this->set_host(hst);
    eps_s=this->set_host(sol);
    
    A = new std::complex<double>*[3];
    for(int i = 0; i <= 2; i++){
        A[i] = new std::complex<double>[3];
	}

	vec = new std::complex<double>[4];
    nv  = new double[4];
    odv = new std::complex<double> [3];
    
    fro=this->frohlich(2., 3.4, eps_b, eps_s, rho);
    // normalized variables
    nv = this->normalized_variables ();
    tau1 = nv[0];
    tau2 = nv[1];
    gamd = nv[2];
    ome0 = nv[3];
    GG = this->set_GG(this->G, tau2);
    
    ome  = omeeV/this->Ome_p;

    odv  = this->set_ome_dep_vrbls(ome, ome0, tau2, gamd);
    OmeH = odv[0];
    OmeP = odv[1];
    GamP = odv[2];
    
    p0 = pcfc0(this->eps_inf, eps_b, eps_s, rho);
    p1 = pcfc1(this->eps_inf, eps_b, eps_s, rho);
    p2 = pcfc2(this->eps_inf, eps_b, eps_s, rho);
    p3 = pcfc3(this->eps_inf, eps_b, eps_s, rho);

    A   = coefficients(OmeH, OmeP, GG, GamP, p0, p1, p2);

    gam1 = gamma1(A, Ome);
    gam2 = gamma2(A, Ome);
    xi   = gimmeXi(p0, gam1, gam2, tau1);

    
    if (ome<ome1 || ome>ome2){
        cq0s.real(0);
        cq0s.imag(0);
        } else cq0s = 1./xi.imag()+tau1*(OmeH-img*Ome)/(GG*xi*xi.imag());
    
    if (fabs(cq0s.imag())>5.e-5 && fabs(this->G) > fabs(fro[1])){
        std::cout.precision(20);
        std::cout<<"> Warning! The error on the imaginary part of |q1|² is larger than 5.e-5"
        <<std::endl<<"> err = "<<cq0s.imag()<<std::endl;
        std::cout<<"> ome ="<<omeeV<<", Ome = "<<Ome<<std::endl;
        } else if (fabs(cq0s.imag())>1. && fabs(this->G) > fabs(fro[1])){
            std::cout.precision(20);
            std::cout<<"> Fatal! The error on the imaginary part of |q1|² is larger than 1"
            <<std::endl<<"> err = "<<cq0s.imag()<<std::endl;
            std::cout<<"> ome ="<<omeeV<<", Ome = "<<Ome<<std::endl;
            exit(-67); 
            } else q0s = cq0s.real();
    

	std::complex<double> q0, q1, q2, p3mod;
	q0 = std::sqrt(q0s);
	q1 = gam1*q0;
	q2 = gam2*q0;

	vec[0] = {1,0};
    vec[1] = gam1;
    vec[2] = gam2;
    vec[3] = {0,0};

	p3s = norm(numerical_output(E0, vec, p3))*q0s;
	p3mod = std::sqrt(p3s);

    auto* out = new std::complex<double>[4]{ q0, q1, q2, p3mod };
    return out;
    }

//Aqui termina la definicion de las funciones ISS

void set_metal(const char* mtl, const char* mdl, int sel=0){
	  /** Ome_p [THz]	Plasma frequency (gold)		**/
	  /** Gam_d [THz]	Drude damping coefficient (gold)**/
	  double *reps, *ieps, la, Dr, Di;
	  char *jc;
	  rows=0;
	  jc = new char [35];
	  if(strcmp(mtl,"gold")==0){
	      Ome_p=8.91;
	      Gam_d=0.0759;
	      eps_inf=9.0685;
	      if(strcmp(mdl,"drude")==0){
		} else if(strcmp(mdl,"spline")==0){
		spln = 1;
		strcpy(jc,"../data/input/goldJCeV.dat");
			std::ifstream inp(jc);
			if (!inp) {
				std::cerr << "Error: could not open " << jc << "\n";
				exit(1);
			}

			std::vector<double> omem_vec, reps_vec, ieps_vec;
			double om, re, im;

			std::string line;
			while (std::getline(inp, line)) {
				std::istringstream iss(line);
				if (iss >> la >> om >> re >> im >> Dr >> Di) {
					omem_vec.push_back(om);
					reps_vec.push_back(re);
					ieps_vec.push_back(im);
				}
			}

			inp.close();

			rows = omem_vec.size();
			omem = new double[rows];
			reps = new double[rows];
			ieps = new double[rows];

			for (int i = 0; i < rows; ++i) {
				omem[i] = omem_vec[i];
				reps[i] = reps_vec[i];
				ieps[i] = ieps_vec[i];
			}

			acc   = gsl_interp_accel_alloc();
			reeps = gsl_spline_alloc(gsl_interp_cspline, rows);
			imeps = gsl_spline_alloc(gsl_interp_cspline, rows);

			gsl_spline_init(reeps, omem, reps, rows);
			gsl_spline_init(imeps, omem, ieps, rows);
		}
	      else {
		std::cout<<"you chosed "<<mdl<<" while option for gold are: drude, drude-lorentz, spline"<<std::endl;
		exit(1);
		}
	    } else  if(strcmp(mtl,"silver")==0){
	      Ome_p=9.6;
	      Gam_d=0.0228;
	      eps_inf=5.3;
		if(strcmp(mdl,"drude")==0){
		} else if (strcmp(mdl, "spline") == 0) {
		    spln = 1;
			strcpy(jc, "../data/input/silverJCeV.dat");
			std::ifstream inp(jc);
			if (!inp) {
				std::cerr << "Error: could not open " << jc << "\n";
				exit(1);
			}

			std::vector<double> omem_vec, reps_vec, ieps_vec;
			double la, Dr, Di, om, re, im;

			std::string line;
			while (std::getline(inp, line)) {
				std::istringstream iss(line);
				if (iss >> la >> om >> re >> im >> Dr >> Di) {
					omem_vec.push_back(om);
					reps_vec.push_back(re);
					ieps_vec.push_back(im);
				}
			}

			inp.close();

			rows = omem_vec.size();
			omem = new double[rows];
			reps = new double[rows];
			ieps = new double[rows];

			for (int i = 0; i < rows; ++i) {
				omem[i] = omem_vec[i];
				reps[i] = reps_vec[i];
				ieps[i] = ieps_vec[i];
			}

			acc   = gsl_interp_accel_alloc();
			reeps = gsl_spline_alloc(gsl_interp_cspline, rows);
			imeps = gsl_spline_alloc(gsl_interp_cspline, rows);

			gsl_spline_init(reeps, omem, reps, rows);
			gsl_spline_init(imeps, omem, ieps, rows);
		}
	      else {
		std::cout<<"you chosed "<<mdl<<" while option for silver are: spline"<<std::endl;
		exit(1);
		}
	    }
	    else{
		std::cout<<mtl<<" not found, option are: gold, silver"<<std::endl;
		exit(1);
		}
	  ome_p=Ome_p; //corretto
	  ome_p2=ome_p*ome_p;
	  }

	double set_host(const char* hst){
	  double eps_0;
	  if(strcmp(hst,"silica")==0) eps_0 = 2.1316; // Si02
	  else if(strcmp(hst,"solgel")==0) eps_0 = 2.1331;//3.8
	  else if(strcmp(hst,"water")==0) eps_0 = 1.7689;
	  else if(strcmp(hst,"glass")==0) eps_0 = 1.2247;
	  else if(strcmp(hst,"PMMA")==0) eps_0 = 2.2201;
	  else if(strcmp(hst,"ethanol")==0) eps_0 = 1.8496;//1.1662;
	  else{
		std::cout<<"option are: silica, solgel, water, glass, PMMA, ethanol"<<std::endl;
		exit(1);
		}
	  return eps_0;
	  }
	  
	void set_active(const char* mod){
	  if(strcmp(mod,"lorentz")==0) act=1.;
	    else if(strcmp(mod,"flat")==0) act=0.;
	    else{
	      std::cout<<"option are: flat, lorentz"<<std::endl;
	      exit(1);
	      }
	  }

std::complex<double> metal(double ome) {
    std::complex<double> eps;

    if (spln == 1) {
        double x_min = omem[0];
        double x_max = omem[rows - 1];

        if (std::real(ome) < x_min || std::real(ome) > x_max) {
            std::cerr << "⚠️ GSL interpolation error: omega = " << std::real(ome)
                      << " is outside the interpolation range [" << x_min << ", " << x_max << "]\n";
            std::cerr << "Aborting safely to avoid crash.\n";
            std::exit(EXIT_FAILURE);
        }

        eps.real(gsl_spline_eval(reeps, ome, acc));
        eps.imag(gsl_spline_eval(imeps, ome, acc));
    } else {
        eps = eps_inf - ome_p2 / (ome * (ome + img * Gam_d));
    }

    ceps_inf = eps + ome_p2 / (ome * (ome + img * Gam_d));
    return eps;
}

	  
	std::complex<double> confinement(double ome){
	  std::complex<double> dlt;
	  dlt=ome_p2/(ome*(ome+img*Gam_d))-ome_p2/(ome*(ome+img*(Gam_d+Gam)));
	  return dlt;
	  }


	std::complex<double> active(double ome, double epsh){
	  std::complex<double> eps;
      double nG = fabs(G);
      eps = std::complex<double>(epsh, 0.) + act*nG*Dome/(2.*(ome - ome_0) + img*Dome) - (1 - act)*nG*img;
      return eps;
    }
	  

    double *frohlich(double a0, double b0, double eps_b, double eps_s=0, double rap=0){
	static double reso[2]={};
	double ome_sp, Gth;
	double eqn9_A, eqn9_D, a, b, d, e;    //declare some needed variables
	double rap3=rap*rap*rap;
	std::complex<double> eps1_A, eps1_B, eps1_D, F_A, F_B, F_D;
	e = 1.e-10;

	double ome_sm=666.;
	double ome21;//, eps_im
	a=a0;
	b=b0;
	ome_sp =.5*(b0+a0);
	ome_sm=666.;
	while (abs(ome_sp-ome_sm)>1.e-20){ /* ome_sp dipende da rho e da ome21, riposiziono ome21 su ome_sp in modo che sia solo funzione di rho*/
	    ome21=ome_sp;
	    a=a0;
	    b=b0;
	    while (fabs(a-b)>=e){   /*if the mod of a-b is greater than the accuracy desired keep bisecting the interval*/
		
		d=0.5*(a+b);        //bisect the interval and find the value of c

		eps1_A = metal(a);
		F_A = F(eps1_A, eps_s, rap3);
		eqn9_A = eqn9(a, ome21, F_A, eps_b, Dome);

		eps1_D = metal(d);
		F_D = F(eps1_D, eps_s, rap3);
		eqn9_D = eqn9(d, ome21, F_D, eps_b, Dome);        

		if (eqn9_D*eqn9_A<=0.){    
		    b=d;    
		    }
		else{
		    a=d;    
		    }
		}
	ome_sm=ome_sp;
	ome_sp = d;   // spasing freq
	}
	ome_sp   = ome_sp;
	Gth = -Dome*imag(F_D)/Dome;
	
	reso[0]=ome_sp;
	reso[1]=Gth;
	
	return reso;
	}

    std::complex<double> set_GG(double G, double tau2){
		std::complex<double> GG;
		GG   = -img*fabs(G)/tau2;
		return GG;
		}
    
	double *normalized_variables (){
		double *nv, tau1, tau2, gamd, ome0;
		nv = new double[4];
				
		tau2 = 2.*Ome_p/Dome;
		tau1 = 5.*tau2;
		
		gamd = .5*Gam_d/Ome_p;
		ome0 = ome_0/Ome_p;
		
		nv[0] = tau1;
		nv[1] = tau2;
		nv[2] = gamd;
		nv[3] = ome0;
		
		return nv;
		}

	double Rabi_frec(double G, double tau2, double E0){
		double omeR = 2.*sqrt(fabs(G)/tau2)*E0*1.e7;
		return omeR;
		}
		
	 std::complex<double> *set_ome_dep_vrbls(double ome, double ome0, double tau2, double gamd){
		std::complex<double> OmeH, GamP, OmeP;
		std::complex<double> *odv;
		odv = new std::complex<double> [3];
		
		OmeH = img*(ome-ome0)-1./tau2;
		GamP = 1./(2.*(gamd-img*ome));
		OmeP = (ome*(ome+2.*img*gamd))*GamP; 
		
		odv[0] = OmeH;
		odv[1] = OmeP;
		odv[2] = GamP;
		
		return odv;
		}
		
	std::vector<std::pair<double,std::complex<double>>>  steady_state(char* mdl , char* mtl ,char* hst ,double omemi, double omema, int omeN=10000, char* sol=0, double rho=0){
	    int omi, ns=0;
	    double eps_b, eps_s, ome, omeeV, dome  = (omema-omemi)/omeN;

	    if ((sol != NULL)) {
		eps_s=set_host(sol);
		ns=1;
		}
		else eps_s=0;
	    eps_b=set_host(hst);
	    std::complex<double> eps1, eps2, alph;

	    std::vector<std::pair<double,std::complex<double>>> valph;
	    
	    double tau2, gamd, ome0, eps3, pG = fabs(G), *nv;
	    std::complex<double> GG, OmeH, GamP, OmeP, **A, detA, *p0, *p1, *p2, *kap, *odv;
	    
	    A = new std::complex<double>*[3];
	    for(int i = 0; i <= 2; i++)
		A[i] = new std::complex<double>[3];
	    
	    kap = new std::complex<double>[3];
	    
		nv = new double[4];

		odv = new std::complex<double> [3];

	    std::fstream comp, stat, egva;
	    stat.open("../data/output/stationary.dat", std::ios::out);
	    comp.open("../data/output/compounds.dat", std::ios::out);
	    egva.open("../data/output/eigenvalues.dat", std::ios::out);
	    
	    stat<<"# * STATIC * "<<std::endl;
	    stat<<"#"<<std::endl;
	    stat<<"# file generated by ** cup.H ** "<<std::endl;
	    stat<<"#"<<std::endl;
	    stat<<"# PARAMETERS: "<<std::endl;
	    stat<<"# G = "<<pG<<std::endl;
	    stat<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
	    stat<<"# ------------- "<<std::endl;
	    stat<<"# ome (eV)\treal(alph)\timag(alph)"<<std::endl;
	    stat<<"#"<<std::endl;


	    comp<<"# * COMPONENTS * "<<std::endl;
	    comp<<"#"<<std::endl;
	    comp<<"# file generated by ** cup.H ** "<<std::endl;
	    comp<<"#"<<std::endl;
	    comp<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
	    comp<<"# ------------- "<<std::endl;
	    comp<<"# ome (eV)\treal(eps1)\timag(eps1)\treal(eps2)\timag(eps2)"<<std::endl;
	    comp<<"#"<<std::endl;
	    
	    egva<<"# * EIGENVALUES * "<<std::endl;
	    egva<<"#"<<std::endl;  
	    egva<<"# file generated by ** cup.H ** "<<std::endl;
	    egva<<"#"<<std::endl;
	    egva<<"# PARAMETERS: "<<std::endl;
	    egva<<"# G = "<<pG<<std::endl;
	    egva<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
	    egva<<"# ------------- "<<std::endl;
	    egva<<"# ome (eV)\treal(k1)\timag(k1)\treal(k2)\timag(k2)\treal(k3)\timag(k3)"<<std::endl;
	    egva<<"#"<<std::endl;

	    /** STEADY STATE **/
	
	    // normalized variables
	    nv = normalized_variables ();
		tau2 = nv[1];
		gamd = nv[2];
		ome0 = nv[3];
		GG = set_GG(G, tau2);

	    for (omi=0; omi<=omeN; omi++){
		omeeV = omemi + omi*dome;
		if (ns == 0){
		    eps1 = metal(omeeV);
		    eps2 = active(omeeV,eps_b);
		    }
		if (ns == 1){
		    eps1 = active(omeeV,eps_b);
		    eps2 = metal(omeeV);
		    eps3 = eps_s;
		    }
		alph = polarizability(eps1, eps2, eps3, rho);
		
		valph.push_back(std::make_pair(omeeV, alph));
		
		stat<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<omeeV<<     // 1 ome
		      "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<real(alph)<< // 2 Re(alpha_R)
		      "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<imag(alph)<< // 3 Im(alpha_R)
		      std::endl;

		comp<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<omeeV<<     // 1 ome
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<real(eps1)<< // 2 Re(eps1)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<imag(eps1)<< // 3 Im(eps1)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<real(eps2)<< // 4 Re(eps2)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<imag(eps2)<< // 5 Im(eps2)
		    std::endl;

		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		
		ome = omeeV/Ome_p;
		
		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);

		OmeH = odv[0];
		OmeP = odv[1];
		GamP = odv[2];

		A    = coefficients(OmeH, OmeP, GG, GamP, p0, p1, p2);
		detA = determinant(A);
		kap  = eigenvalues(A);
		egva<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<omeeV<<              // 1 ome
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*real(kap[0])<<  // 2 Re(kap1)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*imag(kap[0])<<  // 3 Im(kap1)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*real(kap[1])<<  // 4 Re(kap2)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*imag(kap[1])<<  // 5 Im(kap2)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*real(kap[2])<<  // 6 Re(kap3)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*imag(kap[2])<<  // 7 Im(kap3)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<real(detA)<<	  // 8 real(detA)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<imag(detA)<<	  // 9 imag(detA)
		std::endl;
		}
	    return valph;
	    }

	double nfindOme(int nfft, std::vector<std::complex<double>> wave, double dt){
	    double Ome, rOme, iOme, iMax=0, rMax=0, nom;
	    double data[2*nfft];
		std::vector<std::pair<double,std::complex<double>>> vfour;
		gsl_fft_complex_wavetable * wavetable;
		gsl_fft_complex_workspace * workspace;

		wavetable = gsl_fft_complex_wavetable_alloc (nfft);
		workspace = gsl_fft_complex_workspace_alloc (nfft);
	    std::ofstream  wve("../data/output/wave.dat");
	    std::ofstream  frr("../data/output/fourier.dat");
    	    std::ofstream  flg("../data/output/fourier.log");
	    flg<<nfft<<" "<<dt<<" "<<Ome_p<<std::endl;
		for (int i=0; i<nfft; i++) {
			REAL(data,i) = wave[i].real();
			IMAG(data,i) = wave[i].imag();
			wve<<REAL(data,i)<<" "<<IMAG(data,i)<<std::endl;
			}
		gsl_fft_complex_forward(data, 1, nfft, wavetable, workspace);

		for (int i=0; i<nfft; i++) {
			nom=2.*pi*i/(nfft*dt);
			if (nom>1.)  nom=nom-2*pi;
			if (nom<-1.) nom=nom+2*pi;
			vfour.push_back(std::make_pair(nom*Ome_p,std::complex<double>(REAL(data,i), IMAG(data,i))));
// 			frr<<nom<<" "<<REAL(data,i)<<" "<<IMAG(data,i)<<std::endl;
			if (rMax<fabs(REAL(data,i))){
				rMax=fabs(REAL(data,i));
				rOme=i;
				}
			if (iMax<fabs(IMAG(data,i))){
				iMax=fabs(IMAG(data,i));
				iOme=i;
				}
			}
		Ome=pi*(rOme+iOme)/(nfft*dt);
		if (Ome>1.)  Ome=Ome-2*pi;
		if (Ome<-1.) Ome=Ome+2*pi;
		Ome=Ome_p*Ome;
		gsl_fft_complex_wavetable_free (wavetable);
		gsl_fft_complex_workspace_free (workspace);
		sort(vfour.begin(), vfour.end(), 
              [](const std::pair<double, std::complex<double>>& a, const std::pair<double, std::complex<double>>& b) {
                  return a.first < b.first;
              });
		for (int i=0; i<nfft; i++)
			frr<<vfour[i].first<<" "<<real(vfour[i].second)<<" "<<imag(vfour[i].second)<<std::endl;
	    return Ome;
	    }

// 	std::complex<double> omegatest(char* mdl , char* mtl ,char* hst, double nEsat, double omeeV, double tpump, double E0, std::fstream& file, char* sol=0, double rho=0){

// 		//Variables that once defined do not change
// 		double eps_b, eps_s, ome, dt;
// 		// std::vector<std::complex<double>> wave;
// 		int ns=0;
// 		//int steste=1; //It does change, but it is not required to be on the loop

// 		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
// 		double omep, t, tau1, tau2, gamd, ome0, *nv, tildeN, pG = fabs(G);
// 		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP, **A, *B, *p0, *p1, *p2, *p3, *p, dip, steady_dip, po, *odv;
		
// 		//double omeR;

// 			if ((sol != NULL)) {
// 		    eps_s=set_host(sol);
// 		    ns=1;
// 		    }
// 		    else eps_s=0;
// 		eps_b=set_host(hst);

// 		A = new std::complex<double>*[3];
// 		for(int i = 0; i <= 2; i++)
// 		    A[i] = new std::complex<double>[3];
		
// 		B = new std::complex<double>[3];


// 		nv = new double[4];

// 		odv = new std::complex<double> [3];

// 		// std::fstream dyna;
// 		// dyna.open("../data/output/probetime.dat", std::ios::out);
		
		
// 		// dyna<<"# * PROBETEST * "<<std::endl;
// 		// dyna<<"#"<<std::endl;
// 		// dyna<<"# file generated by ** cup.H ** "<<std::endl;
// 		// dyna<<"#"<<std::endl;
// 		// dyna<<"# PARAMETERS: "<<std::endl;
// 		// dyna<<"# G = "<<pG<<std::endl;
// 		// dyna<<"# ome = "<<omeeV<<" eV"<<std::endl;
// 		// dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
// 		// dyna<<"# ------------- "<<std::endl;
// 		// dyna<<"# t (ps)\t E0"<<std::endl;
// 		// dyna<<"#"<<std::endl;

// 		omep = eV2j*Ome_p*2.*M_PI/h;
// 		// normalized variables
// 		nv = normalized_variables ();
// 		tau1 = nv[0];
// 		tau2 = nv[1];
// 		gamd = nv[2];
// 		ome0 = nv[3];
// 		GG = set_GG(G, tau2);

		
// 		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
// 		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
// 		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
// 		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
// 		if (ns == 1) p=p3;
// 		    else p=p0;

// 		ome = omeeV/Ome_p;

// 		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
		
// 		OmeH = odv[0];
// 		OmeP = odv[1];
// 		GamP = odv[2];

// 		chiH = -GG/OmeH;
// 		chiP = -GamP/OmeP;
		
// 		eps1e = eps_inf + chiP;
// 		eps2e = eps_b   + chiH;

// 		double T=100000*omep*1.e-12; // in ome_p

// 		// Initial conditions
// 		dt  = 1.;
// 		//double t_stop = 0.0;

// 		int Nt=T/dt;
		
// 		std::complex<double>* amplitudes = simple_intensity_steady_state(mdl, mtl, hst, omeeV, sol, rho, E0);

// 		//Variables that must be initialized once the for loop starts
// 		// double dnormem=1.e20, dnavr=0, dnavrm=1.e20;
// 		// int wb=0.;
// 		// double dpTOT=0, t0=0.0;
// 		// int AVE=0;

// 		double N=0;
// 		std::complex<double>* q = new std::complex<double>[3]{0.0, 0.0, 0.0};
// 		std::complex<double>* Q = new std::complex<double>[3]{amplitudes[0], amplitudes[1], amplitudes[2]};

// 		steady_dip = numerical_output(E0, Q, p);
// 		// std::cout << std::endl;
// 		// std::cout << "E0 = " << E0 << " Steady Dip = "<< norm(steady_dip);
// 		// std::cout << std::endl;
// 		// std::cout << "E0 = " << E0 << " P3 Dip = "<< norm(amplitudes[3]);
// 		// std::cout << std::endl;

// 		// let's begin without pump
// 		tildeN = 0;
		
// 		//omeR = Rabi_frec(G, tau2, E0);
// 		int i=0, pip=0;;
// 		//double Tperiod = 50;
// 		//double duty = 0.5;
// 		double t_ps;


// 		while (i<=Nt){// && (fabs(norm(dip)-dnormem)> 1.e-20 || t*1.e+12/omep<12*tpump)){
// 			t=i*dt;
// 			i++;

// 			//Introduccion de una onda cuadrada para el bombeo
// 			t_ps    = t * 1.e+12 / omep;  // tiempo actual en ps
// 			dflt(&tildeN, t_ps, tpump, &pip);
// 			//squarewave(&tildeN, Tperiod, duty, t_ps,tpump, &GG, G, tau2);
// 			//sinusoidalwave(&tildeN, Tperiod, t_ps, &GG, G, tau2);

// 			if (ns == 1) {
// 				po  = numerical_output(E0, q, p0);
// 				N   = Runge_Kutta_mono_4 (N, -1./tau1, tildeN/tau1+imag(q[0]*conj(po)), dt);
// 				} else N = 1;
// 			A   = coefficients(OmeH, OmeP, GG*N, GamP, p0, p1, p2); 
// 			B   = inhomogeneous(GG*N, GamP, E0, p0, p1, p2);
// 			q   = Runge_Kutta_4 (q, A, B, dt);
			
// 			for(int ix = 0; ix <= 2; ix++) delete [] A[ix];
// 			delete [] A;
// 			delete [] B;
			
// 			//dnormem = steste*norm(dip);
// 			dip = numerical_output(E0, q, p);
// 			alph = dip/E0;

// 			if (t*1.e+12/omep>2.*tpump){
// 				if(norm(dip)>=0.9*norm(steady_dip))
// 				{
// 					// std::cout << std::endl;
// 					// std::cout << "E0 = " << E0 << " Dip = " << norm(dip) <<" t = " << t_ps;
// 					// std::cout << std::endl;
// 					break;
// 				}
// // 				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-30 && steste==1){ 
// // 					//t_stop = t_ps;
// // 					std::cout << std::endl;
// // 					//std::cout << " 1st " << "E0 = " << E0 << " " << "t = " << t_ps;
// // 					std::cout << " 1st " << "E0 = " << E0 << " Dip = " << norm(dip) <<" t = " << t_ps;
// // 					std::cout << std::endl;
// // 					break;
// // 				}
// // 				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-10) wb=1;
// // // 				if (wb==1 && i%iprint==0){
// // 				if (wb==1){
// // 					wave.push_back(q[0]);
// // 					}
// 				// if (AVE==0){
// 				// 	t0=t;
// 				// 	AVE=1;
// 				// 	}
// 				// dpTOT = dpTOT+norm(dip);
// 				// if (t - t0 >= 1000/omeR) {
// 				// 	dnavrm = dnavr;
// 				// 	dnavr  = dpTOT / (t - t0);
// 				// 	AVE = 0; dpTOT = 0.0;

// 					// double denom2 = std::max(dnavr, 1e-300);
// 					// if (fabs((dnavr - dnavrm) / denom2) < 1e-7 && steste==1) {
// 					// 	std::cout << std::endl;
// 					// 	//std::cout << " 2ndst " << "E0 = " << E0 << " " << "t = " << t_ps;
// 					// 	std::cout << " 2nd " << "E0 = " << E0 << " Dip = " << norm(dip) <<" t = " << t_ps;
// 					// 	std::cout << std::endl;
// 					// 	break;
// 					// }
// 				// }
// 				}
// 			//if (wave.size()>=150000) wave.erase(wave.begin(), wave.begin() +  50000);
// 			}

// 		file<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<< t_ps << // 1 time (ps)
// 			"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<omeeV<<     // 2 E0
// 			std::endl;

// 		// if (ns == 1){
// 		// 	if (wave.size()>0) nOme=nfindOme(wave.size(), wave, dt);
// 		// 		else nOme=0.;
// 		// 	}

// 		// wave.clear();
// 		delete [] q;

		
// // 		std::cout<<"nOme = "<<nOme<<std::endl;
// 		//wave.clear();
// 		delete [] p0;
// 		delete [] p1;
// 		delete [] p2;
// 		delete [] p3;

// 		return dip;
// 		}	

    // std::complex<double> enhancednumerical(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double tpump, char* sol=0, double rho=0)
	// {

	// 	double eps_b, eps_s, ome, dt;
	// 	int ns = 0;
	// 	int steste = 0;   // para el criterio viejo

	// 	std::complex<double> eps1e, eps2e, eps1, eps2, alph;
	// 	double omep, t, tau1, tau2, gamd, ome0, *nv, tildeN;
	// 	std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP;
	// 	std::complex<double> **A, *B, *p0, *p1, *p2, *p3, *p, dip, po, *odv;

	// 	if (sol != NULL) {
	// 		eps_s = set_host(sol);
	// 		ns    = 1;
	// 	} else {
	// 		eps_s = 0;
	// 	}
	// 	eps_b = set_host(hst);

	// 	omep = eV2j * Ome_p * 2. * M_PI / h;

	// 	nv   = normalized_variables();
	// 	tau1 = nv[0];
	// 	tau2 = nv[1];
	// 	gamd = nv[2];
	// 	ome0 = nv[3];
	// 	GG   = set_GG(G, tau2);

	// 	p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
	// 	p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
	// 	p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
	// 	p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
	// 	p  = (ns == 1 ? p3 : p0);

    // 	std::fstream dyna;

    //     dyna.open("../data/output/numtime.dat", std::ios::out);
		
	// 	dyna<<"# * NUMERICAL * "<<std::endl;
	// 	dyna<<"#"<<std::endl;
	// 	dyna<<"# file generated by ** cup.H ** "<<std::endl;
	// 	dyna<<"#"<<std::endl;
	// 	dyna<<"# PARAMETERS: "<<std::endl;
	// 	dyna<<"# G = "<<G<<std::endl;
	// 	dyna<<"# ome = "<<omeeV<<" eV"<<std::endl;
	// 	// dyna<<"# BigOme = " << BigOme << " ,ome_em = " << omeeV - BigOme*Ome_p << " ,omeeV/Ome_p = " << omeeV/Ome_p <<std::endl;
	// 	dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
	// 	dyna<<"# ------------- "<<std::endl;
	// 	dyna<<"# t (ps)\treal(alph)\timag(alph)\treal(dip)\timag(dip)"<<std::endl;
	// 	dyna<<"#"<<std::endl;


	// 	ome  = omeeV / Ome_p;
	// 	odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
	// 	OmeH = odv[0];
	// 	OmeP = odv[1];
	// 	GamP = odv[2];

	// 	chiH  = -GG / OmeH;
	// 	chiP  = -GamP / OmeP;
	// 	eps1e = eps_inf + chiP;
	// 	eps2e = eps_b   + chiH;

	// 	double T = 147000 * omep * 1.e-12; // en unidades de ome_p
	// 	dt  = 1.0;
	// 	int Nt = (int)std::ceil(T / dt);

	// 	int Nprint = (Nt > 50000 ? 50000 : Nt);
	// 	if (Nprint < 1) Nprint = 1;
	// 	int iprint = Nt / Nprint;
	// 	if (iprint < 1) iprint = 1;


	// 	double N = 0.0;
	// 	std::complex<double>* q = new std::complex<double>[3]{0.0, 0.0, 0.0};

	// 	// Para el criterio de convergencia viejo
	// 	double dnormem = 1.e20, dnavr = 0.0, dnavrm = 1.e20;

	// 	double dpTOT   = 0.0, t0 = 0.0;
	// 	int    AVE     = 0;
	// 	double omeR    = Rabi_frec(G, tau2, E0);

	// 	tildeN = 0.0;
	// 	int i   = 0;
	// 	int pip = 0;
	// 	double t_ps  = 0.0;
	// 	double dip_ss = 0.0;   // aquí guardaremos el dipolo estacionario
	// 	double tol = 1e-30;   // o 1e-12 si quieres más estricto
	// 	double newdt = 2*dt;

	// 	while (i <= Nt) {
	// 		t = i*newdt;
	// 		i++;

	// 		t_ps = t * 1.e+12 / omep;
	// 		dflt(&tildeN, t_ps, tpump, &pip);

	// 		if (ns == 1) {
	// 			po = numerical_output(E0, q, p0);
	// 			N  = Runge_Kutta_mono_4(
	// 					N,
	// 					-1./tau1,
	// 					tildeN/tau1 + imag(q[0]*conj(po)),
	// 					newdt);
	// 		} else {
	// 			N = 1;
	// 		}

	// 		A = coefficients(OmeH, OmeP, GG*N, GamP, p0, p1, p2);
	// 		B = inhomogeneous(GG*N, GamP, E0, p0, p1, p2);
	// 		q = Runge_Kutta_4(q, A, B, newdt);

	// 		for (int ix = 0; ix < 3; ++ix)
	// 			delete [] A[ix];
	// 		delete [] A;
	// 		delete [] B;

	// 		// --- criterio viejo: usar variación relativa de |dip| ---
	// 		dnormem = steste * std::norm(dip);  // |dip|^2 previo
	// 		dip     = numerical_output(E0, q, p);
	// 		alph    = dip / E0;
	// 		dip_ss = std::sqrt(std::norm(dip));
	// 		if(i%iprint==0){
	// 			dyna<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<< // 1 time (ps)
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(alph)<<     // 2 Re(alpha)
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(alph)<<     // 3 Im(alpha)
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(dip)<<     // 4 Re(dip)
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(dip)<<     // 5 Im(dip)
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<dip_ss<<             // 6 dip_ss
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<N<<             // 6 N
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(po)<<      // 7 Re(p0)
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(po)<<      // 8 Im(p0) 
	// 				std::endl;
	// 			}

	// 		if (t * 1.e+12 / omep > 2. * tpump) {

	// 			double nd = std::norm(dip);
	// 			if (nd > 0.0) {
	// 				if (std::fabs((nd - dnormem) / nd) < tol && steste==1) break;
	// 			}

	// 			if (AVE == 0) {
	// 				t0  = t;
	// 				AVE = 1;
	// 			}

	// 			dpTOT += norm(dip);

	// 			if (t - t0 >= 1000 / omeR) {
	// 				dnavrm = dnavr;
	// 				dnavr  = dpTOT / (t - t0);
	// 				AVE    = 0;
	// 				dpTOT  = 0.0;

	// 				if ( fabs((dnavr - dnavrm) / dnavr) < 1.e-7
	// 					&& steste == 1 ) {
	// 					// std::cout << std::endl;
	// 					// std::cout << " 2nd omeeV = " << omeeV
	// 					// 		<< " Dip = " << dip_ss
	// 					// 		<< " t = "   << t_ps;
	// 					// std::cout << std::endl;
	// 					break;
	// 				}
	// 			}

	// 			// if (i >= Nt_guard) {
	// 			// 	std::cout << "WARN: no convergio para omeeV=" << omeeV
	// 			// 			<< " en T=" << t_ps << " ps\n";
	// 			// 	std::cout << std::endl;
	// 			// 		std::cout << " WARN omeeV = " << omeeV
	// 			// 				<< " Dip = " << dip_ss
	// 			// 				<< " t = "   << t_ps;
	// 			// 		std::cout << std::endl;
	// 			// 	break;
	// 			// }

	// 		}
	// 	} // fin while

	// 	delete [] q;
	// 	delete [] odv;
	// 	delete [] p0;
	// 	delete [] p1;
	// 	delete [] p2;
	// 	delete [] p3;
	// 	delete [] nv;


	// 	return dip;
	// }

	std::complex<double> omegatest(
		char* mdl , char* mtl ,char* hst,
		double nEsat, double tpump,
		double E0, double ome1, double ome2,
		char* sol, double rho)
	{

		double eps_b, eps_s, ome, dt;
		int ns = 0;
		int steste = 1;   // para el criterio viejo

		std::complex<double> eps1e, eps2e, eps1, eps2;
		double omep, t, tau1, tau2, gamd, ome0, *nv, tildeN;
		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP;
		std::complex<double> **A, *B, *p0, *p1, *p2, *p3, *p, dip, po, *odv;

		if (sol != NULL) {
			eps_s = set_host(sol);
			ns    = 1;
		} else {
			eps_s = 0;
		}
		eps_b = set_host(hst);

		omep = eV2j * Ome_p * 2. * M_PI / h;

		nv   = normalized_variables();
		tau1 = nv[0];
		tau2 = nv[1];
		gamd = nv[2];
		ome0 = nv[3];
		GG   = set_GG(G, tau2);

		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
		p  = (ns == 1 ? p3 : p0);

		// Número de frecuencias en el barrido:
		int Nfreq = 6;
		double omeeV = 0;

		// Paso en frecuencia (eV)
		double center = ome_0;         // centro deseado
        //std::cout << "omemi = " << ome1*Ome_p << ", omemax = " << ome2*Ome_p << std::endl;

		//double fix = ome1*Ome_p + 0.04;
		//double space = fabs(center-fix)/3;
		double omemin = center - 0.01;
		double omemax = center + 0.01;


		double dome = (omemax - omemin) / (Nfreq - 1);


    	std::fstream dyna;

        dyna.open("../data/output/numtime.dat", std::ios::out);
		

		dyna<<"# * FREQUENCY TEST * "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# file generated by ** cup.H ** "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# PARAMETERS: "<<std::endl;
		dyna<<"# G = "<< fabs(G) <<std::endl;
		dyna<<"# E0 = "<<E0<<std::endl;
        dyna<<"# Esat = "<<nEsat<<std::endl;
		dyna<<"# omemin =" << ome1*Ome_p<<std::endl;
		dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		dyna<<"# ------------- "<<std::endl;
		dyna<<"# omeeV\t Dip_max"<<std::endl;
		dyna<<"#"<<std::endl;

		for (int k = 0; k < Nfreq; ++k) {
        	omeeV = omemin + k * dome;
			//omeeV = omemax - k * dome;

			ome  = omeeV / Ome_p;
			odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
			OmeH = odv[0];
			OmeP = odv[1];
			GamP = odv[2];

			chiH  = -GG / OmeH;
			chiP  = -GamP / OmeP;
			eps1e = eps_inf + chiP;
			eps2e = eps_b   + chiH;

			double T = 10000 * omep * 1.e-12; // en unidades de ome_p
			dt  = 1.;
			int Nt = T / dt;

			double N = 0.0;
			std::complex<double>* q = new std::complex<double>[3]{0.0, 0.0, 0.0};

			// Para el criterio de convergencia viejo
			double dnormem = 1.e20, dnavr = 0.0, dnavrm = 1.e20;

			double dpTOT   = 0.0, t0 = 0.0;
			int    AVE     = 0;
			double omeR    = Rabi_frec(G, tau2, E0);

			tildeN = 0.0;
			int i   = 0;
			int pip = 0;
			double t_ps  = 0.0;
			double dip_ss = 0.0;   // aquí guardaremos el dipolo estacionario

			const int Nt_guard = (int)(1*Nt); // ejemplo: 20% extra, tú decides

			while (i <= Nt) {
				t = i * dt;
                i++;

				t_ps = t * 1.e+12 / omep;
				dflt(&tildeN, t_ps, tpump, &pip);

				if (ns == 1) {
					po = numerical_output(E0, q, p0);
					N  = Runge_Kutta_mono_4(
							N,
							-1./tau1,
							tildeN/tau1 + imag(q[0]*conj(po)),
							dt);
				} else {
					N = 1;
				}

				A = coefficients(OmeH, OmeP, GG*N, GamP, p0, p1, p2);
				B = inhomogeneous(GG*N, GamP, E0, p0, p1, p2);
				q = Runge_Kutta_4(q, A, B, dt);

				for (int ix = 0; ix < 3; ++ix)
					delete [] A[ix];
				delete [] A;
				delete [] B;


				dnormem = steste * std::norm(dip);  // |dip|^2 previo
				dip     = numerical_output(E0, q, p);
				// --- criterio viejo: usar variación relativa de |dip| ---
				dip_ss = std::sqrt(std::norm(dip));

				if (t * 1.e+12 / omep > 2. * tpump) {

					// 1) cambio relativo muy pequeño: estado estacionario alcanzado
					if ( fabs((norm(dip) - dnormem) / norm(dip)) < 1.e-30
						&& steste == 1 ) {
						std::cout << std::endl;
						std::cout << " 1st  omeeV = " << omeeV
								<< " Dip = " << dip_ss
								<< " t = "   << t_ps;
						std::cout << std::endl;
						break;
					}

					if (AVE == 0) {
						t0  = t;
						AVE = 1;
					}

					dpTOT += norm(dip);

					if (t - t0 >= 1000 / omeR) {
						dnavrm = dnavr;
						dnavr  = dpTOT / (t - t0);
						AVE    = 0;
						dpTOT  = 0.0;

						if ( fabs((dnavr - dnavrm) / dnavr) < 1.e-7
							&& steste == 1 ) {
							std::cout << std::endl;
							std::cout << " 2nd omeeV = " << omeeV
									<< " Dip = " << dip_ss
									<< " t = "   << t_ps;
							std::cout << std::endl;
							break;
						}
					}

					if (i >= Nt_guard) {
						std::cout << "WARN: no convergio para omeeV=" << omeeV
								<< " en T=" << t_ps << " ps\n";
						std::cout << std::endl;
							std::cout << " WARN omeeV = " << omeeV
									<< " Dip = " << dip_ss
									<< " t = "   << t_ps;
							std::cout << std::endl;
						break;
					}

				}
			} // fin while
			dyna << "  " << std::setw(13) << std::setiosflags(std::ios::left) << omeeV
			<< "\t" << std::setw(13) << std::setiosflags(std::ios::left) << dip_ss
			<< std::endl;
			delete [] q;
			delete [] odv;
		}

		delete [] p0;
		delete [] p1;
		delete [] p2;
		delete [] p3;
		delete [] nv;


		return dip;
	}

		std::complex<double> omegatest_N(
		char* mdl , char* mtl ,char* hst,
		double nEsat, double tpump,
		double E0, double ome1, double ome2,
		char* sol, double rho)
	{
		double eps_b, eps_s, ome, dt;
		int ns = 0;

		std::complex<double> eps1e, eps2e, eps1, eps2;
		double Bigome, ome_em, omep, t, tau1, tau2, gamd, ome0, *nv, tildeN;
		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP;
		std::complex<double> **A, *B, *p0, *p1, *p2, *p3, *p, dip, po, *odv;

		if (sol != NULL) {
			eps_s = set_host(sol);
			ns = 1;
		} else {
			eps_s = 0;
		}
		eps_b = set_host(hst);

		omep = eV2j*Ome_p*2.*M_PI/h;

		nv   = normalized_variables();
		tau1 = nv[0];
		tau2 = nv[1];
		gamd = nv[2];
		ome0 = nv[3];
		GG   = set_GG(G, tau2);

		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
		p  = (ns == 1 ? p3 : p0);

		double omeeV =0, pss=0;
		std::vector<double> omeeV_list, pss_list;

		std::ifstream steady("../data/output/numtime.dat");
		std::string line;
		while (std::getline(steady, line)) {
			if (line.size() == 0 || line[0] == '#') continue; // saltar comentarios
			std::stringstream ss(line);
			ss >> omeeV >> pss;
			omeeV_list.push_back(omeeV);
			pss_list.push_back(pss);
		}
		steady.close();

		std::fstream dyna;
        dyna.open("../data/output/probetime.dat", std::ios::out);
		

		dyna<<"# * FREQUENCY TEST * "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# file generated by ** cup.H ** "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# PARAMETERS: "<<std::endl;
		dyna<<"# G = "<< fabs(G) <<std::endl;
		dyna<<"# E0 = "<<E0<<std::endl;
        dyna<<"# Esat = "<<nEsat<<std::endl;
		dyna<<"# omemin = " << ome1*Ome_p << ", omemax =" << ome2*Ome_p<<std::endl;
		dyna<<"# ------------- "<<std::endl;;
		dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		dyna<<"# ------------- "<<std::endl;
		dyna<<"# t (ps)\t omeeV \t omeeV/Ome_p (norm) \t Bigome (norm) \t ome_em (eV) \t Dip_max"<<std::endl;
		dyna<<"#"<<std::endl;

		for (size_t k = 0; k < omeeV_list.size(); ++k) {

			omeeV = omeeV_list[k];
			pss = pss_list[k];

			ome = omeeV/Ome_p;

			odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
			OmeH = odv[0];
			OmeP = odv[1];
			GamP = odv[2];

			chiH  = -GG/OmeH;
			chiP  = -GamP/OmeP;
			eps1e = eps_inf + chiP;
			eps2e = eps_b   + chiH;

			double T = 10000 * omep * 1.e-12; // en unidades de ome_p
			dt  = 1.;
			int Nt = T/dt;

			double N = 0;
			std::complex<double>* q = new std::complex<double>[3]{0.0, 0.0, 0.0};

			tildeN = 0.0;
			int i   = 0;
			int pip = 0;
			double t_ps  = 0.0;

			Bigome = find_Omega(omeeV, hst, sol, rho);
			ome_em = omeeV - Bigome*Ome_p;

			const int Nt_guard = (int)(1*Nt); // ejemplo: 20% extra, tú decides
			
			while (i <= Nt) {
				t = i * dt;
				i++;

				t_ps = t * 1.e+12 / omep;
				dflt(&tildeN, t_ps, tpump, &pip);

				if (ns == 1) {
					po = numerical_output(E0, q, p0);
					N  = Runge_Kutta_mono_4(N, -1./tau1, tildeN/tau1 + imag(q[0]*conj(po)), dt);
				} else {
					N = 1;
				}

				A = coefficients(OmeH, OmeP, GG*N, GamP, p0, p1, p2);
				B = inhomogeneous(GG*N, GamP, E0, p0, p1, p2);
				q = Runge_Kutta_4(q, A, B, dt);

				for(int ix = 0; ix < 3; ++ix) delete [] A[ix];
				delete [] A;
				delete [] B;

				dip  = numerical_output(E0, q, p);
				double mod_dip = std::sqrt(std::norm(dip));

				//Numero complejo
				
				if (t*1.e+12/omep > 2.*tpump) {
					if (mod_dip >= 0.9*pss) {
						std::cout << std::endl;
						std::cout << "omeeV = " << omeeV << " Dip = " << pss << " 0.9Dip = "<< mod_dip <<" t = " << t_ps;
						std::cout << std::endl;
						break;
					}
				}

				if (i >= Nt_guard) {
						std::cout << "WARN: no convergio para omeeV=" << omeeV
								<< " en T=" << t_ps << " ps\n";
						std::cout << std::endl;
							std::cout << " WARN omeeV = " << omeeV
									<< " Error = " << fabs(0.9*pss-mod_dip)/0.9*pss
									<< " t = "   << t_ps;
							std::cout << std::endl;
						break;
				}
			}

			dyna << "  " << std::setw(8)  << std::setiosflags(std::ios::left) << t_ps
			<< "\t" << std::setw(13) << std::setiosflags(std::ios::left) << omeeV
			<< "\t" << std::setw(13) << std::setiosflags(std::ios::left) << ome
			<< "\t" << std::setw(13) << std::setiosflags(std::ios::left) << Bigome
			<< "\t" << std::setw(13) << std::setiosflags(std::ios::left) << ome_em
			<< "\t" << std::setw(13) << std::setiosflags(std::ios::left) << pss
			<< std::endl;

			delete [] q;
			delete [] odv;
		}

		delete [] p0;
		delete [] p1;
		delete [] p2;
		delete [] p3;
		delete [] nv;


		return dip;
	}


//      std::complex<double> probetest(char* mdl , char* mtl ,char* hst, double nEsat, double omeeV, double tpump, char* sol=0, double rho=0){

// 		//Variables that once defined do not change
// 		double eps_b, eps_s, ome, dt;
// 		// std::vector<std::complex<double>> wave;
// 		int ns=0;
// 		//int steste=1; //It does change, but it is not required to be on the loop

// 		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
// 		double omep, t, tau1, tau2, gamd, ome0, *nv, tildeN, pG = fabs(G);
// 		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP, **A, *B, *p0, *p1, *p2, *p3, *p, dip, steady_dip, po, *odv;
		
// 		//double omeR;

// 			if ((sol != NULL)) {
// 		    eps_s=set_host(sol);
// 		    ns=1;
// 		    }
// 		    else eps_s=0;
// 		eps_b=set_host(hst);

// 		A = new std::complex<double>*[3];
// 		for(int i = 0; i <= 2; i++)
// 		    A[i] = new std::complex<double>[3];
		
// 		B = new std::complex<double>[3];


// 		nv = new double[4];

// 		odv = new std::complex<double> [3];

// 		std::fstream dyna;
// 		dyna.open("../data/output/probetime.dat", std::ios::out);
		
		
// 		dyna<<"# * PROBETEST * "<<std::endl;
// 		dyna<<"#"<<std::endl;
// 		dyna<<"# file generated by ** cup.H ** "<<std::endl;
// 		dyna<<"#"<<std::endl;
// 		dyna<<"# PARAMETERS: "<<std::endl;
// 		dyna<<"# G = "<<pG<<std::endl;
// 		dyna<<"# ome = "<<omeeV<<" eV"<<std::endl;
// 		dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
// 		dyna<<"# ------------- "<<std::endl;
// 		dyna<<"# t (ps)\t E0"<<std::endl;
// 		dyna<<"#"<<std::endl;

// 		omep = eV2j*Ome_p*2.*M_PI/h;
// 		// normalized variables
// 		nv = normalized_variables ();
// 		tau1 = nv[0];
// 		tau2 = nv[1];
// 		gamd = nv[2];
// 		ome0 = nv[3];
// 		GG = set_GG(G, tau2);

		
// 		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
// 		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
// 		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
// 		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
// 		if (ns == 1) p=p3;
// 		    else p=p0;

// 		ome = omeeV/Ome_p;

// 		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
		
// 		OmeH = odv[0];
// 		OmeP = odv[1];
// 		GamP = odv[2];

// 		chiH = -GG/OmeH;
// 		chiP = -GamP/OmeP;
		
// 		eps1e = eps_inf + chiP;
// 		eps2e = eps_b   + chiH;

// 		double T=100000*omep*1.e-12; // in ome_p

// 		// Initial conditions
// 		dt  = 1.;
// 		//double t_stop = 0.0;

// 		int Nt=T/dt;

// 		//Introduction of the for-loop
// 		const int Npts = 100;
// 		const double Emin = 1e-15 * nEsat;
// 		const double Emax = 1e-2  * nEsat;
// 		const double ratio = Emax / Emin;

// 		for (int k = 0; k < Npts; ++k) {
// 		    double frac = double(k) / (Npts - 1);
// 		    double E0 = Emin * std::pow(ratio, frac);

// 			std::complex<double>* amplitudes = simple_intensity_steady_state(mdl, mtl, hst, omeeV, sol, rho, E0);

// 			//Variables that must be initialized once the for loop starts
// 			//double dnormem=1.e20, dnavr=0, dnavrm=1.e20;
// 			//int wb=0.;
// 			//double dpTOT=0, t0=0.0;
// 			//int AVE=0;

// 			double N=0;
// 			std::complex<double>* q = new std::complex<double>[3]{0.0, 0.0, 0.0};
// 			std::complex<double>* Q = new std::complex<double>[3]{amplitudes[0], amplitudes[1], amplitudes[2]};

// 			steady_dip = numerical_output(E0, Q, p);
// 			std::cout << std::endl;
// 			std::cout << "E0 = " << E0 << " Steady Dip = "<< norm(steady_dip);
// 			std::cout << std::endl;
// 			std::cout << "E0 = " << E0 << " P3 Dip = "<< norm(amplitudes[3]);
// 			std::cout << std::endl;

// 			// let's begin without pump
// 			tildeN = 0;
			
// 			//omeR = Rabi_frec(G, tau2, E0);
// 			int i=0, pip=0;;
// 			//double Tperiod = 50;
// 			//double duty = 0.5;
//             double t_ps;


// 			while (i<=Nt){// && (fabs(norm(dip)-dnormem)> 1.e-20 || t*1.e+12/omep<12*tpump)){
// 				t=i*dt;
// 				i++;

// 				//Introduccion de una onda cuadrada para el bombeo
// 				t_ps    = t * 1.e+12 / omep;  // tiempo actual en ps
// 				dflt(&tildeN, t_ps, tpump, &pip);
// 				//squarewave(&tildeN, Tperiod, duty, t_ps,tpump, &GG, G, tau2);
// 				//sinusoidalwave(&tildeN, Tperiod, t_ps, &GG, G, tau2);

// 				if (ns == 1) {
// 					po  = numerical_output(E0, q, p0);
// 					N   = Runge_Kutta_mono_4 (N, -1./tau1, tildeN/tau1+imag(q[0]*conj(po)), dt);
// 					} else N = 1;
// 				A   = coefficients(OmeH, OmeP, GG*N, GamP, p0, p1, p2); 
// 				B   = inhomogeneous(GG*N, GamP, E0, p0, p1, p2);
// 				q   = Runge_Kutta_4 (q, A, B, dt);
				
// 				for(int ix = 0; ix <= 2; ix++) delete [] A[ix];
// 				delete [] A;
// 				delete [] B;
				
// 				//dnormem = steste*norm(dip);
// 				dip = numerical_output(E0, q, p);
// 				alph = dip/E0;

// 				if (t*1.e+12/omep>2.*tpump){
// 					if(norm(dip)>=0.9*norm(steady_dip))
// 					{
// 						std::cout << std::endl;
// 						std::cout << "E0 = " << E0 << " Dip = " << norm(dip) <<" t = " << t_ps;
// 						std::cout << std::endl;
// 						break;
// 					}
// 	// 				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-30 && steste==1){ 
//   	// 					//t_stop = t_ps;
// 	// 					std::cout << std::endl;
// 	// 					//std::cout << " 1st " << "E0 = " << E0 << " " << "t = " << t_ps;
// 	// 					std::cout << " 1st " << "E0 = " << E0 << " Dip = " << norm(dip) <<" t = " << t_ps;
// 	// 					std::cout << std::endl;
// 	// 					break;
//  	// 				}
// 	// 				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-10) wb=1;
// 	// // 				if (wb==1 && i%iprint==0){
// 	// 				if (wb==1){
// 	// 					wave.push_back(q[0]);
// 	// 					}
// 					// if (AVE==0){
// 					// 	t0=t;
// 					// 	AVE=1;
// 					// 	}
// 					// dpTOT = dpTOT+norm(dip);
// 					// if (t - t0 >= 1000/omeR) {
// 					// 	dnavrm = dnavr;
// 					// 	dnavr  = dpTOT / (t - t0);
// 					// 	AVE = 0; dpTOT = 0.0;

// 						// double denom2 = std::max(dnavr, 1e-300);
// 						// if (fabs((dnavr - dnavrm) / denom2) < 1e-7 && steste==1) {
// 						// 	std::cout << std::endl;
// 						// 	//std::cout << " 2ndst " << "E0 = " << E0 << " " << "t = " << t_ps;
// 						// 	std::cout << " 2nd " << "E0 = " << E0 << " Dip = " << norm(dip) <<" t = " << t_ps;
// 						// 	std::cout << std::endl;
// 						// 	break;
// 						// }
// 					// }
// 				  }
// 				//if (wave.size()>=150000) wave.erase(wave.begin(), wave.begin() +  50000);
// 				}

// 			dyna<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<< t_ps << // 1 time (ps)
// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<E0<<     // 2 E0
// 				std::endl;

// 			// if (ns == 1){
// 			// 	if (wave.size()>0) nOme=nfindOme(wave.size(), wave, dt);
// 			// 		else nOme=0.;
// 			// 	}

// 			// wave.clear();
// 			delete [] q;
// 		}
		
// // 		std::cout<<"nOme = "<<nOme<<std::endl;
// 		//wave.clear();
// 		delete [] p0;
// 		delete [] p1;
// 		delete [] p2;
// 		delete [] p3;

// 		return dip;
// 		}



	std::complex<double> probetest(char* mdl , char* mtl ,char* hst,
								double nEsat, double omeeV, double tpump,
								char* sol, double rho)
	{
		// Variables que no cambian dentro del bucle en E0
		double eps_b, eps_s, ome, dt;
		std::vector<std::complex<double>> wave;
		int ns = 0;
		int steste = 1;  // necesario para el criterio viejo

		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
		double omep, dip_ss, t, tau1, tau2, gamd, ome0, *nv, tildeN, pG = fabs(G);
		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP;
		std::complex<double> **A, *B, *p0, *p1, *p2, *p3, *p, dip, po, *odv;

		// Host, metal, etc.
		if (sol != NULL) {
			eps_s = set_host(sol);
			ns    = 1;
		} else {
			eps_s = 0;
		}
		eps_b = set_host(hst);

		A = new std::complex<double>*[3];
		for (int i = 0; i <= 2; ++i)
			A[i] = new std::complex<double>[3];

		B   = new std::complex<double>[3];
		nv  = new double[4];
		odv = new std::complex<double>[3];

		// Archivo de salida global del probetest
		std::fstream dyna;
		dyna.open("../data/output/numtime.dat", std::ios::out);

		dyna << "# * PROBETEST * " << std::endl;
		dyna << "#" << std::endl;
		dyna << "# file generated by ** cup.H ** " << std::endl;
		dyna << "#" << std::endl;
		dyna << "# PARAMETERS: " << std::endl;
		dyna << "# G = " << pG << std::endl;
		dyna << "# ome = " << omeeV << " eV" << std::endl;
		dyna << "# a " << mdl << " model was used for " << mtl << std::endl;
		dyna << "# ------------- " << std::endl;
		dyna << "# E0\t Max Dip\t t" << std::endl;
		dyna << "#" << std::endl;

		// Normalización de frecuencias, etc.
		omep = eV2j * Ome_p * 2. * M_PI / h;
		nv   = normalized_variables();
		tau1 = nv[0];
		tau2 = nv[1];
		gamd = nv[2];
		ome0 = nv[3];
		GG   = set_GG(G, tau2);

		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
		if (ns == 1) p = p3; else p = p0;

		ome = omeeV / Ome_p;
		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
		OmeH = odv[0];
		OmeP = odv[1];
		GamP = odv[2];

		chiH  = -GG / OmeH;
		chiP  = -GamP / OmeP;
		eps1e = eps_inf + chiP;
		eps2e = eps_b   + chiH;

		double T = 100000 * omep * 1.e-12; // en unidades de ome_p
		dt = 1.;
		int Nt = T / dt;

		// ---- Barrido en E0 (logarítmico) ----
		const int    Npts  = 500;
		const double Emin  = 1e-15 * nEsat;
		const double Emax  = nEsat;
		const double ratio = Emax / Emin;

		for (int k = 0; k < Npts; ++k) {

			double frac = double(k) / (Npts - 1);
			double E0   = Emin * std::pow(ratio, frac);

			// Variables del criterio de convergencia (viejo) — se reinician por cada E0
			double dnormem = 1.e20, dnavr = 0.0, dnavrm = 1.e20;
			int    wb      = 0;
			double dpTOT   = 0.0, t0 = 0.0;
			int    AVE     = 0;
			double omeR    = Rabi_frec(G, tau2, E0);

			double N = 0.0;
			std::complex<double>* q = new std::complex<double>[3]{0.0, 0.0, 0.0};

			// empezamos sin bombeo
			tildeN = 0.0;
			int    i   = 0;
			int    pip = 0;
			double t_ps = 0.0;

			while (i <= Nt) {
				t = i * dt;
				i++;

				// tiempo en ps
				t_ps = t * 1.e+12 / omep;

				// bombeo
				dflt(&tildeN, t_ps, tpump, &pip);

				if (ns == 1) {
					po = numerical_output(E0, q, p0);
					N  = Runge_Kutta_mono_4(
							N,
							-1./tau1,
							tildeN/tau1 + imag(q[0]*conj(po)),
							dt);
				} else {
					N = 1;
				}

				A = coefficients(OmeH, OmeP, GG*N, GamP, p0, p1, p2);
				B = inhomogeneous(GG*N, GamP, E0, p0, p1, p2);
				q = Runge_Kutta_4(q, A, B, dt);

				for (int ix = 0; ix <= 2; ++ix)
					delete [] A[ix];
				delete [] A;
				delete [] B;

				// ---- criterio viejo: se basa en el cambio relativo de |dip| ----
				dnormem = steste * norm(dip);
				dip     = numerical_output(E0, q, p);
				dip_ss = std::sqrt(std::norm(dip));
				alph    = dip / E0;
				

				if (t * 1.e+12 / omep > 2. * tpump) {

					// 1) cambio relativo muy pequeño: estado estacionario alcanzado
					if ( fabs((norm(dip) - dnormem) / norm(dip)) < 1.e-30
						&& steste == 1 ) {
						std::cout << std::endl;
						std::cout << " 1st  E0 = " << E0
								<< " Dip = " << dip_ss
								<< " t = "   << t_ps;
						std::cout << std::endl;
						break;
					}

					// 2) cuando el cambio relativo cae por debajo de 1e-10 empezamos a acumular para promedio
					if ( fabs((norm(dip) - dnormem) / norm(dip)) < 1.e-10 )
						wb = 1;

					if (wb == 1) {
						wave.push_back(q[0]);
					}

					if (AVE == 0) {
						t0  = t;
						AVE = 1;
					}

					dpTOT += norm(dip);

					if (t - t0 >= 1000 / omeR) {
						dnavrm = dnavr;
						dnavr  = dpTOT / (t - t0);
						AVE    = 0;
						dpTOT  = 0.0;

						if ( fabs((dnavr - dnavrm) / dnavr) < 1.e-7
							&& steste == 1 ) {
							std::cout << std::endl;
							std::cout << " 2nd E0 = " << E0
									<< " Dip = " << dip_ss
									<< " t = "   << t_ps;
							std::cout << std::endl;
							break;
						}
					}
				}

				if (wave.size() >= 150000)
					wave.erase(wave.begin(), wave.begin() + 50000);
			} // fin while sobre t
			

			// escribir (t_final, E0) para este E0
			// << "  " << std::setw(8)  << std::setiosflags(std::ios::left) << t_ps
				dyna << " " << std::setw(13) << std::setiosflags(std::ios::left) << E0
				<< "\t" << std::setw(13) << std::setiosflags(std::ios::left) << dip_ss
				<< "\t" << std::setw(8)  << std::setiosflags(std::ios::left) << t_ps
				<< std::endl;

			wave.clear();
			delete [] q;
		} // fin for sobre E0

		delete [] p0;
		delete [] p1;
		delete [] p2;
		delete [] p3;

		dyna.close();

		return dip;  // valor de dip para el último E0 simulado
	}

			std::complex<double> probetest_N(char* mdl , char* mtl ,char* hst,
								double nEsat, double omeeV, double tpump,
								char* sol, double rho)
	{
		// Variables que no cambian dentro del bucle en E0
		double eps_b, eps_s, ome, dt;
		std::vector<std::complex<double>> wave;
		int ns = 0;

		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
		double Eratio, omep, t, tau1, tau2, gamd, ome0, *nv, tildeN, pG = fabs(G);
		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP;
		std::complex<double> **A, *B, *p0, *p1, *p2, *p3, *p, dip, po, *odv;

		// Host, metal, etc.
		if (sol != NULL) {
			eps_s = set_host(sol);
			ns    = 1;
		} else {
			eps_s = 0;
		}
		eps_b = set_host(hst);

		A = new std::complex<double>*[3];
		for (int i = 0; i <= 2; ++i)
			A[i] = new std::complex<double>[3];

		B   = new std::complex<double>[3];
		nv  = new double[4];
		odv = new std::complex<double>[3];

		//Abrimos el archivo de probetest

		std::vector<double> E0_list, pss_list;

		std::ifstream steady("../data/output/numtime.dat");
		std::string line;
		while (std::getline(steady, line)) {
			if (line.size() == 0 || line[0] == '#') continue; // saltar comentarios
			std::stringstream ss(line);
			double E0, pss;
			ss >> E0 >> pss;
			E0_list.push_back(E0);
			pss_list.push_back(pss);
		}
		steady.close();

		// Archivo de salida global del probetest
		std::fstream dyna;
		dyna.open("../data/output/probetime.dat", std::ios::out);

		dyna << "# * PROBETEST * " << std::endl;
		dyna << "#" << std::endl;
		dyna << "# file generated by ** cup.H ** " << std::endl;
		dyna << "#" << std::endl;
		dyna << "# PARAMETERS: " << std::endl;
		dyna << "# G = " << pG << std::endl;
		dyna << "# ome = " << omeeV << " eV" << std::endl;
		dyna<<  "# Esat = "<<nEsat<<std::endl;
		dyna << "# a " << mdl << " model was used for " << mtl << std::endl;
		dyna << "# ------------- " << std::endl;
		dyna << "# t (ps)"<< std::setw(8) << "E0" << std::setw(13) << "E0/Esat" << std::setw(13) << "Dipmax" << std::endl;
		dyna << "#" << std::endl;


		// Normalización de frecuencias, etc.
		omep = eV2j * Ome_p * 2. * M_PI / h;
		nv   = normalized_variables();
		tau1 = nv[0];
		tau2 = nv[1];
		gamd = nv[2];
		ome0 = nv[3];
		GG   = set_GG(G, tau2);

		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
		if (ns == 1) p = p3; else p = p0;

		ome = omeeV / Ome_p;
		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
		OmeH = odv[0];
		OmeP = odv[1];
		GamP = odv[2];

		chiH  = -GG / OmeH;
		chiP  = -GamP / OmeP;
		eps1e = eps_inf + chiP;
		eps2e = eps_b   + chiH;

		double T = 100000 * omep * 1.e-12; // en unidades de ome_p
		dt = 1.;
		int Nt = T / dt;


		for (size_t k = 0; k < E0_list.size(); ++k) {

			double E0 = E0_list[k];
			double pss = pss_list[k];
			double target = 0.9 * pss;

			double N = 0.0;
			std::complex<double>* q = new std::complex<double>[3]{0.0, 0.0, 0.0};

			// empezamos sin bombeo
			tildeN = 0.0;
			int    i   = 0;
			int    pip = 0;
			double t_ps = 0.0;

			while (i <= Nt) {
				t = i * dt;
				i++;

				// tiempo en ps
				t_ps = t * 1.e+12 / omep;

				// bombeo
				dflt(&tildeN, t_ps, tpump, &pip);

				if (ns == 1) {
					po = numerical_output(E0, q, p0);
					N  = Runge_Kutta_mono_4(
							N,
							-1./tau1,
							tildeN/tau1 + imag(q[0]*conj(po)),
							dt);
				} else {
					N = 1;
				}

				A = coefficients(OmeH, OmeP, GG*N, GamP, p0, p1, p2);
				B = inhomogeneous(GG*N, GamP, E0, p0, p1, p2);
				q = Runge_Kutta_4(q, A, B, dt);

				for (int ix = 0; ix <= 2; ++ix)
					delete [] A[ix];
				delete [] A;
				delete [] B;

				dip = numerical_output(E0, q, p);
				double mod_dip = std::sqrt(std::norm(dip));
				Eratio = E0 / nEsat;
				alph    = dip / E0;

				if (t*1.e+12/omep > 2.0*tpump) {
					if (mod_dip >= target) {
						std::cout << std::endl;
						std::cout << "E0 = " << E0 << " Dip = " << pss << " 0.9Dip = "<< mod_dip <<" t = " << t_ps;
						std::cout << std::endl;
						break;
					}
				}

				if (wave.size() >= 150000)
					wave.erase(wave.begin(), wave.begin() + 50000);
			} // fin while sobre t

			// escribir (t_final, E0) para este E0
			dyna << "  " << std::setw(8)  << std::setiosflags(std::ios::left) << t_ps
				<< "\t" << std::setw(13) << std::setiosflags(std::ios::left) << E0
				<< "\t" << std::setw(13) << std::setiosflags(std::ios::left) << Eratio
				<< "\t" << std::setw(13) << std::setiosflags(std::ios::left) << pss
				<< std::endl;

			wave.clear();
			delete [] q;
		} // fin for sobre E0

		delete [] p0;
		delete [] p1;
		delete [] p2;
		delete [] p3;

		dyna.close();

		return dip;  // valor de dip para el último E0 simulado
	}

		std::complex<double> taucharacterization(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double T, double tpump, char* sol=0, double rho=0){
		double eps_b, eps_s, ome, dt, dnormem=1.e20, dnavr=0, dnavrm=1.e20;
		std::vector<std::complex<double>> wave;
		int ns=0, wb=0.;
		if ((sol != NULL)) {
		    eps_s=set_host(sol);
		    ns=1;
		    }
		    else eps_s=0;
		eps_b=set_host(hst);
		int steste=0;
		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
		double omep, t, tau1, tau2, gamd, ome0, pG = fabs(G), *nv, tildeN;
		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP, **A, *B, *p0, *p1, *p2, *p3, *p, dip, po, *odv;

		double dpTOT=0, t0;
		double omeR;
		int AVE=0;
		
		A = new std::complex<double>*[3];
		for(int i = 0; i <= 2; i++)
		    A[i] = new std::complex<double>[3];
		
		B = new std::complex<double>[3];


		nv = new double[4];

		odv = new std::complex<double> [3];
		
		if (T==0){ 
		    T=100000;
		    steste=1;
		    }

		std::fstream dyna, fnct;
		dyna.open("../data/output/numtime.dat", std::ios::out);
		fnct.open("../data/output/numfunc.dat", std::ios::out);
		//prb.open("../data/output/probetime.dat", std::ios::out);

		// Convert Δ [eV] to τ [ps]: τ = 2ħ / Δ
	    double ntau2 = 2.0 * hbar * j2eV / Dome; // in seconds
	    double ntau2_ps = ntau2 * ps;
		double ntau1_ps = 5.0 * ntau2_ps;

		
		dyna<<"# * NUMERICAL * "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# file generated by ** cup.H ** "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# PARAMETERS: "<<std::endl;
		dyna<<"# G = "<<pG<<std::endl;
		dyna<<"# ome = "<<omeeV<<" eV"<<std::endl;
		dyna<<"# linewidth Δ = "<<Dome<<" eV"<<std::endl;
		dyna<<"# tau 1 = "<<ntau1_ps<<" ps"<<std::endl;
		dyna<<"# tau 2 = "<<ntau2_ps<<" ps"<<std::endl;
		dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		dyna<<"# ------------- "<<std::endl;
		dyna<<"# t (ps)\treal(alph)\timag(alph)\treal(dip)\timag(dip)"<<std::endl;
		dyna<<"#"<<std::endl;

		fnct<<"# * NUMERICAL * "<<std::endl;
		fnct<<"#"<<std::endl;
		fnct<<"# file generated by ** cup.H ** "<<std::endl;
		fnct<<"#"<<std::endl;
		fnct<<"# PARAMETERS: "<<std::endl;
		fnct<<"# G = "<<pG<<std::endl;
		fnct<<"# ome = "<<omeeV<<" eV"<<std::endl;
		fnct<<"# linewidth Δ = "<<Dome<<" eV"<<std::endl;
		fnct<<"# tau 1 = "<<ntau1_ps<<" ps"<<std::endl;
		fnct<<"# tau 2 = "<<ntau2_ps<<" ps"<<std::endl;
		fnct<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		fnct<<"# ------------- "<<std::endl;
		fnct<<"# t (ps)\treal(q0)\timag(q0)\treal(q1)\timag(q1)\treal(q2)\timag(q2)"<<std::endl;
		fnct<<"#"<<std::endl;

		omep = eV2j*Ome_p*2.*M_PI/h;
		// normalized variables
		nv = normalized_variables ();
		tau1 = nv[0];
		tau2 = nv[1];
		gamd = nv[2];
		ome0 = nv[3];
		GG = set_GG(G, tau2);

		
		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
		if (ns == 1) p=p3;
		    else p=p0;

		ome = omeeV/Ome_p;
		
		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
		
		OmeH = odv[0];
		OmeP = odv[1];
		GamP = odv[2];

		chiH = -GG/OmeH;
		chiP = -GamP/OmeP;

		eps1e = eps_inf + chiP;
		eps2e = eps_b   + chiH;


		// NUMERICAL SOLUTIONS
		T=T*omep*1.e-12; // in ome_p
		// Initial conditions
		dt  = 1.;
		
		double N=0;
		std::complex<double> *q ;
		q = new std::complex<double>[3] {0,0};

		int Nt=T/dt;
		int Nprint=Nt;

		if(Nt > 50000) Nprint=50000;    
// 		if(Nt > 500000) Nprint=500000;    
		int iprint=Nt/Nprint;
		GG = set_GG(G, tau2);
		// let's begin without pump
		tildeN = 0;
	    
		omeR = Rabi_frec(G, tau2, E0);
		int i=0;

		// These variables just work on specific functions
		//int pip = 0;
		double Tperiod = 50;
		//double Eperiod;
		//double duty = 0.5;
		//double E0_temp = E0;

		while (i<=Nt){// && (fabs(norm(dip)-dnormem)> 1.e-20 || t*1.e+12/omep<12*tpump)){
		    t=i*dt;
		    i++;

		//Introduccion de una onda cuadrada para el bombeo
		double t_ps    = t * 1.e+12 / omep;  // tiempo actual en ps

		//Introduction of a square wave as the electric field
		//double local_t = fmod(t_ps, Eperiod);
			    //if (local_t < duty * Eperiod) {
			        //E0 = 0;
			    //} else {
			        //E0 = E0_temp;
			    //}
		
		//dflt(&tildeN, t_ps, tpump, &pip);
		//squarewave(&tildeN, Tperiod, duty, t_ps,tpump, &GG, G, tau2);
		sinusoidalwave(&tildeN, Tperiod, t_ps, &GG, G, tau2);

		//Aqui termina el cambio

		    if (ns == 1) {
				po  = numerical_output(E0, q, p0);
				N   = Runge_Kutta_mono_4 (N, -1./tau1, tildeN/tau1+imag(q[0]*conj(po)), dt);
				} else N = 1;
		    A   = coefficients(OmeH, OmeP, GG*N, GamP, p0, p1, p2); 
		    B   = inhomogeneous(GG*N, GamP, E0, p0, p1, p2);
		    q   = Runge_Kutta_4 (q, A, B, dt);
		    
		    for(int ix = 0; ix <= 2; ix++) delete [] A[ix];
		    delete [] A;
		    delete [] B;
			
		    dnormem = steste*norm(dip);
		    dip = numerical_output(E0, q, p);
		    alph = dip/E0;
				if(i%iprint==0){
				dyna<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<< // 1 time (ps)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(alph)<<     // 2 Re(alpha)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(alph)<<     // 3 Im(alpha)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(dip)<<     // 4 Re(dip)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(dip)<<     // 5 Im(dip)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<N<<             // 6 N
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(po)<<      // 7 Re(p0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(po)<<      // 8 Im(p0) 
					std::endl;
				fnct<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<< // 1 time (ps)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[0])<<   // 2 Re(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[0])<<   // 3 Im(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[1])<<   // 4 Re(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[1])<<   // 5 Im(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[2])<<   // 6 Re(q2)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[2])<<   // 7 Im(q2)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<N<<    // 8 N
					std::endl;
				}
			if (t*1.e+12/omep>2.*tpump){
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-30 && steste==1){
					//std::cout << std::endl;
					//std::cout << " 1st " << "E0 = " << E0 << " " << "t = " << t_ps;
					//std::cout << std::endl;
					break;
				}
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-10) wb=1;
// 				if (wb==1 && i%iprint==0){
				if (wb==1){
				    wave.push_back(q[0]);
				    }
				if (fabs((dnavr-dnavrm)/dnavr) < 1.e-7 && steste==1) {
					//std::cout << std::endl;
					//std::cout << " 2nd " << "E0 = " << E0 << " " << "t = " << t_ps;
					//std::cout << std::endl;
					break;}
				if (AVE==0){
					t0=t;
					AVE=1;
					}
				dpTOT = dpTOT+norm(dip);
				if (t-t0 >= 1000/omeR){
					dnavrm=dnavr;
					dnavr=dpTOT/(t-t0);
					AVE=0;
					dpTOT =0;
					}
				}
		    if (wave.size()>=150000) wave.erase(wave.begin(), wave.begin() +  50000);
		    }
		if (ns == 1){
			if (wave.size()>0) nOme=nfindOme(wave.size(), wave, dt);
				else nOme=0.;
			}
// 		std::cout<<"nOme = "<<nOme<<std::endl;
		wave.clear();
		delete [] p0;
		delete [] p1;
		delete [] p2;
		delete [] p3;
		delete [] q;
		return dip;
		}

		 std::complex<double> numericalOEF(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double T, double tpump, char* sol=0, double rho=0){
		double eps_b, eps_s, ome, dt, dnormem=1.e20, dnavr=0, dnavrm=1.e20;
		std::vector<std::complex<double>> wave;
		int ns=0, wb=0.;
		if ((sol != NULL)) {
		    eps_s=set_host(sol);
		    ns=1;
		    }
		    else eps_s=0;
		eps_b=set_host(hst);
		int steste=0;
		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
		double omep, t, tau1, tau2, gamd, ome0, pG = fabs(G), *nv, tildeN;
		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP, **A, *B, *p0, *p1, *p2, *p3, *p, dip, po, *odv;

		double dpTOT=0, t0;
		double omeR;
		int AVE=0;
		
		A = new std::complex<double>*[3];
		for(int i = 0; i <= 2; i++)
		    A[i] = new std::complex<double>[3];
		
		B = new std::complex<double>[3];


		nv = new double[4];

		odv = new std::complex<double> [3];
		
		if (T==0){ 
		    T=100000;
		    steste=1;
		    }

		std::fstream dyna, fnct;
		dyna.open("../data/output/numtime.dat", std::ios::out);
		fnct.open("../data/output/numfunc.dat", std::ios::out);
		//prb.open("../data/output/probetime.dat", std::ios::out);
		
		dyna<<"# * NUMERICAL * "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# file generated by ** cup.H ** "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# PARAMETERS: "<<std::endl;
		dyna<<"# G = "<<pG<<std::endl;
		dyna<<"# ome = "<<omeeV<<" eV"<<std::endl;
		dyna<<"# Electric field amplitude = "<<E0<<std::endl;
		dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		dyna<<"# ------------- "<<std::endl;
		dyna<<"# t (ps)\treal(alph)\timag(alph)\treal(dip)\timag(dip)"<<std::endl;
		dyna<<"#"<<std::endl;

		fnct<<"# * NUMERICAL * "<<std::endl;
		fnct<<"#"<<std::endl;
		fnct<<"# file generated by ** cup.H ** "<<std::endl;
		fnct<<"#"<<std::endl;
		fnct<<"# PARAMETERS: "<<std::endl;
		fnct<<"# G = "<<pG<<std::endl;
		fnct<<"# ome = "<<omeeV<<" eV"<<std::endl;
		fnct<<"# Electric field amplitude = "<<E0<<std::endl;
		fnct<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		fnct<<"# ------------- "<<std::endl;
		fnct<<"# t (ps)\treal(q0)\timag(q0)\treal(q1)\timag(q1)\treal(q2)\timag(q2)"<<std::endl;
		fnct<<"#"<<std::endl;

		omep = eV2j*Ome_p*2.*M_PI/h;
		// normalized variables
		nv = normalized_variables ();
		tau1 = nv[0];
		tau2 = nv[1];
		gamd = nv[2];
		ome0 = nv[3];
		GG = set_GG(G, tau2);
		
		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
		if (ns == 1) p=p3;
		    else p=p0;

		ome = omeeV/Ome_p;
		
		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
		
		OmeH = odv[0];
		OmeP = odv[1];
		GamP = odv[2];

		chiH = -GG/OmeH;
		chiP = -GamP/OmeP;

		eps1e = eps_inf + chiP;
		eps2e = eps_b   + chiH;


		// NUMERICAL SOLUTIONS
		T=T*omep*1.e-12; // in ome_p
		// Initial conditions
		dt  = 1.;
		
		double N=0;
		std::complex<double> *q ;
		q = new std::complex<double>[3] {0,0};

		int Nt=T/dt;
		int Nprint=Nt;

		if(Nt > 50000) Nprint=50000;    
// 		if(Nt > 500000) Nprint=500000;    
		int iprint=Nt/Nprint;
		GG = set_GG(G, tau2);
		// let's begin without pump
		tildeN = 0;
	    
		int i=0;

		// These variables just work on specific functions
		//int pip = 0;
		double Tperiod = 50;
		double Eperiod = 20;
		double duty = 0.5;
		double E0_temp = E0;

		while (i<=Nt){// && (fabs(norm(dip)-dnormem)> 1.e-20 || t*1.e+12/omep<12*tpump)){
		    t=i*dt;
		    i++;

		//Introduccion de una onda cuadrada para el bombeo
		double t_ps    = t * 1.e+12 / omep;  // tiempo actual en ps

		//Introduction of a square wave as the electric field
		double local_t = fmod(t_ps, Eperiod);
			    if (local_t < duty * Eperiod) {
			        E0 = 0;
			    } else {
			        E0 = E0_temp;
			    }
		omeR = Rabi_frec(G, tau2, E0);
		
		//dflt(&tildeN, t_ps, tpump, &pip);
		//squarewave(&tildeN, Tperiod, duty, t_ps,tpump, &GG, G, tau2);
		sinusoidalwave(&tildeN, Tperiod, t_ps, &GG, G, tau2);

		//Aqui termina el cambio

		    if (ns == 1) {
				po  = numerical_output(E0, q, p0);
				N   = Runge_Kutta_mono_4 (N, -1./tau1, tildeN/tau1+imag(q[0]*conj(po)), dt);
				} else N = 1;
		    A   = coefficients(OmeH, OmeP, GG*N, GamP, p0, p1, p2); 
		    B   = inhomogeneous(GG*N, GamP, E0, p0, p1, p2);
		    q   = Runge_Kutta_4 (q, A, B, dt);
		    
		    for(int ix = 0; ix <= 2; ix++) delete [] A[ix];
		    delete [] A;
		    delete [] B;
			
		    dnormem = steste*norm(dip);
		    dip = numerical_output(E0, q, p);
		    //alph = dip/E0;
				if(i%iprint==0){
				dyna<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<< // 1 time (ps)
					//"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(alph)<<     // 2 Re(alpha)
					//"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(alph)<<     // 3 Im(alpha)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(dip)<<     // 4 Re(dip)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(dip)<<     // 5 Im(dip)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<N<<             // 6 N
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(po)<<      // 7 Re(p0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(po)<<      // 8 Im(p0) 
					std::endl;
				fnct<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<< // 1 time (ps)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[0])<<   // 2 Re(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[0])<<   // 3 Im(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[1])<<   // 4 Re(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[1])<<   // 5 Im(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[2])<<   // 6 Re(q2)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[2])<<   // 7 Im(q2)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<N<<    // 8 N
					std::endl;
				}
			if (t*1.e+12/omep>2.*tpump){
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-30 && steste==1){
					//std::cout << std::endl;
					//std::cout << " 1st " << "E0 = " << E0 << " " << "t = " << t_ps;
					//std::cout << std::endl;
					break;
				}
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-10) wb=1;
// 				if (wb==1 && i%iprint==0){
				if (wb==1){
				    wave.push_back(q[0]);
				    }
				if (fabs((dnavr-dnavrm)/dnavr) < 1.e-7 && steste==1) {
					//std::cout << std::endl;
					//std::cout << " 2nd " << "E0 = " << E0 << " " << "t = " << t_ps;
					//std::cout << std::endl;
					break;}
				if (AVE==0){
					t0=t;
					AVE=1;
					}
				dpTOT = dpTOT+norm(dip);
				if (t-t0 >= 1000/omeR){
					dnavrm=dnavr;
					dnavr=dpTOT/(t-t0);
					AVE=0;
					dpTOT =0;
					}
				}
		    if (wave.size()>=150000) wave.erase(wave.begin(), wave.begin() +  50000);
		    }
		if (ns == 1){
			if (wave.size()>0) nOme=nfindOme(wave.size(), wave, dt);
				else nOme=0.;
			}
// 		std::cout<<"nOme = "<<nOme<<std::endl;
		wave.clear();
		delete [] p0;
		delete [] p1;
		delete [] p2;
		delete [] p3;
		delete [] q;
		return dip;
		}

	std::complex<double> numerical(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double T, double tpump, char* sol=0, double rho=0){
		double eps_b, eps_s, ome, dt, dnormem=1.e20, dnavr=0, dnavrm=1.e20;
		std::vector<std::complex<double>> wave;
		int ns=0;
		int wb=0.;
		if ((sol != NULL)) {
		    eps_s=set_host(sol);
		    ns=1;
		    }
		    else eps_s=0;
		eps_b=set_host(hst);
		int steste=0;
		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
		double omep, t, tau1, tau2, gamd, ome0, pG = fabs(G), *nv, tildeN;
		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP, **A, *B, *p0, *p1, *p2, *p3, *p, dip, po, *odv;
		double dip_ss;

		double dpTOT=0, t0;
		double omeR;
		int AVE=0;
		
		A = new std::complex<double>*[3];
		for(int i = 0; i <= 2; i++)
		    A[i] = new std::complex<double>[3];
		
		B = new std::complex<double>[3];


		nv = new double[4];

		odv = new std::complex<double> [3];
		
		if (T==0){ 
		    T=100000;
		    steste=1;
		    }

		std::fstream dyna, fnct;
		dyna.open("../data/output/numtime.dat", std::ios::out);
		fnct.open("../data/output/numfunc.dat", std::ios::out);
		//prb.open("../data/output/probetime.dat", std::ios::out);
		
		dyna<<"# * NUMERICAL * "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# file generated by ** cup.H ** "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# PARAMETERS: "<<std::endl;
		dyna<<"# G = "<<pG<<std::endl;
		dyna<<"# ome = "<<omeeV<<" eV"<<std::endl;
		dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		dyna<<"# ------------- "<<std::endl;
		dyna<<"# t (ps)\treal(alph)\timag(alph)\treal(dip)\timag(dip)\t dip_modulus\t N"<<std::endl;
		dyna<<"#"<<std::endl;

		fnct<<"# * NUMERICAL * "<<std::endl;
		fnct<<"#"<<std::endl;
		fnct<<"# file generated by ** cup.H ** "<<std::endl;
		fnct<<"#"<<std::endl;
		fnct<<"# PARAMETERS: "<<std::endl;
		fnct<<"# G = "<<pG<<std::endl;
		fnct<<"# ome = "<<omeeV<<" eV"<<std::endl;
		fnct<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		fnct<<"# ------------- "<<std::endl;
		fnct<<"# t (ps)\treal(q0)\timag(q0)\treal(q1)\timag(q1)\treal(q2)\timag(q2)"<<std::endl;
		fnct<<"#"<<std::endl;

		omep = eV2j*Ome_p*2.*M_PI/h;
		// normalized variables
		nv = normalized_variables ();
		tau1 = nv[0];
		tau2 = nv[1];
		gamd = nv[2];
		ome0 = nv[3];
		GG = set_GG(G, tau2);

		
		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
		if (ns == 1) p=p3;
		    else p=p0;

		ome = omeeV/Ome_p;
		
		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
		
		OmeH = odv[0];
		OmeP = odv[1];
		GamP = odv[2];

		chiH = -GG/OmeH;
		chiP = -GamP/OmeP;

		eps1e = eps_inf + chiP;
		eps2e = eps_b   + chiH;


		// NUMERICAL SOLUTIONS
		T=T*omep*1.e-12; // in ome_p
		// Initial conditions
		dt  = 1.;
		
		double N=0;
		std::complex<double> *q ;
		dip = 0.0;
		q = new std::complex<double>[3]{0.0,0.0,0.0};


		int Nt=T/dt;
		int Nprint=Nt;

		//if(Nt > 50000) Nprint=50000;    
 		if(Nt > 500000) Nprint=500000;    
		int iprint=Nt/Nprint;
		GG = set_GG(G, tau2);
		// let's begin without pump
		tildeN = 0;
	    
		omeR = Rabi_frec(G, tau2, E0);
		int i=0;

		// These variables just work on specific functions
		//int pip = 0;
		double Tperiod = 50;
		//double Eperiod;
		//double duty = 0.5;
		//double E0_temp = E0;

		while (i<=Nt){// && (fabs(norm(dip)-dnormem)> 1.e-20 || t*1.e+12/omep<12*tpump)){
		    t=i*dt;
			i++;

		//Introduccion de una onda cuadrada para el bombeo
		double t_ps    = t * 1.e+12 / omep;  // tiempo actual en ps

		//Introduction of a square wave as the electric field
		//double local_t = fmod(t_ps, Eperiod);
			    //if (local_t < duty * Eperiod) {
			        //E0 = 0;
			    //} else {
			        //E0 = E0_temp;
			    //}
		
		//dflt(&tildeN, t_ps, tpump, &pip);
		//squarewave(&tildeN, Tperiod, duty, t_ps,tpump, &GG, G, tau2);
		sinusoidalwave(&tildeN, Tperiod, t_ps, &GG, G, tau2);

		//Aqui termina el cambio

		    if (ns == 1) {
				po  = numerical_output(E0, q, p0);
				N   = Runge_Kutta_mono_4 (N, -1./tau1, tildeN/tau1+imag(q[0]*conj(po)), dt);
				} else N = 1;
		    A   = coefficients(OmeH, OmeP, GG*N, GamP, p0, p1, p2);
		    B   = inhomogeneous(GG*N, GamP, E0, p0, p1, p2);
		    q   = Runge_Kutta_4 (q, A, B, dt);
		    
		    for(int ix = 0; ix <= 2; ix++) delete [] A[ix];
		    delete [] A;
		    delete [] B;

			std::complex<double> dip_prev = dip;
			dip = numerical_output(E0, q, p);
			dnormem = steste * std::norm(dip_prev);
			dip_ss = std::sqrt(std::norm(dip));
			
		    // dnormem = steste*norm(dip);
		    // dip = numerical_output(E0, q, p);
			// dip_ss = std::sqrt(std::norm(dip));
		    alph = dip/E0;

				if(i%iprint==0){
				dyna<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<< // 1 time (ps)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(alph)<<     // 2 Re(alpha)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(alph)<<     // 3 Im(alpha)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(dip)<<     // 4 Re(dip)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(dip)<<     // 5 Im(dip)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<dip_ss<<      // dip modulus
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<N<<             // 6 N
					// "\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(po)<<      // 7 Re(p0)
					// "\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(po)<<      // 8 Im(p0)
					std::endl;
				fnct<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<< // 1 time (ps)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[0])<<   // 2 Re(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[0])<<   // 3 Im(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[1])<<   // 4 Re(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[1])<<   // 5 Im(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[2])<<   // 6 Re(q2)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[2])<<   // 7 Im(q2)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<N<<    // 8 N
					std::endl;
				}
			//if (t*1.e+12/omep>2.*tpump){
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-30 && steste==1){
					//std::cout << std::endl;
					//std::cout << " 1st " << "E0 = " << E0 << " " << "t = " << t_ps;
					//std::cout << std::endl;
					break;
				}
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-10) wb=1;
// 				if (wb==1 && i%iprint==0){
				if (wb==1){
				    wave.push_back(q[0]);
				    }
				if (fabs((dnavr-dnavrm)/dnavr) < 1.e-7 && steste==1) {
					//std::cout << std::endl;
					//std::cout << " 2nd " << "E0 = " << E0 << " " << "t = " << t_ps;
					//std::cout << std::endl;
					break;}
				if (AVE==0){
					t0=t;
					AVE=1;
					}
				dpTOT = dpTOT+norm(dip);
				if (t-t0 >= 1000/omeR){
					dnavrm=dnavr;
					dnavr=dpTOT/(t-t0);
					AVE=0;
					dpTOT =0;
					}
				//}
		    if (wave.size()>=150000) wave.erase(wave.begin(), wave.begin() +  50000);
		    }
		if (ns == 1){
			if (wave.size()>0) nOme=nfindOme(wave.size(), wave, dt);
				else nOme=0.;
			}
// 		std::cout<<"nOme = "<<nOme<<std::endl;
		wave.clear();
		delete [] p0;
		delete [] p1;
		delete [] p2;
		delete [] p3;
		delete [] q;
		return dip;
		}

	// 	std::complex<double> omegatest(
	// 	char* mdl , char* mtl ,char* hst,
	// 	double nEsat, double tpump,
	// 	double E0, double ome1, double ome2,
	// 	char* sol, double rho)
	// {

	// 	int Nfreq = 1;
	// 	double omeeV = 0.0;
	// 	double T = 0;
	// 	//double omemin = 2.52;
	// 	//double omemax = ome2*Ome_p;
	// 	//double dome = (omemax - omemin) / (Nfreq - 1);

    // 	std::fstream dyna;

    //     dyna.open("../data/output/numtime.dat", std::ios::out);
		

	// 	// dyna<<"# * FREQUENCY TEST * "<<std::endl;
	// 	// dyna<<"#"<<std::endl;
	// 	// dyna<<"# file generated by ** cup.H ** "<<std::endl;
	// 	// dyna<<"#"<<std::endl;
	// 	// dyna<<"# PARAMETERS: "<<std::endl;
	// 	// dyna<<"# G = "<< fabs(G) <<std::endl;
	// 	// dyna<<"# E0 = "<<E0<<std::endl;
    //     // dyna<<"# Esat = "<<nEsat<<std::endl;
	// 	// dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
	// 	// dyna<<"# ------------- "<<std::endl;
	// 	// dyna<<"# t\t omeeV\t Dip_max"<<std::endl;
	// 	// dyna<<"#"<<std::endl;

	// 	OmegaPoint r{};

	// 	for (int k = 0; k < Nfreq; ++k) {
    //     	//omeeV = omemin + k * dome;
	// 		//omeeV = omemax - k * dome;
	// 		omeeV = 2.52;

	// 		r = numerical1(mdl, mtl, hst, E0, omeeV, T, tpump, sol, rho);
	// 			dyna<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<r.t_ps<< // 1 time (ps)
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<omeeV<<      // dip_ss
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<r.dip_ss<<      // dip modulus
	// 				std::endl;
	// 		}
			
	// 	 return {0.0, 0.0}; // o lo que quieras devolver
	// 	}


	// 	OmegaPoint numerical1(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double T, double tpump, char* sol=0, double rho=0){
		
	// 	double eps_b, eps_s, ome, dt;
	// 	int ns=0;
	// 	if ((sol != NULL)) {
	// 	    eps_s=set_host(sol);
	// 	    ns=1;
	// 	    }
	// 	    else eps_s=0;
	// 	eps_b=set_host(hst);
	// 	int steste=1;
	// 	T=147000;

	// 	std::complex<double> eps1e, eps2e, eps1, eps2, alph;
	// 	double omep, t, tau1, tau2, gamd, ome0, *nv, tildeN;
	// 	std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP, **A, *B, *p0, *p1, *p2, *p3, *p, dip, po, *odv;
	// 	double dip_ss = 0.0, t_ps = 0.0;

	// 	double dpTOT=0, t0;
	// 	double omeR;
	// 	int AVE=0;
		
	// 	A = new std::complex<double>*[3];
	// 	for(int i = 0; i <= 2; i++)
	// 	    A[i] = new std::complex<double>[3];
		
	// 	B = new std::complex<double>[3];

	// 	nv = new double[4];

	// 	odv = new std::complex<double> [3];

    // 	std::fstream dyna;

    //     dyna.open("../data/output/probetime.dat", std::ios::out);
		

	// 	dyna<<"# * FREQUENCY TEST * "<<std::endl;
	// 	dyna<<"#"<<std::endl;
	// 	dyna<<"# file generated by ** cup.H ** "<<std::endl;
	// 	dyna<<"#"<<std::endl;
	// 	dyna<<"# PARAMETERS: "<<std::endl;
	// 	dyna<<"# G = "<< fabs(G) <<std::endl;
	// 	dyna<<"# E0 = "<<E0<<std::endl;
    //     //dyna<<"# Esat = "<<nEsat<<std::endl;
	// 	dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
	// 	dyna<<"# ------------- "<<std::endl;
	// 	dyna<<"# t\t omeeV\t Dip_max"<<std::endl;
	// 	dyna<<"#"<<std::endl;



	// 	omep = eV2j*Ome_p*2.*M_PI/h;
	// 	// normalized variables
	// 	nv = normalized_variables ();
	// 	tau1 = nv[0];
	// 	tau2 = nv[1];
	// 	gamd = nv[2];
	// 	ome0 = nv[3];
	// 	GG = set_GG(G, tau2);

	// 	p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
	// 	p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
	// 	p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
	// 	p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
	// 	if (ns == 1) p=p3;
	// 	    else p=p0;

	// 	ome = omeeV/Ome_p;
		
	// 	odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
		
	// 	OmeH = odv[0];
	// 	OmeP = odv[1];
	// 	GamP = odv[2];

	// 	chiH = -GG/OmeH;
	// 	chiP = -GamP/OmeP;

	// 	eps1e = eps_inf + chiP;
	// 	eps2e = eps_b   + chiH;


	// 	// NUMERICAL SOLUTIONS
	// 	T=T*omep*1.e-12; // in ome_p
	// 	// Initial conditions
	// 	dt  = 1.;
		
	// 	double N=0;
	// 	std::complex<double> *q ;
	// 	dip = 0.0;
	// 	q = new std::complex<double>[3]{0.0,0.0,0.0};

	// 	int Nt=T/dt;
	// 	int Nprint=Nt;

	// 	//if(Nt > 50000) Nprint=50000;    
 	// 	if(Nt > 500000) Nprint=500000;    
	// 	int iprint=Nt/Nprint;

	// 	GG = set_GG(G, tau2);
	// 	// let's begin without pump
	// 	tildeN = 0;
		
	// 	double dnormem=1.e20, dnavr=0, dnavrm=1.e20;
	    
	// 	omeR = Rabi_frec(G, tau2, E0);
	// 	int i=0;

	// 	// These variables just work on specific functions
	// 	int pip = 0;


	// 		// std::cout << std::setprecision(17);
	// 		// std::cout << "\n===== DEBUG numerical1 =====\n";
	// 		// std::cout << "eps_b   = " << eps_b << "\n";
	// 		// std::cout << "eps_s   = " << eps_s << "\n";
	// 		// std::cout << "tau1    = " << tau1 << "\n";
	// 		// std::cout << "tau2    = " << tau2 << "\n";
	// 		// std::cout << "gamd    = " << gamd << "\n";
	// 		// std::cout << "ome0    = " << ome0 << "\n";
	// 		// std::cout << "GG      = " << GG << "\n";
	// 		// std::cout << "OmeH    = " << OmeH << "\n";
	// 		// std::cout << "OmeP    = " << OmeP << "\n";
	// 		// std::cout << "GamP    = " << GamP << "\n";
	// 		// std::cout << "omeeV   = " << omeeV << "\n";
	// 		// std::cout << "============================\n";

	// 	while (i<=Nt){// && (fabs(norm(dip)-dnormem)> 1.e-20 || t*1.e+12/omep<12*tpump)){
	// 	    t=i*dt;
	// 		i++;

	// 	//Introduccion de una onda cuadrada para el bombeo
	// 	t_ps    = t * 1.e+12 / omep;  // tiempo actual en ps
		
	// 	dflt(&tildeN, t_ps, tpump, &pip);
	// 	//squarewave(&tildeN, Tperiod, duty, t_ps,tpump, &GG, G, tau2);
	// 	//sinusoidalwave(&tildeN, Tperiod, t_ps, &GG, G, tau2);

	// 	//Aqui termina el cambio

	// 	    if (ns == 1) {
	// 			po  = numerical_output(E0, q, p0);
	// 			N   = Runge_Kutta_mono_4 (N, -1./tau1, tildeN/tau1+imag(q[0]*conj(po)), dt);
	// 			} else N = 1;
	// 	    A   = coefficients(OmeH, OmeP, GG*N, GamP, p0, p1, p2);
	// 	    B   = inhomogeneous(GG*N, GamP, E0, p0, p1, p2);
	// 	    q   = Runge_Kutta_4 (q, A, B, dt);
		    
	// 	    for(int ix = 0; ix <= 2; ix++) delete [] A[ix];
	// 	    delete [] A;
	// 	    delete [] B;

	// 		std::complex<double> dip_prev = dip;
	// 		dip = numerical_output(E0, q, p);
	// 		dnormem = steste * std::norm(dip_prev);
	// 		dip_ss = std::sqrt(std::norm(dip));

	// 		if (i <= 3) {
	// 		 std::cout << std::setprecision(17)
	// 		<< "i="<<i
	// 		<< " t_ps="<<t_ps
	// 		<< " tildeN="<<tildeN
	// 		<< " N="<<N
	// 		<< " q0="<<q[0] << " q1="<<q[1] << " q2="<<q[2]
	// 		<< " dip="<<dip
	// 		<< "\n";
	// 		}

	// 		if(i%iprint==0){
	// 			dyna<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<< // 1 time (ps)
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(dip)<<     // 4 Re(dip)
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(dip)<<     // 5 Im(dip)
	// 				"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<dip_ss<<      // dip modulus

	// 				std::endl;
	// 			}

	// 		//if (t*1.e+12/omep>2.*tpump){
	// 			if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-30 && steste==1){
	// 				break;
	// 			}
	// 			if (fabs((dnavr-dnavrm)/dnavr) < 1.e-7 && steste==1) {
	// 				break;}
	// 			if (AVE==0){
	// 				t0=t;
	// 				AVE=1;
	// 				}
	// 			dpTOT = dpTOT+norm(dip);
	// 			if (t-t0 >= 1000/omeR){
	// 				dnavrm=dnavr;
	// 				dnavr=dpTOT/(t-t0);
	// 				AVE=0;
	// 				dpTOT =0;
	// 				}
	// 			//}
	// 	    }
	// 	delete [] p0;
	// 	delete [] p1;
	// 	delete [] p2;
	// 	delete [] p3;
	// 	delete [] q;
		
	// 	return OmegaPoint{t_ps, dip_ss};
	// 	}
    

	    std::complex<double> analytical(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double T, double tpump, char* sol=0, double rho=0){
		double eps_b, eps_s, ome, dt, dnormem=1.e20;
    		if ((sol != NULL)) eps_s=set_host(sol);
		    else eps_s=0;
		eps_b=set_host(hst);
		int steste=0, amp=0;
		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
		double omep, t, tau1, tau2, gamd, ome0, pG = fabs(G), *nv, N, tildeN;
		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP, **A, *B, *p0, *p1, *p2, *p3, *p, dip, *odv;
				
		A = new std::complex<double>*[3];
		for(int i = 0; i <= 2; i++)
		    A[i] = new std::complex<double>[3];
		
		B = new std::complex<double>[3];
		
		nv = new double[4];

		odv = new std::complex<double> [3];
		
		if (T==0){ 
		    T=100000;
		    steste=1;
		    }
		
		std::fstream dyna, fnct;
		dyna.open("../data/output/anltime.dat", std::ios::out);
		fnct.open("../data/output/anlfunc.dat", std::ios::out);
		
		dyna<<"# * ANALYTICAL * "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# file generated by ** cup.H ** "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# PARAMETERS: "<<std::endl;
		dyna<<"# G = "<<pG<<std::endl;
		dyna<<"# ome = "<<omeeV<<" eV"<<std::endl;
		dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		dyna<<"# ------------- "<<std::endl;
		dyna<<"# t (ps)\treal(alph)\timag(alph)\treal(dip)\timag(dip)"<<std::endl;
		dyna<<"#"<<std::endl;

		fnct<<"# * ANALYTICAL * "<<std::endl;
		fnct<<"#"<<std::endl;
		fnct<<"# file generated by ** cup.H ** "<<std::endl;
		fnct<<"#"<<std::endl;
		fnct<<"# PARAMETERS: "<<std::endl;
		fnct<<"# G = "<<pG<<std::endl;
		fnct<<"# ome = "<<omeeV<<" eV"<<std::endl;
		fnct<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		fnct<<"# ------------- "<<std::endl;
		fnct<<"# t (ps)\treal(q0)\timag(q0)\treal(q1)\timag(q1)\treal(q2)\timag(q2)"<<std::endl;
		fnct<<"#"<<std::endl;

		omep = eV2j*Ome_p*2.*M_PI/h;// omep = eV2j*Ome_p/h;  //converto in Hz
		// normalized variables
		nv = normalized_variables ();
		tau1 = nv[0];
		tau2 = nv[1];
		gamd = nv[2];
		ome0 = nv[3];
		GG = set_GG(G, tau2);

		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
		if ((sol != NULL)) p=p3;
		    else p=p0;
		ome = omeeV/Ome_p;
		
			
		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
		
		OmeH = odv[0];
		OmeP = odv[1];
		GamP = odv[2];

		chiH = -GG/OmeH;
		chiP = -GamP/OmeP;

		eps1e = eps_inf + chiP;
		eps2e = eps_b   + chiH;

		// coefficients of the differential system of equations
// 		A  = coefficients(OmeH, OmeP, GG, GamP, p0, p1, p2);  
// 		B  = inhomogeneous(GG, GamP, E0, p0, p1, p2);

		// ANALYTICAL SOLUTIONS
		T=T*omep*1.e-12; // in ome_p
		// Initial conditions
		dt  = 1.;
		std::complex<double> *q, *qss, *kap, **EVE, *C;
		
		q   = new std::complex<double>[3];
		for(int i = 0; i <= 1; i++) q[i] = std::complex<double> (0., 0.);
		
		qss = new std::complex<double>[3];
		kap = new std::complex<double>[3];		
		EVE = new std::complex<double>*[3];
		for(int i = 0; i <= 2; i++)
		    EVE[i] = new std::complex<double>[3];
		C   = new std::complex<double>[3];

		int Nt=T/dt;
		int Nprint=Nt;

		if(Nt > 500000) Nprint=500000;    
		int iprint=Nt/Nprint;

		// let's begin without pump
		tildeN = 0;
		N      = 0;
		GG  = set_GG(G*N, tau2);; 
		A   = coefficients(OmeH, OmeP, GG, GamP, p0, p1, p2); 
		B   = inhomogeneous(GG, GamP, E0, p0, p1, p2);
		qss = steady_state_solution(A, B);
		kap = eigenvalues(A);
		EVE = eigenvectors(A);
		C   = constants(A, B, q);
		double t0 = 0.;
		int i=0, pip=0;
		while (i<=Nt){
		    t=i*dt;
		    i++;
		    if(t*1.e+12/omep >= tpump && pip ==0) {
			// the pump is switched on at t = tpump
			tildeN = 1.;
			t0  = t;
			pip = 1.;
			GG  = set_GG(G*tildeN, tau2);
			A   = coefficients(OmeH, OmeP, GG, GamP, p0, p1, p2);
			kap = eigenvalues(A);
			EVE = eigenvectors(A);
			B   = inhomogeneous(GG, GamP, E0, p0, p1, p2);
			qss = steady_state_solution(A, B);
			C   = constants(A, B, q);
			if (steste == 1){
				for (int ei =0; ei<=2; ei++)
					if (kap[ei].real() > 0) amp = 1;
						else amp =0;
				}
			}
		    N   = tildeN*(1.-exp(-(t-t0)/(tau1)));
		    for(int ii = 0; ii <= 2; ii++)
				q[ii]   = qss[ii] + C[0]*EVE[ii][0]*exp(kap[0]*(t-t0))
								  + C[1]*EVE[ii][1]*exp(kap[1]*(t-t0))
								  + C[2]*EVE[ii][2]*exp(kap[2]*(t-t0));
    
		    dnormem = steste*norm(dip);
		    dip = numerical_output(E0, q, p);
		    alph = dip/E0;
		    
		    if(i%iprint==0){
				dyna<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<<   // 1 time (ps)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(alph)<<   // 2 Re(alpha)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(alph)<<   // 3 Im(alpha)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(dip)<<    // 4 Re(dip)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(dip)<<    // 5 Im(dip)
					std::endl;
				fnct<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<<   // 1 time (ps)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[0])<<   // 2 Re(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[0])<<   // 3 Im(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[1])<<   // 4 Re(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[1])<<   // 5 Im(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[2])<<   // 6 Re(q2)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[2])<<   // 7 Im(q2)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<N<<            // 8 Im(q2)
					std::endl;
				}
			if (amp == 1 && norm(dip)/E0 > 5.e26) break;
			if (amp == 0 && t*1.e+12/omep>2.*tpump) 
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-20) break;
		    }
		delete [] A;
		delete [] B;
		delete [] C;
		delete [] qss;
		delete [] EVE;
		return dip;
		}
	};
