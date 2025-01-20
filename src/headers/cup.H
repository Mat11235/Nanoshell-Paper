/*
 * This file is part of the Nano-Shell Simulation Project.
 * 
 * Copyright (C) 2025 Alessandro Veltri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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

double eV2nm(double ome){
  double lam;
  lam=cc/(ome*eV2j/h)*1.e9;
  return lam;
  }

class nanosphere{
    private:
	gsl_interp_accel *acc;
	gsl_interp_accel *ach;
	gsl_spline *reeps;
	gsl_spline *imeps;
	gsl_spline *reeph;
	gsl_spline *imeph;
	char *wrk;
    public:
	std::complex<double> img;
	double lrntz, act;
	double pi, vF, enha;
	double A, r1, r2, rap;
	double G, nOme;
	double ome_p, ome_L, tau_d, tau_L, Dome, lam_0, eps_0;
	double ome_p2,ome_L2, tau_d2, tau_L2, ome_0, gamd, kap;
	double ome_o2, ome, f, n, enne, omeB;
	double Ome_p, Gam_d, Gam, eps_inf, Ome_L, Gam_L, Deps, DELTAOme;
	std::complex<double> eps3;
	void init(){
	  /** all values are in MKS **/
	  img = std::complex<double> (0., 1.);
	  pi = acos(-1.);
	  }

	void set_metal(const char* mtl, const char* mdl, int sel=0){
	  /** Ome_p [THz]	Plasma frequency (gold)		**/
	  /** Ome_L [THz]	Lorentz frequency (gold)	**/
	  /** Gam_d [THz]	Drude damping coefficient (gold)**/
	  /** Gam_L [THz]	Lorentz spectral width (gold)	**/
	  /** Deps weighting factor.			**/
	  double *omem, *reps, *ieps, la, Dr, Di, lc;
	  char line [900], *jc;
	  int rows=0;
	  jc = new char [35];
	  if(strcmp(mtl,"gold")==0){
	      vF=1.39e6;
	      /** m/s **/
// 	      lc=5e-9; // va bene per lo spettro
// 	      lc=1.5e-8; // va bene per le misure sul setup senza pompa
	      lc=1.e-8; // va bene per le misure sul setup
// 	      Gam=h*j2eV*vF/(r1*2.*pi);
	      Gam=h*j2eV*vF/(lc*2.*pi);
// 	      std::cout<<r1<<" Gam = "<<Gam<<std::endl;
	      Ome_p=8.91;
	      Gam_d=0.0759;
// 	      std::cout<<"Gam_d = "<<Gam_d<<std::endl;
	      eps_inf=9.0685;
	      if(strcmp(mdl,"drude")==0){
		lrntz=0;
		} else if(strcmp(mdl,"drude-lorentz")==0){
		  if (sel==1){
		    Ome_p=8.91;//2113.6;
		    Gam_d=0.0759;//15.92;
		    eps_inf=5.9673;
		    Deps=1.09;
		    Ome_L=1.90724;//wait!!
		    Gam_L=11.8238;//wait!!
		    }
		  else if (sel==2){
		    Ome_p=2106.78;
		    Gam_d=2.46;
		    eps_inf=6.617;
		    Deps=1.73;
		    Ome_L=730.56;
		    Gam_L=222.0;
		}else exit(1);
		lrntz=1;
	      } else if(strcmp(mdl,"spline")==0){
		if (sel==1) strcpy(jc,"../data/input/goldJCeV.dat");
		else if (sel==2) strcpy(jc,"../data/input/goldWeV.dat");
		else exit(-1);
		std::fstream inp;
		inp.open(jc,std::ios::in);
		while(!inp.eof()){
		  inp.getline (line,900);
		  rows++;
		  }
		inp.close();
		omem=new double [rows];
		reps=new double [rows];
		ieps=new double [rows];
		inp.open(jc, std::ios::in);
		for (int i=0; i<rows; i++)
		  inp>>la>>omem[i]>>reps[i]>>ieps[i]>>Dr>>Di;
		acc = gsl_interp_accel_alloc ();
		reeps = gsl_spline_alloc (gsl_interp_cspline, rows);
		imeps = gsl_spline_alloc (gsl_interp_cspline, rows);
		gsl_spline_init (reeps, omem, reps, rows);
		gsl_spline_init (imeps, omem, ieps, rows);
		lrntz=-1;
		}
	      else {
		std::cout<<"you chosed "<<mdl<<" while option for gold are: drude, drude-lorentz, spline"<<std::endl;
		exit(1);
		}
	    } else  if(strcmp(mtl,"silver")==0){
	      vF=1.39e6;
	      /** m/s **/
	      lc=1.e-8; // 
	      Gam=h*j2eV*vF/(lc*2.*pi);
	      Ome_p=9.6;
	      Gam_d=0.0228;
	      eps_inf=5.3;
		if(strcmp(mdl,"drude")==0){
		lrntz=0;
		} else if(strcmp(mdl,"spline")==0){
		strcpy(jc,"../data/input/silverJCeV.dat");
		std::fstream inp;
		inp.open(jc,std::ios::in);
		while(!inp.eof()){
		  inp.getline (line,900);
		  rows++;
		  }
		inp.close();
		omem=new double [rows];
		reps=new double [rows];
		ieps=new double [rows];

		inp.open(jc, std::ios::in);
		for (int i=0; i<rows; i++)
		  inp>>la>>omem[i]>>reps[i]>>ieps[i]>>Dr>>Di;
		acc = gsl_interp_accel_alloc ();
		reeps = gsl_spline_alloc (gsl_interp_cspline, rows);
		imeps = gsl_spline_alloc (gsl_interp_cspline, rows);
		gsl_spline_init (reeps, omem, reps, rows);
		gsl_spline_init (imeps, omem, ieps, rows);
		lrntz=-1;
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
	  ome_L=Ome_L*1.e12;
	  tau_d=1./(Gam_d*1.e12);
	  tau_L=1./(Gam_L*1e12);
	  ome_p2=ome_p*ome_p;
	  ome_L2=ome_L*ome_L;
	  ome_o2=ome_0*ome_0;
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
	  char line [900], *mc;
	  double *omem, *reph, *ieph;
	  int rows=0;
	  mc = new char [25];
	  if(strcmp(mod,"lorentz")==0) act=1.;
	    else if(strcmp(mod,"flat")==0) act=0.;
	    else if(strcmp(mod,"interpole")==0){
		std::fstream inp;
		strcpy(mc,"../data/input/mc.dat");
		inp.open("../data/input/mc.dat",std::ios::in);
		while(!inp.eof()){
		  inp.getline (line,900);
		  rows++;
		  }
		rows--;
		inp.close();
		omem=new double [rows];
		reph=new double [rows];
		ieph=new double [rows];
		inp.open("../data/input/mc.dat", std::ios::in);
		for (int i=0; i<rows; i++){
		  inp>>omem[i]>>reph[i]>>ieph[i];
		  reph[i]=reph[i]-2.220;
		  ieph[i]=ieph[i]+0.12;
		  }
		ach = gsl_interp_accel_alloc ();
		reeph = gsl_spline_alloc (gsl_interp_cspline, rows);
		imeph = gsl_spline_alloc (gsl_interp_cspline, rows);
		gsl_spline_init (reeph, omem, reph, rows);
		gsl_spline_init (imeph, omem, ieph, rows);
		act=66.;
	      }
	    else{
	      std::cout<<"option are: flat, lorentz, interpole"<<std::endl;
	      exit(1);
	      }
	  }

	std::complex<double> metal(double ome){
	  std::complex<double> eps;
	  if (lrntz == -1){
	     eps.real(gsl_spline_eval (reeps, ome, acc));
	     eps.imag(gsl_spline_eval (imeps, ome, acc));
	    } else eps= eps_inf - ome_p2/(ome*(ome+img*Gam_d));// - lrntz*Deps*ome_L2/((ome*ome-ome_L2)+img*ome/tau_L);
	  return eps;
	  }
	  
	std::complex<double> confinement(double ome){
	  std::complex<double> dlt;
	  dlt=ome_p2/(ome*(ome+img*Gam_d))-ome_p2/(ome*(ome+img*(Gam_d+Gam)));
	  return dlt;
	  }


	std::complex<double> active(double ome, double epsh){
	  std::complex<double> eps, epsH;
	  if (act!=66.){
		  epsH = std::complex<double> (epsh,0.);
	      eps=epsH-act*G*Dome/(2.*(ome-ome_0)+img*Dome)+(1-act)*G*img;
	      } else {
		eps.real(epsh+enha*gsl_spline_eval (reeph, ome, ach));
		eps.imag(enha*gsl_spline_eval (imeph, ome, ach));
		}
	  kap =-.5*imag(eps)/(sqrt(real(eps)));//1/sqrt(2) // sqrt(.5*(-real(eps)+abs(eps)));
	  gamd = 4.*pi*kap*ome/(j2eV*h*cc);
	  enne = gamd/(.95*8.e-20);//.95 quantum yeld 8.e-20 cross section
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
	Gth = Dome*imag(F_D)/Dome;
	
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
	    
	    double tau2, gamd, ome0, eps3, pG = fabs(G), *p0, *p1, *p2, *nv;
	    std::complex<double> GG, OmeH, GamP, OmeP, **A, detA, *kap, *odv;
	    
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

		p0 = pcfc0(eps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(eps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(eps_inf, eps_b, eps_s, rho);
		
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
	    
	    std::complex<double> numerical(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double T, double tpump, char* sol=0, double rho=0){
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
		double omep, t, tau1, tau2, gamd, ome0, pG = fabs(G), *p0, *p1, *p2, *p3, *p, *nv, tildeN;
		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP, **A, *B, dip, po, *odv;

		double dpTOT=0, t0;
		double omeR;
		int AVE=0;
		
		A = new std::complex<double>*[3];
		for(int i = 0; i <= 2; i++)
		    A[i] = new std::complex<double>[3];
		
		B = new std::complex<double>[3];


		nv = new double[4];

		odv = new std::complex<double> [3];
		
		p0 = new double[3];
		p1 = new double[3];
		p2 = new double[3];
		p3 = new double[3];
		
		if (T==0){ 
		    T=100000;
		    steste=1;
		    }

		std::fstream dyna, fnct;
		dyna.open("../data/output/numtime.dat", std::ios::out);
		fnct.open("../data/output/numfunc.dat", std::ios::out);
		
		dyna<<"# * NUMERICAL * "<<std::endl;
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

		omep = eV2j*Ome_p*2.*M_PI/h;// omep = eV2j*Ome_p/h;  //converto in Hz
		// normalized variables
		nv = normalized_variables ();
		tau1 = nv[0];
		tau2 = nv[1];
		gamd = nv[2];
		ome0 = nv[3];
// 		std::cout<<"tau1 = "<<tau1/omep*1.e+12<<" ps"<<std::endl;
// 		exit(0);
		GG = set_GG(G, tau2);
// 		std::cout<<"|Esat|² = 1./("<<fabs(G)<<"*"<<tau1<<") = "<<1./(fabs(G)*tau1)<<std::endl;
// 		std::cout<<" ../bin/Esat   = "<<sqrt(1./(fabs(G)*tau1))<<std::endl;
		
		p0 = pcfc0(eps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(eps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(eps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(eps_inf, eps_b, eps_s, rho);
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
		int i=0, pip=0;;
		while (i<=Nt){// && (fabs(norm(dip)-dnormem)> 1.e-20 || t*1.e+12/omep<12*tpump)){
		    t=i*dt;
		    i++;
		    if(t*1.e+12/omep >= tpump && pip ==0) {
				tildeN =1.;
// 				GG = set_GG(G, tau2);
				pip = 1.;
				}
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
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-30 && steste==1) break;
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-10) wb=1;
// 				if (wb==1 && i%iprint==0){
				if (wb==1){
				    wave.push_back(q[0]);
				    }
				if (fabs((dnavr-dnavrm)/dnavr) < 1.e-7 && steste==1) break;
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
    
	    std::complex<double> analytical(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double T, double tpump, char* sol=0, double rho=0){
		double eps_b, eps_s, ome, dt, dnormem=1.e20;
    		if ((sol != NULL)) eps_s=set_host(sol);
		    else eps_s=0;
		eps_b=set_host(hst);
		int steste=0, amp=0;
		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
		double omep, t, tau1, tau2, gamd, ome0, pG = fabs(G), *p0, *p1, *p2, *p3, *p, *nv, N, tildeN;
		std::complex<double> GG, OmeH, chiH, GamP, OmeP, chiP, **A, *B, dip, *odv;
				
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

		p0 = pcfc0(eps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(eps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(eps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(eps_inf, eps_b, eps_s, rho);
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
