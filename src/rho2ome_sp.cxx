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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <armadillo>
#include "headers/math33.H"
#include "headers/single.H"
#include "headers/cup.H"
#include "headers/ns_ISS.H"

using namespace std;

/** Compila con: 
g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/rho2ome_sp.cxx -o ../bin/rho2ome_sp -lgsl -lgslcblas -lm -larmadillo
**/

// complex<double> img = complex<double> (0., 1.);
double pi=acos(-1.);

// complex<double> F(complex<double> eps2, double eps3, double rap3){
//     complex<double> eff;
//     eff=eps2*(2.*eps2*(rap3-1.)-2.*eps3*(rap3+2.))/(eps2*(2.*rap3+1.)+2.*eps3*(1.-rap3));
//     return eff;
//     }
// double eqn9(double ome, double ome0, complex<double> F, double eps_b, double Dome){
//     double eq9;
//     eq9=Dome*(eps_b-real(F))-2.*(ome - ome0)*imag(F);
//     return eq9;
//     }
int main(int argc, char ** argv){
    if (argc==1){
        cout<<endl<<" Usage: "<<argv[0]<<" <rho>"<<endl<<endl;
        exit(1);
        }
  double rap, rho=atof(argv[1]);
  double omemi, omema, eps_b, eps3, rap3, ome0, omep, ome_sp; // 0.66 para 2.98 para -0.55
  complex<double> eps1, eps2, eps_h, F_ome;
  complex<double> m11, m12, m13, m21, m22, m23, m31, m32, m33;
  
  complex<double> q1, sig, q2, E0;
  double ieh_sp, Del, N0=1.;

  
  char mtl[16], mdl[16], hst[16], sol[16], active[16];
  nanosphere ns;
  fstream nano, xtrm;  
  
  nano.open("../data/input/nanosphere_eV.dat", ios::in);
  xtrm.open("../data/input/extremes.dat", ios::in);
    
  nano>>ns.r1>>ns.Dome>>ns.ome_0>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rap>>hst;
  
  ns.r1=ns.r1*1.e-9;
  ns.r2=ns.r1/ns.rap;

  ns.init();
  ns.set_metal(mtl,mdl,1);
  
  eps_b=ns.set_host(hst);  // poniendo los valores de la permitividad silicon
  eps3=ns.set_host(sol);   // water
  ns.set_active(active);  // material activo con una lorentziana
  
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//bisection method
   
    cout.precision(7);        //set the precision
    cout.setf(ios::fixed);
    double eqn9_A, eqn9_D, a, b, a0, b0, d, e;    //declare some needed variables
    complex<double> eps2_A, eps2_B, eps2_D, F_A, F_B, F_D;
    xtrm>>a0>>b0;
    e = 1.e-10;

   double tau1, tau2, ome_sm=666., p1_st2, p1_sp2;

   tau2=2./ns.Dome;
   tau1=0.1*tau2;
  
   ome0 = eV2erg*ns.ome_0/he;	//converto in Hz
   omep = eV2erg*ns.Ome_p/he;	//converto in Hz
   tau2=2.*he/(eV2erg*ns.Dome);	//converto in s
   tau1=0.1*tau2;
    
   ome0 = ome0/omep;	// 
   tau2=tau2*omep;	//converto in omep ok, y la calculo con la ecuacion que 
   tau1=tau1*omep;	//converto in omep


///////// GAIN THRESHOLD ///////////
double ome21, eps_im;

eps_im = ns.G;
    ns.rap=rho;
    rap3=ns.rap*ns.rap*ns.rap;
    a=a0;
    b=b0;
    ome_sp =.5*(b0+a0);
    ome_sm=666.;
    while (abs(ome_sp-ome_sm)>1.e-20){ /* ome_sp dipende da rho e da ome21, riposiziono ome21 su ome_sp in modo che sia solo funzione di rho*/
        ome21=ome_sp;
        a=a0;
        b=b0;
        while (fabs(a-b)>=e){   /*if the mod of a-b is greater than the accuracy desired keep bisecting the interval*/
            
            d=(a+b)/2.0;        //bisect the interval and find the value of c

            eps2_A = ns.metal(a);
            F_A = F(eps2_A, eps3, rap3);
            eqn9_A = eqn9(a, ome21, F_A, eps_b, ns.Dome);

            eps2_D = ns.metal(d);
            F_D = F(eps2_D, eps3, rap3);
            eqn9_D = eqn9(d, ome21, F_D, eps_b, ns.Dome);        

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
    ome_sp   = eV2erg*ome_sp/he;
    ome_sp   = ome_sp/omep;
    ome0   = eV2erg*ome21/he;	
    ome0   = ome0/omep;

    Del=2./tau2;

//        ieh_sp = (2.*(ome_sp-ome0)*(eps_b-real(F_D))+Del*imag(F_D))/Del;
    ieh_sp = Del*imag(F_D)/Del;     
    p1_st2=N0/(eps_im*tau1);
//        p1_sp2 = p1_st2*(4.*(ome_sp-ome0)*(ome_sp-ome0)+Del*Del)*(ieh_sp-eps_im)/(ieh_sp*Del*Del);
    p1_sp2 = p1_st2*Del*Del*(ieh_sp-eps_im)/(ieh_sp*Del*Del);
    
    ome_sp   = ome_sp*omep;
    ome_sp   = ome_sp*he/eV2erg;            

    cout<<" "<<ome_sp<<" "<<ieh_sp<<" "<<p1_sp2;
  return 0;
  }
