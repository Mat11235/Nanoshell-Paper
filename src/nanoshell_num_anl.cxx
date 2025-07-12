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
#include <armadillo>
#include "headers/math33.H"
#include "headers/nanoshell.H"
#include "headers/cup.H"
#include "headers/ns_ISS.H"

/*
g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/nanoshell_num_anl.cxx -o ../bin/nsh -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv){
    double   omeeV, omemi, omema, eps_b, E0, T, tpump, eps3, rho, *fro, Ome;
    complex<double> eps1, eps2, alph, alph_num, alph_anl, p3;

    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    if (argv[1]==0){
        cout<<endl<<"  Usage: "<<argv[0]<<" <omega in eV>"<<endl<<endl;
        exit(0);
        }
    omeeV=atof(argv[1]);
    
    nanosphere ns;
    ifstream nano("../data/input/nanosphere_eV.dat");
    if (!nano) {
        cerr << "Error: Cannot open input file" << endl;
        return 1;
    }
    ifstream time("../data/input/time.dat");
    if (!time) {
        cerr << "Error: Cannot open input file" << endl;
        return 1;
    }
    ofstream frlc("../data/output/frohlich.dat");
    if (!frlc) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }
    ofstream omga("../data/output/omega.dat");
    if (!omga) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }
    
    nano>>ns.r1>>ns.Dome>>ns.ome_0>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    time>>T>>tpump;    
    
    if (E0==0.) E0=1.e-30; // zero is problematic as a value for E0 
    
    ns.init();
    ns.set_metal(mtl,mdl,1);
    eps3=ns.set_host(sol);
    eps_b=ns.set_host(hst);
    ns.set_active(active);
    
    fro=ns.frohlich(omemi, omema, eps_b, eps3, rho);
    cout.precision(10);
    cout.setf(ios::fixed);
    frlc<<fro[0]<<" "<<fro[1]<<endl;
    

    ns.ome_0=fro[0];

    cout<<fro[0]<<" "<<fro[1]<<endl;
    
    ns.steady_state(mdl, mtl, hst, omemi, omema, 1000, sol, rho);
    alph_num=ns.numerical(mdl, mtl, hst, E0, omeeV, T, tpump, sol, rho)/E0;
    alph_anl=ns.analytical(mdl, mtl, hst, E0, omeeV, T, tpump, sol, rho)/E0;


    
    eps1 = ns.active(omeeV,eps_b);
    eps2 = ns.metal(omeeV);
    
    
    alph = polarizability(eps1, eps2, eps3, rho);
    
    omga<<omeeV<<" "<<real(alph)<<" "<<imag(alph)<<endl;

    cout<<"reference  "<<omeeV<<" "<<real(alph)<<" "<<imag(alph)<<endl;
    cout<<"analytical "<<omeeV<<" "<<real(alph_anl)<<" "<<imag(alph_anl)<<endl;
    cout<<"numerical  "<<omeeV<<" "<<real(alph_num)<<" "<<imag(alph_num)<<endl;

    Ome = find_Omega(ns, omeeV, hst, sol, rho);
    cout<<" Ome = "<<Ome*ns.Ome_p<<endl;
    cout<<"nOme = "<<ns.nOme<<endl;
//     intensity_steady_state(ns, mdl, mtl, hst, omemi, omema, sol, rho);

    return 0;
    }
    
