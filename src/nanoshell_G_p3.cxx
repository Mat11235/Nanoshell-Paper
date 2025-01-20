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
g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/nanoshell_G_p3.cxx -o ../bin/Gap -lgsl -lgslcblas -lm -larmadillo
*/

double interpolate(double x1, double x2, double f1, double f2, double x) {
    // Perform linear interpolation
    return f1 + (f2 - f1) * (x - x1) / (x2 - x1);
}

std::pair<double, double> find_zeros(const std::vector<double>& x, const std::vector<double>& fx, double epsilon = 1e-6, int max_iterations = 100) {
    double zero_1 = 0.0;
    double zero_2 = 0.0;
    double c, fc;
    
    for (int i = 0; i < int(x.size()); i++) {
        double x1 = x[i];
        double x2 = x[i + 1];
        double f1 = fx[i];
        double f2 = fx[i + 1];
        
        
        if (f1 * f2 < 0) {
            double a = x1;
            double b = x2;
            int iterations = 0;
            while (std::abs(b - a) > epsilon && iterations < max_iterations) {
                iterations++;
                c = (a + b) / 2.0;
                fc = interpolate(x1, x2, f1, f2, c);

                if (f1 * fc < 0) {
                    b = c;
                }
                else {
                    a = c;
                }
            }

            if (zero_1 == 0.0) {
                zero_1 = c;
            }
            else {
                zero_2 = c;
                break;
            }
        }
        if (i==int(x.size()) - 1 && zero_2 == 0.0) {
            zero_2 = x[i];
            break;
            }
        
    }
    return std::make_pair(zero_1, zero_2);
}

std::vector<double> extract_ome(std::vector<std::pair<double,std::complex<double>>> vectorOfPairs){
    std::vector<double> extractFirst;
        int i =0;
        for (const auto& pair : vectorOfPairs) {
            i++;
            extractFirst.push_back(pair.first);
            if (i==int(vectorOfPairs.size())-1) break;
        }
        return extractFirst;
        }
    
std::vector<double> extract_ralph(std::vector<std::pair<double,std::complex<double>>> vectorOfPairs){
    std::vector<double> extractSecnd;
    int i =0;
    for (const auto& pair : vectorOfPairs) {
        i++;
        extractSecnd.push_back(real(pair.second));
        if (i==int(vectorOfPairs.size())-1) break;
        }
    return extractSecnd;
    }


std::pair<double, double> fnd_extrm(std::vector<std::pair<double,std::complex<double>>> vkape, double Ome_p){
    std::vector<double> tmp;
    double extr_1, extr_2;
    int never=1;
    for (int ii=0; ii<int(vkape.size()); ii++){
        if (real(vkape[ii].second)>=0) {
            tmp.push_back(vkape[ii].first-imag(vkape[ii].second)*Ome_p);
            never=0;
            }
        }
    sort(tmp.begin(), tmp.end());
    if (never ==0 ){
        extr_1=tmp[0];
        extr_2=tmp[tmp.size()-1];
        } else{
        extr_1=666;
        extr_2=777;
        } 
    return std::make_pair(extr_1, extr_2);
    }
    
std::vector<double> extract_ialph(std::vector<std::pair<double,std::complex<double>>> vectorOfPairs){
    std::vector<double> extractSecnd;
    int i =0;
    for (const auto& pair : vectorOfPairs) {
        i++;
        extractSecnd.push_back(imag(pair.second));
        if (i==int(vectorOfPairs.size())-1) break;
        }
    return extractSecnd;
    }

bool compareSecond(const std::pair<double, double>& a, const std::pair<double, double>& b) {
    return fabs(a.second) > fabs(b.second);
    }

using namespace std;

int main(){
    double omemi, omema, E0, eps_b, eps3, *fro, rho, omeeV, kex1, kex2, dG, wG, Ome, *res, ewh;
    double Isat, tildeN=1, ntau1, ntau2;
    int omeN = 10000, GN=400;
    complex<double> alph, p3nm, eps1, eps2;
    
    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    
    res     = new double[7];
    
    std::vector<std::pair<double,std::complex<double>>> vkape;
    std::pair<double,double> kzero;
    
    nanosphere ns;
    
    fstream nano, emix;
    nano.open("../data/input/nanosphere_eV.dat", ios::in);
    emix.open("../data/output/emission_maximum.dat", ios::out);
    
    nano>>ns.r1>>ns.Dome>>ns.ome_0>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    
    E0=1.e-8;
    ns.init();
    eps3=ns.set_host(sol);
    eps_b=ns.set_host(hst);
    ns.set_metal(mtl,mdl,1);
    ns.set_active(active);
    
    fro=ns.frohlich(omemi, omema, eps_b, eps3, rho);   
//     dG=2.*fro[1]/GN;
    dG=40.*fro[1]/GN;
   
    omeeV    = fro[0];
    ns.ome_0 = omeeV;

    eps2=ns.metal(omeeV);
    
    ntau2 = 2./ns.Dome;
    ntau1 = 5.*ntau2;
    
    for (int jG=0; jG<=GN; jG++){
        ns.G=jG*dG;
        Isat=tildeN/(fabs(ns.G)*ntau1);
        wG=ns.G/fro[1];
        cout<<jG<<" "<<ns.G<<" "<<wG<<endl;

        eps1    = ns.active(omeeV,eps_b);
        alph    = polarizability(eps1, eps2, eps3, rho);
        Ome     = find_Omega1(ns, omeeV, hst, sol, rho);

        vkape   = gimme_emi_kap(ns, mdl, mtl, hst, omemi, omema, omeN, sol, rho);    
        
        kzero   = fnd_extrm(vkape, ns.Ome_p);
        
        kex1    = kzero.first;
        kex2    = kzero.second;

        res     = ISS_results(ns, Ome, kex1/ns.Ome_p, kex2/ns.Ome_p, omeeV, hst, sol, rho);       

        p3nm    = ns.numerical(mdl, mtl, hst, E0, omeeV, 0, 0, sol, rho);
        
        if (abs(kex2-kex1)>10) ewh = 0;
            else ewh = kex2-kex1;
        emix<<ns.G<<" "<<wG<<" "<<norm(alph*E0)<<" "<<res[5]<<" "<<norm(p3nm)<<" "<<ewh<<" "<<kex1<<" "<<Isat<<endl;
        
        vkape.clear();
        }

    return 0;
    }
    
