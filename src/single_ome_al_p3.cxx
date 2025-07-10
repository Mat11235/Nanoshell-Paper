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
#include "headers/single.H"
#include "headers/cup.H"
#include "headers/ns_ISS.H"

/*
g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/single_ome_al_p3.cxx -o ../bin/sap -lgsl -lgslcblas -lm -larmadillo
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
        
using namespace std;

int main(){
    double   omemi, omema, E0, eps_b, *fro, ome1, ome2, kex1, kex2;
    complex<double> eps1, eps2, alph, alph_num, alph_anl;
    int omeN = 10000;
    
    char mtl[16], mdl[16], sol[16], active[16];
    std::vector<std::pair<double,std::complex<double>>> valph;
    std::vector<std::pair<double,std::complex<double>>> vkape;
    std::vector<std::pair<double,double>> p3;
    std::vector<double> ralph;
    std::vector<double> ialph;
    std::vector<double> vome;
    
    std::pair<double,double> rzero;
    std::pair<double,double> izero;
    std::pair<double,double> kzero;
    
    nanosphere ns;
    ns.init();
    
    ifstream nano("../data/input/nanosphere_eV.dat");
    if (!nano) {
        cerr << "Error: Cannot open input file" << endl;
        return 1;
    }

    nano>>ns.r1>>ns.Dome>>ns.ome_0>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0;
    
    if (E0==0.) E0=1.e-30; // zero is problematic as a value for E0 
    
    ns.init();
    eps_b=ns.set_host(sol);
    ns.set_metal(mtl,mdl,1);
    ns.set_active(active);

    fro=ns.frohlich(omemi, omema, eps_b);
    
    valph = ns.steady_state(mdl, mtl, sol, omemi, omema, omeN);

    ralph = extract_ralph(valph);
    ialph = extract_ialph(valph);
    vome  = extract_ome(valph);

    rzero = find_zeros(vome, ralph);
    izero = find_zeros(vome, ialph);

    vkape = gimme_emi_kap(ns, mdl, mtl, sol, omemi, omema, omeN);    

    kzero = fnd_extrm(vkape, ns.Ome_p);
    
    kex1 = kzero.first;
    kex2 = kzero.second;
    
    if (izero.first==0 && fabs(izero.second-omema)<1.e4) {
        ome1 = 666*omema;
        ome2 = 777*omema;
        } else {
        ome1 = (rzero.first > izero.first) ? rzero.first : izero.first;
        ome2 = (rzero.second < izero.second) ? rzero.second : izero.second;
        }

    ofstream cop3("../data/output/oGp/ome_p3.dat");
    if (!cop3) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }
    ofstream coal("../data/output/oGp/ome_al.dat");
    if (!coal) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }
    ofstream cext("../data/output/oGp/ome_ex.dat");
    if (!cext) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }


    cext<<ome1<<" "<<ome2<<" "<<fro[0]<<" "<<kex1<<" "<<kex2<<endl;

    for (int ii=0; ii<int(vkape.size()); ii++)
        coal<<" "<<vkape[ii].first<<" "<<real(valph[ii].second)<<" "<<imag(valph[ii].second)<<
                                    " "<<real(vkape[ii].second)<<" "<<imag(vkape[ii].second)*ns.Ome_p<<endl;
    coal.close();
    
    return 0;
    }
    
