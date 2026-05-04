//This function makes a sweep over the frequency, mantaining fixed the gain and the electric field magnitude.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <string>
#include "../src/headers/math33.H"
#include "../src/headers/nanoshell.H"
#include "../src/headers/cup.H"

/*
g++ -Wall -I/usr/include/ -L/usr/local/lib electricfieldtest.cxx -o tim -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv) {
    // Initialize variables for input parameters, and create an instance of the nanosphere class
    double omemi, omema, E0, rho, *fro, tpump, T, eps3, eps_b;
    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    // You do not have to input anything
    nanosphere  simulation;
    simulation.init();
    
    fstream nano, time;

    nano.open("../data/input/nanosphere_eV.dat", ios::in);
    time.open("../data/input/time.dat", ios::in);
    
    // omemi and omema define an initial arbitrary frequency interval.
    // The border-check routine reduces this range to the smaller interval, where emission is detected.
    nano>>simulation.r1>>simulation.Dome>>simulation.ome_0>>simulation.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    time>>T>>tpump;  
        
    simulation.set_metal(mtl,mdl,1);
    simulation.set_active(active);
    eps3=simulation.set_host(sol);
    eps_b=simulation.set_host(hst);
    
    fro=simulation.frohlich(omemi, omema, eps_b, eps3, rho);

    // Calculate the saturation electric field Esat
    double ntau1, ntau2;
    ntau2 = 2./simulation.Dome;
    ntau1 = 5.*ntau2;

    double nEsat = sqrt(1/(fabs(simulation.G)*ntau1));

    // Set E0 such that the normalized electric field Enorm = E0 / Esat = 1e-9
    E0 = nEsat*(1e-9);
    //double nEsatup2 = 1/(fabs(simulation.G)*ntau1);

    // Border check
//----------------------------------
    int omeN=1000;
    int omi, ch=1, o1=0, o2=0, oB=0;
    double  omeeV, dome  = (omema-omemi)/omeN, Ome, Omem=0., ome;
    double ome1, ome2;
    
        dome  = (omema-omemi)/100000;
        for (omi=0; omi<=100000; omi++){
            omeeV = omemi + omi*dome;
            Omem = Ome;
            Ome = simulation.find_Omega(omeeV, hst, sol, rho);
            ome  = omeeV/simulation.Ome_p;
            if (o1==1 && Omem != 0) ch  = Omem*Ome/fabs(Omem*Ome);
            if (fabs(Ome)>1.e-30 && o1==0){
                ome1=ome;
                o1=1;
                }
            if (fabs(Ome)<1.e-30 && o2==0 && o1==1 && oB==1){
                ome2=ome-dome/simulation.Ome_p;
                o2=1;
                }
            if (ch < 0 && o1==1 && oB==0) {
                simulation.omeB     = ome;
                oB=1;
                }
            }
        if (fabs(simulation.G)<fabs(fro[1])) {
            ome1=simulation.omeB/simulation.Ome_p+dome;
            ome2=simulation.omeB/simulation.Ome_p-dome;
            }
        if (ome1*simulation.Ome_p==omemi || ome2 ==0){
            std::cout<<"Fatal! the frequency range is too small!"<<std::endl;
            exit(78);
            }
//----------------------------------

    // Inform the user about the test
    cout << "Calculating the time_behavior up to "<<T<<" ps\n";
    cout << "                switching the pump on at "<<tpump<<"ps ...\n\n";
    cout << "* Beware if  "<<simulation.G<<" is greater than "<<fro[1]<<"\n";
    cout << "* you can get exponential behaviors in the\n";
    cout << "* analytical results for some frequencies\n\n";

    cout << "Parameters:\n";
    cout << "  Metal model: " << mdl << "\n";
    cout << "  Metal type: " << mtl << "\n";
    cout << "  Core material: " << hst << "\n";
    cout << "  Solvent: " << sol << "\n";
    cout << "  Radius ratio: " << rho << "\n\n";

    // The first loop is to save the maximun value of the dipole magnitude
    cout << "Running sweep over the frequency...\n";
    simulation.freqtest_first(mdl, mtl, hst, nEsat, tpump, E0, ome1, ome2, sol, rho);

    // The second loop is to save the time it takes the dipole magnitude to reach 0.9 of the maximun (Switch-on time)
    cout << "Running final loop test...\n";
    simulation.freqtest_second(mdl, mtl, hst, nEsat, tpump, E0, ome1, ome2, sol, rho);
   
    // Save the results to a file, I did not really used this, so I did not change the file name
    ofstream output("results/time_behavior.log");
    if (output.is_open()) {
        output << "Steady-state polarizability calculation results:\n";
	output << "  Metal model: " << mdl << "\n";
        output << "  Metal type: " << mtl << "\n";
        output << "  Core material: " << hst << "\n";
        output << "  Time range: [0, " << T << "] ps\n";
        output << "  The pump is switched on at "<<tpump<<" ps\n";
        output << "  Solvent: " << sol << "\n";
        output << "  Radius ratio: " << rho << "\n";
        output << "\nCheck simulation output for detailed results.\n";
        output.close();
        cout.precision(10);
        cout.setf(ios::fixed);
        cout <<"> Output saved in the following files:"<<endl;
	    cout<<">   ../data/output/freqtime.dat"<<endl;
        cout<<"> log saved in results/time_behavior.log"<<endl;
    } else {
        cerr << "Error: Could not open file for writing results.\n";
    }

    return 0;
}
