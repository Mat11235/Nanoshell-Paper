//This function makes a sweep over the electric field (probefield) amplitude, mantaining fixed the gain and the frequency.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <string>
#include "../src/headers/math33.H"
#include "../src/headers/nanoshell.H"
#include "../src/headers/cup.H"

/*
g++ -Wall -I../src/headers -L../src/lib electric_field_test.cxx -o electric_field_test -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv) {
    // Initialize variables for input parameters, and create an instance of the nanosphere class
    double   omeeV, omemi, omema, E0, rho, *fro, tpump, T, eps3, eps_b;
    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    // You have to input the frequency (this is fixed)
    if (argv[1]==0){
        cout<<endl<<"  Usage: "<<argv[0]<<" <omega in eV>"<<endl<<endl;
        exit(0);
        }
    omeeV=atof(argv[1]);
    nanosphere  simulation;
    simulation.init();
    
    fstream nano, time;
    
    nano.open("../data/input/nanosphere_eV.dat", ios::in);
    time.open("../data/input/time.dat", ios::in);
    // Even if the variable E0 is initialized, this value is not used, then is changed during the sweep
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
    //double nEsatup2 = 1/(fabs(simulation.G)*ntau1);

	
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
    cout << "Running sweep over the probefield amplitude...\n";
    simulation.probetest_first(mdl, mtl, hst, nEsat, omeeV, tpump, sol, rho);

    // The second loop is to save the time it takes the dipole magnitude to reach 0.9 of the maximun (Switch-on time)
    cout << "Running final loop test...\n";
    simulation.probetest_second(mdl, mtl, hst, nEsat, omeeV, tpump, sol, rho);

    // Save the results to a file, I did not really used this, so I did not change the file name
    ofstream output("results/time_behavior.log");
    if (output.is_open()) {
        output << "Steady-state polarizability calculation results:\n";
	output << "  Frequency: " << omeeV << " eV\n";
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
	    cout<<">   ../data/output/probetime.dat"<<endl;
        cout<<"> log saved in results/time_behavior.log"<<endl;
    } else {
        cerr << "Error: Could not open file for writing results.\n";
    }

    return 0;
}
