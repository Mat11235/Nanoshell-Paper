/*
 * Driver: time_behavior
 * Goal : Configure a nanoshell/nanosphere system and run a time-domain
 *        simulation at a single frequency (omeeV), with optional E0 override.
 *
 * Usage:
 *   ./tim <omega_eV> [E0_override]
 *     - omega_eV   : driving/analysis frequency in eV (required)
 *     - E0_override: sets E0 from CLI (same units as file E0); if omitted,
 *                    uses the value from nanosphere_eV.dat
 *
 * Inputs:
 *   ../data/input/nanosphere_eV.dat
 *     r1 Dome ome_0 G omemi omema mtl mdl active sol E0 rho hst
 *   ../data/input/time.dat
 *     T tpump
 *
 * Outputs (summary + whatever numericalOEF writes internally):
 *   results/time_behavior.log
 *   (Project-convention files: ../data/output/anltime.dat, anlfunc.dat,
 *    numtime.dat, numfunc.dat — adjust message if your routine writes other names.)
 */

// --- CLI args --------------------------------------------------------------
// argv[1] : omeeV (eV) — required driving frequency
// argv[2] : EfA (E0 override) — optional; if absent, keep E0 from file
// NOTE: for production, prefer `if (argc < 2)` over `if (argv[1]==0)` to guard input.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <string>
#include "../src/headers/math33.H"
#include "../src/headers/nanoshell.H"
#include "../src/headers/cup.H"

/*
g++ -Wall -I/usr/include/ -L/usr/local/lib time_behavior.cxx -o tim -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv) {
    // Initialize variables for input parameters
    // Create an instance of the nanosphere class
    double   omeeV, EfA, omemi, omema, E0, rho, *fro, tpump, T, eps3, eps_b;
    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    if (argv[1]==0){
        cout<<endl<<"  Usage: "<<argv[0]<<" <omega in eV>"<<endl<<endl;
        exit(0);
        }
    omeeV=atof(argv[1]);
	EfA = 1;                  // valor “no seteado” por defecto (Electric field amplitude)
	if (argc >= 3) EfA = atof(argv[2]);  // segundo argumento opcional in eV 
    nanosphere  simulation;
    simulation.init();
    
    fstream nano, time;

    nano.open("../data/input/nanosphere_eV.dat", ios::in);
    time.open("../data/input/time.dat", ios::in);
    
    nano>>simulation.r1>>simulation.Dome>>simulation.ome_0>>simulation.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    time>>T>>tpump;  

	E0=EfA;
        
    simulation.set_metal(mtl,mdl,1);
    simulation.set_active(active);
    eps3=simulation.set_host(sol);
    eps_b=simulation.set_host(hst);
    
    // Perform the time_behavior calculation

 
    fro=simulation.frohlich(omemi, omema, eps_b, eps3, rho);
    // Calculate the saturation electric field Esat
    //double ntau1, ntau2;
    //ntau2 = 2./simulation.Dome;
    //ntau1 = 5.*ntau2;

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


    cout << "Running calculation...\n";
    simulation.numericalOEF(mdl, mtl, hst, E0, omeeV, T, tpump, sol, rho);
    // Output the results

    // Save the results to a file
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
	    cout<<">   ../data/output/anltime.dat"<<endl;
	    cout<<">   ../data/output/anlfunc.dat"<<endl;
	    cout<<">   ../data/output/numtime.dat"<<endl;
	    cout<<">   ../data/output/numfunc.dat"<<endl;
        cout<<"> log saved in results/time_behavior.log"<<endl;
    } else {
        cerr << "Error: Could not open file for writing results.\n";
    }

    return 0;
}
