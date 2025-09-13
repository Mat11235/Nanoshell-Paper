/*
 * Program: time_behavior
 * Purpose: Configure a nanosphere/nanoshell simulation and run a
 *          tau-characterization at a single frequency (omeeV).
 *
 * CLI usage:
 *   ./tim <omega_eV> [linewidth_eV]
 *     - omega_eV     : excitation/analysis frequency in electron-volts (required).
 *     - linewidth_eV : optional damping/linewidth to override the file value.
 *
 * Inputs (text files):
 *   ../data/input/nanosphere_eV.dat
 *     -> r1  Dome  ome_0  G  omemi  omema  mtl  mdl  active  sol  E0  rho  hst
 *   ../data/input/time.dat
 *     -> T  tpump
 *
 * Outputs (summary + whatever taucharacterization writes):
 *   results/time_behavior.log
 *   (Common output conventions in this project include:
 *     ../data/output/anltime.dat, ../data/output/anlfunc.dat,
 *     ../data/output/numtime.dat,  ../data/output/numfunc.dat)
 *
 * Notes:
 *   - The optional second argument lets you override "Dome" (linewidth).
 *   - If not provided, the Dome read from nanosphere_eV.dat is used.
 */

// Variable legend:
//   omeeV     : analysis frequency (eV), from CLI
//   linewidth : optional CLI linewidth to override Dome (eV); -1 means "no override"
//   omemi/omema: frequency bounds read from file (used by frohlich here)
//   E0        : field amplitude read from file (used by taucharacterization)
//   rho       : radius ratio / filling factor
//   T         : total time window (ps)
//   tpump     : pump-on time (ps)
//   eps3/eps_b: dielectric constants for solvent/core returned by set_host
//   Dome      : linewidth/damping parameter inside simulation (overridden if CLI is provided)
//   G         : gain/coupling parameter; compared against fro[1] for stability warning

// --- Parse command-line arguments ------------------------------------------------
// argv[1]: required frequency in eV (omeeV)
// argv[2]: optional linewidth in eV (to override the file's Dome)
// Format: ./tim 2.8122 0.05 the last parameter is the linewidth

#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <string>
#include "../src/headers/math33.H"
#include "../src/headers/nanoshell.H"
#include "../src/headers/cup.H"

using namespace std;

int main(int argc, char** argv) {
    // Initialize variables for input parameters
    // Create an instance of the nanosphere class
	
    double   omeeV, linewidth, omemi, omema, E0, rho, *fro, tpump, T, eps3, eps_b;
    char mtl[16], mdl[16], hst[16], sol[16], active[16];

	
    if (argv[1]==0){
        cout<<endl<<"  Usage: "<<argv[0]<<" <omega in eV>"<<endl<<endl;
        exit(0);
        }
	
    omeeV=atof(argv[1]);
	linewidth = -1.0;                  // if not given, a default value is taken
	if (argc >= 3) linewidth = atof(argv[2]);  // Second argument optional in eV
    nanosphere  simulation;
    simulation.init();
    
    fstream nano, time;

    nano.open("../data/input/nanosphere_eV.dat", ios::in);
    time.open("../data/input/time.dat", ios::in);
    
    nano>>simulation.r1>>simulation.Dome>>simulation.ome_0>>simulation.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    time>>T>>tpump;  

	simulation.Dome = linewidth;
        
    simulation.set_metal(mtl,mdl,1);
    simulation.set_active(active);
    eps3=simulation.set_host(sol);
    eps_b=simulation.set_host(hst);
    
    // Perform the time_behavior calculation
 
    fro=simulation.frohlich(omemi, omema, eps_b, eps3, rho);
    // Calculate of tau, but this is done inside the function tau_carachterization on cup.h
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

    cout << "Running tau charachterization calculation...\n";
    simulation.taucharacterization(mdl, mtl, hst, E0, omeeV, T, tpump, sol, rho);
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
