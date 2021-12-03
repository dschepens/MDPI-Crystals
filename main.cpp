#include <iostream>
#include "Functions.h"

using std::cout;
using std::endl;

int main(){
    
    // set grid cube size
    // smaller cube sizes give more accurate results and take longer to run
    double cube = 0.1;

    // read molecule files
    Molecule m1 = read_molecule("/file_path_to_molecule1.pdb");
    Molecule m2 = read_molecule("/file_path_to_molecule2.pdb");

    // calculate volumes
    double m1_m2_overlap_vol = get_overlap_volume(m1, m2, cube); // calculate overlap volume of 2 molecules
    double m1_vol = get_overlap_volume(m1, m1, cube); // calculate volume of molecule 1 by overlapping with itself
    double m2_vol = get_overlap_volume(m2, m2, cube); // calculate volume of molecule 2 by overlapping with itself
    
    // output individual volumes
    cout << "Overlap vol: " << m1_m2_overlap_vol << endl;
    cout << "Molecule 1 vol: " << m1_vol << endl;
    cout << "Molecule 2 vol: " << m2_vol << endl;
    
    // output overlap percentage
    cout << "Overlap vol % relative to Molecule 1: " << m1_m2_overlap_vol/m1_vol << endl;
    cout << "Overlap vol % relative to Molecule 1: " << m1_m2_overlap_vol/m2_vol << endl;

    
    return 0;
}
