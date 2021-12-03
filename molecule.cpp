// Molecule object 

#include <stdio.h>
#include "molecule.h"

// Empty Constructor
Molecule::Molecule() {}

// Constructor with predefined vector of Atoms
Molecule::Molecule(std::vector<Atom> a){
    atoms = a;
}

// get function for the vector of Atoms
std::vector<Atom> Molecule::get_atoms(){
    return atoms;
}

// get size of the vector of Atoms
int Molecule::size(){
    return atoms.size();
}

// add an Atom object to the vector of Atoms
void Molecule::add_atom(Atom atom){
    atoms.push_back(atom);
}


// erase an Atom object from the vector using the index
void Molecule::erase_atom(int i){
    atoms.erase(atoms.begin()+i-1);
}
