// Molecule object

// Characteristics:
// - vector of Atom objects

// Functions:
// - Standard get function for vector of Atom objects
// - Molecule size function which returns size of the vector
// - Add Atom object to the vector
// - Erase Atom from the vector knowing the index of the Atom object in the vector

#ifndef molecule_h
#define molecule_h

#include "atom.h"

class Molecule{
private:
    std::vector<Atom> atoms;
public:
    // Constructors
    Molecule();
    Molecule(std::vector<Atom> a);
    
    // Get functions
    std::vector<Atom> get_atoms();
    int size();
    
    void add_atom(Atom atom); // adds an Atom object to the vector of Atoms
    void erase_atom(int i); // erases an Atom object from the vector using the index
};

#endif /* molecule_h */
