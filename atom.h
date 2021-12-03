// Atom object

// Characteristics:
// - string name
// - double x, y, z coordinates in 3D space
// - double r radius

// Functions:
// - standard set of get functions


#ifndef atom_h
#define atom_h

#include <vector>
#include <cmath>
#include <tuple>
#include <set>
#include <fstream>

#include <stdio.h>

class Atom{
private:
    std::string name;
    double x;
    double y;
    double z;
    double r;
public:
    // Constructor
    Atom(std::string n, double x_c, double y_c, double z_c, double rad);
    
    // get functions
    std::string get_name();
    double get_x();
    double get_y();
    double get_z();
    double get_r();
};;

#endif
