// Atom object functions

#include "atom.h"

// Constructor
Atom::Atom(std::string n, double x_c, double y_c, double z_c, double rad){
    name = n;
    x = x_c;
    y = y_c;
    z = z_c;
    r = rad;
}

// get function for name
std::string Atom::get_name(){
    return name;
}

// get function for x coordinate
double Atom::get_x(){
    return x;
}

// get function for y coordinate
double Atom::get_y(){
    return y;
}

// get function for z coordinate
double Atom::get_z(){
    return z;
}

// get function for radius
double Atom::get_r(){
    return r;
}
