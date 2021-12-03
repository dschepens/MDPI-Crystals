// All functions needed to calculate the overlap volumes

#ifndef Functions_h
#define Functions_h

#include "atom.h"
#include "molecule.h"



// ------------------------------------------------------------- //
// functions needed for the get_overlapping_points_atoms function

bool atoms_contained(double x, double y, double z, std::vector<Atom> atoms){
    for (int i = 0; i < atoms.size(); i++) {
        double d = sqrt(pow(x - atoms[i].get_x(), 2) + pow(y - atoms[i].get_y(), 2) + pow(z - atoms[i].get_z(), 2));
        if (d > atoms[i].get_r()) {
            return false;
        }
    }
    return true;
}

std::vector<double> minmax(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double cube_size){
    std::vector<double> minmax;
    
    minmax.push_back(floor(xmin) + cube_size/2);
    minmax.push_back(ceil(xmax));
    minmax.push_back(floor(ymin) + cube_size/2);
    minmax.push_back(ceil(ymax));
    minmax.push_back(floor(zmin) + cube_size/2);
    minmax.push_back(ceil(zmax));

    return minmax;
}

// ------------------------------------------------------------- //





// ------------------------------------------------------------- //

// finding all overlapping points/cubes in 2 Atoms that overlap
// - two conditions exist for Atom overlaps:

// -- first: a smaller Atom is inside another bigger atom
// --- in this condition, the grid is formed around the smaller Atom
// --- each point that is inside the smaller Atom is added to the list

// -- second: there is a general overlap of 2 Atoms
// --- in this condition, the grid is formed around the overlap area
// --- the overlap area is found through calculating Xmin, Xmax, Ymin, Ymax, Zmin, Zmax
// --- simple geometry is used to calculate the points for the grid
// --- the links for references are provided down below
std::set<std::tuple<double, double, double>> get_overlapping_points_atoms(Atom atom1, Atom atom2, double cube_size){
    
    std::set<std::tuple<double, double, double>> points;
    
    double x1 = atom1.get_x();
    double x2 = atom2.get_x();
    
    double y1 = atom1.get_y();
    double y2 = atom2.get_y();
    
    double z1 = atom1.get_z();
    double z2 = atom2.get_z();
    
    double r1 = atom1.get_r();
    double r2 = atom2.get_r();
    
    double d = sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
    
    
    
    // first condition
    // - if one sphere is inside another sphere
    // - the grid is formed around the smaller sphere
    // - all points from smaller sphere are added to the set of overlapping points
    if(d <= std::max(r1,r2) - std::min(r1,r2) || d < std::max(r1,r2)){
        
        // if sphere 1 is smaller than sphere 2
        if(r1<=r2){
            std::vector<double> min_max = minmax(x1-r1, x1+r1, y1-r1, y1+r1, z1-r1, z1+r1, cube_size);
            
            for(double n = min_max[0]; n <= min_max[1]; n += cube_size){
                for(double m = min_max[2]; m <= min_max[3]; m += cube_size){
                    for(double l = min_max[4]; l <= min_max[5]; l += cube_size){
                        n = round(n*100)/100;
                        m = round(m*100)/100;
                        l = round(l*100)/100;
                        if(pow(n - x1,2) + pow(m - y1,2) + pow(l - z1,2) < pow(r1,2) && pow(n - x2, 2) + pow(m - y2, 2) + pow(l - z2,2) < pow(r2, 2)){
                            std::tuple<double, double, double> o1 (n,m,l);
                            points.insert(o1);
                        }
                    }
                }
            }
        }
        // if sphere 2 is smaller than sphere 1
        else{
            std::vector<double> min_max = minmax(x2-r2, x2+r2, y2-r2, y2+r2, z2-r2, z2+r2, cube_size);
            
            for(double n = min_max[0]; n <= min_max[1]; n += cube_size){
                for(double m = min_max[2]; m <= min_max[3]; m += cube_size){
                    for(double l = min_max[4]; l <= min_max[5]; l += cube_size){
                        n = round(n*100)/100;
                        m = round(m*100)/100;
                        l = round(l*100)/100;
                        if(pow(n - x1,2) + pow(m - y1,2) + pow(l - z1,2) < pow(r1,2) && pow(n - x2, 2) + pow(m - y2, 2) + pow(l - z2,2) < pow(r2, 2)){
                            std::tuple<double, double, double> o1 (n,m,l);
                            points.insert(o1);
                        }
                    }
                }
            }
        }
    }
    
    
    // second condition
    // - two spheres overlap
    // - the grid is formed around Xmin, Xmax, Ymin, Ymax, Zmin, Zmax points of overlap
    // - the formulas and code for finding Xmin, Xmax, Ymin, Ymax have been borrowed and modified from
    // -- https://mathworld.wolfram.com/Circle-CircleIntersection.html
    // -- https://github.com/benfred/bens-blog-code/blob/master/circle-intersection/circle-intersection.js
    // - the formulas for finding Zmin and Zmax have been borrowed from
    // -- http://www.ambrsoft.com/TrigoCalc/Sphere/TwoSpheres/Intersection.htm
    // please refer to the original links for a detailed explanation
    else
    {
        std::vector<Atom> atoms;
        atoms.push_back(atom1);
        atoms.push_back(atom2);
        
        // finding Xmax, Xmin, Ymin, Ymax
        double rad = (r1 * r1 - r2 * r2 + d * d) / (2 * d);
        double h = sqrt(r1 * r1 - rad * rad);
        double x0 = x1 + (rad) * (x2 - x1) / d;
        double y0 = y1 + (rad) * (y2 - y1) / d;
        double rx = -(x2 - x1) * (h / d);
        double ry = -(y2 - y1) * (h / d);
        
        double xmin = std::min(x0+rx, x0-rx);
        double xmax = std::max(x0+rx, x0-rx);
        double ymin = std::min(y0+ry, y0-ry);
        double ymax = std::max(y0+ry, y0-ry);
        
        // making sure that Xmax, Xmin, Ymin, Ymax are contained in the spheres
        for(int i = 0; i < 2; i++){
            Atom a = atoms[i];
            
            if ((a.get_x() - a.get_r() < xmin) && atoms_contained(a.get_x()-a.get_r(), a.get_y(), a.get_z(), atoms)) {
                xmin = a.get_x() - a.get_r();
            }
            if ((a.get_x() + a.get_r() > xmax) && atoms_contained(a.get_x()+a.get_r(), a.get_y(), a.get_z(), atoms)) {
                xmax = a.get_x() + a.get_r();
            }
            if ((a.get_y() - a.get_r() < ymin) && atoms_contained(a.get_x(), a.get_y()-a.get_r(), a.get_z(), atoms)) {
                ymin = a.get_y() - a.get_r();
            }
            if ((a.get_y() + a.get_r() > xmax) && atoms_contained(a.get_x(), a.get_y()+a.get_r(), a.get_z(), atoms)) {
                ymax = a.get_y() + a.get_r();
            }
        }
        
        
        
        // finding Zmin, Zmax
        double A = 2 * (x2-x1);
        double B = 2 * (y2-y1);
        double C = 2 * (z2-z1);
        double De = x1*x1 - x2*x2 + y1*y1 - y2*y2 + z1*z1 - z2*z2 - r1*r1 +r2*r2;
        double t = (x1*A + y1*B + z1*C + De)/(A*(x1-x2) + B*(y1-y2) + C*(z1-z2));
        double z_center = z1 + t*(z2-z1);
        
        double zmin = z_center - std::min(r1,r2);
        double zmax = z_center + std::min(r1,r2);
                
        std::vector<double> min_max = minmax(xmin, xmax, ymin, ymax, zmin, zmax, cube_size);
        
        // iterating through the grid and finding points of overlap
        // a point is added to the set only if the distances from the point to the centers of the spheres
        // are smaller than the radii
        for(double n = min_max[0]; n <= min_max[1]; n += cube_size){
            for(double m = min_max[2]; m <= min_max[3]; m += cube_size){
                for(double l = min_max[4]; l <= min_max[5]; l += cube_size){
                    n = round(n*100)/100;
                    m = round(m*100)/100;
                    l = round(l*100)/100;
                    if(pow(n - x1,2) + pow(m - y1,2) + pow(l - z1,2) < pow(r1,2) && pow(n - x2, 2) + pow(m - y2, 2) + pow(l - z2,2) < pow(r2, 2)){
                        std::tuple<double, double, double> o1 (n,m,l);
                        points.insert(o1);
                    }
                }
            }
        }
    }
    
    
    return points;
}

// ------------------------------------------------------------- //





// ------------------------------------------------------------- //

// find which Atoms overlap in 2 Molecules

// - considers each Atom in Molecule 1 with all the Atoms in Molecule 2
// - overlapping Atoms are stored in pairs in a vector
// - calculates the distance between the centers of 2 Atoms that are being compared
// - if the distance between centers are smalled than 2 radii combined, then the atoms overlap
// - returns the vector of pairs of Atoms that overlap
std::vector<std::pair<int,int>> get_overlapping_atoms(Molecule molecule1, Molecule molecule2){
    std::vector<std::pair<int,int>> overlaps; // create vector of overlapping Atoms
    
    // get all atoms from 2 molecules
    std::vector<Atom> m1 = molecule1.get_atoms();
    std::vector<Atom> m2 = molecule2.get_atoms();
    
    // compare each Atom in both Molecules with each other
    for(int i = 0; i < m1.size(); i++){
        for(int j = 0; j < m2.size(); j++){
            // calculate the distance between the centers of 2 Atoms
            double d = sqrt(pow(m1[i].get_x()-m2[j].get_x(),2) + pow(m1[i].get_y()-m2[j].get_y(),2) + pow(m1[i].get_z()-m2[j].get_z(),2));
            
            // - if the distance between centers are smalled than 2 radii combined, then the atoms overlap
            if(d < m1[i].get_r() + m2[j].get_r()){
                std::pair<int,int> o1 (i,j);
                overlaps.push_back(o1);
                 
            }
        }
    }
    
    return overlaps;
}
// ------------------------------------------------------------- //



// ------------------------------------------------------------- //

// finding overlapping points in 2 Molecule objects

// - uses get_overlapping_points_atoms function to find overlapping points in pairs of Atoms
// - the pairs of Atoms that overlap are known
// - once overlapping points in pairs of Atoms are found, the points are added to the larger list of overlapping points in 2 Molecules
std::set<std::tuple<double,double,double>> get_overlapping_points(Molecule molecule1, Molecule molecule2, std::vector<std::pair<int,int>> molecule_ovlps, double cube_size){
    std::set<std::tuple<double, double, double>> points;
    std::vector<Atom> m1 = molecule1.get_atoms();
    std::vector<Atom> m2 = molecule2.get_atoms();
    
    for(int i = 0; i < molecule_ovlps.size(); i++){
        Atom atom1 = m1[std::get<0>(molecule_ovlps[i])];
        Atom atom2 = m2[std::get<1>(molecule_ovlps[i])];
        
        // after finding all overlapping points in 2 Atom objects, add the points to the general set of all overlapping points in 2 Molecule Objects
        std::set<std::tuple<double,double,double>> ov = get_overlapping_points_atoms(atom1, atom2, cube_size);
        points.insert(ov.begin(), ov.end());
    }
    return points;
}

// ------------------------------------------------------------- //






// ------------------------------------------------------------- //

// get overlap volume from 2 Molecule objects and defined grid cube size

// - uses get_overlapping_atoms function to find all atoms that overlap
// - uses molecule_points function to find all points in the grid that are overlapping in the 2 Molecule objects
// - all overlapping points are stored in a set which does not add duplicates
// - hence, if an overlapping point was added before, it will not be added again
// - the volume is then calculated by calculating the volumes of all cubes/points that are overlapping
double get_overlap_volume(Molecule molecule1, Molecule molecule2, double cube){
    std::vector<std::pair<int,int>> atom_overlaps = get_overlapping_atoms(molecule1, molecule2); //identify which atoms overlap
    std::set<std::tuple<double,double,double>> molecule_points = get_overlapping_points(molecule1, molecule2, atom_overlaps, cube); // find all overlapping points
    double molecule_volume = molecule_points.size() * pow(cube,3); // calculate volume
    
    return molecule_volume;
}
// ------------------------------------------------------------- //






// ------------------------------------------------------------- //

// read molecule from pdb file

// - reads pdb file as a text file and store each line into a string
// - gets each atom name, x, y, z coordinates from each line
// - radii are predefined based on atom name
Molecule read_molecule(std::string path){
    // initialize empty molecule
    Molecule m;
    
    // open pdb file and store each line into a string then close file
    std::ifstream f (path);
    std::string line;
    std::string t;
    while (getline(f,line))
    {
        t+=line;
    }
    f.close();
    
    
    // define starting positions for numbers, name, x, y, z
    int num_pos = 10;
    int n_pos = 13;
    int x_pos = 32;
    int y_pos = 40;
    int z_pos = 47;
    
    // define string for name, x, y, and z
    // string x, y, z will be converted to double when initializing Atom objects
    std::string n;
    std::string x;
    std::string y;
    std::string z;
    double rad = 0.0;
    
    // initialize increments
    int inc = 67;
    if(t[81] == 'A')
        inc = 81;
    
    
    int i = 13;
    
    // the while loop runs till getting to Atom connections in PDB files
    // Atom connections are not needed
    while(t[num_pos] != ' '){
        if(i == n_pos){
            n_pos += inc;
            n = t[i];
            n += t[i+1];
            i+=19;
        }
        if(i == x_pos){
            for(int j = i; j < i+7; j++){
                x += t[j];
            }
            x_pos += inc;
            i+=8;
        }
        if(i == y_pos){
            for(int j = i; j < i+7; j++){
                y += t[j];
            }
            y_pos += inc;
            i+=7;
        }
        if(i == z_pos){
            for(int j = i; j < i+7; j++){
                z += t[j];
            }
            z_pos += inc;
            
            if(n=="C ")
                rad = 1.7;
            else if(n=="F ")
                rad = 1.47;
            else if(n=="O ")
                rad = 1.52;
            else if(n=="H ")
                rad = 1.09;
            else if(n=="I ")
                rad = 1.98;
            else if(n=="Cl")
                rad = 1.75;
            else if(n=="Br")
                rad = 1.85;
            else if(n=="N ")
                rad = 1.55;
            
            m.add_atom(Atom(n, std::stod(x), std::stod(y), std::stod(z), rad));
            x = "";
            y = "";
            z = "";
            rad = 0;
            i+=inc-34;
        }
        num_pos += inc;
    }
    
    
    return m;
}
// ------------------------------------------------------------- //

#endif
