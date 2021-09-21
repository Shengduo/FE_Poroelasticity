/** @file Node.hh
 * Declaration of class Node and class CohesiveNode
 */
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
using namespace std;

//------------------------------------------------------------------------------------------------------------------------------------------------------
/** Class Node, declares node
 * contains vector XYZ (2 or 3 dim)
 * node ID
 * nodal degree of freedom
 * and its match to the global one
 */
class Node {
// PRIVATE MEMBERS
protected:
    // Space dim
    int _spaceDim;

    // ID number of the node
    int _ID;

    // Nodal coordinates
    vector<double> _nodalXYZ;

    // Nodal variables existence, false OR true, 
    // [displacement, velocity, pressure, trace_strain]
    // vector<string> *_varNames;

    // Nodal degree of freedom, 0 free, 1 fixed
    // displacement(spaceDim), velocity(spaceDim), pressure(1) 
    vector<int> _nodalDOF;

    // Nodal values (0 - mass density; 1,2 - body force (force per unit volume); ...)
    vector<double> _nodalProperties;
    
// PUBLIC MEMBERS
public:
    /** TO DO 
     * Change this nodalBodyForce to private, currently for debugging purpose..
     */
    vector<double> _nodalBodyForce;
    
    // Default Constructor
    Node();

    // Constructor 1
    Node(int ID, int spaceDim = 2);

    // Constructor 2
    Node(int ID, 
         const vector<double> & XYZ, 
         const vector<int> & DOF, 
         int spaceDim = 2, 
         double density = 1., 
         const vector<double> *bodyForce = NULL);

    // Destructor
    ~Node();

    // Set spaceDim
    void setSpaceDim(int spaceDim);

    // Set ID
    void setID(int ID);

    // Set nodal coordinates XYZ
    void setXYZ(const vector<double> & XYZ);

    // Set DOF 1
    void setDOF(const vector<int> & DOF);

    // Set DOF 2
    void setDOF(int index, int DOF);

    // Set mass density
    void setMassDensity(double density);

    // Set body force
    void setBodyForce(const vector<double> *bodyForce);

    // Get spaceDim
    int getSpaceDim() const;

    // Get ID
    int getID() const;

    // Get XYZ
    const vector<double> & getXYZ() const;

    // Get DOF 1
    const vector<int> & getDOF() const;

    // Get DOF 2
    int getDOF(int index) const;

    // Get mass density
    double getMassDensity() const;

    // Get body force
    vector<double> getBodyForce() const;

    // Output nodal information to a file
    void outputInfo(ofstream & myFile, bool outputElse = false) const;

// NOT IMPLEMENTED
private:
    Node(const Node &); ///< Copy constructor, not implemented
    
    const Node& operator=(const Node &); ///< Move-assign operator, not implemented
};

//------------------------------------------------------------------------------------------------------------------------------------------------------
/** Class CohesiveNode, 
 * deals with the cohesive nodes initialization and calculation.
 * has DOF: {lambda(spaceDim), pressureFault}
 */
class CohesiveNode : public Node {
// PUBLIC MEMBERS
public:
    // Constructor 1
    CohesiveNode(int ID, int spaceDim = 2);

    // Constructor 2
    CohesiveNode(int ID, const vector<double> & XYZ, const vector<int> & DOF, int spaceDim = 2);

// NOT IMPLEMENTED
private:
    // Copy constructor
    CohesiveNode(const CohesiveNode &);

    // Move-assign operator
    const CohesiveNode & operator=(const CohesiveNode &);
};