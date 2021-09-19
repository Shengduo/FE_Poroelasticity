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

// PUBLIC MEMBERS
public:
    // Default Constructor
    Node();

    // Constructor 1
    Node(int ID, int spaceDim = 2);

    // Constructor 2
    Node(int ID, const vector<double> & XYZ, const vector<int> & DOF, int spaceDim = 2);

    // Destructor
    ~Node();

    // Set spaceDim
    void setSpaceDim(int spaceDim);

    // Set ID
    void setID(int ID);

    // Set nodal coordinates XYZ
    void setXYZ(const vector<double> & XYZ);

    // Set DOF
    void setDOF(const vector<int> & DOF);

    // Get spaceDim
    int getSpaceDim() const;

    // Get ID
    int getID() const;

    // Get XYZ
    const vector<double> & getXYZ() const;

    // Get DOF
    const vector<int> & getDOF() const;

    // Get the variables;
    // vector<string> getVarNames() const;

    // Output nodal information to a file
    void outputInfo(ofstream & myFile, bool outputDOF = false) const;

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