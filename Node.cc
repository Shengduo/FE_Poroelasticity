#include "Node.hh"
using namespace std;

/** ------------------------------------------------------------------------------------------------------------------------------------------------------
 * Methods of class Node
 * With default constructor
 */

// Default Constructor
Node::Node() {};

// Constructor 1
Node::Node(int ID, int spaceDim) {
    _ID = ID;
    _spaceDim = spaceDim;
    _nodalXYZ.resize(_spaceDim, 0.);

    // Default explicit, without trace_strain
    _nodalDOF.resize(2 * _spaceDim + 2, 0);
};

// Constructor 2
Node::Node(int ID, const vector<double> & XYZ, const vector<int> & DOF, int spaceDim) {
    _ID = ID;
    _spaceDim = spaceDim;
    _nodalXYZ.resize(_spaceDim, 0.);

    // Default explicit, without trace_strain
    _nodalDOF.resize(2 * _spaceDim + 2, 0);

    // Set nodal coordinates and degree of freedom
    this->setXYZ(XYZ);
    this->setDOF(DOF);
};

// Destructor
Node::~Node(){};

// Set spaceDim
void Node::setSpaceDim(int spaceDim) {
    _spaceDim = spaceDim;
};

// Get spaceDim
int Node::getSpaceDim() const {
    return _spaceDim;
}

// Set ID
void Node::setID(int ID) {
    _ID = ID;
};

// Get ID
int Node::getID() const {
    return _ID;
};

// Set nodal coordinates
void Node::setXYZ(const vector<double> & XYZ) {
    if (XYZ.size() != _spaceDim) {
        throw("Input nodal coordinates are not compatible with space dimension!");
    }
    _nodalXYZ.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; i++) _nodalXYZ[i] = XYZ[i];
};

// Get the coordinates of a node
const vector<double> & Node::getXYZ() const {
    return _nodalXYZ;
};

// Set DOF
void Node::setDOF(const vector<int> & DOF) {
    if (DOF.size() != _nodalDOF.size()) {
        throw("Input DOF incompatible with nodal DOF!");
    }
    for(int i = 0; i < _nodalDOF.size(); i++) _nodalDOF[i] = DOF[i];
};

// Get DOF
const vector<int> & Node::getDOF() const {
    return _nodalDOF;
};

// Output nodal information to a file
void Node::outputInfo(ofstream & myFile, bool outputDOF) const {
    myFile << setw(10) << _ID << " ";
    for (int i = 0; i < _nodalXYZ.size(); i++) myFile << setw(12) << _nodalXYZ[i] << " ";
    if (outputDOF) {
        for (int i = 0; i < _nodalDOF.size(); i++) myFile << setw(10) << _nodalDOF[i] << " ";
    }
    myFile << "\n";
};


/** ------------------------------------------------------------------------------------------------------------------------------------------------------
 * Methods of class CohesiveNode
 * No default constructor
 */

// Constructor 1 for CohesiveNode
CohesiveNode::CohesiveNode(int ID, int spaceDim) {
    setID(ID);
    setSpaceDim(spaceDim);
    _nodalXYZ.resize(_spaceDim, 0.);

    // {Lagrange multiplier, fault_pressure}
    _nodalDOF.resize(_spaceDim + 1, 0);
};

// Constructor 2
CohesiveNode::CohesiveNode(int ID, const vector<double> & XYZ, const vector<int> & DOF, int spaceDim) {
    setSpaceDim(spaceDim);
    setID(ID);
    setSpaceDim(spaceDim);
    setXYZ(XYZ);
    _nodalDOF.resize(spaceDim + 1);
    setDOF(DOF);
};
