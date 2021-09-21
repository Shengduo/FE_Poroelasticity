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
    
    // Nodal body force reset
    _nodalBodyForce.resize(_spaceDim, 0.);
};

// Constructor 2
Node::Node(int ID, 
           const vector<double> & XYZ, 
           const vector<int> & DOF, 
           int spaceDim, 
           double density, 
           const vector<double>* bodyForce) {
    _ID = ID;
    _spaceDim = spaceDim;
    _nodalXYZ.resize(_spaceDim, 0.);

    // Default explicit, without trace_strain
    _nodalDOF.resize(2 * _spaceDim + 2, 0);

    // Nodal body force reset
    _nodalBodyForce.resize(_spaceDim, 0.);

    // Set nodal coordinates and degree of freedom
    this->setXYZ(XYZ);
    this->setDOF(DOF);
    this->setMassDensity(density);
    this->setBodyForce(bodyForce);
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

// Set DOF 1
void Node::setDOF(const vector<int> & DOF) {
    if (DOF.size() != _nodalDOF.size()) {
        throw("Input DOF incompatible with nodal DOF!");
    }
    for(int i = 0; i < _nodalDOF.size(); i++) _nodalDOF[i] = DOF[i];
};

// Set DOF 2
void Node::setDOF(int index, int DOF) {
    if (index >= _nodalDOF.size()) throw("Cannot run setDOF2 due to overflow!");
    _nodalDOF[index] = DOF;
};

// Get DOF 1
const vector<int> & Node::getDOF() const {
    return _nodalDOF;
};

// Get DOF 2
int Node::getDOF(int index) const {
    return _nodalDOF[index];
};

// Set mass density
void Node::setMassDensity(double density) {
    if (_nodalProperties.size() < 1) _nodalProperties.resize(1);
    _nodalProperties[0] = density;
};

// Get mass density
double Node::getMassDensity() const {
    if (_nodalProperties.size() < 1) throw ("Mass density not initialized!");
    return _nodalProperties[0];
}

// Set body force
void Node::setBodyForce(const vector<double> *bodyForce) {
    if (_nodalProperties.size() < 1 + _spaceDim) _nodalProperties.resize(1 + _spaceDim);
    if (bodyForce) {
        for (int i = 0; i < _spaceDim; i++) {
            _nodalProperties[1 + i] = (*bodyForce)[i];
        }
    }
    else {
        for (int i = 0; i < _spaceDim; i++) {
            _nodalProperties[1 + i] = 0.;
        }
    }
};

// Get body force
vector<double> Node::getBodyForce() const {
    if (_nodalProperties.size() < 1 + _spaceDim) throw("Body force not initialized!");
    vector<double> res(_spaceDim, 0.);
    for (int i = 0; i < _spaceDim; i++) {
        res[i] = _nodalProperties[1 + i];
    }
    return res;
};

// Output nodal information to a file
void Node::outputInfo(ofstream & myFile, bool outputElse) const {
    myFile << setw(10) << _ID << " ";
    for (int i = 0; i < _nodalXYZ.size(); i++) myFile << setw(12) << _nodalXYZ[i] << " ";
    if (outputElse) {
        for (int i = 0; i < _nodalDOF.size(); i++) myFile << setw(10) << _nodalDOF[i] << " ";
        for (int i = 0; i < _nodalProperties.size(); i++) 
            myFile << setw(12) << _nodalProperties[i] << " ";
        for (int i = 0; i < _nodalBodyForce.size(); i++)
            myFile << setw(12) << _nodalBodyForce[i] << " ";
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

    // Set nodal body force to 0.
    _nodalBodyForce.resize(_spaceDim, 0.);
};

// Constructor 2
CohesiveNode::CohesiveNode(int ID, const vector<double> & XYZ, const vector<int> & DOF, int spaceDim) {
    setSpaceDim(spaceDim);
    setID(ID);
    setSpaceDim(spaceDim);
    setXYZ(XYZ);
    _nodalDOF.resize(spaceDim + 1);
    setDOF(DOF);
    _nodalBodyForce.resize(_spaceDim, 0.);
};
