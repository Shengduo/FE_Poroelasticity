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
    
    // Resize solution s and s_t
    s.resize(_nodalDOF.size(), 0.);
    s_t.resize(_nodalDOF.size(), 0.);

    // Nodal body force reset
    _nodalBodyForce.resize(_spaceDim, 0.);
    
};

// Constructor 2
Node::Node(int ID,
           const vector<double> &XYZ,
           const vector<int> &DOF,
           int spaceDim,
           double density,
           const vector<double> *bodyForce,
           double lambda,
           double shearModulus,
           double biotAlpha,
           double biotMp,
           double fluidMobility,
           double fluidViscosity,
           double fluidDensity,
           double porosity,
           const vector<double> *fluidBodyForce,
           double source, 
           const vector<double> *traction) {
    _ID = ID;
    _spaceDim = spaceDim;
    _nodalXYZ.resize(_spaceDim, 0.);

    // Default explicit, without trace_strain
    _nodalDOF.resize(2 * _spaceDim + 2, 0);

    // Nodal body force reset
    _nodalBodyForce.resize(_spaceDim, 0.);

    // Reset _nodalProperties.
    _nodalProperties.resize(14, 0.);
    
    // Resize s and s_t
    s.resize(_nodalDOF.size(), 0.);
    s_t.resize(_nodalDOF.size(), 0.);

    // Set nodal coordinates and degree of freedom
    this->setXYZ(XYZ);
    this->setDOF(DOF);
    this->setMassDensity(density);
    this->setBodyForce(bodyForce);
    this->setLambda(lambda);
    this->setShearModulus(shearModulus);
    this->setBiotAlpha(biotAlpha);
    this->setBiotMp(biotMp);
    this->setFluidMobility(fluidMobility);
    this->setFluidViscosity(fluidViscosity);
    this->setFluidDensity(fluidDensity);
    this->setPorosity(porosity);
    this->setFluidBodyForce(fluidBodyForce);
    this->setSource(source);
    this->setTraction(traction);
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
        throw "Input nodal coordinates are not compatible with space dimension!";
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
    if (_nodalDOF.size() < 2 * _spaceDim + 2) _nodalDOF.resize(2 * _spaceDim + 2);
    if (DOF.size() != _nodalDOF.size()) {
        throw "Input DOF incompatible with nodal DOF!";
    }
    for(int i = 0; i < _nodalDOF.size(); i++) _nodalDOF[i] = DOF[i];
};

// Set DOF 2
void Node::setDOF(int index, int DOF) {
    if (index >= _nodalDOF.size()) throw "Cannot run setDOF2 due to overflow!";
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
    if (_nodalProperties.size() < 1) throw "Mass density not initialized!";
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
    if (_nodalProperties.size() < 1 + _spaceDim) throw "Body force not initialized!";
    vector<double> res(_spaceDim, 0.);
    for (int i = 0; i < _spaceDim; i++) {
        res[i] = _nodalProperties[1 + i];
    }
    return res;
};

/** Initialize s = [displacement (-+), velocity (-+), pressure (-+), trace_strain(-+), lambda, pressure_fault, theta] */
void Node::initializeS(const vector<double> & initialS) {
    if (initialS.size() != 2 * _spaceDim + 2) throw "Initial s vector not compatible with nodal DOFs!";
    
    s.resize(initialS.size());
    for (int i = 0; i < initialS.size(); i++) s[i] = initialS[i];
};

/** Push initial s to the global vector s */
void Node::pushS(Vec & globalS) const {
    for (int i = 0; i < _nodalDOF.size(); i++) {
        if (_nodalDOF[i] != -1) {
            VecSetValue(globalS, _nodalDOF[i], s[i], INSERT_VALUES);
        }
    }
};

/** Get current s from the global vector s */
void Node::fetchS(const Vec & globalS) {
    if (s.size() < _nodalDOF.size()) s.resize(_nodalDOF.size());
    for (int i = 0; i < _nodalDOF.size(); i++) {
        if (_nodalDOF[i] != -1) {
            VecGetValues(globalS, 1, &(_nodalDOF[i]), &(s[i]));
        }
    }
};

/** Get current s_t from the global vector s_t */
void Node::fetchS_t(const Vec & globalS_t) {
    if (s_t.size() < _nodalDOF.size()) s_t.resize(_nodalDOF.size());
    for (int i = 0; i < _nodalDOF.size(); i++) {
        if (_nodalDOF[i] != -1) {
            VecGetValues(globalS_t, 1, &(_nodalDOF[i]), &(s_t[i]));
        }
    }
};

/** Output nodal information to a file */
void Node::outputInfo(ofstream & myFile, bool outputElse) const {
    myFile <<"Node ID: " <<  setw(10) << _ID << " " << "Node XYZ: ";
    for (int i = 0; i < _nodalXYZ.size(); i++) myFile << setw(12) << _nodalXYZ[i] << " ";
    if (outputElse) {
        myFile << "Node DOF: ";
        for (int i = 0; i < _nodalDOF.size(); i++) myFile << setw(10) << _nodalDOF[i] << " ";
        myFile << "Node Properties: ";
        for (int i = 0; i < _nodalProperties.size(); i++) 
            myFile << setw(12) << _nodalProperties[i] << " ";
        myFile << "Node Forces: ";
        for (int i = 0; i < _nodalBodyForce.size(); i++)
            myFile << setw(12) << _nodalBodyForce[i] << " ";
        myFile << "S: ";
        for (int i = 0; i < s.size(); i++)
            myFile << setw(12) << s[i] << " "; 
        myFile << "S_t: ";
        for (int i = 0; i < s_t.size(); i++)
            myFile << setw(12) << s_t[i] << " "; 
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

    // {lambda, fault_pressure, theta}
    _nodalDOF.resize(2 * (2 * _spaceDim + 2) + _spaceDim + 3, 0);

    // Set nodal body force to 0.
    _nodalBodyForce.resize(_spaceDim, 0.);
};

// Constructor 2
CohesiveNode::CohesiveNode(int ID,
                           const vector<double> &XYZ,
                           const vector<int> &DOF,
                           const vector<Node*> &lowerUpperNodes, 
                           int spaceDim,
                           int density,
                           const vector<double> *bodyForce,
                           double fluidMobility_x, 
                           double fluidMobility_z,
                           double fluidViscosity,
                           double porosity,
                           double thickness,
                           double beta_p,
                           double beta_sigma,
                           double rateStateA,
                           double rateStateB,
                           double DRateState,
                           const vector<double> *fluidBodyForce,
                           double source, 
                           const vector<double> *initialSlip,
                           double V_reference, 
                           double f_reference) {
    setSpaceDim(spaceDim);
    setID(ID);
    setSpaceDim(spaceDim);
    setXYZ(XYZ);
    _nodalDOF.resize(2 * (2 * _spaceDim + 2) + spaceDim + 3);
    setDOF(DOF);
    setLowerUpperNodes(lowerUpperNodes);
    _nodalBodyForce.resize(_spaceDim, 0.);
    _nodalProperties.resize(20, 0.);
    
    // Resize solution s and s_t
    s.resize(_nodalDOF.size(), 0.);
    s_t.resize(_nodalDOF.size(), 0.);

    setMassDensity(density);
    setBodyForce(bodyForce);
    setFluidMobility_x(fluidMobility_x);
    setFluidMobility_z(fluidMobility_z);
    setFluidViscosity(fluidViscosity);
    setPorosity(porosity);
    setThickness(thickness);
    setBeta_p(beta_p);
    setBeta_sigma(beta_sigma);
    setRateStateA(rateStateA);
    setRateStateB(rateStateB);
    setDRateState(DRateState);
    setFluidBodyForce(fluidBodyForce);
    setSource(source);
    setSlip(initialSlip);
    setVReference(V_reference);
    setFReference(f_reference);
};

/** Initialize s = [u(-+), v(-+), p(-+), trace_strain(-+), lambda, fault_pressure, theta, slip_rate] */
void CohesiveNode::initializeS(const vector<double> & initialS) {
    if (initialS.size() != _spaceDim + 3) throw "Initial s vector not compatible with cohesive nodal DOFs!";
    s.resize(2 * (2 * _spaceDim + 2) + initialS.size());
    for (int i = 0; i < initialS.size(); i++) s[i + 2 * (2 * _spaceDim + 2)] = initialS[i];

    // For upper and lower nodes
    // u(-+)
    s[0] = _lowerUpperNodes[0]->s[0]; 
    s[1] = _lowerUpperNodes[0]->s[1]; 
    s[2] = _lowerUpperNodes[1]->s[0]; 
    s[3] = _lowerUpperNodes[1]->s[1]; 

    // v(-+)
    s[4] = _lowerUpperNodes[0]->s[2]; 
    s[5] = _lowerUpperNodes[0]->s[3]; 
    s[6] = _lowerUpperNodes[1]->s[2]; 
    s[7] = _lowerUpperNodes[1]->s[3]; 

    // p(-+)
    s[8] = _lowerUpperNodes[0]->s[4]; 
    s[9] = _lowerUpperNodes[1]->s[4]; 

    // trace_strain(-+)
    s[10] = _lowerUpperNodes[0]->s[5]; 
    s[11] = _lowerUpperNodes[1]->s[5]; 
};

/** Push initial s to the global vector s */
void CohesiveNode::pushS(Vec & globalS) const {
    for (int i = 2 * (2 * _spaceDim + 2); i < _nodalDOF.size(); i++) {
        if (_nodalDOF[i] != -1) {
            VecSetValue(globalS, _nodalDOF[i], s[i], INSERT_VALUES);
        }
    }
};