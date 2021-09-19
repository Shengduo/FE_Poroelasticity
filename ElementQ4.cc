/** @file ElementQ4.cc
 * Definition of class ElementQ4
 */
#include "ElementQ4.hh"

//------------------------------------------------------------------------------------------------------------------------------------------------------
/** Methods for class ElementQ4
 * Default Constructor implemented
 */
// Default Constructor
ElementQ4::ElementQ4() {};

// Constructor
ElementQ4::ElementQ4(int ID, const vector<Node*> & NID) {
    this->_ID = ID;
    this->_NID.resize(NID.size());
    for (int i = 0; i < NID.size(); i++) {
        this->_NID[i] = NID[i];
    }
};

// Destructor
ElementQ4::~ElementQ4(){
    // Delete vector<Node*> _NID
    for (int i = 0; i < _NID.size(); i++) {
        delete _NID[i];
    }
};

// Shape function N at (ksi, eta)
vector<double> ElementQ4::N(double ksi, double eta) {
    vector<double> N(4, 0.);
    N[0] = (1. - ksi) * (1. - eta) / 4.;
    N[1] = (1. + ksi) * (1. - eta) / 4.;
    N[2] = (1. + ksi) * (1. + eta) / 4.; 
    N[3] = (1. - ksi) * (1. + eta) / 4.;
    return N;
};

// Gradient of shape function B at (ksi, eta)
vector<double> ElementQ4::B(double ksi, double eta) {
    vector<double> B(8, 0.); 
    // [\partial N[0] / \partial \ksi, \partial N[0] / \partial \eta, \partial N[1] / \partial \ksi...]
    B[0] = - (1. - eta) / 4.;
    B[1] = - (1. - ksi) / 4.;
    B[2] = (1. - eta) / 4.;
    B[3] = - (1. + ksi) / 4.;
    B[4] = (1. + eta) / 4.;
    B[5] = (1. + ksi) / 4.;  
    B[6] = - (1. + eta) / 4.;
    B[7] = (1. - ksi) / 4.;
    return B;
};

// Jacobian at any given location in base space (ksi, eta), J = det(\partial (x,y) / \partial (ksi, eta))
double ElementQ4::J(double ksi, double eta) const {
    double matrix [4] = {0., 0., 0., 0.};    
    for (int i = 0; i < 4; i++) {
        matrix[0] += _NID[i]->getXYZ()[0] * B(ksi, eta)[2 * i];
        matrix[1] += _NID[i]->getXYZ()[0] * B(ksi, eta)[2 * i + 1];
        matrix[2] += _NID[i]->getXYZ()[1] * B(ksi, eta)[2 * i];
        matrix[3] += _NID[i]->getXYZ()[1] * B(ksi, eta)[2 * i + 1];
    }
    return matrix[0] * matrix[3] - matrix[1] * matrix[2]; 
}

// Set Element ID
void ElementQ4::setID(int ID) {
    _ID = ID;
};

// Get Element ID
int ElementQ4::getID() const {
    return _ID;
};

// Set Element NID
void ElementQ4::setNID(const vector<Node*> & NID) {
    _NID.resize(NID.size());
    for (int i = 0; i < NID.size(); i++) {
        _NID[i] = NID[i];
    }
};

// Get Element NID
const vector<Node*> & ElementQ4::getNID() const {
    return _NID;
};

// Output information
void ElementQ4::outputInfo(ofstream & myFile) const {
    myFile << setw(10) << getID() << " ";
    for (int i = 0; i < getNID().size(); i++) {
        myFile << setw(10) << getNID()[i]->getID() << " ";
    }
    myFile << "\n";
};

//------------------------------------------------------------------------------------------------------------------------------------------------------
/** Methods of class ElementQ4Cohesive
 * No default constructor implemented
 */
// Constructor 1 for ElementQ4Cohesive
ElementQ4Cohesive::ElementQ4Cohesive(int ID, const vector<CohesiveNode*> & NID) {
    setID(ID);
    setNID(NID);
};

// Destructor
ElementQ4Cohesive::~ElementQ4Cohesive() {
    for (int i = 0; i < _NID.size(); i++) {
        delete _NID[i];
    }
};

// Set element NID
void ElementQ4Cohesive::setNID(const vector<CohesiveNode*> & NID) {
    _NID.resize(NID.size());
    for (int i = 0; i < NID.size(); i++) {
        _NID[i] = NID[i];
    }
};

// Get element NID
const vector<CohesiveNode*> & ElementQ4Cohesive::getNID() const {
    return _NID;
};

// Output element info
void ElementQ4Cohesive::outputInfo(ofstream & myFile) const {
    myFile << setw(10) << getID() << " ";
    for (int i = 0; i < getNID().size(); i++) {
        myFile << setw(10) << getNID()[i]->getID() << " ";
    }
    myFile << "\n"; 
};
