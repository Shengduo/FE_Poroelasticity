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
    // DEBUG LINES
    // cout << "J= " << matrix[0] * matrix[3] - matrix[1] * matrix[2] << "\n";
    return matrix[0] * matrix[3] - matrix[1] * matrix[2]; 
};

// Constants for 2-point gaussian integral
const vector<double> ElementQ4::IntPos = {- 1. / sqrt(3.), 1. / sqrt(3.)};
const vector<double> ElementQ4::IntWs = {1., 1.};

// IntegratorNf, integrates a vector input inside an element, both sides using shape function.
// first-dim: vector of nodes, second-dim: values (vector)
void ElementQ4::IntegratorNf(vector<vector<double>> & res, 
                             const vector<vector<double>> & NodeValues) const {
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    if (NodeValues.size() != nOfNodes) throw("Not all nodal values are provided for ElementQ4 Integrator!");
    // DEBUG LINES
    /**
    for (int i = 0; i < NodeValues.size(); i++) {
        for (int j = 0; j < NodeValues[i].size(); j++) {
            cout << NodeValues[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
    */
    // Set res first to all 0.
    for (int i = 0; i < res.size(); i++) {
        for (int j = 0; j < res[i].size(); j++) res[i][j] = 0.;
    }

    // First Calculate intN^T N, four point Gaussian integral
    vector<double> IntNTN(nOfNodes * nOfNodes, 0.);

    for (int i = 0; i < nOfNodes; i++) {
        for (int j = 0; j < nOfNodes; j++) {
            for (int k = 0; k < nOfIntPts; k++) {
                for (int l = 0; l < nOfIntPts; l++) {
                    // \int_{\Omega} f = \sigma _{i,j} w_i * w_j * J(i, j) * f(i,j)
                    IntNTN[i * 4 + j] += N(IntPos[k], IntPos[l])[i] * N(IntPos[k], IntPos[l])[j] * IntWs[k] * IntWs[l] * J(IntPos[k], IntPos[l]);
                }
            }
        }
    } 

    // DEBUG LINES
    /**
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << IntNTN[4 * i + j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
    */ 

    // First loop through all fields
    for (int j = 0; j < NodeValues[0].size(); j++) {
        // Then loop through all points
        for (int i = 0; i < NodeValues.size(); i++) {
            // Calculate IntNTN(i-th row) * f[j]
            for (int k = 0; k < 4; k++) {
                res[i][j] += IntNTN[i * 4 + k] * NodeValues[k][j];
            }
        }
    }
    // DEBUG LINES
    /**
    for (int i = 0; i < NodeValues.size(); i++) {
        for (int j = 0; j < NodeValues[i].size(); j++) {
            cout << res[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
    */
};

// IntegratorNfN, integrates a vector input inside an element, both sides using shape function
// first-dim: vector of nodes^2, second-dim: values (vector)
void ElementQ4::IntegratorNfN(vector<vector<double>> & res, const vector<vector<double>> & NodeValues) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    if (res.size() != nOfNodes * nOfNodes) res.resize(nOfNodes * nOfNodes);
    if (NodeValues.size() != nOfNodes) throw("Not all nodal values are provided for ElementQ4 Integrator!");
    int nOfFields = NodeValues[0].size();
    double pointValue = 0.;
    // Set res to all 0;
    for (int i = 0; i < res.size(); i++) {
        for (int j = 0; j < res[i].size(); j++) res[i][j] = 0.;
    }

    for (int f = 0; f < nOfFields; f++) {
        // Calculate res[i,j][f]
        for (int i = 0; i < nOfNodes; i++) {
            for (int j = 0; j < nOfNodes; j++) {
                for (int k = 0; k < nOfIntPts; k++) {
                    for (int l = 0; l < nOfIntPts; l++) {
                        pointValue = 0.;
                        for (int p = 0; p < nOfNodes; p++) pointValue += NodeValues[p][f] * N(IntPos[k], IntPos[l])[p];
                        res[nOfNodes * i + j][f] += 
                            N(IntPos[k], IntPos[l])[i] * N(IntPos[k], IntPos[l])[j] * IntWs[k] * IntWs[l] * J(IntPos[k], IntPos[l]) * pointValue;
                    }
                }                
            }
        }

        // DEBUG LINES
        /**
        for (int i = 0; i < nOfNodes; i++) {
            for (int j = 0; j < nOfNodes; j++) {
                cout << setw(12) << res[nOfNodes * i + j][f] << " ";
            }
            cout << "\n";
        }
        cout << "\n";
        */
    }
};


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

// Shape function N, vector of 2 on 2 nodes, evaluated in base space (ksi)
vector<double> ElementQ4Cohesive::N(double ksi) {
    return vector<double> {(1. - ksi) / 2., (1. + ksi) / 2.};
};

// Gradient of shape function B, vector of 2 * 1 on 2 nodes, evaluated in base space (ksi)
vector<double> ElementQ4Cohesive::B(double ksi) {
    return vector<double> {- 1. / 2., 1. / 2.};
};

// Jacobian at any given location in base space (ksi, eta), J = d|x| / d ksi
double ElementQ4Cohesive::J(double ksi) const {
    return sqrt(pow(- _NID[0]->getXYZ()[0] + _NID[1]->getXYZ()[0], 2) + 
                pow(- _NID[0]->getXYZ()[1] + _NID[1]->getXYZ()[1], 2)) / 2.; 
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

// IntegratorNf, integrates a vector input inside an element, both sides using shape function.
// first-dim: vector of nodes, second-dim: values (vector)
void ElementQ4Cohesive::IntegratorNf(vector<vector<double>> & res, 
                                     const vector<vector<double>> & NodeValues) const {
    if (NodeValues.size() != 2) throw("Not all nodal values are provided for ElementQ4Cohesive Integrator!");
    for (int i = 0; i < res.size(); i++) {
        for (int j = 0; j < res[i].size(); j++) res[i][j] = 0.;
    }

    // First Calculate intN^T N, four point Gaussian integral
    vector<double> IntNTN(2 * 2, 0.);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                // \int_{\Omega} f = \sigma _{i} w_i * J(i) * f(i)
                IntNTN[i * 2 + j] += N(IntPos[k])[i] * N(IntPos[k])[j] * IntWs[k] * J(IntPos[k]);
            }
        }
    } 

    // First loop through all fields
    for (int j = 0; j < NodeValues[0].size(); j++) {
        // Then loop through all points
        for (int i = 0; i < NodeValues.size(); i++) {
            // Calculate IntNTN(i-th row) * f[j]
            for (int k = 0; k < 2; k++) {
                res[i][j] += IntNTN[i * 2 + k] * NodeValues[k][j];
            }
        }
    }
};

// IntegratorNfN, integrates a vector input inside an element, both sides using shape function
// first-dim: vector of nodes^2, second-dim: values (vector)
void ElementQ4Cohesive::IntegratorNfN(vector<vector<double>> & res, const vector<vector<double>> & NodeValues) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    if (res.size() != nOfNodes * nOfNodes) res.resize(nOfNodes * nOfNodes);
    if (NodeValues.size() != nOfNodes) throw("Not all nodal values are provided for ElementQ4Cohesive Integrator!");
    int nOfFields = NodeValues[0].size();
    double pointValue = 0.;
    // Set res to all 0;
    for (int i = 0; i < res.size(); i++) {
        for (int j = 0; j < res[i].size(); j++) res[i][j] = 0;
    }

    for (int f = 0; f < nOfFields; f++) {
        // Calculate res[i,j][f]
        for (int i = 0; i < nOfNodes; i++) {
            for (int j = 0; j < nOfNodes; j++) {
                for (int k = 0; k < nOfIntPts; k++) {
                    // Compute point value.
                    pointValue = 0.;
                    for (int p = 0; p < nOfNodes; p++) pointValue += NodeValues[p][f] * N(IntPos[k])[p];
                    res[nOfNodes * i + j][f] += 
                        N(IntPos[k])[i] * N(IntPos[k])[j] * IntWs[k] * J(IntPos[k]) * pointValue;
                }                
            }
        }
    }
};

// Output element info
void ElementQ4Cohesive::outputInfo(ofstream & myFile) const {
    myFile << setw(10) << getID() << " ";
    for (int i = 0; i < getNID().size(); i++) {
        myFile << setw(10) << getNID()[i]->getID() << " ";
    }
    myFile << "\n";
};
