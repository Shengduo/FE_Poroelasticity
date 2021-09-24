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
};

// Constants for 2-point gaussian integral
const vector<double> ElementQ4::IntPos = {- 1. / sqrt(3.), 1. / sqrt(3.)};
const vector<double> ElementQ4::IntWs = {1., 1.};

/** calculate inverse of Jacobian at any given location in base space (ksi, eta), 
 * invJ = (\partial (x,y) / \partial (ksi, eta)) ^ (-1)
 * returns false if J is singular
 */
bool ElementQ4::InvJ(vector<double> & res, double ksi, double eta) const {
    // Reset all res to 0.
    res.resize(4);
    fill(res.begin(), res.end(), 0.0);
    
    double matrix [4] = {0., 0., 0., 0.};
    // Calculate J_matrix
    for (int i = 0; i < 4; i++) {
        res[0] += _NID[i]->getXYZ()[0] * B(ksi, eta)[2 * i];
        res[1] += _NID[i]->getXYZ()[0] * B(ksi, eta)[2 * i + 1];
        res[2] += _NID[i]->getXYZ()[1] * B(ksi, eta)[2 * i];
        res[3] += _NID[i]->getXYZ()[1] * B(ksi, eta)[2 * i + 1];
    }

    // Calculate J
    double J = res[0] * res[3] - res[1] * res[2];
    
    // If J_matrix is singular
    if (fabs(J) <= 1e-15) {
        fill(res.begin(), res.end(), 0.);
        return false;
    }

    // Calculate invJ, using J as a temp
    double temp = res[3];
    res[3] = res[0] / J;
    res[0] = temp / J;
    res[1] = - res[1] / J;
    res[2] = - res[2] / J;
    return true;
};

/** Evaluate a function F at (ksi, eta) */
void ElementQ4::evaluateF(vector<double> & res, double ksi, double eta, 
                          const vector<vector<double>> & NodeValues) const {
    // Make sure sizes match
    if (NodeValues.size() != _NID.size()) throw "Not all nodal values are provided for evaluate";
    
    // Initialize result based on input
    res.resize(NodeValues[0].size(), 0.);
    for (int i = 0; i < res.size(); i++) res[i] = 0.;

    // Loop through all fields
    for (int f = 0; f < NodeValues[0].size(); f++) {
        // Loop through all nodes
        for (int n = 0; n < NodeValues.size(); n++) {
            res[f] += N(ksi, eta)[n] * NodeValues[n][f];
        }
    }
};

/** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at (ksi, eta) 
 * in LOGICAL space with given nodal values.
 * Calculated by using shape function to map
 */
void ElementQ4::evaluateF_x(vector<double> & res, double ksi, double eta, 
                            const vector<vector<double>> & NodeValues) const {
    // Number of nodes
    int nOfNodes = this->getNID().size();
    if (NodeValues.size() != nOfNodes) 
        throw "In evaluateF_x, nodeValues do not match number of nodes!";
    
    // Space dim is 2 for Q4
    int spaceDim = 2;

    // Set res to 0.;
    res.resize(spaceDim * NodeValues[0].size());
    fill(res.begin(), res.end(), 0.0);

    vector<double> invJ(spaceDim * spaceDim, 0.);
    if (!InvJ(invJ, ksi, eta)) throw "J is singular for evaluateF_x!";

    vector<double> res_s(res.size(), 0.);    
    // Loop through fields
    for(int f = 0; f < NodeValues[0].size(); f++) {
        // Loop through nodes
        for (int n = 0; n < NodeValues.size(); n++) {
            // Loop through spaceDim
            for (int d = 0; d < spaceDim; d++) {
                // Calculate \partial F / \partial \ksi
                res_s[d + f * spaceDim] += B(ksi, eta)[n * spaceDim + d] * NodeValues[n][f];
            }
        }
    }

    /** Apply invJ = (\partial x \partial \ksi) ^ -1 */
    // All fields
    for (int f = 0; f < NodeValues[0].size(); f++) {
        // All spaceDims
        for (int d = 0; d < spaceDim; d++) {
            res[f * spaceDim + d] = invJ[0 * spaceDim + d] * res_s[f * spaceDim + 0] 
                                    + invJ[1 * spaceDim + d] * res_s[f * spaceDim + 1];
        }
    }
};

// IntegratorNf, integrates a vector input inside an element, both sides using shape function.
// first-dim: vector of nodes, second-dim: values (vector)
void ElementQ4::IntegratorNf(vector<vector<double>> & res, 
                             const vector<vector<double>> & NodeValues) const {
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    if (NodeValues.size() != nOfNodes) throw "Not all nodal values are provided for ElementQ4 IntegratorNf!";
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
};

/** IntegratorNfN, integrates a vector input inside an element, both sides using shape function
 * first-dim: vector of nodes^2, second-dim: values (vector)
 */
void ElementQ4::IntegratorNfN(vector<vector<double>> & res, const vector<vector<double>> & NodeValues) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    if (res.size() != nOfNodes * nOfNodes) res.resize(nOfNodes * nOfNodes);
    if (NodeValues.size() != nOfNodes) throw "Not all nodal values are provided for ElementQ4 IntegratorNfN!";
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
    }
};

/** IntegratorBfB, integrates a vector input inside an element, 
 * both sides using gradient of shape function
 * RES:
 * first-dim: vector of nodes^2, 
 * NODEVALUES:
 * first-dim vector of nodes, 
 * second-dim 2 * 2 matrix.
 */
void ElementQ4::IntegratorBfB(vector<double> & res,
                              const vector<vector<double>> & NodeValues) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    if (res.size() != nOfNodes * nOfNodes) res.resize(nOfNodes * nOfNodes);
    if (NodeValues.size() != nOfNodes) throw "Not all nodal values are provided for ElementQ4 IntegratorBfB!";
    double pointValue = 0.;
    vector<double> pointD(getNID()[0]->getSpaceDim()^2, 0.);
    // Set res to all 0;
    for (int i = 0; i < res.size(); i++) {
        res[i] = 0.;
    }

    // Calculate res[i,j]
    // In the integral (B^T D B)_{i,j} = B_pi D_pq B qj
    for (int i = 0; i < nOfNodes; i++) {
        for (int j = 0; j < nOfNodes; j++) {
            for (int k = 0; k < nOfIntPts; k++) {
                for (int l = 0; l < nOfIntPts; l++) {
                    // Constants for the integral
                    pointValue = IntWs[k] * IntWs[l] * J(IntPos[k], IntPos[l]);
                    evaluateF(pointD, IntPos[k], IntPos[l], NodeValues); 
                    for (int p = 0; p < 2; p++) {
                        for (int q = 0; q < 2; q++) {
                            res[i * nOfNodes + j] += pointValue 
                                                     * B(IntPos[k], IntPos[l])[p + 2 * i] 
                                                     * pointD[2 * p + q] 
                                                     * B(IntPos[k], IntPos[l])[q + 2 * j]; 
                        }
                    }                   
                }
            }                
        }
    }
};

/** IntegratorBfN, integrates a vector input inside an element, 
 * both sides using gradient of shape function
 * first-dim: vector of nodes^2, second-dim: values (vector)
 * NODEVALUES, first dim: nodes, second dim: spaceDim by 1 matrix, stored as a vector
 */
void ElementQ4::IntegratorBfN(vector<double> & res,
                              const vector<vector<double>> & NodeValues) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    int spaceDim = this->getNID()[0]->getSpaceDim(); 
    if (res.size() != nOfNodes * nOfNodes) res.resize(nOfNodes * nOfNodes);
    if (NodeValues.size() != nOfNodes) throw "Not all nodal values are provided for ElementQ4 IntegratorBfN!";
    if (NodeValues[0].size() != spaceDim * 1) throw "Nodal matrix provided not compatible with IntegratorBfN!";
    
    // Temp values and vectors
    double pointValue = 0.;
    vector<double> pointVector(NodeValues[0].size(), 0.);

    // Set res to all 0;
    for (int i = 0; i < res.size(); i++) {
        res[i] = 0.;
    }

    // Calculate res[i,j]
    // In the integral (B^T F N)_{i,j} = B_pi F_pq N qj
    for (int i = 0; i < nOfNodes; i++) {
        for (int j = 0; j < nOfNodes; j++) {
            for (int k = 0; k < nOfIntPts; k++) {
                for (int l = 0; l < nOfIntPts; l++) {
                    // Constants for the integral
                    pointValue = IntWs[k] * IntWs[l] * J(IntPos[k], IntPos[l]);
                    evaluateF(pointVector, IntPos[k], IntPos[l], NodeValues); 
                    for (int p = 0; p < spaceDim; p++) {
                        for (int q = 0; q < 1; q++) {
                            res[i * nOfNodes + j] += pointValue 
                                                     * B(IntPos[k], IntPos[l])[p + spaceDim * i] 
                                                     * pointVector[1 * p + q] 
                                                     * N(IntPos[k], IntPos[l])[q + j]; 
                        }
                    }                   
                }
            }                
        }
    }
};

/** IntegratorNfB, integrates a vector input inside an element, 
 * both sides using gradient of shape function
 * first-dim: vector of nodes^2, second-dim: values (vector)
 * NODEVALUES, first dim: nodes, second dim: 1 by spaceDim matrix, stored as a vector)
 */
void ElementQ4::IntegratorNfB(vector<double> & res,
                   const vector<vector<double>> & NodeValues) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    int spaceDim = this->getNID()[0]->getSpaceDim(); 
    if (res.size() != nOfNodes * nOfNodes) res.resize(nOfNodes * nOfNodes);
    if (NodeValues.size() != nOfNodes) throw "Not all nodal values are provided for ElementQ4 IntegratorNfB!";
    if (NodeValues[0].size() != 1 * spaceDim) throw "Nodal matrix provided not compatible with IntegratorNfB!";
    
    // Temp values and vectors
    double pointValue = 0.;
    vector<double> pointVector(NodeValues[0].size(), 0.);

    // Set res to all 0;
    for (int i = 0; i < res.size(); i++) {
        res[i] = 0.;
    }

    // Calculate res[i,j]
    // In the integral (N^T F B)_{i,j} = N_pi F_pq B qj
    for (int i = 0; i < nOfNodes; i++) {
        for (int j = 0; j < nOfNodes; j++) {
            for (int k = 0; k < nOfIntPts; k++) {
                for (int l = 0; l < nOfIntPts; l++) {
                    // Constants for the integral
                    pointValue = IntWs[k] * IntWs[l] * J(IntPos[k], IntPos[l]);
                    evaluateF(pointVector, IntPos[k], IntPos[l], NodeValues); 
                    for (int p = 0; p < 1; p++) {
                        for (int q = 0; q < spaceDim; q++) {
                            res[i * nOfNodes + j] += pointValue 
                                                     * N(IntPos[k], IntPos[l])[p + i] 
                                                     * pointVector[p + q * 1] 
                                                     * B(IntPos[k], IntPos[l])[q + spaceDim * j]; 
                        }
                    }                   
                }
            }                
        }
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