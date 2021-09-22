/** @file ElementQ4Cohesive.cc
 * source file for class ElementQ4Cohesive
 */

#include "ElementQ4Cohesive.hh"

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

/** Evaluate a vector at ksi in logical space with given nodal values.
 * Calculated by using shape function to map.
 */
void ElementQ4Cohesive::evaluateF(vector<double> & res, double ksi,
                                  const vector<vector<double>> & NodeValues) const {
    // Make sure sizes match
    if (NodeValues.size() != _NID.size()) throw("Not all nodal values are provided for evaluateF");
    
    // Initialize result based on input
    res.resize(NodeValues[0].size(), 0.);
    for (int i = 0; i < res.size(); i++) res[i] = 0.;

    // Loop through all fields
    for (int f = 0; f < NodeValues[0].size(); f++) {
        // Loop through all nodes
        for (int n = 0; n < NodeValues.size(); n++) {
            res[f] += N(ksi)[n] * NodeValues[n][f];
        }
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

/** IntegratorNfN, integrates a vector input inside an element, both sides using shape function
 * first-dim: vector of nodes^2, second-dim: values (vector)
 */
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

/** IntegratorBfB, integrates a vector input inside an element, 
 * both sides using gradient of shape function
 * RES:
 * first-dim: vector of nodes^2, 
 * NODEVALUES:
 * first-dim vector of nodes, 
 * second-dim 2 * 2 matrix.
 */
void ElementQ4Cohesive::IntegratorBfB(vector<double> & res,
                                      const vector<vector<double>> & NodeValues) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    if (res.size() != nOfNodes * nOfNodes) res.resize(nOfNodes * nOfNodes);
    if (NodeValues.size() != nOfNodes) throw("Not all nodal values are provided for ElementQ4 IntegratorBfB!");
    double pointValue = 0.;
    vector<double> pointD(NodeValues[0].size(), 0.);
    // Set res to all 0;
    for (int i = 0; i < res.size(); i++) {
        res[i] = 0.;
    }

    // Calculate res[i,j]
    // In the integral (B^T D B)_{i,j} = B_pi D_pq B qj
    for (int i = 0; i < nOfNodes; i++) {
        for (int j = 0; j < nOfNodes; j++) {
            for (int k = 0; k < nOfIntPts; k++) {
                // Constants for the integral
                pointValue = IntWs[k] * J(IntPos[k]);
                evaluateF(pointD, IntPos[k], NodeValues);
                for (int p = 0; p < 1; p++) {
                    for (int q = 0; q < 1; q++) {
                        res[i * nOfNodes + j] += pointValue 
                                * B(IntPos[k])[p + i] 
                                * pointD[p + q] 
                                * B(IntPos[k])[q + j]; 
                    }
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