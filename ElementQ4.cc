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
    elementDOF = 0;
    this->_ID = ID;
    this->_NID.resize(NID.size());
    for (int i = 0; i < NID.size(); i++) {
        this->_NID[i] = NID[i];
        elementDOF += NID[i]->getDOF().size();
    }
    
    int nOfIntPts = IntPos.size();
    
    // Calculate Bvector, Nvector, pointValue 
    localGlobalIndices = new PetscInt [elementDOF];
    Bvector.resize(pow(nOfIntPts, 2));
    Nvector.resize(pow(nOfIntPts, 2));
    pointValue.resize(pow(nOfIntPts, 2));

    for (int i = 0; i < nOfIntPts; i++) {
        for (int j = 0; j < nOfIntPts; j++) {
            Bvector[i * nOfIntPts + j] = B_x(IntPos[i], IntPos[j]);
            Nvector[i * nOfIntPts + j] = N(IntPos[i], IntPos[j]);
            pointValue[i * nOfIntPts + j] = J(IntPos[i], IntPos[j]) * IntWs[i] * IntWs[j];
        }
    }
    int index = 0;
    // Calculate localGlobalIndices
    for (int n = 0; n < 4; n++) {
        for (int i = 0; i < NID[n]->getDOF().size(); i++) {
            localGlobalIndices[index] = NID[n]->getDOF()[i];
            index += 1;
        }
    }

};

// Destructor
ElementQ4::~ElementQ4(){
    // Delete vector<Node*> _NID
    for (int i = 0; i < _NID.size(); i++) {
        delete _NID[i];
    }

    // Delete elementDOF
    delete [] localGlobalIndices;
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

/** Shape function N, vector of 4 on 4 nodes, 
 * evaluated in base space (ksi, eta)
 * at i, j of elemental shape matrix
 */
double ElementQ4::N(const vector<double> & Nvector, int i, int j) const {
    // Number of Nodes
    int nOfDofs = _NID[0]->getDOF().size();
    double res = 0.;
    
    if (i == j % nOfDofs) {
        res = Nvector[j / nOfDofs];
    }

    return res;
};

// Gradient of shape function N at (ksi, eta)
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

// Gradient_x of shape function B at (ksi, eta)
vector<double> ElementQ4::B_x(double ksi, double eta) const{
    // Space dimension is 2 for Q4
    int spaceDim = 2;
    int nOfNodes = _NID.size();
    vector<double> B_eta = B(ksi, eta); 
    // [\partial N[0] / \partial \ksi, \partial N[0] / \partial \eta, \partial N[1] / \partial \ksi...]
    vector<double> invJ(spaceDim * spaceDim, 0.);
    InvJ(invJ, ksi, eta);
    vector<double> B_x(B_eta.size(), 0.);

    for (int n = 0; n < nOfNodes; n++) {
        for (int i = 0; i < spaceDim; i++) {
            for (int j = 0; j < spaceDim; j++) {
                B_x[i + n * spaceDim] += B_eta[j + n * spaceDim] * invJ[j * spaceDim + i];
            }            
        }
    }
    return B_x;
};

// Gradient of shape function B_x[i, j] at (ksi, eta)
double ElementQ4::B_x(const vector<double> & Bvector, int i, int j) const {
    // [\partial N[0] / \partial \ksi, \partial N[0] / \partial \eta, \partial N[1] / \partial \ksi...]
    // Some constants
    int spaceDim = 2;
    int nOfDofs = 2 * spaceDim + 2;
    double res = 0.;
    if (i / spaceDim == j % nOfDofs) {
        res = Bvector[j / nOfDofs * spaceDim + i % spaceDim];
    }
    return res;
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
    
    // double matrix [4] = {0., 0., 0., 0.};
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

// ================ Integrators of ElementQ4 ========================================
/** IntegratorNf, integrates a vector input inside an element, 
 * left side using shape function
 * RES: nOfNodes * nOfDofs
 * NODEVALUES: dim 1, nOfNodes; dim 2, nOfDofs
 * FLAG: 0 - NodeValues are given at nodes,
 *       1 - NodeValues are given at integration points
 */
void ElementQ4::IntegratorNf(double *res, 
                             int resSize, 
                             const vector<vector<double>> & NodeValues, 
                             int flag) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfDofs = this->getNID()[0]->getDOF().size();
    int spaceDim = this->getNID()[0]->getSpaceDim();
    int nOfIntPts = IntPos.size();
    if (NodeValues.size() != nOfNodes) throw "Not all nodes are provided for ElementQ4 IntegratorNf!";
    if (NodeValues[0].size() != nOfDofs) throw "Not all nodal DOFs are provided for ElementQ4 IntegratorNf!";
    if (resSize != nOfNodes * nOfDofs) throw "ResSize error for ElementQ4 IntegratorNf!";

    // Pre-Calculate values of N and f at integration points
    // vector<double> pointValue(pow(nOfIntPts,spaceDim));
    // vector<vector<double>> Nvector(pow(nOfIntPts,spaceDim));
    vector<vector<double>> fvector(pow(nOfIntPts, spaceDim));
    for (int i = 0; i < nOfIntPts; i++) {
        for (int j = 0; j < nOfIntPts; j++) {
            // Nvector[i * nOfIntPts + j] = N(IntPos[i], IntPos[j]);
            // pointValue[i * nOfIntPts + j] = J(IntPos[i], IntPos[j]) * IntWs[i] * IntWs[j];
            if (flag == 1) {
                fvector[i * nOfIntPts + j] = NodeValues[i * nOfIntPts + j];
            }
            else {
                evaluateF(fvector[i * nOfIntPts + j], IntPos[i], IntPos[j], NodeValues);
            }            
        }
    }

    // Computing the integral
    // res_i = N_{p, i} f_p
    int IntPtIndex;
    for (int i = 0; i < nOfNodes * nOfDofs; i++) {
        for (int p = 0; p < nOfDofs; p++) {
            // All integration points
            for (int k = 0; k < nOfIntPts; k++) {
                for (int l = 0; l < nOfIntPts; l++) {
                    IntPtIndex = k * nOfIntPts + l;
                    res[i] += N(Nvector[IntPtIndex], p, i) * pointValue[IntPtIndex] * fvector[IntPtIndex][p];
                }
            }
        }
    }
};

/** IntegratorBf, integrates a vector input inside an element, 
 * left side using shape function
 * RES: nOfNodes * nOfDofs
 * NODEVALUES: dim 1, nOfNodes; dim 2, nOfDofs * spaceDim
 * FLAG: 0 - nodevalues are given at nodes, 
 *       1 - nodevalues are given at integration points
 */
void ElementQ4::IntegratorBf(double *res, 
                             int resSize, 
                             const vector<vector<double>> & NodeValues, 
                             int flag) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfDofs = this->getNID()[0]->getDOF().size();
    int spaceDim = this->getNID()[0]->getSpaceDim();
    int nOfIntPts = IntPos.size();
    if (NodeValues.size() != nOfNodes) throw "Not all nodes are provided for ElementQ4 IntegratorBf!";
    if (NodeValues[0].size() != nOfDofs * spaceDim) throw "Not all nodal DOFs are provided for ElementQ4 IntegratorBf!";
    if (resSize != nOfDofs * nOfNodes) throw "resSize error for ElementQ4 IntegratorBf!";

    // Pre-Calculate values of N and f at integration points
    // vector<double> pointValue(pow(nOfIntPts,spaceDim));
    // vector<vector<double>> Bvector(pow(nOfIntPts,spaceDim));
    vector<vector<double>> fvector(pow(nOfIntPts, spaceDim));
    for (int i = 0; i < nOfIntPts; i++) {
        for (int j = 0; j < nOfIntPts; j++) {
            // Bvector[i * nOfIntPts + j] = B_x(IntPos[i], IntPos[j]);
            // pointValue[i * nOfIntPts + j] = J(IntPos[i], IntPos[j]) * IntWs[i] * IntWs[j];
            if (flag == 1) {
                fvector[i * nOfIntPts + j] = NodeValues[i * nOfIntPts + j];
            }
            else {
                evaluateF(fvector[i * nOfIntPts + j], IntPos[i], IntPos[j], NodeValues);
            }            
        }
    }

    // Computing the integral
    // res_i = B_{p, i} f_p
    int IntPtIndex;
    for (int i = 0; i < nOfNodes * nOfDofs; i++) {
        for (int p = 0; p < nOfDofs * spaceDim; p++) {
            // All integration points
            for (int k = 0; k < nOfIntPts; k++) {
                for (int l = 0; l < nOfIntPts; l++) {
                    IntPtIndex = k * nOfIntPts + l;
                    res[i] += B_x(Bvector[IntPtIndex], p, i) * pointValue[IntPtIndex] * fvector[IntPtIndex][p];
                }
            }
        }
    }
};

/** IntegratorNfN, integrates a vector input inside an element, 
 * both sides using shape function (nOfDofs, nOfDofs * nOfNodes)
 * first-dim: vector of nodes, second-dim: values (vector, (nOfDofs, nOfDofs))
 * FLAG: 0 - nodevalues are given at nodes,
 *       1 - nodevalues are given at integration points
 */
void ElementQ4::IntegratorNfN(double *res, 
                              int resSize, 
                              const vector<vector<double>> & NodeValues, 
                              int flag) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    int nOfDofs = getNID()[0]->getDOF().size();
    int nOfColsRes = nOfNodes * nOfDofs;

    if (resSize != pow(nOfColsRes, 2)) throw "resSize error for ElementQ4 IntegratorNfN!";
    if (NodeValues.size() != nOfNodes) throw "Not all nodal values are provided for ElementQ4 IntegratorNfN!";
    if (NodeValues[0].size() != nOfDofs * nOfDofs) throw "Mass matrix size not compatible with Element Q4 IntegratorNfN!";
    
    // Pre-calculate and store some values
    // vector<double> pointValue(nOfIntPts * nOfIntPts, 0.);
    vector<vector<double>> pointM(nOfIntPts * nOfIntPts);
    // vector<vector<double>> Nvector(nOfIntPts * nOfIntPts);
    int intPtIndex;
    for (int i = 0; i < nOfIntPts; i++) {
        for (int j = 0; j < nOfIntPts; j++) {
            intPtIndex = nOfIntPts * i + j;
            // pointValue[intPtIndex] = IntWs[i] * IntWs[j] * J(IntPos[i], IntPos[j]);
            if (flag == 1) {
                pointM[intPtIndex] = NodeValues[intPtIndex];
            }
            else {
                evaluateF(pointM[intPtIndex], IntPos[i], IntPos[j], NodeValues);
            }
            // Nvector[intPtIndex] = N(IntPos[i], IntPos[j]);
        }
    }

    int resIJindex;
    // Calculate res[i,j] = N_{pi} M_{pq} N_{qj}
    for (int i = 0; i < nOfColsRes; i++) {
        for (int j = 0; j < nOfColsRes; j++) {
            resIJindex = nOfColsRes * i + j;
            for (int k = 0; k < nOfIntPts; k++) {
                for (int l = 0; l < nOfIntPts; l++) {
                    intPtIndex = nOfIntPts * k + l;
                    for (int p = 0; p < nOfDofs; p++) {
                        for (int q = 0; q < nOfDofs; q++) {
                            res[resIJindex] += N(Nvector[intPtIndex], p, i) 
                                               * pointM[intPtIndex][p * nOfDofs + q] 
                                               * N(Nvector[intPtIndex], q, j)
                                               * pointValue[intPtIndex];
                        }
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
 * second-dim (spaceDim * nDof) ^ 2 matrix.
 * FLAG: 0 - nodevalues are given at nodes
 *       1 - nodevalues are given at integration points
 */
void ElementQ4::IntegratorBfB(double *res,
                              int resSize, 
                              const vector<vector<double>> & NodeValues, 
                              int flag) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfDofs = this->getNID()[0]->getDOF().size();
    int nOfIntPts = IntPos.size();
    int spaceDim = this->getNID()[0]->getSpaceDim();
    int nOfColsF = nOfDofs * spaceDim;
    int nOfColsRes = nOfDofs * nOfNodes;

    if (resSize != nOfColsRes * nOfColsRes) throw "resSize error in ElementQ4 IntegratorBfB!";
    
    if (NodeValues.size() != nOfNodes) 
        throw "Not all nodal values are provided for ElementQ4 IntegratorBfB!";
    
    if (NodeValues[0].size() != nOfColsF * nOfColsF) 
        throw "Input f is not compatible with nodal dofs in ElementQ4 IntegratorBfB!";

    // Pre-store integration constants at different points
    // vector<double> pointValue(nOfIntPts * nOfIntPts);
    vector<vector<double>> pointD(nOfIntPts * nOfIntPts, vector<double>(NodeValues[0].size(), 0.));

    // Stores the B_x vectors
    // vector<vector<double>> Bvector(nOfIntPts * nOfIntPts);

    // Calculate Bvector, pointValue and pointD
    for (int i = 0; i < nOfIntPts; i++) {
        for (int j = 0; j < nOfIntPts; j++) {
            // pointValue[i * nOfIntPts + j] = IntWs[i] * IntWs[j] * J(IntPos[i], IntPos[j]);
            // Bvector[i * nOfIntPts + j] = B_x(IntPos[i], IntPos[j]);
            if (flag == 1) {
                pointD[i * nOfIntPts + j] = NodeValues[i * nOfIntPts + j];
            }
            else {
                evaluateF(pointD[i * nOfIntPts + j], IntPos[i], IntPos[j], NodeValues);
            }                         
        }
    }
    
    // Calculate res[i,j]
    // In the integral (B^T D B)_{i,j} = B_pi D_pq B qj
    int intPtIndex = 0;
    for (int i = 0; i < nOfColsRes; i++) {
        for (int j = 0; j < nOfColsRes; j++) {
            for (int k = 0; k < nOfIntPts; k++) {
                for (int l = 0; l < nOfIntPts; l++) {
                    for (int p = 0; p < nOfColsF; p++) {
                        for (int q = 0; q < nOfColsF; q++) {
                            intPtIndex = k * nOfIntPts + l;
                            res[i * nOfColsRes + j] += pointValue[intPtIndex]
                                                     * B_x(Bvector[intPtIndex], p, i)
                                                     // Debug line
                                                     * pointD[intPtIndex][nOfColsF * p + q] 
                                                     * B_x(Bvector[intPtIndex], q, j); 
                        }
                    }                   
                }
            }                
        }
    }
};

/** IntegratorBfN, integrates a vector input inside an element, 
 * left gradient of shape function, 
 * right shape function
 * RES: first-dim: vector of nOfDof^2
 * NODEVALUES: first dim: spaceDim * nOfDof, 
 * second dim: nOfDof
 * FLAG: 0 - nodevalues given at nodes, 
 *       1 - at integration points
 */
void ElementQ4::IntegratorBfN(double *res,
                              int resSize, 
                              const vector<vector<double>> & NodeValues, 
                              int flag) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    int spaceDim = this->getNID()[0]->getSpaceDim();
    int nOfDofs = this->getNID()[0]->getDOF().size();
    int nRowsF = spaceDim * nOfDofs;
    int nColsF = nOfDofs;
    int nColsRes = nOfDofs * nOfNodes; 
    if (resSize != nColsRes * nColsRes) throw "resSize error in ElementQ4 IntegratorBfN!";
    if (NodeValues.size() != nOfNodes) throw "Not all nodal values are provided for ElementQ4 IntegratorBfN!";
    if (NodeValues[0].size() != nRowsF * nColsF) throw "Nodal matrix provided not compatible with IntegratorBfN!";

    // Pre-calculate and store some values
    // vector<double> pointValue(nOfIntPts * nOfIntPts, 0.);
    vector<vector<double>> pointF(nOfIntPts * nOfIntPts);
    // vector<vector<double>> Nvector(nOfIntPts * nOfIntPts);
    // vector<vector<double>> Bvector(nOfIntPts * nOfIntPts);
    int intPtIndex;
    for (int i = 0; i < nOfIntPts; i++) {
        for (int j = 0; j < nOfIntPts; j++) {
            intPtIndex = nOfIntPts * i + j;
            // pointValue[intPtIndex] = IntWs[i] * IntWs[j] * J(IntPos[i], IntPos[j]);
            
            if (flag == 1) {
                pointF[intPtIndex] = NodeValues[intPtIndex];
            }
            else {
                evaluateF(pointF[intPtIndex], IntPos[i], IntPos[j], NodeValues);
            }

            // Nvector[intPtIndex] = N(IntPos[i], IntPos[j]);
            // Bvector[intPtIndex] = B_x(IntPos[i], IntPos[j]);
        }
    }

    // Stores IJ index in the vector
    int resIJindex;
    // Calculate res[i,j]
    // In the integral (B^T F N)_{i,j} = B_pi F_pq N qj
    for (int i = 0; i < nColsRes; i++) {
        for (int j = 0; j < nColsRes; j++) {
            resIJindex = i * nColsRes + j;
            // Loop through integration points
            for (int k = 0; k < nOfIntPts; k++) {
                for (int l = 0; l < nOfIntPts; l++) {
                    // Constants for the integral
                    intPtIndex = nOfIntPts * k + l;
                    for (int p = 0; p < nRowsF; p++) {
                        for (int q = 0; q < nColsF; q++) {
                            res[resIJindex] += pointValue[intPtIndex] 
                                               * B_x(Bvector[intPtIndex], p, i)
                                               * pointF[intPtIndex][p * nColsF + q]
                                               * N(Nvector[intPtIndex], q, j);
                        }
                    }                   
                }
            }                
        }
    }
};

/** IntegratorNfB, integrates a vector input inside an element, 
 * left side shape function, right side gradient of shape function, 
 * RES: first-dim: vector of nodes^2, second-dim: values (vector)
 * NODEVALUES: first dim: nodes, second dim: 1 by spaceDim matrix, stored as a vector)
 * FLAG: 0 - nodevalues are given at nodes, 
 *       1 - nodevalues are given at integration points
 */
void ElementQ4::IntegratorNfB(double *res,
                              int resSize, 
                              const vector<vector<double>> & NodeValues, 
                              int flag) const {
    // Some constants
    int nOfNodes = this->getNID().size();
    int nOfIntPts = IntPos.size();
    int spaceDim = this->getNID()[0]->getSpaceDim();
    int nOfDofs = this->getNID()[0]->getDOF().size();
    int nRowsF = nOfDofs;
    int nColsF = spaceDim * nOfDofs;
    int nColsRes = nOfDofs * nOfNodes; 
    if (resSize != nColsRes * nColsRes) throw "resSize error for ElementQ4 IntegratorNfB!";
    if (NodeValues.size() != nOfNodes) throw "Not all nodal values are provided for ElementQ4 IntegratorNfB!";
    if (NodeValues[0].size() != nRowsF * nColsF) throw "Nodal matrix provided not compatible with IntegratorNfB!";

    // Pre-calculate and store some values
    // vector<double> pointValue(nOfIntPts * nOfIntPts, 0.);
    vector<vector<double>> pointF(nOfIntPts * nOfIntPts);
    // vector<vector<double>> Nvector(nOfIntPts * nOfIntPts);
    // vector<vector<double>> Bvector(nOfIntPts * nOfIntPts);
    int intPtIndex;
    for (int i = 0; i < nOfIntPts; i++) {
        for (int j = 0; j < nOfIntPts; j++) {
            intPtIndex = nOfIntPts * i + j;
            // pointValue[intPtIndex] = IntWs[i] * IntWs[j] * J(IntPos[i], IntPos[j]);
            if (flag == 1) {
                pointF[intPtIndex] = NodeValues[intPtIndex];
            }
            else {
                evaluateF(pointF[intPtIndex], IntPos[i], IntPos[j], NodeValues);
            }
            // Nvector[intPtIndex] = N(IntPos[i], IntPos[j]);
            // Bvector[intPtIndex] = B_x(IntPos[i], IntPos[j]);
        }
    }

    // Stores IJ index in the vector
    int resIJindex;
    // Calculate res[i,j]
    // In the integral (N^T F B)_{i,j} = N_pi F_pq B qj
    for (int i = 0; i < nColsRes; i++) {
        for (int j = 0; j < nColsRes; j++) {
            resIJindex = i * nColsRes + j;
            // Loop through integration points
            for (int k = 0; k < nOfIntPts; k++) {
                for (int l = 0; l < nOfIntPts; l++) {
                    // Constants for the integral
                    intPtIndex = nOfIntPts * k + l;
                    for (int p = 0; p < nRowsF; p++) {
                        for (int q = 0; q < nColsF; q++) {
                            res[resIJindex] += pointValue[intPtIndex] 
                                                * N(Nvector[intPtIndex], p, i)
                                                * pointF[intPtIndex][p * nColsF + q]
                                                * B_x(Bvector[intPtIndex], q, j);
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

// Get ElementDOF
int ElementQ4::getElementDOF() const {
    return elementDOF;
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


//============ Element Jacobians and residuals =====================================================

/** Calculate element jacobian Jf */
void ElementQ4::JF(Mat & globalJF, double *localJF, int localJFSize, int Kernel, double s_tshift) const {
    // Zero localJF
    for (int i = 0; i < localJFSize; i++) localJF[i] = 0.;

    // Switch kernel
    switch (Kernel) {
        // 2D Quasi-static linear elasticity
        case 0: {
            int nOfNodes = this->getNID().size();
            int spaceDim = this->getNID()[0]->getSpaceDim();
            int nOfDofs = this->getNID()[0]->getDOF().size();
            int nCols = spaceDim * nOfDofs;
            // int nColsJF = nOfDofs * nOfNodes;
            vector<vector<double>> NodeJFs(this->getNID().size(), vector<double>(nCols * nCols, 0.));

            // Call pointwise Jacobian
            for (int i = 0; i < nOfNodes; i++) {
                ElasticKernel::Jf3(NodeJFs[i], spaceDim, this->getNID()[i]->getNodalProperties());
            }

            // Calculate B^T D B
            IntegratorBfB(localJF, localJFSize, NodeJFs);
            
            // Assemble to globalJF
            JFPush(globalJF, localJF, localJFSize);      

            /** DEBUG LINES
            ofstream myFile;
            myFile.open("NodeJFs.txt");
            // Nodal D matrix
            for (int n = 0; n < nOfNodes; n++) {
                myFile << "Node ID: " << this->getNID()[n]->getID() << "\n";
                for (int i = 0; i < nCols; i++) {
                    for (int j = 0; j < nCols; j++) {
                        myFile << setw(12) << NodeJFs[n][i * nCols + j] << " ";
                    }
                    myFile << "\n";
                }
                myFile << "\n";
            }    
            */  
        }
        
        // 2D Isotropic linear poroelasticity
        case 1: {
            int nOfNodes = this->getNID().size();
            int spaceDim = this->getNID()[0]->getSpaceDim();
            int nOfDofs = this->getNID()[0]->getDOF().size();
            // int nCols = spaceDim * nOfDofs;
            int nColsJF = nOfDofs * nOfNodes;
            int nOfIntPts = IntPos.size();

            // Initialize pointwise Jfs
            vector<vector<double>> Jf0s(pow(nOfIntPts, spaceDim));
            vector<vector<double>> Jf1s(pow(nOfIntPts, spaceDim));
            vector<vector<double>> Jf2s(pow(nOfIntPts, spaceDim));
            vector<vector<double>> Jf3s(pow(nOfIntPts, spaceDim));

            // Jfs after integration
            vector<double> Jf0(pow(nColsJF, spaceDim), 0.);
            vector<double> Jf1(pow(nColsJF, spaceDim), 0.);
            vector<double> Jf2(pow(nColsJF, spaceDim), 0.);
            vector<double> Jf3(pow(nColsJF, spaceDim), 0.);
            vector<double> Jf(pow(nColsJF, spaceDim), 0.);
            
            // For each integration point
            vector<double> ss(nOfDofs, 0.);
            vector<double> s_xs(nOfDofs * spaceDim, 0.);
            vector<double> s_ts(nOfDofs, 0.);
            vector<double> as;
            
            // Nodal values
            vector<vector<double>> nodalSs(nOfNodes);
            vector<vector<double>> nodalS_ts(nOfNodes);
            vector<vector<double>> nodalAs(nOfNodes);
            for (int n = 0; n < nOfNodes; n++) {
                nodalSs[n] = _NID[n]->s;
                nodalS_ts[n] = _NID[n]->s_t;
                nodalAs[n] = _NID[n]->getNodalProperties();
            }
            
            // Calculate point values at integration points
            for (int i = 0; i < nOfIntPts; i++) {
                for (int j = 0; j < nOfIntPts; j++) {
                    evaluateF(ss, IntPos[i], IntPos[j], nodalSs);
                    evaluateF_x(s_xs, IntPos[i], IntPos[j], nodalSs);
                    evaluateF(as, IntPos[i], IntPos[j], nodalAs);
                    evaluateF(s_ts, IntPos[i], IntPos[j], nodalS_ts);

                    PoroelasticKernel::Jf0(Jf0s[i * nOfIntPts + j],
                                           spaceDim,
                                           ss,
                                           s_xs,
                                           s_ts,
                                           s_tshift,
                                           as,
                                           as,
                                           as);
                    
                    PoroelasticKernel::Jf1(Jf1s[i * nOfIntPts + j],
                                           spaceDim,
                                           ss,
                                           s_xs,
                                           s_ts,
                                           s_tshift,
                                           as,
                                           as,
                                           as);

                    PoroelasticKernel::Jf2(Jf2s[i * nOfIntPts + j],
                                           spaceDim,
                                           ss,
                                           s_xs,
                                           s_ts,
                                           s_tshift,
                                           as,
                                           as,
                                           as);

                    PoroelasticKernel::Jf3(Jf3s[i * nOfIntPts + j],
                                           spaceDim,
                                           ss,
                                           s_xs,
                                           s_ts,
                                           s_tshift,
                                           as,
                                           as,
                                           as);
                }
            }

            // Integration
            IntegratorNfN(localJF, localJFSize, Jf0s, 1);
            IntegratorNfB(localJF, localJFSize, Jf1s, 1);
            IntegratorBfN(localJF, localJFSize, Jf2s, 1);
            IntegratorBfB(localJF, localJFSize, Jf3s, 1);
            
            // Push to the global JF
            JFPush(globalJF, localJF, localJFSize);
        }
        default:
            break;
    }
};


/** Calculate element residual F, 
 * and then push ot globalF
 */
void ElementQ4::elementF(Vec & globalF, double *localF, int localFSize, int Kernel, double s_tshift) const {
    // Zero local F
    for (int i = 0; i < localFSize; i++) localF[i] = 0.;
    
    // Switch kernel 0 - elastic (not implemented), 1 - poroelastic
    switch(Kernel) {
        // Linear poroelastic kernels
        case 1: {
            // Some constants
            int spaceDim = _NID[0]->getSpaceDim();
            int nOfNodes = _NID.size();
            int nOfDofs = _NID[0]->getDOF().size();
            int nOfIntPts = IntPos.size();

            // For each integration point
            vector<vector<double>> F0s(pow(nOfIntPts, spaceDim), vector<double>(nOfDofs, 0.));
            vector<vector<double>> F1s(pow(nOfIntPts, spaceDim), vector<double>(nOfDofs * spaceDim, 0.));
            vector<double> ss(nOfDofs, 0.);
            vector<double> s_xs(nOfDofs * spaceDim, 0.);
            vector<double> s_ts(nOfDofs, 0.);
            vector<double> as;
            
            // Nodal values
            vector<vector<double>> nodalSs(nOfNodes);
            vector<vector<double>> nodalS_ts(nOfNodes);
            vector<vector<double>> nodalAs(nOfNodes);
            
            // Results after integration
            vector<double> F0;
            vector<double> F1;

            for (int n = 0; n < nOfNodes; n++) {
                nodalSs[n] = _NID[n]->s;
                nodalS_ts[n] = _NID[n]->s_t;
                nodalAs[n] = _NID[n]->getNodalProperties();
            }
            // Calculate point values at integration points
            for (int i = 0; i < nOfIntPts; i++) {
                for (int j = 0; j < nOfIntPts; j++) {
                    evaluateF(ss, IntPos[i], IntPos[j], nodalSs);
                    evaluateF_x(s_xs, IntPos[i], IntPos[j], nodalSs);
                    evaluateF(as, IntPos[i], IntPos[j], nodalAs);
                    evaluateF(s_ts, IntPos[i], IntPos[j], nodalS_ts);

                    // DEBUG LINES
                    /**
                    if (_ID == 1) {
                        cout << "Element F_X: \n";
                        cout << "(i, j) = " << i << " " << j << "\n";
                        cout << "nodal s: " << nodalSs[0][4] << " " << nodalSs[1][4] << " " << nodalSs[2][4] << " " << nodalSs[3][4] << "\n";
                        cout << "s, s_x, s_t: " << ss[4] << " " << s_xs[8] << " " << s_xs[9]  << " " << s_ts[4] << "\n"; 
                    }
                    */
                    PoroelasticKernel::F0(F0s[i * nOfIntPts + j],
                                          spaceDim,
                                          ss,
                                          s_xs,
                                          s_ts,
                                          s_tshift,
                                          as,
                                          as,
                                          as);
                    
                    PoroelasticKernel::F1(F1s[i * nOfIntPts + j],
                                          spaceDim,
                                          ss,
                                          s_xs,
                                          s_ts,
                                          s_tshift,
                                          as,
                                          as,
                                          as);
                    // DEBUG LINES
                    /**
                    if (_ID == 1) {
                        cout << "f0p: " << F0s[i * nOfIntPts + j][4] << "\n";
                        cout << "f1p: " << F1s[i * nOfIntPts + j][8] << " " << F1s[i * nOfIntPts + j][9] << "\n";
                    }
                    */               
                }
            }
            // Integrate
            IntegratorNf(localF, localFSize, F0s, 1);
            

            IntegratorBf(localF, localFSize, F1s, 1);

            // DEBUG LINES
            /**
            if (_ID == 1) {
                cout << "F0: ";
                for (auto shit : F0) cout << shit << " ";
                cout << "\n";

                cout << "F1: ";
                for (auto shit : F1) cout << shit << " ";
                cout << "\n";
                //cout << "F0: " << F0[4] << " " << F0[12] << " " << F0[20] << " " << F0[28] << "\n";
                //cout << "F1: " << F1[4] << " " << F1[12] << " " << F1[20] << " " << F1[28] << "\n";
            } 
            */

            // Push to globalF
            elementFPush(globalF, localF, localFSize);
        }
        default :
            break;
    }

};

/** Push local Jf to global Jf */
void ElementQ4::JFPush(Mat & globalJF, double *elementJF, int elementJFSize) const {
    // Push to global JF
    MatSetValues(globalJF, elementDOF, localGlobalIndices, elementDOF, localGlobalIndices, elementJF, ADD_VALUES);
    /**
    // Stores global and local index
    int I, J, localI, localJ;
    double this_value;
    // Assemble to globalJF
    // First node
    for (int n1 = 0; n1 < nOfNodes; n1++)
    {
        // cout << "\n";
        // Then row of local JF3
        for (int i = 0; i < nOfDofs; i++)
        {
            I = this->getNID()[n1]->getDOF()[i];
            if (I == -1)
                continue;

            localI = n1 * nOfDofs + i;
            // Then node
            for (int n2 = 0; n2 < nOfNodes; n2++)
            {
                // Then col of local JF3
                for (int j = 0; j < nOfDofs; j++)
                {
                    J = this->getNID()[n2]->getDOF()[j];
                    if (J == -1)
                        continue;
                    localJ = n2 * nOfDofs + j;
                    this_value = JF[localI * nColsJF + localJ];
                    if (abs(this_value) < 1e-15) continue;
                    MatSetValues(globalJF, 1, &(I), 1, &(J), &(JF[localI * nColsJF + localJ]), ADD_VALUES);
                }
            }
        }
    }

    // Debug Lines
    ofstream myFile;
    myFile.open("JF.txt");

    // Element stiffness matrix
    myFile << "Element ID : " << this->getID() <<"\n";
    for (int n1 = 0; n1 < nOfNodes; n1++) {
        for (int n2 = 0; n2 < nOfNodes; n2++) {
            myFile << setw(12) << JF[n1 * nOfDofs * nColsJF+ n2 * nOfDofs] << " " 
                << setw(12) << JF[n1 * nOfDofs * nColsJF + n2 * nOfDofs + 1] << " ";
        }
        myFile << "\n"; 
        for (int n2 = 0; n2 < nOfNodes; n2++) {
            myFile << setw(12) << JF[(n1 * nOfDofs + 1) * nColsJF+ n2 * nOfDofs] << " " 
                << setw(12) << JF[(n1 * nOfDofs + 1) * nColsJF + n2 * nOfDofs + 1] << " ";
        }
        myFile << "\n"; 
    }

    // Element stiffness matrix including all
    myFile << "\n" << "Global JF :" << "\n";
    for (int i = 0; i < nColsJF; i++) {
        for (int j = 0; j < nColsJF; j++) {
            myFile << setw(12) << JF[i * nColsJF + j] - JF[j * nColsJF + i] << " ";
        }
        myFile << "\n";
    }
    myFile.close();
    */
};

/** Push elementF to globalF */
void ElementQ4::elementFPush(Vec & globalF, double *elementF, int elementFSize) const {

    // Push to global vector
    VecSetValues(globalF, elementFSize, localGlobalIndices, elementF, ADD_VALUES);
};