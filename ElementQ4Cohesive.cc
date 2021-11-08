/** @file ElementQ4Cohesive.cc
 * source file for class ElementQ4Cohesive
 */

#include "ElementQ4Cohesive.hh"

//------------------------------------------------------------------------------------------------------------------------------------------------------
/** Methods of class ElementQ4Cohesive
 * No default constructor implemented
 */

/** Constructor 1 for ElementQ4Cohesive */
ElementQ4Cohesive::ElementQ4Cohesive(int ID, const vector<CohesiveNode*> & NID, vector<clock_t> *clocks, vector<double> *timeConsumed) {
    elementDOF = 0;
    setID(ID);
    setNID(NID);

    // Initialize Timers
    this->clocks = clocks;
    this->timeConsumed = timeConsumed;

    // Initialize constants
    spaceDim = _NID[0]->getSpaceDim();
    nOfIntPts = IntPos.size();
    nOfNodes = _NID.size();
    nOfDofs = _NID[0]->getDOF().size();

    // Initialize vectors
    ss.resize(nOfDofs, 0.);
    s_xs.resize(nOfDofs * spaceDim, 0.);
    s_ts.resize(nOfDofs, 0.);
    as.resize(_NID[0]->_nodalProperties.size(), 0.);

    // Collecting nodal data, pre-allocate
    nodalSs.resize(nOfNodes, NULL);
    nodalS_ts.resize(nOfNodes, NULL);
    nodalAs.resize(nOfNodes, NULL);

    // Jacobians and residuals
    isJfAssembled = PETSC_FALSE;
    Jf0s.resize(pow(nOfIntPts, spaceDim - 1));
    Jf1s.resize(pow(nOfIntPts, spaceDim - 1));
    Jf2s.resize(pow(nOfIntPts, spaceDim - 1));
    Jf3s.resize(pow(nOfIntPts, spaceDim - 1));

    F0s.resize(pow(nOfIntPts, spaceDim - 1));
    F1s.resize(pow(nOfIntPts, spaceDim - 1));

    
    // Calculate Bvector, Nvector, pointValue 
    localGlobalIndices = new PetscInt [elementDOF];
    Bvector.resize(pow(nOfIntPts, spaceDim - 1));
    Nvector.resize(pow(nOfIntPts, spaceDim - 1));
    pointValue.resize(pow(nOfIntPts, spaceDim - 1));

    // TO DO: CHANGE HERE
    for (int i = 0; i < nOfIntPts; i++) {
        Bvector[i] = B_x(IntPos[i]);
        Nvector[i] = N(IntPos[i]);
        pointValue[i] = J(IntPos[i]) * IntWs[i];
    }
    int index = 0;
    // Calculate localGlobalIndices
    for (int n = 0; n < nOfNodes; n++) {
        for (int i = 0; i < NID[n]->getDOF().size(); i++) {
            localGlobalIndices[index] = NID[n]->getDOF()[i];
            index += 1;
        }
    }
};

/** Destructor */
ElementQ4Cohesive::~ElementQ4Cohesive() {
    // Delete vector<Node*> _NID
    for (int i = 0; i < _NID.size(); i++) {
        delete _NID[i];
    }

    // Delete elementDOF
    delete [] localGlobalIndices;
};

/** Set element NID */
void ElementQ4Cohesive::setNID(const vector<CohesiveNode*> & NID) {
    _NID.resize(NID.size());
    for (int i = 0; i < NID.size(); i++) {
        _NID[i] = NID[i];
        elementDOF += NID[i]->getDOF().size();
    }
};

/** Get element NID */
const vector<CohesiveNode*> & ElementQ4Cohesive::getNID() const {
    return _NID;
};

/** Shape function N, vector of 2 on 2 nodes, evaluated in base space (ksi) */
vector<double> ElementQ4Cohesive::N(double ksi) {
    return vector<double> {(1. - ksi) / 2., (1. + ksi) / 2.};
};

/** Shape function N, vector of 2 on 2 nodes, 
 * evaluated in base space (ksi)
 * at i, j of elemental shape matrix
 */
double ElementQ4Cohesive::N(const vector<double> & Nvector, int i, int j) const {
    // Some constants    
    if (i == j % nOfDofs) {
        return Nvector[j / nOfDofs];
    }
    return 0.0;
};

/** Gradient of shape function B, vector of 2 * 1 on 2 nodes, evaluated in base space (ksi) */
vector<double> ElementQ4Cohesive::B(double ksi) {
    return vector<double> {- 1. / 2., 1. / 2.};
};

/** Gradient of shape function B, vector of 2 * 2 on 2 nodes, 2 spaceDims, 
 * evaluated in physical space 
 */
vector<double> ElementQ4Cohesive::B_x(double ksi) const {
    double px_pksi = (_NID[1]->getXYZ()[0] - _NID[0]->getXYZ()[0]) / 2.; 
    double py_pksi = (_NID[1]->getXYZ()[1] - _NID[0]->getXYZ()[1]) / 2.;
    vector<double> res(4, 0.);
    if (abs(px_pksi) > 1e-15) {
        res[0] = - 1. / 2. / px_pksi;
        res[2] = 1. / 2. / px_pksi;
    }
    if (abs(py_pksi) > 1e-15) {
        res[1] = - 1. / 2. / py_pksi;
        res[3] = 1. / 2. / py_pksi;
    }
    return res;
};

/** Gradient of shape function B_x[i, j] at (ksi, eta) */
double ElementQ4Cohesive::B_x(const vector<double> & Bvector, int i, int j) const {
    // [\partial N[0] / \partial \ksi, \partial N[0] / \partial \eta, \partial N[1] / \partial \ksi...]
    // Some constants
    if (i / spaceDim == j % nOfDofs) {
        return  Bvector[j / nOfDofs * spaceDim + i % spaceDim];
    }
    return 0.0;
};

/** Jacobian at any given location in base space (ksi, eta), J = d|x| / d ksi */
double ElementQ4Cohesive::J(double ksi) const {
    return sqrt(pow(- _NID[0]->getXYZ()[0] + _NID[1]->getXYZ()[0], 2) + 
                pow(- _NID[0]->getXYZ()[1] + _NID[1]->getXYZ()[1], 2)) / 2.; 
};


/** Evaluate a function F at (i, j)th integration point */
void ElementQ4Cohesive::evaluateF(vector<double> & res, int i,
                                  const vector<vector<double> *> & NodeValues) const {
    // Make sure sizes match
    if (NodeValues.size() != _NID.size()) throw "Not all nodal values are provided for evaluate";
    
    // Initialize result based on input
    if (res.size() != NodeValues[0]->size()) res.resize(NodeValues[0]->size(), 0.);
    for (int i = 0; i < res.size(); i++) res[i] = 0.;

    // Loop through all fields
    for (int f = 0; f < NodeValues[0]->size(); f++) {
        // Loop through all nodes
        for (int n = 0; n < NodeValues.size(); n++) {
            res[f] += Nvector[i][n] * (*(NodeValues[n]))[f];
        }
    }
};


/** Evaluate a function F at (ksi, eta) */
void ElementQ4Cohesive::evaluateF(vector<double> & res, double ksi,
                                  const vector<vector<double> *> & NodeValues) const {
    // Make sure sizes match
    if (NodeValues.size() != _NID.size()) throw "Not all nodal values are provided for evaluate";
    
    // Initialize result based on input
    res.resize(NodeValues[0]->size(), 0.);
    for (int i = 0; i < res.size(); i++) res[i] = 0.;

    // Loop through all fields
    for (int f = 0; f < NodeValues[0]->size(); f++) {
        // Loop through all nodes
        for (int n = 0; n < NodeValues.size(); n++) {
            res[f] += N(ksi)[n] * (*(NodeValues[n]))[f];
        }
    }
};

/** Evaluate a vector at ksi in logical space with given nodal values.
 * Calculated by using shape function to map.
 */
void ElementQ4Cohesive::evaluateF(vector<double> & res, double ksi,
                                  const vector<vector<double>> & NodeValues) const {
    // Make sure sizes match
    if (NodeValues.size() != _NID.size()) throw "Not all nodal values are provided for evaluateF.";
    
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

/** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at (ksi, eta) 
 * in LOGICAL space with given nodal values.
 * Calculated by using shape function to map
 */
void ElementQ4Cohesive::evaluateF_x(vector<double> & res, int i, 
                                    const vector<vector<double> *> & NodeValues) const {
    // Number of nodes
    // int nOfNodes = this->getNID().size();
    if (NodeValues.size() != nOfNodes) 
        throw "In evaluateF_x, nodeValues do not match number of nodes!";
    
    // Space dim is 2 for Q4
    // int spaceDim = 2;

    // Set res to 0.;
    if (res.size() != spaceDim * NodeValues[0]->size()) res.resize(spaceDim * NodeValues[0]->size());
    fill(res.begin(), res.end(), 0.0);

    // Loop through fields
    for(int f = 0; f < NodeValues[0]->size(); f++) {
        // Loop through nodes
        for (int n = 0; n < NodeValues.size(); n++) {
            // Loop through spaceDim
            for (int d = 0; d < spaceDim; d++) {
                // Calculate \partial F / \partial \ksi
                res[d + f * spaceDim] += Bvector[i][n * spaceDim + d] * (*(NodeValues[n]))[f];
            }
        }
    }
};

/** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at ksi
 * in LOGICAL space with given nodal values.
 * Calculated by using shape function to map
 */
void ElementQ4Cohesive::evaluateF_x(vector<double> & res, double ksi, 
                                    const vector<vector<double> *> & NodeValues) const {
    // Number of nodes
    // int nOfNodes = this->getNID().size();
    if (NodeValues.size() != nOfNodes) 
        throw "In evaluateF_x, nodeValues do not match number of nodes!";

    // Set res to 0.;
    res.resize(spaceDim * NodeValues[0]->size());
    fill(res.begin(), res.end(), 0.0);

    // Loop through fields
    for(int f = 0; f < NodeValues[0]->size(); f++) {
        // Loop through nodes
        for (int n = 0; n < NodeValues.size(); n++) {
            // Loop through spaceDim
            for (int d = 0; d < spaceDim; d++) {
                // Calculate \partial F / \partial x
                res[d + f * spaceDim] += B_x(ksi)[n * spaceDim + d] * (*(NodeValues[n]))[f];
            }
        }
    }
};

/** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at (ksi) 
 * in LOGICAL space with given nodal values.
 * Calculated by using shape function to map
 */
void ElementQ4Cohesive::evaluateF_x(vector<double> & res, double ksi,
                        const vector<vector<double>> & NodeValues) const {
    // Number of nodes
    int nOfNodes = this->getNID().size();
    if (NodeValues.size() != nOfNodes) 
        throw "In evaluateF_x, nodeValues do not match number of nodes!";
    
    // Space dim is 2 for Q4Cohesive
    int spaceDim = 2;

    // Set res to 0.;
    res.resize(spaceDim * NodeValues[0].size());
    fill(res.begin(), res.end(), 0.0);

    // Loop through fields
    for(int f = 0; f < NodeValues[0].size(); f++) {
        // Loop through nodes
        for (int n = 0; n < NodeValues.size(); n++) {
            for (int d = 0; d < spaceDim; d++) {
                // Calculate \partial F / \partial x
                res[f] += B_x(ksi)[n * spaceDim + d] * NodeValues[n][f];
            }
        }
    }
};

// ================ Integrators of ElementQ4Cohesive ========================================
/** IntegratorNf, integrates a vector input inside an element, 
 * left side using shape function
 * RES: nOfNodes * nOfDofs
 * NODEVALUES: dim 1, nOfNodes; dim 2, nOfDofs
 * FLAG: 0 - NodeValues are given at nodes,
 *       1 - NodeValues are given at integration points
 */
void ElementQ4Cohesive::IntegratorNf(double *res, 
                                     int resSize, 
                                     const vector<vector<double>> & NodeValues, 
                                     int flag) const {
    // Some constants
    if (NodeValues.size() != nOfNodes) throw "Not all nodes are provided for ElementQ4 IntegratorNf!";
    if (NodeValues[0].size() != nOfDofs) throw "Not all nodal DOFs are provided for ElementQ4 IntegratorNf!";
    if (resSize != nOfNodes * nOfDofs) throw "ResSize error for ElementQ4 IntegratorNf!";

    // Integration point index
    int IntPtIndex;

    // Pre-Calculate values of N and f at integration points
    if (flag == 1) {
        // Computing the integral
        // res_i = N_{p, i} f_p        
        for (int i = 0; i < nOfNodes * nOfDofs; i++) {
            for (int p = 0; p < nOfDofs; p++) {
                if (p != i % nOfDofs) continue;
                // All integration points
                for (int k = 0; k < nOfIntPts; k++) {
                    IntPtIndex = k;
                    res[i] += Nvector[IntPtIndex][i / nOfDofs]
                                * pointValue[IntPtIndex] 
                                * NodeValues[IntPtIndex][p];
                }
            }
        }
    }
    else {
        vector<vector<double>> fvector(pow(nOfIntPts, spaceDim - 1));
        for (int i = 0; i < nOfIntPts; i++) {
            evaluateF(fvector[i], IntPos[i], NodeValues);                
        }        

        // Computing the integral
        // res_i = N_{p, i} f_p        
        for (int i = 0; i < nOfNodes * nOfDofs; i++) {
            for (int p = 0; p < nOfDofs; p++) {
                if (p != i % nOfDofs) continue;
                // All integration points
                for (int k = 0; k < nOfIntPts; k++) {
                    IntPtIndex = k;
                    res[i] += Nvector[IntPtIndex][i / nOfDofs]
                                * pointValue[IntPtIndex] 
                                * fvector[IntPtIndex][p];
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
    if (NodeValues.size() != nOfNodes) throw "Not all nodes are provided for ElementQ4 IntegratorBf!";
    if (NodeValues[0].size() != nOfDofs * spaceDim) throw "Not all nodal DOFs are provided for ElementQ4 IntegratorBf!";
    if (resSize != nOfDofs * nOfNodes) throw "resSize error for ElementQ4 IntegratorBf!";

    // Integration point index
    int IntPtIndex;

    if (flag == 1) {
        // Computing the integral
        // res_i = B_{p, i} f_p        
        for (int i = 0; i < nOfNodes * nOfDofs; i++) {
            for (int p = 0; p < nOfDofs * spaceDim; p++) {
                if (p / spaceDim != i % nOfDofs) continue;
                // All integration points
                for (int k = 0; k < nOfIntPts; k++) {
                    IntPtIndex = k;
                    res[i] += Bvector[IntPtIndex][i / nOfDofs * spaceDim + p % spaceDim] 
                            * pointValue[IntPtIndex] 
                            * NodeValues[IntPtIndex][p];
                }
            }
        }
    }
    else {
        // Pre-Calculate values of N and f at integration points
        vector<vector<double>> fvector(pow(nOfIntPts, spaceDim));
        for (int i = 0; i < nOfIntPts; i++) {
            for (int j = 0; j < nOfIntPts; j++) {
                evaluateF(fvector[i * nOfIntPts + j], IntPos[i], IntPos[j], NodeValues);                   
            }
        }

        // Computing the integral
        // res_i = B_{p, i} f_p
        
        for (int i = 0; i < nOfNodes * nOfDofs; i++) {
            for (int p = 0; p < nOfDofs * spaceDim; p++) {
                if (p / spaceDim != i % nOfDofs) continue;
                // All integration points
                for (int k = 0; k < nOfIntPts; k++) {
                    IntPtIndex = k;
                    res[i] += Bvector[IntPtIndex][i / nOfDofs * spaceDim + p % spaceDim] 
                            * pointValue[IntPtIndex] 
                            * fvector[IntPtIndex][p];
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
void ElementQ4Cohesive::IntegratorNfN(double *res,
                                      int resSize,
                                      const vector<vector<double>> &NodeValues,
                                      const vector<int> &f_is,
                                      const vector<int> &f_js,
                                      int flag) const {
    // Some constants
    // int nOfNodes = this->getNID().size();
    // int nOfIntPts = IntPos.size();
    // int nOfDofs = getNID()[0]->getDOF().size();
    int nOfColsRes = nOfNodes * nOfDofs;

    if (resSize != pow(nOfColsRes, 2)) throw "resSize error for ElementQ4 IntegratorNfN!";
    if (NodeValues.size() != nOfNodes) throw "Not all nodal values are provided for ElementQ4 IntegratorNfN!";
    if (NodeValues[0].size() != nOfDofs * nOfDofs) throw "Mass matrix size not compatible with Element Q4 IntegratorNfN!";
    
    int intPtIndex;
    int resIJindex;
    int p, q;
    if (flag == 1) {        
        // Calculate res[i,j] = N_{pi} M_{pq} N_{qj}
        for (int i = 0; i < nOfColsRes; i++) {
            for (int j = 0; j < nOfColsRes; j++) {
                resIJindex = nOfColsRes * i + j;
                for (int k = 0; k < nOfIntPts; k++) {
                    intPtIndex = k;
                    for (int temp = 0; temp < f_is.size(); temp++) {
                        p = f_is[temp];
                        q = f_js[temp];
                        if (p != i % nOfDofs || q != j % nOfDofs) continue;
                        res[resIJindex] += Nvector[intPtIndex][i / nOfDofs] 
                                            * NodeValues[intPtIndex][p * nOfDofs + q] 
                                            * Nvector[intPtIndex][j / nOfDofs]
                                            * pointValue[intPtIndex];
                    }
                }                
            }
        }
    }
    else {
        // Pre-calculate and store some values
        vector<vector<double>> pointM(nOfIntPts);
        for (int i = 0; i < nOfIntPts; i++) {
            intPtIndex = i;
            evaluateF(pointM[intPtIndex], IntPos[i], NodeValues);
        }

        // Calculate res[i,j] = N_{pi} M_{pq} N_{qj}
        for (int i = 0; i < nOfColsRes; i++) {
            for (int j = 0; j < nOfColsRes; j++) {
                resIJindex = nOfColsRes * i + j;
                for (int k = 0; k < nOfIntPts; k++) {
                    intPtIndex = k;
                    for (int temp = 0; temp < f_is.size(); temp++) {
                        p = f_is[temp];
                        q = f_js[temp];
                        if (p != i % nOfDofs || q != j % nOfDofs) continue;
                        res[resIJindex] += Nvector[intPtIndex][i / nOfDofs] 
                                            * pointM[intPtIndex][p * nOfDofs + q] 
                                            * Nvector[intPtIndex][j / nOfDofs]
                                            * pointValue[intPtIndex];
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
void ElementQ4Cohesive::IntegratorBfB(double *res,
                                      int resSize,
                                      const vector<vector<double>> &NodeValues,
                                      const vector<int> &f_is,
                                      const vector<int> &f_js,
                                      int flag) const {
    // Some constants
    int nOfColsF = nOfDofs * spaceDim;
    int nOfColsRes = nOfDofs * nOfNodes;

    if (resSize != nOfColsRes * nOfColsRes) throw "resSize error in ElementQ4 IntegratorBfB!";
    
    if (NodeValues.size() != nOfNodes) 
        throw "Not all nodal values are provided for ElementQ4 IntegratorBfB!";
    
    if (NodeValues[0].size() != nOfColsF * nOfColsF) 
        throw "Input f is not compatible with nodal dofs in ElementQ4 IntegratorBfB!";
    
    // Some indices
    int intPtIndex = 0;
    int resIJindex = 0;
    int p, q;
    // Pre-store integration constants at different points
    if (flag == 1) {
        // Calculate res[i,j]
        // In the integral (B^T D B)_{i,j} = B_pi D_pq B qj
        for (int i = 0; i < nOfColsRes; i++) {
            for (int j = 0; j < nOfColsRes; j++) {
                resIJindex = i * nOfColsRes + j;
                for (int k = 0; k < nOfIntPts; k++) {
                    intPtIndex = k;
                    for (int temp = 0; temp < f_is.size(); temp++) {
                        p = f_is[temp];
                        q = f_js[temp];
                        if (p / spaceDim != i % nOfDofs || q / spaceDim != j % nOfDofs) continue;
                        res[resIJindex] += Bvector[intPtIndex][i / nOfDofs * spaceDim + p % spaceDim] 
                                            * NodeValues[intPtIndex][p * nOfColsF + q] 
                                            * Bvector[intPtIndex][j / nOfDofs * spaceDim + q % spaceDim]
                                            * pointValue[intPtIndex];
                    }
                }                
            }
        }
    }
        
    else {
        vector<vector<double>> pointD(nOfIntPts, vector<double>(NodeValues[0].size(), 0.));
        // Calculate Bvector, pointValue and pointD
        for (int i = 0; i < nOfIntPts; i++) {
            evaluateF(pointD[i], IntPos[i], NodeValues);
        }

        // Calculate res[i,j]
        // In the integral (B^T D B)_{i,j} = B_pi D_pq B qj
        for (int i = 0; i < nOfColsRes; i++) {
            for (int j = 0; j < nOfColsRes; j++) {
                resIJindex = i * nOfColsRes + j;
                for (int k = 0; k < nOfIntPts; k++) {
                    intPtIndex = k;
                    for (int temp = 0; temp < f_is.size(); temp++) {
                        p = f_is[temp];
                        q = f_js[temp];
                        if (p / spaceDim != i % nOfDofs || q / spaceDim != j % nOfDofs) continue;
                        res[resIJindex] += Bvector[intPtIndex][i / nOfDofs * spaceDim + p % spaceDim] 
                                            * pointD[intPtIndex][p * nOfColsF + q] 
                                            * Bvector[intPtIndex][j / nOfDofs * spaceDim + q % spaceDim]
                                            * pointValue[intPtIndex];
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
void ElementQ4Cohesive::IntegratorBfN(double *res,
                                      int resSize,
                                      const vector<vector<double>> &NodeValues,
                                      const vector<int> &f_is,
                                      const vector<int> &f_js,
                                      int flag) const {
    // Some constants
    int nRowsF = spaceDim * nOfDofs;
    int nColsF = nOfDofs;
    int nColsRes = nOfDofs * nOfNodes; 
    if (resSize != nColsRes * nColsRes) throw "resSize error in ElementQ4Cohesive IntegratorBfN!";
    if (NodeValues.size() != nOfNodes) throw "Not all nodal values are provided for ElementQ4Cohesive IntegratorBfN!";
    if (NodeValues[0].size() != nRowsF * nColsF) throw "Nodal matrix provided for ElementQ4Cohesive not compatible with IntegratorBfN!";
    
    // Integration point index
    int intPtIndex;
    // Stores IJ index in the vector
    int resIJindex;

    int p, q;

    if (flag == 1) {
        // Calculate res[i,j]
        // In the integral (B^T F N)_{i,j} = B_pi F_pq N qj
        for (int i = 0; i < nColsRes; i++) {
            for (int j = 0; j < nColsRes; j++) {
                resIJindex = i * nColsRes + j;
                // Loop through integration points
                for (int k = 0; k < nOfIntPts; k++) {
                    // Constants for the integral
                    intPtIndex = k;
                    for (int temp = 0; temp < f_is.size(); temp++) {
                        p = f_is[temp];
                        q = f_js[temp];
                        if (p / spaceDim != i % nOfDofs || q != j % nOfDofs) continue;
                        res[resIJindex] += pointValue[intPtIndex] 
                                            * Bvector[intPtIndex][i / nOfDofs * spaceDim + p % spaceDim]
                                            * NodeValues[intPtIndex][p * nColsF + q]
                                            * Nvector[intPtIndex][j / nOfDofs];
                    }         
                }                
            }
        }
    }
    else {
        // Pre-calculate and store some values
        vector<vector<double>> pointF(nOfIntPts);
        for (int i = 0; i < nOfIntPts; i++) {
            intPtIndex = i;            
            evaluateF(pointF[intPtIndex], IntPos[i], NodeValues);
        }

        // Calculate res[i,j]
        // In the integral (B^T F N)_{i,j} = B_pi F_pq N qj
        for (int i = 0; i < nColsRes; i++) {
            for (int j = 0; j < nColsRes; j++) {
                resIJindex = i * nColsRes + j;
                // Loop through integration points
                for (int k = 0; k < nOfIntPts; k++) {
                    // Constants for the integral
                    intPtIndex = k;
                    for (int temp = 0; temp < f_is.size(); temp++) {
                        p = f_is[temp];
                        q = f_js[temp];
                        if (p / spaceDim != i % nOfDofs || q != j % nOfDofs) continue;
                        res[resIJindex] += pointValue[intPtIndex] 
                                            * Bvector[intPtIndex][i / nOfDofs * spaceDim + p % spaceDim]
                                            * pointF[intPtIndex][p * nColsF + q]
                                            * Nvector[intPtIndex][j / nOfDofs];
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
void ElementQ4Cohesive::IntegratorNfB(double *res,
                                      int resSize,
                                      const vector<vector<double>> &NodeValues,
                                      const vector<int> &f_is,
                                      const vector<int> &f_js,
                                      int flag) const {
    // Some constants
    int nRowsF = nOfDofs;
    int nColsF = spaceDim * nOfDofs;
    int nColsRes = nOfDofs * nOfNodes; 
    if (resSize != nColsRes * nColsRes) throw "resSize error for ElementQ4 IntegratorNfB!";
    if (NodeValues.size() != nOfNodes) throw "Not all nodal values are provided for ElementQ4 IntegratorNfB!";
    if (NodeValues[0].size() != nRowsF * nColsF) throw "Nodal matrix provided not compatible with IntegratorNfB!";

    // Integration points index
    int intPtIndex;

    // Stores IJ index in the vector
    int resIJindex;
    int p, q;

    if (flag == 1) {
        // Calculate res[i,j]
        // In the integral (N^T F B)_{i,j} = N_pi F_pq B qj
        for (int i = 0; i < nColsRes; i++) {
            for (int j = 0; j < nColsRes; j++) {
                resIJindex = i * nColsRes + j;
                // Loop through integration points
                for (int k = 0; k < nOfIntPts; k++) {
                    // Constants for the integral
                    intPtIndex = k;
                    for (int temp = 0; temp < f_is.size(); temp++) {
                        p = f_is[temp];
                        q = f_js[temp];
                        if (p != i % nOfDofs || q / spaceDim != j % nOfDofs) continue;
                        res[resIJindex] += pointValue[intPtIndex] 
                                            * Nvector[intPtIndex][i / nOfDofs]
                                            * NodeValues[intPtIndex][p * nColsF + q]
                                            * Bvector[intPtIndex][j / nOfDofs * spaceDim + q % spaceDim];
                    }
                }                
            }
        }
    } 
    else {
        // Pre-calculate and store some values
        vector<vector<double>> pointF(nOfIntPts);
        for (int i = 0; i < nOfIntPts; i++) {
            intPtIndex = i;
            evaluateF(pointF[intPtIndex], IntPos[i], NodeValues);
        }

        // Calculate res[i,j]
        // In the integral (N^T F B)_{i,j} = N_pi F_pq B qj
        for (int i = 0; i < nColsRes; i++) {
            for (int j = 0; j < nColsRes; j++) {
                resIJindex = i * nColsRes + j;
                // Loop through integration points
                for (int k = 0; k < nOfIntPts; k++) {
                    // Constants for the integral
                    intPtIndex = k;
                    for (int temp = 0; temp < f_is.size(); temp++) {
                        p = f_is[temp];
                        q = f_js[temp];
                        if (p != i % nOfDofs || q / spaceDim != j % nOfDofs) continue;
                        res[resIJindex] += pointValue[intPtIndex] 
                                            * Nvector[intPtIndex][i / nOfDofs]
                                            * pointF[intPtIndex][p * nColsF + q]
                                            * Bvector[intPtIndex][j / nOfDofs * spaceDim + q % spaceDim];
                    }
                }                
            }
        }
    }
};

/** Output element info */
void ElementQ4Cohesive::outputInfo(ofstream & myFile) const {
    myFile << setw(10) << getID() << " ";
    for (int i = 0; i < getNID().size(); i++) {
        myFile << setw(10) << getNID()[i]->getID() << " ";
    }
    myFile << "\n";
};