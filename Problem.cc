/** @file Problem.cc
 * Source file of class Problem.
 * Initialize, solve and write the results.
 */
#include "Problem.hh"
// Constructor
Problem::Problem(int spaceDim) {
    _spaceDim = spaceDim;
    nodeTime = 0.;
    slip.resize(2, 0.0);
};

// Destructor
Problem::~Problem() {
    if (myGeometry) delete myGeometry;

    // Release Nodes
    deleteNodes();

    // Release Elements
    deleteElements();

    // Delete pre-allocated element residual and jacobian array.
    if (localF) delete [] localF;
    if (localJF) delete [] localJF;

    // Destroy global vectors and matrices
    VecDestroy(&globalF);
    VecDestroy(&globalS);
    VecDestroy(&globalS_t);

    // Mat destroy
    MatDestroy(&globalJF);

    // Delete rows
    if (globalRows) delete [] globalRows;
};

// General initialization
void Problem::initialize(const vector<double> & xRanges, const vector<int> & edgeNums) {
    
    // Initialize geometry2D
    if (_spaceDim == 2) initializeGeometry2D(xRanges, edgeNums);
    
    // Initialize nodes
    initializeNodes();
    
    // Assign Nodal DOFs
    assignNodalDOF();
    
    // Test push global F
    testPushGlobalF();
    
    // Test fetch global F
    testFetchGlobalF();
    
    // Initialize elements
    initializeElements();
    
    // Test integrators
    testIntegrators();
};

// Initialize Geometry2D
void Problem::initializeGeometry2D(const vector<double> & xRanges, const vector<int> & edgeNums) {
    myGeometry = new Geometry2D(xRanges[0], xRanges[1], edgeNums[0], edgeNums[1]);
};

// Initialize nodes, !!currently only 2-dim nodes are available
void Problem::initializeNodes() {
    /**
    // Some input geometrical values
    int spaceDim = _spaceDim;

    // Number of elements and nodes
    // Parameters for bulk Nodes
    double lambda = 1.;
    double shearModulus = 1.;
    double biotAlpha = 0.5;
    double biotMp = 1.;
    double fluidMobility = 1.;
    double fluidViscosity = 1.;
    // Element size
    vector<double> edgeSize (spaceDim);
    edgeSize[0] = myGeometry->xRange / myGeometry->xEdgeNum;
    edgeSize[1] = myGeometry->yRange / myGeometry->yEdgeNum;
    // DOF of fixed velocity and pressure
    vector<int> DOF_vp_fixed (2 * spaceDim + 2, 0); 
    for (int i = spaceDim; i < 2 * spaceDim; i++) DOF_vp_fixed[i] = 1;
    DOF_vp_fixed[2 * spaceDim] = 1; 
    DOF_vp_fixed[2 *spaceDim + 1] = 1;

    // Default DOF
    vector<int> DOF_default (2 * spaceDim + 2, 0);
    DOF_default[2 * spaceDim + 1] = 1;
    
    // Mass Density and Body force.
    double massDensity = 1.0;
    vector<double> bodyForce  = {0, -10.0};
    
    // Upper subzone nodes
    upperNodes.resize(myGeometry->nOfNodes);
    int nodeID = 0;
    int nodeID_in_set = 0;
    vector<double> thisXYZ(spaceDim);
    vector<double> initialS(spaceDim * 2 + 2, 1.);
    // First y
    for (int j = 0; j < myGeometry->yNodeNum; j++) {
        // Then x
        for (int i = 0; i < myGeometry->xNodeNum; i++) {
            // Reset the coordinates
            thisXYZ[0] = i * edgeSize[0];
            thisXYZ[1] = j * edgeSize[1];
            massDensity = thisXYZ[0] + thisXYZ[1];
            bodyForce  = {-thisXYZ[0], -thisXYZ[1]};
            // Initialize a node
            upperNodes[nodeID_in_set] = 
                new Node(nodeID, thisXYZ, DOF_default, spaceDim, 
                         massDensity, 
                         &bodyForce, 
                         lambda, 
                         shearModulus, 
                         biotAlpha, 
                         biotMp, 
                         fluidMobility, 
                         fluidViscosity);
            if (j == myGeometry->yNodeNum - 1) // Upper surface
                upperNodes[nodeID_in_set]->setDOF(DOF_vp_fixed);
            upperNodes[nodeID_in_set]->initializeS(initialS);
            nodeID += 1;
            nodeID_in_set += 1;
            
            
            
        }
    }

    fill(initialS.begin(), initialS.end(), 2.0);
    // Lower subzone nodes
    lowerNodes.resize(myGeometry->nOfNodes);
    // Reset in-set node ID
    nodeID_in_set = 0;
    // First y
    for (int j = 0; j < myGeometry->yNodeNum; j++) {
        // Then x
        for (int i = 0; i < myGeometry->xNodeNum; i++) {
            // Reset the coordinates
            thisXYZ[0] = i * edgeSize[0];
            thisXYZ[1] = - j * edgeSize[1];
            massDensity = thisXYZ[0] + thisXYZ[1];
            bodyForce  = {-thisXYZ[0], -thisXYZ[1]};
            // Initialize a node
            lowerNodes[nodeID_in_set] = 
                new Node(nodeID, thisXYZ, DOF_default, spaceDim, 
                         massDensity, 
                         &bodyForce, 
                         lambda, 
                         shearModulus, 
                         biotAlpha, 
                         biotMp, 
                         fluidMobility, 
                         fluidViscosity);
            
            if (j == myGeometry->yNodeNum - 1) // Lower surface
                lowerNodes[nodeID_in_set]->setDOF(DOF_vp_fixed);
            lowerNodes[nodeID_in_set]->initializeS(initialS);
            nodeID += 1;
            nodeID_in_set += 1;
            
        }
    }
    nodeID_in_set = 0;
    
    // Cohesive Nodes;
    cohesiveNodes.resize(myGeometry->xNodeNum);
    initialS.resize(spaceDim + 2);
    fill(initialS.begin(), initialS.end(), 3.0);
    
    // Properly set the DOFs and lowerUpperNodes
    vector<int> DOF_cohesive_default(2 * (2 * spaceDim + 2) + spaceDim + 2, 2);
    for (int i = 2 * (2 * spaceDim + 2); i < 2 * (2 * spaceDim + 2) + spaceDim + 2; i++) DOF_cohesive_default[i] = 0;
    vector<Node*> lowerUpperNodes(2, NULL);
    
    for (int i = 0; i < myGeometry->xNodeNum; i++) {
        // lowerUpperNodes
        lowerUpperNodes[0] = lowerNodes[i];
        lowerUpperNodes[1] = upperNodes[i];

        thisXYZ[0] = i * edgeSize[0];
        thisXYZ[1] = 0.;
        massDensity = thisXYZ[0] + thisXYZ[1];
        bodyForce  = {-thisXYZ[0], -thisXYZ[1]};
        cohesiveNodes[nodeID_in_set] = new CohesiveNode(nodeID, thisXYZ, DOF_cohesive_default, lowerUpperNodes, spaceDim);
        cohesiveNodes[nodeID_in_set]->setMassDensity(massDensity);
        cohesiveNodes[nodeID_in_set]->setBodyForce(&bodyForce);
        cohesiveNodes[nodeID_in_set]->initializeS(initialS);
        nodeID += 1;
        nodeID_in_set += 1;
        
    };
    _totalNofNodes = nodeID;
    ofstream myFile;
    myFile.open("Testlog.txt");
    myFile << "=================== NodeInfoBefore ======================================" << "\n";
    for (Node* node : upperNodes) node->outputInfo(myFile, true);
    for (Node* node : lowerNodes) node->outputInfo(myFile, true);
    for (CohesiveNode* node : cohesiveNodes) node->outputInfo(myFile, true);
    myFile.close();
    */
};

// Assign Nodal DOFs;
void Problem::assignNodalDOF() {
    // Pointer to currentDOF
    _totalDOF = 0;
    
    // Assign upperzone
    for (int i = 0; i < upperNodes.size(); i++) {
        for (int j = 0; j < upperNodes[i]->getDOF().size(); j++) {
            if (upperNodes[i]->getDOF(j) == 0) {
                upperNodes[i]->setDOF(j, _totalDOF);
                _totalDOF += 1;
            }
            else {
                upperNodes[i]->setDOF(j, -1); 
            }
        }
    }

    // Assign lowerzone
    for (int i = 0; i < lowerNodes.size(); i++) {
        for (int j = 0; j < lowerNodes[i]->getDOF().size(); j++) {
            if (lowerNodes[i]->getDOF(j) == 0) {
                lowerNodes[i]->setDOF(j, _totalDOF);
                _totalDOF += 1;
            }
            else {
                lowerNodes[i]->setDOF(j, -1); 
            }
        }
    }

    // Assign cohesivezone
    for (int i = 0; i < cohesiveNodes.size(); i++) {
        for (int j = 0; j < cohesiveNodes[i]->getDOF().size(); j++) {
            if (cohesiveNodes[i]->getDOF(j) == 0) {
                cohesiveNodes[i]->setDOF(j, _totalDOF);
                _totalDOF += 1;
            }
            else if (cohesiveNodes[i]->getDOF(j) == 1){
                cohesiveNodes[i]->setDOF(j, -1); 
            }
        }
        cohesiveNodes[i]->updateDOF();
    }
};

// TEST Try push to globalF
void Problem::testPushGlobalF() {
    // Set size of global F
    VecCreate(PETSC_COMM_WORLD, &globalF);
    VecSetSizes(globalF, PETSC_DECIDE, _totalDOF);
    VecSetFromOptions(globalF);

    // Upper zone
    for (Node* node : upperNodes) {
        node->pushS(globalF);
    }
    // lower zone
    for (Node* node : lowerNodes) {
        node->pushS(globalF);
    }
    // Cohesive zone
    for (CohesiveNode* node : cohesiveNodes) {
        node->pushS(globalF);
    }
    // Debug lines
    // Output the global F
    //cout << "TEST: Global F\n";
    // VecView(globalF, PETSC_VIEWER_STDOUT_SELF);
};

// TEST Try Fetching from globalF
void Problem::testFetchGlobalF() {
    // Upper zone
    for (Node* node : upperNodes) {
        node->fetchS_t(globalF);
    }
    // lower zone
    for (Node* node : lowerNodes) {
        node->fetchS_t(globalF);
    }
    // Cohesive zone
    for (CohesiveNode* node : cohesiveNodes) {
        node->fetchS_t(globalF);
    }
};

// Initialize Elements, !! CURRENTLY only ElementQ4s are available
void Problem::initializeElements() {
    // Using the NID to assign values to each element
    vector<Node*> NID(4, NULL);
    
    // Subzone ElementQ4s
    upperElements.resize(myGeometry->nOfElements, NULL);
    lowerElements.resize(myGeometry->nOfElements, NULL);

    for (int i = 0; i < myGeometry->xEdgeNum; i++) {
        for (int j = 0; j < myGeometry->yEdgeNum; j++) {
            // Setting NID for upper subzone
            NID = {upperNodes[j * myGeometry->xNodeNum + i], 
                   upperNodes[j * myGeometry->xNodeNum + i + 1], 
                   upperNodes[(j + 1) * myGeometry->xNodeNum + i + 1], 
                   upperNodes[(j + 1) * myGeometry->xNodeNum + i]};

            upperElements[j * myGeometry->xEdgeNum + i] = 
                new ElementQ4(j * myGeometry->xEdgeNum + i, NID);

            // Setting NID for lower subzone
            NID = {lowerNodes[(j + 1) * myGeometry->xNodeNum + i], 
                   lowerNodes[(j + 1) * myGeometry->xNodeNum + i + 1], 
                   lowerNodes[j * myGeometry->xNodeNum + i + 1], 
                   lowerNodes[j * myGeometry->xNodeNum + i]};

            lowerElements[j * myGeometry->xEdgeNum + i] = 
                new ElementQ4(j * myGeometry->xEdgeNum + i + myGeometry->nOfElements, NID);
            
        }
    }
    
    // Cohesive elements
    int cohesiveElementST = 2 * myGeometry->nOfElements;
    vector<CohesiveNode*> NID_cohesive(2);
    cohesiveElements.resize(myGeometry->xEdgeNum);
    for (int i = 0; i < myGeometry->xEdgeNum; i++) {
        NID_cohesive = {cohesiveNodes[i], 
               cohesiveNodes[i + 1]};
        cohesiveElements[i] = 
            new ElementQ4Cohesive(cohesiveElementST + i, NID_cohesive); 
    }

    ofstream myFile;
    myFile.open("Testlog.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== ElementInfo ======================================" << "\n";
    for (ElementQ4* thisElement : upperElements) thisElement->outputInfo(myFile);
    for (ElementQ4* thisElement : lowerElements) thisElement->outputInfo(myFile);
    for (ElementQ4Cohesive* thisElement : cohesiveElements) thisElement->outputInfo(myFile);
    myFile.close();
};

// Delete Nodes;
void Problem::deleteNodes() {
    // Release upper node pointers
    for (int i = 0; i < upperNodes.size(); i++) {
        delete upperNodes[i];
    }

    // Release lower node pointers
    for (int i = 0; i < lowerNodes.size(); i++) {
        delete lowerNodes[i];
    }

    // Release cohesive node pointers
    for (int i = 0; i < cohesiveNodes.size(); i++) {
        delete cohesiveNodes[i];
    }
}

// Delete Elements
void Problem::deleteElements() {
    // Release upper element pointers
    for (int i = 0; i < upperElements.size(); i++) {
        delete upperElements[i];
    }

    // Release lower element pointers
    for (int i = 0; i < lowerElements.size(); i++) {
        delete lowerElements[i];
    }

    // Release cohesive element pointers
    for (int i = 0; i < cohesiveElements.size(); i++) {
        delete cohesiveElements[i];
    }
};

// Calculate Nodal bodyForces, NOW for TESTING INTEGRATORS
void Problem::computeBodyForces() {
    /** 
    vector<vector<double>> eleBodyForces(4, vector<double> (_spaceDim, 0.));
    vector<vector<double>> nodalValues(4, vector<double> (_spaceDim, 0.));

    // Loop through upperZone elements
    for (ElementQ4* element : upperElements) {
        // Set the function-to-be-integrated
        for (int i = 0; i < element->getNID().size(); i++) {
            nodalValues[i] = {element->getNID()[i]->getMassDensity() 
                                * element->getNID()[i]->getBodyForce()[0],
                              element->getNID()[i]->getMassDensity() 
                                * element->getNID()[i]->getBodyForce()[1]};
        }
        // Call integrator
        element->IntegratorNf(eleBodyForces, nodalValues);
        
        // Add to the nodal bodyforces
        for (int i = 0; i < element->getNID().size(); i++) {
            for (int j = 0; j < _spaceDim; j++) {
                element->getNID()[i]->_nodalBodyForce[j] += eleBodyForces[i][j];
            }
        }
    }

    // Loop through lowerZone elements
    for (ElementQ4* element : lowerElements) {
        // Set the function-to-be-integrated
        for (int i = 0; i < element->getNID().size(); i++) {
            nodalValues[i] = {element->getNID()[i]->getMassDensity() 
                                * element->getNID()[i]->getBodyForce()[0],
                              element->getNID()[i]->getMassDensity() 
                                * element->getNID()[i]->getBodyForce()[1]};
        }
        // Call integrator
        element->IntegratorNf(eleBodyForces, nodalValues);
        
        // Add to the nodal bodyforces
        for (int i = 0; i < element->getNID().size(); i++) {
            for (int j = 0; j < _spaceDim; j++) {
                element->getNID()[i]->_nodalBodyForce[j] += eleBodyForces[i][j];
            }
        }
    }

    // Loop through cohesive nodes
    eleBodyForces.resize(2);
    nodalValues.resize(2);

    for (ElementQ4Cohesive* element : cohesiveElements) {
        // Set the function-to-be-integrated
        for (int i = 0; i < element->getNID().size(); i++) {
            nodalValues[i] = {element->getNID()[i]->getMassDensity() 
                                * element->getNID()[i]->getBodyForce()[0],
                              element->getNID()[i]->getMassDensity() 
                                * element->getNID()[i]->getBodyForce()[1]};
        }
        // Call integrator
        element->IntegratorNf(eleBodyForces, nodalValues);
        
        // Add to the nodal bodyforces
        for (int i = 0; i < element->getNID().size(); i++) {
            for (int j = 0; j < _spaceDim; j++) {
                element->getNID()[i]->_nodalBodyForce[j] += eleBodyForces[i][j];
            }
        }
    }
    // Output to log file
    ofstream myFile;
    myFile.open("Testlog.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== NodeInfoAfter ======================================" << "\n";

    for (Node* node : upperNodes) node->outputInfo(myFile, true);
    for (Node* node : lowerNodes) node->outputInfo(myFile, true);
    for (CohesiveNode* node : cohesiveNodes) node->outputInfo(myFile, true);
    myFile.close();
    */
};

// Test all integrators
void:: Problem::testIntegrators() {
    // Test: compute body forces
    computeBodyForces();

    // Test: test gradient function
    // testEvaluateF_x();

    // Test: global mass matrix
    testIntegratorNfN();

    // Test: global stiffness matrix
    testIntegratorBfB();

    // Test: BfN
    testIntegratorBfN();

    // Test: NfB
    testIntegratorNfB();
};

/** Test gradient operator */
void Problem::testEvaluateF_x() const {
    vector<vector<double>> NodeValues(4, vector<double>(3, 0.));
    vector<double> res(_spaceDim * 3, 0.);
    vector<double> IntPos = {- 1./sqrt(3.), 1./sqrt(3.)};
    
    // Loop through upper elements
    for (ElementQ4 *element : upperElements) {
        // cout <<"Element No." << setw(10) << element->getID() << " " << "\n";
        for (int i = 0; i < element->getNID().size(); i++) {
            NodeValues[i][0] = element->getNID()[i]->getMassDensity();
            NodeValues[i][1] = element->getNID()[i]->getBodyForce()[0];
            NodeValues[i][2] = element->getNID()[i]->getBodyForce()[1];
        }
        
        for (int i = 0; i < IntPos.size(); i++) {
            for (int j = 0; j < IntPos.size(); j++) {
                element->evaluateF_x(res, IntPos[i], IntPos[j], NodeValues);
                cout << "Integration point: " << i * _spaceDim + j;
                for (int k = 0; k < res.size(); k++) {
                    cout << setw(12) << res[k];
                }
                cout << "\n";
            }
        }
        cout << "\n";
    }

    // Loop through upper elements
    for (ElementQ4 *element : lowerElements) {
        // cout <<"Element No." << setw(10) << element->getID() << " " << "\n";
        for (int i = 0; i < element->getNID().size(); i++) {
            NodeValues[i][0] = element->getNID()[i]->getMassDensity();
            NodeValues[i][1] = element->getNID()[i]->getBodyForce()[0];
            NodeValues[i][2] = element->getNID()[i]->getBodyForce()[1];
        }

        // Debug lines
        for (int i = 0; i < IntPos.size(); i++) {
            for (int j = 0; j < IntPos.size(); j++) {
                element->evaluateF_x(res, IntPos[i], IntPos[j], NodeValues);
                cout << "Integration point: " << i * _spaceDim + j;
                for (int k = 0; k < res.size(); k++) {
                    cout << setw(12) << res[k];
                }
                cout << "\n";
            }
        }
        cout << "\n";
    }

    NodeValues.resize(2);
    // Loop through cohesive elements
    for (ElementQ4Cohesive *element : cohesiveElements) {
        cout <<"Element No." << setw(10) << element->getID() << " " << "\n";
        for (int i = 0; i < element->getNID().size(); i++) {
            NodeValues[i][0] = element->getNID()[i]->getMassDensity();
            NodeValues[i][1] = element->getNID()[i]->getBodyForce()[0];
            NodeValues[i][2] = element->getNID()[i]->getBodyForce()[1];
        }

        for (int i = 0; i < IntPos.size(); i++) {
            element->evaluateF_x(res, IntPos[i], NodeValues);
            cout << "Integration point: " << i;
            for (int k = 0; k < res.size(); k++) {
                cout << setw(12) << res[k];
            }
            cout << "\n";
        }
        cout << "\n";
    }
};

/** Test integratorNfN */
void Problem::testIntegratorNfN() const {
    // Set mass density
    /** 
    int nOfDofs = upperElements[0]->getNID()[0]->getDOF().size(); 
    int nOfNodes = upperElements[0]->getNID().size();
    int eleMassCols = nOfDofs * nOfNodes;
    int &NVCols = nOfDofs;
    vector<double> massDensity(pow(NVCols, 2), 0.);
    for (int i = 0; i < NVCols; i++) {
        massDensity[i * NVCols + i] = 1.;
    }

    // Initialize some results
    vector<double> globalMassMatrix(_totalDOF * _totalDOF, 0.0);
    vector<double> eleMassMatrix;
    vector<vector<double>> nodalValues(nOfNodes, massDensity);

    // The global i, j indices
    int I, J;
    // Loop through upper Elements
    for (ElementQ4 *element : upperElements) {
        // Set nodal values to nodal density
       
        // Calculate the nodal values
        element->IntegratorNfN(eleMassMatrix, nodalValues);
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalMassMatrix[_totalDOF * I + J] 
                            += eleMassMatrix[(n1 * nOfDofs + k) * eleMassCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        }        
    }
    
    // Loop through lower Elements
    for (ElementQ4 *element : lowerElements) {
        // Calculate the nodal values
        element->IntegratorNfN(eleMassMatrix, nodalValues);        
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalMassMatrix[_totalDOF * I + J] 
                            += eleMassMatrix[(n1 * nOfDofs + k) * eleMassCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        } 
    }
    
    // Set mass density
    nOfDofs = cohesiveElements[0]->getNID()[0]->getDOF().size(); 
    nOfNodes = cohesiveElements[0]->getNID().size();
    eleMassCols = nOfDofs * nOfNodes;
    // int &NVCols = nOfDofs;
    massDensity.resize(pow(NVCols, 2), 0.);
    fill(massDensity.begin(), massDensity.end(), 0.);
    for (int i = 0; i < NVCols; i++) {
        massDensity[i * NVCols + i] = 1.;
    }
    // Loop through cohesive Elements
    eleMassMatrix.resize(pow(eleMassCols, 2));
    nodalValues.resize(nOfNodes);
    fill(nodalValues.begin(), nodalValues.end(), massDensity);
    // Loop through cohesive Elements
    for (ElementQ4Cohesive *element : cohesiveElements) {
        // Calculate the nodal values
        element->IntegratorNfN(eleMassMatrix, nodalValues);        
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalMassMatrix[_totalDOF * I + J] 
                            += eleMassMatrix[(n1 * nOfDofs + k) * eleMassCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        } 
    }

    // Printout the matrix
    ofstream myFile;
    myFile.open("Testlog.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== GlobalNfNMatrix ======================================" << "\n";
    printMatrix(myFile, globalMassMatrix, _totalDOF, _totalDOF);
    myFile.close();
    */
};

/** Test integratorBfB */
void Problem::testIntegratorBfB() const {
    /** 
    // Set Nodal stiffness matrix
    int spaceDim = upperElements[0]->getNID()[0]->getSpaceDim();
    int nOfDofs = upperElements[0]->getNID()[0]->getDOF().size(); 
    int nOfNodes = upperElements[0]->getNID().size();
    int eleStiffCols = nOfDofs * nOfNodes;
    int NVCols = nOfDofs * spaceDim;

    vector<double> pointStiff(pow(NVCols, 2), 0.);
    for (int i = 0; i < NVCols; i++) {
        pointStiff[i * NVCols + i] = 1.;
    }
    // Initialize some results
    vector<double> globalStiffMatrix(_totalDOF * _totalDOF, 0.0);
    vector<double> eleStiffMatrix(pow(eleStiffCols, 2), 0.);

    // Set nodal values to Identity
    vector<vector<double>> nodalValues(nOfNodes, pointStiff);

    // The global i, j indices
    int I, J;
    // Loop through upper Elements
    for (ElementQ4 *element : upperElements) {
        // Set nodal values to nodal density
       
        // Calculate the nodal values
        element->IntegratorBfB(eleStiffMatrix, nodalValues);
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalStiffMatrix[_totalDOF * I + J] 
                            += eleStiffMatrix[(n1 * nOfDofs + k) * eleStiffCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        }        
    }
    
    // Loop through lower Elements
    for (ElementQ4 *element : lowerElements) {
        // Set nodal values to nodal density
       
        // Calculate the nodal values
        element->IntegratorBfB(eleStiffMatrix, nodalValues);
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalStiffMatrix[_totalDOF * I + J] 
                            += eleStiffMatrix[(n1 * nOfDofs + k) * eleStiffCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        }        
    }
    
    // Loop through cohesive Elements
    // Set mass density
    nOfDofs = cohesiveElements[0]->getNID()[0]->getDOF().size(); 
    nOfNodes = cohesiveElements[0]->getNID().size();
    eleStiffCols = nOfDofs * nOfNodes;
    NVCols = nOfDofs * spaceDim;
    pointStiff.resize(pow(NVCols, 2), 0.);
    fill(pointStiff.begin(), pointStiff.end(), 0.);
    for (int i = 0; i < NVCols; i++) {
        pointStiff[i * NVCols + i] = 1.;
    }
    // Loop through cohesive Elements
    eleStiffMatrix.resize(pow(eleStiffCols, 2));
    nodalValues.resize(nOfNodes);
    fill(nodalValues.begin(), nodalValues.end(), pointStiff);
    // Loop through cohesive Elements
    for (ElementQ4Cohesive *element : cohesiveElements) {
        element->IntegratorBfB(eleStiffMatrix, nodalValues);
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalStiffMatrix[_totalDOF * I + J] 
                            += eleStiffMatrix[(n1 * nOfDofs + k) * eleStiffCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        }        
    }

    // Printout the matrix
    ofstream myFile;
    myFile.open("Testlog.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== GlobalBfBMatrix ======================================" << "\n";
    printMatrix(myFile, globalStiffMatrix, _totalDOF, _totalDOF); 
    myFile.close();
    */
};

/** Test integratorBfN */
void Problem::testIntegratorBfN() const {
    /** 
    // Set Nodal stiffness matrix
    int spaceDim = upperElements[0]->getNID()[0]->getSpaceDim();
    int nOfDofs = upperElements[0]->getNID()[0]->getDOF().size(); 
    int nOfNodes = upperElements[0]->getNID().size();
    int eleBfNCols = nOfDofs * nOfNodes;
    int fCols = nOfDofs;
    int fRows = nOfDofs * spaceDim;

    vector<double> pointF(fRows * fCols, 1.);
    // Initialize some results
    vector<double> globalBfNMatrix(pow(_totalDOF, 2), 0.);
    vector<double> eleBfNMatrix(pow(eleBfNCols, 2), 0.);

    // Set nodal values to Identity
    vector<vector<double>> nodalValues(nOfNodes, pointF);

    // The global i, j indices
    int I, J;
    // Loop through upper Elements
    for (ElementQ4 *element : upperElements) {
        // Set nodal values to nodal density
       
        // Calculate the nodal values
        element->IntegratorBfN(eleBfNMatrix, nodalValues);
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalBfNMatrix[_totalDOF * I + J] 
                            += eleBfNMatrix[(n1 * nOfDofs + k) * eleBfNCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        }        
    }
    
    // Loop through lower Elements
    for (ElementQ4 *element : lowerElements) {
        // Set nodal values to nodal density
        // Calculate the nodal values
        element->IntegratorBfN(eleBfNMatrix, nodalValues);
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalBfNMatrix[_totalDOF * I + J] 
                            += eleBfNMatrix[(n1 * nOfDofs + k) * eleBfNCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        }        
    }
    
    // Loop through cohesive Elements
    nOfDofs = cohesiveElements[0]->getNID()[0]->getDOF().size(); 
    nOfNodes = cohesiveElements[0]->getNID().size();
    eleBfNCols = nOfDofs * nOfNodes;
    fRows = nOfDofs * spaceDim;
    fCols = nOfDofs;
    pointF.resize(fCols * fRows, 0.);
    fill(pointF.begin(), pointF.end(), 1.);

    // Loop through cohesive Elements
    eleBfNMatrix.resize(pow(eleBfNCols, 2));
    nodalValues.resize(nOfNodes);
    fill(nodalValues.begin(), nodalValues.end(), pointF);
    // Loop through cohesive Elements
    for (ElementQ4Cohesive *element : cohesiveElements) {
        // Set nodal values to nodal density
        // Calculate the nodal values
        element->IntegratorBfN(eleBfNMatrix, nodalValues);
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalBfNMatrix[_totalDOF * I + J] 
                            += eleBfNMatrix[(n1 * nOfDofs + k) * eleBfNCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        }        
    }

    // Printout the matrix
    ofstream myFile;
    myFile.open("Testlog.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== GlobalBfNMatrix ======================================" << "\n";
    printMatrix(myFile, globalBfNMatrix, _totalDOF, _totalDOF); 
    myFile.close();
    */
};

/** Test integratorNfB */
void Problem::testIntegratorNfB() const {
    /** 
    // Set Nodal stiffness matrix
    int spaceDim = upperElements[0]->getNID()[0]->getSpaceDim();
    int nOfDofs = upperElements[0]->getNID()[0]->getDOF().size(); 
    int nOfNodes = upperElements[0]->getNID().size();
    int eleNfBCols = nOfDofs * nOfNodes;
    int fRows = nOfDofs;
    int fCols = nOfDofs * spaceDim;

    vector<double> pointF(fRows * fCols, 1.);
    // Initialize some results
    vector<double> globalNfBMatrix(pow(_totalDOF, 2), 0.);
    vector<double> eleNfBMatrix(pow(eleNfBCols, 2), 0.);

    // Set nodal values to Identity
    vector<vector<double>> nodalValues(nOfNodes, pointF);

    // The global i, j indices
    int I, J;
    // Loop through upper Elements
    for (ElementQ4 *element : upperElements) {
        // Set nodal values to nodal density
       
        // Calculate the nodal values
        element->IntegratorNfB(eleNfBMatrix, nodalValues);
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalNfBMatrix[_totalDOF * I + J] 
                            += eleNfBMatrix[(n1 * nOfDofs + k) * eleNfBCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        }        
    }
    
    // Loop through lower Elements
    for (ElementQ4 *element : lowerElements) {
        // Set nodal values to nodal density
        // Calculate the nodal values
        element->IntegratorNfB(eleNfBMatrix, nodalValues);
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalNfBMatrix[_totalDOF * I + J] 
                            += eleNfBMatrix[(n1 * nOfDofs + k) * eleNfBCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        }        
    }
    
    // Loop through cohesive Elements
    // Set mass density
    nOfDofs = cohesiveElements[0]->getNID()[0]->getDOF().size(); 
    nOfNodes = cohesiveElements[0]->getNID().size();
    eleNfBCols = nOfDofs * nOfNodes;
    fRows = nOfDofs;
    fCols = nOfDofs * spaceDim;
    pointF.resize(fCols * fRows, 0.);
    fill(pointF.begin(), pointF.end(), 1.);

    // Loop through cohesive Elements
    eleNfBMatrix.resize(pow(eleNfBCols, 2));
    nodalValues.resize(nOfNodes);
    fill(nodalValues.begin(), nodalValues.end(), pointF);
    // Loop through cohesive Elements
    for (ElementQ4Cohesive *element : cohesiveElements) {
        // Set nodal values to nodal density
        // Calculate the nodal values
        element->IntegratorNfB(eleNfBMatrix, nodalValues);
        for (int n1 = 0; n1 < nOfNodes; n1++) {
            for (int n2 = 0; n2 < nOfNodes; n2++) {
                for (int k = 0; k < nOfDofs; k++) {
                    for (int l = 0; l < nOfDofs; l++) {
                        I = element->getNID()[n1]->getDOF(k);
                        if (I == -1) continue;
                        J = element->getNID()[n2]->getDOF(l);
                        if (J == -1) continue;
                        globalNfBMatrix[_totalDOF * I + J] 
                            += eleNfBMatrix[(n1 * nOfDofs + k) * eleNfBCols + (n2 * nOfDofs + l)];
                    }
                }
            }
        }        
    }

    // Printout the matrix
    ofstream myFile;
    myFile.open("Testlog.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== GlobalNfBMatrix ======================================" << "\n";
    printMatrix(myFile, globalNfBMatrix, _totalDOF, _totalDOF); 
    myFile.close();
    */
};

/** Non-zeros per row in global matrices */
void Problem::getNNZPerRow(PetscInt *nnz) {
    // First clear nnz
    for (int i = 0; i < _totalDOF; i++) {
        nnz[i] = 0;
    }

    // Store non-zero entries
    vector<int> nonZeros(pow(_totalDOF, _spaceDim), 0);
    globalRows = new PetscInt [_totalDOF];

    // Some constants
    int nOfNodes;
    int nOfDofs;
    
    // Stores global and local index
    int I, J;
    
    // Upper elements non-zeros
    for (ElementQ4 *element : upperElements) {
        // First node
        for (int n1 = 0; n1 < nOfNodes; n1++)
        {
            nOfNodes = element->getNID().size();
            nOfDofs = element->getNID()[0]->getDOF().size();
            // cout << "\n";
            // Then row of local JF3
            for (int i = 0; i < nOfDofs; i++)
            {
                I = element->getNID()[n1]->getDOF()[i];
                if (I == -1)
                    continue;
                // Then node
                for (int n2 = 0; n2 < nOfNodes; n2++)
                {
                    // Then col of local JF3
                    for (int j = 0; j < nOfDofs; j++)
                    {
                        J = element->getNID()[n2]->getDOF()[j];
                        if (J == -1)
                            continue;  
                        nonZeros[I * _totalDOF + J] = 1;                      
                    }
                }
            }
        }
    }

    // Upper elements non-zeros
    for (ElementQ4 *element : lowerElements) {
        // First node
        for (int n1 = 0; n1 < nOfNodes; n1++)
        {
            nOfNodes = element->getNID().size();
            nOfDofs = element->getNID()[0]->getDOF().size();
            // cout << "\n";
            // Then row of local JF3
            for (int i = 0; i < nOfDofs; i++)
            {
                I = element->getNID()[n1]->getDOF()[i];
                if (I == -1)
                    continue;
                // Then node
                for (int n2 = 0; n2 < nOfNodes; n2++)
                {
                    // Then col of local JF3
                    for (int j = 0; j < nOfDofs; j++)
                    {
                        J = element->getNID()[n2]->getDOF()[j];
                        if (J == -1)
                            continue;  
                        nonZeros[I * _totalDOF + J] = 1;                      
                    }
                }
            }
        }
    }

    // Upper elements non-zeros
    for (ElementQ4Cohesive *element : cohesiveElements) {
        // First node
        for (int n1 = 0; n1 < nOfNodes; n1++)
        {
            nOfNodes = element->getNID().size();
            nOfDofs = element->getNID()[0]->getDOF().size();
            // cout << "\n";
            // Then row of local JF3
            for (int i = 0; i < nOfDofs; i++)
            {
                I = element->getNID()[n1]->getDOF()[i];
                if (I == -1)
                    continue;
                // Then node
                for (int n2 = 0; n2 < nOfNodes; n2++)
                {
                    // Then col of local JF3
                    for (int j = 0; j < nOfDofs; j++)
                    {
                        J = element->getNID()[n2]->getDOF()[j];
                        if (J == -1)
                            continue;  
                        nonZeros[I * _totalDOF + J] = 1;                      
                    }
                }
            }
        }
    }
    
    // Calculate all non-zeros
    for (int i = 0; i < _totalDOF; i++) {
        globalRows[i] = i;
        for (int j = 0; j < _totalDOF; j++) {
            nnz[i] += nonZeros[i * _totalDOF + j];
        }
    }

    _totalNonZeros = 0;
    for (int i = 0; i < _totalDOF; i++) {
        _totalNonZeros += nnz[i];
    }
    cout << "Total Non-zeros number in JF is: " << _totalNonZeros << "\n";
};

/** Printout a matrix */
void Problem::printMatrix(ofstream & myFile, const vector<double>& Matrix, int nRows, int nCols) const {
    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            myFile << setw(12) << Matrix[i * nCols + j] << " ";
        }
        myFile << "\n";
    }
    myFile << "\n";

    if (nRows == nCols) {
        myFile << "Matrix is Square!" << "\n";
        double err = 0.;
        double sum = 0.;
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                err += abs(Matrix[i * nCols + j] - Matrix[j * nCols + i]);
                sum += Matrix[i * nCols + j]; 
            }
        }
        myFile << "Error from Symmetric: " << setw(12) << err << "\n";
        myFile << "Sum of all entries: " << setw(12) << sum << "\n";
        myFile << "\n";
    }    
};

// ============= Test Elastic Solution ============================================================
/** Only has 1 block upperNodes and upperElements
 * Test Petsc, Mat, Vec, KSP solver, integratorBfB
 */
// Initialization of elastic problem
void Problem::initializeElastic(const vector<double> & xRanges, const vector<int> & edgeNums) {
    // Initialize geometry2D
    if (_spaceDim == 2) initializeGeometry2D(xRanges, edgeNums);
    
    // Initialize nodes
    initializeNodesElastic();
    
    // Assign Nodal DOFs
    assignNodalDOFElastic();
    
    // Test push global F
    testPushGlobalFElastic();
        
    // Initialize elements
    initializeElementsElastic();

    // Test push global JF
    testPushGlobalJFElastic();

    // Linear solver
    solveElastic();

    // Get back result
    testFetchGlobalSElastic();
};

// Initialization of Nodes
void Problem::initializeNodesElastic() {
    // Some input geometrical values
    int spaceDim = _spaceDim;

    // Number of elements and nodes
    // Parameters for bulk Nodes
    double lambda = 1.;
    double shearModulus = 1.;
    double biotAlpha = 1.;
    double biotMp = 1.;
    double fluidMobility = 1.;
    double fluidViscosity = 1.;
    // Element size
    vector<double> edgeSize (spaceDim);
    edgeSize[0] = myGeometry->xRange / myGeometry->xEdgeNum;
    edgeSize[1] = myGeometry->yRange / myGeometry->yEdgeNum;
    
    // DOF of fixed u_x, velocity and pressure
    vector<int> DOF_x_fixed (2 * spaceDim + 2, 1); 
    DOF_x_fixed[1] = 0;

    // DOF of fixed u_y, velocity and pressure
    vector<int> DOF_y_fixed (2 * spaceDim + 2, 1); 
    DOF_y_fixed[0] = 0;

    // DOF of fixed x, y, velocity and pressure
    vector<int> DOF_xy_fixed (2 * spaceDim + 2 ,1);

    // Default DOF
    vector<int> DOF_default (2 * spaceDim + 2, 1);
    DOF_default[0] = 0;
    DOF_default[1] = 0;

    // Mass Density and Body force.
    double massDensity = 1.0;
    vector<double> bodyForce  = {0.0, 0.0};
    
    // Upper subzone nodes
    upperNodes.resize(myGeometry->nOfNodes);
    int nodeID = 0;
    int nodeID_in_set = 0;
    vector<double> thisXYZ(spaceDim);
    vector<double> initialS(spaceDim * 2 + 2, 0.);
    // First y
    for (int j = 0; j < myGeometry->yNodeNum; j++) {
        // Then x
        for (int i = 0; i < myGeometry->xNodeNum; i++) {
            // Reset the coordinates
            thisXYZ[0] = i * edgeSize[0];
            thisXYZ[1] = j * edgeSize[1];
            massDensity = thisXYZ[0] + thisXYZ[1];
            // upper surface, put the nodal forces into initialS;
            if (j == myGeometry->yNodeNum - 1 ) {
                if (i == 0 || i == myGeometry->xNodeNum - 1) {
                    initialS[1] = -1.0;
                }
                else {
                    initialS[1] = -2.0;
                }
            }
            else {
                initialS[1] = 0.0;
            }

            // Initialize a node
            upperNodes[nodeID_in_set] = 
                new Node(nodeID, thisXYZ, DOF_default, spaceDim, 
                         massDensity, 
                         &bodyForce, 
                         lambda, 
                         shearModulus, 
                         biotAlpha, 
                         biotMp, 
                         fluidMobility, 
                         fluidViscosity);
            if (j == 0) // Lower surface
                upperNodes[nodeID_in_set]->setDOF(DOF_y_fixed);
            if (i == 0) // Left surface
                upperNodes[nodeID_in_set]->setDOF(DOF_x_fixed);
            if (i == 0 && j == 0) // Fixed point
                upperNodes[nodeID_in_set]->setDOF(DOF_xy_fixed);
            upperNodes[nodeID_in_set]->initializeS(initialS);
            nodeID += 1;
            nodeID_in_set += 1;
        }
    }
    
    _totalNofNodes = nodeID;
    ofstream myFile;
    myFile.open("Testlog_Elastic.txt");
    myFile << "=================== NodeInfoBefore ======================================" << "\n";
    for (Node* node : upperNodes) node->outputInfo(myFile, true);
    myFile.close();
};

// Assign global ID for each DOF
void Problem::assignNodalDOFElastic() {
    // Pointer to currentDOF
    _totalDOF = 0;
    
    // Assign upperzone
    for (int i = 0; i < upperNodes.size(); i++) {
        for (int j = 0; j < upperNodes[i]->getDOF().size(); j++) {
            if (upperNodes[i]->getDOF(j) == 0) {
                upperNodes[i]->setDOF(j, _totalDOF);
                _totalDOF += 1;
            }
            else {
                upperNodes[i]->setDOF(j, -1); 
            }
        }
    }

    // Output to log file
    ofstream myFile;
    myFile.open("Testlog_Elastic.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== NodeInfoAfter ======================================" << "\n";

    for (Node* node : upperNodes) node->outputInfo(myFile, true);
    myFile.close();
};

// TEST Try push to globalF, load only pushes the upper surface force
void Problem::testPushGlobalFElastic() {
    // Set size of global F
    VecCreate(PETSC_COMM_WORLD, &globalF);
    VecSetSizes(globalF, PETSC_DECIDE, _totalDOF);
    VecSetFromOptions(globalF);

    // Upper zone
    for (Node* node : upperNodes) {
        node->pushS(globalF);
    }
    
    // Output the global F
    /**
    cout << "TEST: Global F\n";
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF,PETSC_VIEWER_ASCII_MATLAB);
    VecView(globalF, PETSC_VIEWER_STDOUT_SELF);
    */
};

// TEST Try push to globalJF, load only pushes the upper surface force
void Problem::testPushGlobalJFElastic() {
    /**
    // Set size of global JF
    MatCreate(PETSC_COMM_WORLD, &globalJF);
    MatSetSizes(globalJF, PETSC_DECIDE, PETSC_DECIDE, _totalDOF, _totalDOF);
    MatSetFromOptions(globalJF);
    MatSetUp(globalJF);

    // Upper zone
    for (ElementQ4* element : upperElements) {
        element->JF(globalJF, 0);
    }

    // Assemble globalJF
    MatAssemblyBegin(globalJF, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalJF, MAT_FINAL_ASSEMBLY);

    // Output the global F
    cout << "TEST: Global JF\n";
    MatView(globalJF, PETSC_VIEWER_STDOUT_SELF);
    */
};

// Initialization of Elastic Elements
void Problem::initializeElementsElastic() {
    // Using the NID to assign values to each element
    vector<Node*> NID(4, NULL);
    
    // Subzone ElementQ4s
    upperElements.resize(myGeometry->nOfElements, NULL);

    for (int i = 0; i < myGeometry->xEdgeNum; i++) {
        for (int j = 0; j < myGeometry->yEdgeNum; j++) {
            // Setting NID for upper subzone
            NID = {upperNodes[j * myGeometry->xNodeNum + i], 
                   upperNodes[j * myGeometry->xNodeNum + i + 1], 
                   upperNodes[(j + 1) * myGeometry->xNodeNum + i + 1], 
                   upperNodes[(j + 1) * myGeometry->xNodeNum + i]};

            upperElements[j * myGeometry->xEdgeNum + i] = 
                new ElementQ4(j * myGeometry->xEdgeNum + i, NID);            
        }
    }
    ofstream myFile;
    myFile.open("Testlog_Elastic.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== ElementInfo ======================================" << "\n";
    for (ElementQ4* thisElement : upperElements) thisElement->outputInfo(myFile);
    for (ElementQ4* thisElement : lowerElements) thisElement->outputInfo(myFile);
    for (ElementQ4Cohesive* thisElement : cohesiveElements) thisElement->outputInfo(myFile);
    myFile.close();
};

/** Linearly solve Ks = F for elastic problems */
void Problem::solveElastic() {
    // Set size of global S
    VecCreate(PETSC_COMM_WORLD, &globalS);
    VecSetSizes(globalS, PETSC_DECIDE, _totalDOF);
    VecSetFromOptions(globalS);

    /** Initialize KSP */
    // Linear solver
    KSP ksp;

    // Preconditioner 
    PC pc;                 
    KSPCreate(PETSC_COMM_WORLD, &ksp);

    // Set operators
    KSPSetOperators(ksp, globalJF, globalJF);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);

    // Solve globalJF * globalS = globalF
    KSPSolve(ksp, globalF, globalS);

    // Destroy KSP
    KSPDestroy(&ksp);
};

/** Get back result from globalS */
void Problem::testFetchGlobalSElastic() {
    // Read results into each s
    for (Node* node : upperNodes) {
        node->fetchS_t(globalS);
    }

    // Output to log file
    ofstream myFile;
    myFile.open("Testlog_Elastic.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== NodeInfoAfterSolving ======================================" << "\n";

    for (Node* node : upperNodes) node->outputInfo(myFile, true);
    myFile.close();
};

// ============= Test PoroelasticElastic Solution ============================================================
/** Only has 1 block upperNodes and upperElements
 * Test Petsc, Mat, Vec, KSP solver, integratorBfB
 */
// Initialization of poroelastic problem
void Problem::initializePoroElastic(const vector<double> & xRanges, const vector<int> & edgeNums, int bulkKernel, int cohesiveKernel, double endingTime, double dt, string outputPrefix) {
    // Initialize geometry2D
    if (_spaceDim == 2) initializeGeometry2D(xRanges, edgeNums);

    // Initialize nodes
    initializeNodesPoroElastic();

    // Assign Nodal DOFs
    assignNodalDOFPoroElastic();

    // Initialize the timer
    timeConsumed.resize(4, 0.);
    clocks.resize(timeConsumed.size() + 1);

    // Initialize elements
    initializeElementsPoroElastic();

    // Initialize Mats and Vecs, Mats and TS
    initializePetsc();

    // Set outputprefix
    this->outputPrefix = outputPrefix;
    
    // Set the kernels
    this->_bulkKernel = bulkKernel;
    this->_cohesiveKernel = cohesiveKernel;

    // Non-Linear solver
    solvePoroElastic(endingTime, dt);

    // Cout timeConsumed
    for (int i = 0; i < timeConsumed.size(); i++) {
        cout << "Time Consumed [" << i << "] is " << timeConsumed[i] << "\n";
    }
    
};

// Initialization of Nodes
void Problem::initializeNodesPoroElastic() {
    // Some input geometrical values
    int spaceDim = _spaceDim;

    // Number of elements and nodes
    // Parameters for bulk Nodes
    // Mass Density and Body force.
    double massDensity = 1.0;
    vector<double> bodyForce  = {0.0, 0.0};
    double lambda = 1.;
    double shearModulus = 1.;
    double biotAlpha = 1.;
    double biotMp = 1.;
    double fluidMobility = 1.;
    double fluidViscosity = 1.;
    double fluidDensity = 1.;
    double porosity = 0.;
    vector<double> fluidBodyForce = {0.0, 0.0};
    double source = 0.;

    // For cohesive nodal properties
    double fluidMobilityX = 0.000005;
    double fluidMobilityZ = 0.000020;
    double faultPorosity = 0.1;
    double thickness = 0.001;
    double betaP = 1.0;
    double betaSigma = 1.0;
    double rateStateA = 0.008;
    double rateStateB = 0.012;
    double DRateState = 1.0;
    double faultSource = 0.0;
    double fReference = 0.6;
    double VReference = 1.0;
    vector<double> initialSlip = {0.0, 0.0};
    // double activesource = 1.0;

    // Element size
    vector<double> edgeSize (spaceDim);
    edgeSize[0] = myGeometry->xRange / myGeometry->xEdgeNum;
    edgeSize[1] = myGeometry->yRange / myGeometry->yEdgeNum;
    
    

    /** First test a case with fixed displacement
     * p, u, trace_strain to vary and see if the diffusion makes sense
     */  
    // Default DOF
    vector<int> DOF_default (2 * spaceDim + 2, 0); 
    DOF_default[2] = 1; DOF_default[3] = 1;  // Fix V

    vector<int> DOFCohesive_default (17, 2);
    for (int i = 12; i < 17; i++) DOFCohesive_default[i] = 0; 
    // DOFCohesive_default[15] = 1;             // Fix psi for now
    
    // DEBUG LINE
    DOF_default[4] = 1; // Fix p
    // DOFCohesive_default[8] = 1; 
    // DOFCohesive_default[9] = 1; 
    // DOFCohesive_default[13] = 1;
    DOFCohesive_default[14] = 1; 
    // DOFCohesive_default[15] = 1; // Fix psi
    // DOFCohesive_default[16] = 1; // Fix V

    // All locked DOF
    vector<int> DOF_allLocked (2 * spaceDim + 2, 1);
    vector<int> DOFCohesive_allLocked (17, 2);
    for (int i = 12; i < 17; i++) DOFCohesive_allLocked[i] = 1; 

    // Upper subzone nodes
    upperNodes.resize(myGeometry->nOfNodes);
    int nodeID = 0;
    int nodeID_in_set = 0;
    vector<double> thisXYZ(spaceDim);
    vector<double> initialS(spaceDim * 2 + 2, 0.);
    // First y
    for (int j = 0; j < myGeometry->yNodeNum; j++) {
        // Then x
        for (int i = 0; i < myGeometry->xNodeNum; i++) {
            // Reset the coordinates
            thisXYZ[0] = i * edgeSize[0];
            thisXYZ[1] = j * edgeSize[1];
            // massDensity = thisXYZ[0] + thisXYZ[1];
            
            // upper surface, put pressure = 1 into initialS;
            if (j == myGeometry->yNodeNum - 1 ) {
                initialS[2 * spaceDim] = 0.0;
            }
            // lower surface, put pressure = 0 into initialS;
            else if (j == 0) {
                initialS[2 * spaceDim] = -0.0;
            }
            else {
                initialS[2 * spaceDim] = 0.0;
            }

            // Initialize a node
            upperNodes[nodeID_in_set] = 
                new Node(nodeID, thisXYZ, DOF_default, spaceDim, 
                         massDensity, 
                         &bodyForce, 
                         lambda, 
                         shearModulus, 
                         biotAlpha, 
                         biotMp, 
                         fluidMobility, 
                         fluidViscosity, 
                         fluidDensity, 
                         porosity,
                         &fluidBodyForce, 
                         source);
            
            // Fix pressure on the lower boundary 
            if (j == 0) {
                // Fix p
                upperNodes[nodeID_in_set]->setDOF(2 * spaceDim, 1);
                // upperNodes[nodeID_in_set]->setDOF(1, 1);
                if (i == 0) upperNodes[nodeID_in_set]->setDOF(0, 1);
                // upperNodes[nodeID_in_set]->setDOF(0, 1);
            }

            // Give upper surface a traction
            if (j == myGeometry->yNodeNum - 1) {
                // Apply pure shear traction (1.0);
                vector<double> traction = {0.0, -0.02};
                upperNodes[nodeID_in_set]->setTraction(&traction);
            }

            // Give source in the middle
            if (j == myGeometry->yNodeNum / 2  && i == myGeometry->xNodeNum / 2) {
                // Fix p, ux, uy
                upperNodes[nodeID_in_set]->setSource(0.0);
            }
            
            upperNodes[nodeID_in_set]->initializeS(initialS);
            nodeID += 1;
            nodeID_in_set += 1;
        }
    }

    // Lower subzone nodes
    lowerNodes.resize(myGeometry->nOfNodes);
    nodeID_in_set = 0;
    // First y
    for (int j = 0; j < myGeometry->yNodeNum; j++) {
        // Then x
        for (int i = 0; i < myGeometry->xNodeNum; i++) {
            // Reset the coordinates
            thisXYZ[0] = i * edgeSize[0];
            thisXYZ[1] = -j * edgeSize[1];
            // massDensity = thisXYZ[0] + thisXYZ[1];
            
            // lower surface, put pressure = 1 into initialS;
            if (j == myGeometry->yNodeNum - 1) {
                initialS[2 * spaceDim] = 0.0;
            }
            // lower surface, put pressure = 0 into initialS;
            else if (j == 0) {
                initialS[2 * spaceDim] = 0.0;
            }
            else {
                initialS[2 * spaceDim] = 0.0;
            }

            // Initialize a node
            lowerNodes[nodeID_in_set] = 
                new Node(nodeID, thisXYZ, DOF_default, spaceDim, 
                         massDensity, 
                         &bodyForce, 
                         lambda, 
                         shearModulus, 
                         biotAlpha, 
                         biotMp, 
                         fluidMobility, 
                         fluidViscosity, 
                         fluidDensity, 
                         porosity,
                         &fluidBodyForce, 
                         source);
            
            // Fix some boundaries on the bottom
            if (j == 0) {
                if (i == 0) lowerNodes[nodeID_in_set]->setDOF(0, 1);
            }
            if (j == myGeometry->yNodeNum - 1) {
                // Fix p, ux, uy
                lowerNodes[nodeID_in_set]->setDOF(2 * spaceDim, 1);
                if (i == 0) lowerNodes[nodeID_in_set]->setDOF(0, 1);
                lowerNodes[nodeID_in_set]->setDOF(1, 1);
            }

            lowerNodes[nodeID_in_set]->initializeS(initialS);
            nodeID += 1;
            nodeID_in_set += 1;
        }
    }

    // Cohesive zone nodes
    cohesiveNodes.resize(myGeometry->xNodeNum);
    nodeID_in_set = 0;
    initialS.resize(spaceDim + 3);
    fill(initialS.begin(), initialS.end(), 0.0);
    initialS[4] = 0;

    // Initialize \psi = f_* + b log (V_* \theta / D_RS)
    // double theta = 1.0e6;

    // initialS[3] is \psi
    // initialS[3] = fReference - rateStateB * log(VReference * rateStateB / DRateState);
    initialS[3] = fReference;
    // initialS[4] = VReference;
    cout << "InitialS[3] is:" << initialS[3] << "\n";
    cout << "InitialS[4] is:" << initialS[4] << "\n";
    vector<Node*> lowerUpperNodes (2, NULL);

    for (int i = 0; i < myGeometry->xNodeNum; i++) {
        // Reset the coordinates
        thisXYZ[0] = i * edgeSize[0];
        thisXYZ[1] = 0.;

        // Get lower and upper nodes
        lowerUpperNodes[0] = lowerNodes[i];
        lowerUpperNodes[1] = upperNodes[i];
        // Initialize a node
        cohesiveNodes[nodeID_in_set] =
            new CohesiveNode(nodeID,
                             thisXYZ,
                             DOFCohesive_default,
                             lowerUpperNodes,
                             spaceDim,
                             massDensity,
                             &bodyForce,
                             fluidMobilityX,
                             fluidMobilityZ,
                             fluidViscosity, 
                             faultPorosity,
                             thickness,
                             betaP,
                             betaSigma,
                             rateStateA,
                             rateStateB,
                             DRateState,
                             &fluidBodyForce,
                             faultSource, 
                             &initialSlip, 
                             VReference, 
                             fReference);
    
        cohesiveNodes[nodeID_in_set]->initializeS(initialS);
        
        nodeID += 1;
        nodeID_in_set += 1;
    }
    
    _totalNofNodes = nodeID;
    ofstream myFile;
    myFile.open("Testlog_PoroElastic1.txt");
    myFile << "=================== NodeInfoBefore ======================================" << "\n";
    for (Node* node : upperNodes) node->outputInfo(myFile, true);
    for (Node* node : lowerNodes) node->outputInfo(myFile, true);
    for (CohesiveNode* node : cohesiveNodes) node->outputInfo(myFile, true);
    myFile.close();
};

// Assign global ID for each DOF
void Problem::assignNodalDOFPoroElastic() {
    // Pointer to currentDOF
    _totalDOF = 0;
    
    // Assign upperzone
    for (int i = 0; i < upperNodes.size(); i++) {
        for (int j = 0; j < upperNodes[i]->getDOF().size(); j++) {
            if (upperNodes[i]->getDOF(j) == 0) {
                upperNodes[i]->setDOF(j, _totalDOF);
                _totalDOF += 1;
            }
            else {
                upperNodes[i]->setDOF(j, -1); 
            }
        }
    }

    // Assign lowerzone
    for (int i = 0; i < lowerNodes.size(); i++) {
        for (int j = 0; j < lowerNodes[i]->getDOF().size(); j++) {
            if (lowerNodes[i]->getDOF(j) == 0) {
                lowerNodes[i]->setDOF(j, _totalDOF);
                _totalDOF += 1;
            }
            else {
                lowerNodes[i]->setDOF(j, -1); 
            }
        }
    }

    // Assign cohesivezone
    for (int i = 0; i < cohesiveNodes.size(); i++) {
        for (int j = 0; j < cohesiveNodes[i]->getDOF().size(); j++) {
            if (cohesiveNodes[i]->getDOF(j) == 0) {
                cohesiveNodes[i]->setDOF(j, _totalDOF);
                _totalDOF += 1;
            }
            // Locked degree of freedom
            else if (cohesiveNodes[i]->getDOF(j) == 1){
                cohesiveNodes[i]->setDOF(j, -1); 
            }
        }
        cohesiveNodes[i]->updateDOF();
    }

    // Output to log file
    ofstream myFile;
    myFile.open("Testlog_PoroElastic1.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== NodeInfoAfter ======================================" << "\n";

    for (Node* node : upperNodes) node->outputInfo(myFile, true);
    for (Node* node : lowerNodes) node->outputInfo(myFile, true);
    for (CohesiveNode* node : cohesiveNodes) node->outputInfo(myFile, true);
    myFile.close();
};

// Initialize Vec F, S, S_t, Mat JF
void Problem::initializePetsc() {
    // Set size of global F
    VecCreate(PETSC_COMM_WORLD, &globalF);
    VecSetSizes(globalF, PETSC_DECIDE, _totalDOF);
    VecSetFromOptions(globalF);
    VecSetOption(globalF, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

    // Initialize globalS
    VecCreate(PETSC_COMM_WORLD, &globalS);
    VecSetSizes(globalS, PETSC_DECIDE, _totalDOF);
    VecSetFromOptions(globalS);
    VecSetOption(globalS, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

    // Initialize globalS_t
    VecCreate(PETSC_COMM_WORLD, &globalS_t);
    VecSetSizes(globalS_t, PETSC_DECIDE, _totalDOF);
    VecSetFromOptions(globalS_t);
    VecSetOption(globalS_t, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    
    // Initialize globalJF
    PetscInt *globalJF_nnz = new PetscInt [_totalDOF];
    getNNZPerRow(globalJF_nnz);

    MatCreateSeqAIJ(PETSC_COMM_WORLD, _totalDOF, _totalDOF, 0, globalJF_nnz, &globalJF);
    // MatSetSizes(globalJF, PETSC_DECIDE, PETSC_DECIDE, _totalDOF, _totalDOF);
    delete [] globalJF_nnz;

    MatSetFromOptions(globalJF);
    MatSetOption(globalJF, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    MatSetUp(globalJF);

    // Zero the system Vecs and Mats
    VecZeroEntries(globalS);
    VecZeroEntries(globalS_t);
    VecZeroEntries(globalF);

    // Set initial conditions by pushing to globalS
    for (Node* node : upperNodes) {
        node->pushS(globalS);
    }

    for (Node* node : lowerNodes) {
        node->pushS(globalS);
    }

    for (CohesiveNode* node : cohesiveNodes) {
        node->pushS(globalS);
    }
}

// Initialization of Elastic Elements
void Problem::initializeElementsPoroElastic() {
    // Using the NID to assign values to each element
    vector<Node*> NID(4, NULL);
    vector<CohesiveNode*> NID_cohesive(2, NULL);
    int GlobalEleID = 0;
    vector<int> upperFace = {2, 3};
    vector<vector<int>> loadFaces(1, upperFace);
    vector<vector<int>> noLoadFaces;
    // Upperzone ElementQ4s
    upperElements.resize(myGeometry->nOfElements, NULL);
    for (int i = 0; i < myGeometry->xEdgeNum; i++) {
        for (int j = 0; j < myGeometry->yEdgeNum; j++) {
            // Setting NID for upper subzone
            NID = {upperNodes[j * myGeometry->xNodeNum + i], 
                   upperNodes[j * myGeometry->xNodeNum + i + 1], 
                   upperNodes[(j + 1) * myGeometry->xNodeNum + i + 1], 
                   upperNodes[(j + 1) * myGeometry->xNodeNum + i]};

            // Apply load only to the upper surface
            if (j == myGeometry->yEdgeNum - 1) {
                upperElements[j * myGeometry->xEdgeNum + i] = 
                    new ElementQ4(GlobalEleID, NID, loadFaces, &clocks, &timeConsumed); 
            }
            else {
                upperElements[j * myGeometry->xEdgeNum + i] = 
                    new ElementQ4(GlobalEleID, NID, noLoadFaces, &clocks, &timeConsumed); 
            }
            
            // Augmented global ID
            GlobalEleID += 1;           
        }
    }

    // No load on lower surface
    loadFaces.clear();

    // Lowerzone ElementQ4s
    lowerElements.resize(myGeometry->nOfElements, NULL);
    for (int i = 0; i < myGeometry->xEdgeNum; i++) {
        for (int j = 0; j < myGeometry->yEdgeNum; j++) {
            // Setting NID for upper subzone
            NID = {lowerNodes[(j + 1) * myGeometry->xNodeNum + i], 
                   lowerNodes[(j + 1) * myGeometry->xNodeNum + i + 1], 
                   lowerNodes[j * myGeometry->xNodeNum + i + 1], 
                   lowerNodes[j * myGeometry->xNodeNum + i]};

            lowerElements[j * myGeometry->xEdgeNum + i] = 
                new ElementQ4(GlobalEleID, NID, noLoadFaces, &clocks, &timeConsumed); 
            // Augment global node ID
            GlobalEleID += 1;           
        }
    }

    // Cohesivezone ElementQ4Cohesives
    cohesiveElements.resize(myGeometry->xEdgeNum, NULL);
    for (int i = 0; i < myGeometry->xEdgeNum; i++) {
        // Setting NID for upper subzone
        NID_cohesive = {cohesiveNodes[i], cohesiveNodes[i + 1]};

        cohesiveElements[i] = 
            new ElementQ4Cohesive(GlobalEleID, NID_cohesive, &clocks, &timeConsumed); 
        // Augment global node ID
        GlobalEleID += 1;           
    }

    ofstream myFile;
    myFile.open("Testlog_PoroElastic1.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== ElementInfo ======================================" << "\n";
    for (ElementQ4* thisElement : upperElements) thisElement->outputInfo(myFile);
    for (ElementQ4* thisElement : lowerElements) thisElement->outputInfo(myFile);
    for (ElementQ4Cohesive* thisElement : cohesiveElements) thisElement->outputInfo(myFile);
    myFile.close();

    // Allocate localF and localJF
    localFSize = upperElements[0]->getElementDOF();
    localJFSize = pow(localFSize, 2);
    localF = new double [localFSize];
    localJF = new double [localJFSize];

    // Allocate localFCohesive and localJFCohesive
    localFCohesiveSize = cohesiveElements[0]->getElementDOF();

    localJFCohesiveSize = pow(localFCohesiveSize, 2);
    localFCohesive = new double [localFCohesiveSize];
    localJFCohesive = new double [localJFCohesiveSize];
};

/** TS (SNES) solver for linear poroelastic problems */
void Problem::solvePoroElastic(double endingTime, double dt) {
    /** Initialize TS */
    // Non-linear time-dependent solver
    TS ts;

    /** Create timestepping solver context */
    TSCreate(PETSC_COMM_WORLD, & ts);    
    TSSetProblemType(ts, TS_NONLINEAR);

    /** Set Implicit Function and Implicit Jacobian */
    TSSetIFunction(ts, globalF, IFunction, this);
    TSSetIJacobian(ts, globalJF, globalJF, IJacobian, this);

    /** Set final time */
    TSSetMaxTime(ts, endingTime);
    
    TSSetTimeStep(ts, dt);
    TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER);

    /** Set the solution */
    TSSetSolution(ts, globalS);

    /** Set TS parameters from user options */
    TSSetFromOptions(ts);
    
    /** Set tolerances */
    PetscReal rtol = 1e-20;
    PetscReal atol = 1e-30;
    TSSetTolerances(ts, atol, NULL, rtol, NULL);
    SNES snes;
    TSGetSNES(ts, &snes);
    KSP ksp;
    SNESGetKSP(snes, &ksp);
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);

    SNESSetTolerances(snes, atol, rtol, rtol, 1000, 1000);     // Default values for the ints
    SNESSetLagJacobian(snes, 1);
    SNESSetLagPreconditioner(snes, 1);
    SNESLineSearch sneslinesearch;
    SNESGetLineSearch(snes, &sneslinesearch);
    SNESLineSearchSetType(sneslinesearch, SNESLINESEARCHBASIC);
    SNESLineSearchSetTolerances(sneslinesearch, PETSC_DEFAULT, 1000, rtol, atol, rtol, 1000); 
    
    KSPSetTolerances(ksp, rtol, atol, 1.1, 10000);   // 

    /** Solve TS */
    TSSetMaxSNESFailures(ts, -1);
    TSSolve(ts, globalS);
    
    TSView(ts, PETSC_VIEWER_STDOUT_SELF);
    SNESView(snes, PETSC_VIEWER_STDOUT_SELF);
    SNESLineSearchView(sneslinesearch, PETSC_VIEWER_STDOUT_SELF);
};

/** Calculating TSIFunction from the current ts and timestep t,
 * given global solution vector s
 * and global solution vector time derivative s_t,
 * stores the result in Vec F
 */
PetscErrorCode Problem::IFunction(TS ts, PetscReal t, Vec s, Vec s_t, Vec F, void *ctx) {
    PetscErrorCode ierr;
    // Convert the pointer
    Problem *myProblem = (Problem*) ctx;
    
    // Clear the vector
    ierr = VecZeroEntries(F);

    // Write vtk file
    if (t > myProblem->nodeTime) {
        // Debug lines
        cout << "IFunction t = " << t << "\n";
        ierr = TSGetStepNumber(ts, &(myProblem->stepNumber));
        myProblem->writeVTU();
        myProblem->nodeTime = t;
        // myProblem->prescribedSlip(t);
    }

    // Fetch all s into the nodes
    for (Node* node : myProblem->upperNodes) {
        node->fetchS(s);
        node->fetchS_t(s_t);
    }

    for (Node* node : myProblem->lowerNodes) {
        node->fetchS(s);
        node->fetchS_t(s_t);
    }
  
    for (CohesiveNode* node : myProblem->cohesiveNodes) {
        node->fetchS(s);
        node->fetchS_t(s_t);
    }    

    // Loop through all elements in upper and lower Elements
    for (ElementQ4 *element : myProblem->upperElements) {
        // Calculate F within the element
        element->elementF(F, myProblem->localF, myProblem->localFSize, myProblem->_bulkKernel, 0.);
    }

    for (ElementQ4 *element : myProblem->lowerElements) {
        // Calculate F within the element
        element->elementF(F, myProblem->localF, myProblem->localFSize, myProblem->_bulkKernel, 0.);
    }

    for (ElementQ4Cohesive *element : myProblem->cohesiveElements) {
        // Calculate F within the element
        element->elementF(F, myProblem->localFCohesive, myProblem->localFCohesiveSize, myProblem->_cohesiveKernel, 0., t);
    }

    // Assemble the global residual function
    ierr = VecAssemblyBegin(F);
    ierr = VecAssemblyEnd(F);
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF,PETSC_VIEWER_ASCII_MATLAB);

    // DEBUG LINES
    
    cout << "Test: s \n";
    VecView(s, PETSC_VIEWER_STDOUT_SELF);
    cout << "\n";
    
    /**
    cout << "s_t: \n";
    VecView(s_t, PETSC_VIEWER_STDOUT_SELF);
    cout << "\n";
    */
    // Output the global F
    cout << "TEST: Global F\n";
    
    // VecView(F, PETSC_VIEWER_STDOUT_SELF);
    double Fnorm;
    VecNorm(F, NORM_2, &Fnorm);
    cout << "Norm of F is: " << Fnorm << "\n"; 
    
    return ierr;
}

/** Calculating TSIJacobian from the current ts and timestep t,
 * given global solution vector s
 * and global solution vector time derivative s_t,
 * stores the result in Mat Amat and Mat Pmat
 */
PetscErrorCode Problem::IJacobian(TS ts, PetscReal t, Vec s, Vec s_t, PetscReal s_tshift, Mat Amat, Mat Pmat, void *ctx) {
    // DEBUG LINES
    cout << "IJacobian t = " << t << "\n";
    // cout << "IJacobian s_tshift = " << s_tshift << "\n";

    // Convert the pointer
    Problem *myProblem = (Problem*) ctx;
    PetscErrorCode ierr;
   
    // Clear the matrix
    PetscBool isAssembled;
    
    // Check if assembled
    ierr = MatAssembled(Pmat, &isAssembled);

    /** No need to fetch here, 
     * since IFunction is always evaluated first 
    // Fetch all s into the nodes
    for (Node* node : myProblem->upperNodes) {
        node->fetchS(s);
        node->fetchS_t(s_t);
    }
    */

    myProblem->nodeTime = t;

    // Loop through all elements in upperElements
    for (ElementQ4 *element : myProblem->upperElements) {
        // Calculate F within the element
        element->JF(Pmat, myProblem->localJF, myProblem->localJFSize, myProblem->_bulkKernel, s_tshift);
    }

    // Loop through all elements in lowerElements
    for (ElementQ4 *element : myProblem->lowerElements) {
        // Calculate F within the element
        element->JF(Pmat, myProblem->localJF, myProblem->localJFSize, myProblem->_bulkKernel, s_tshift);
    }

    // Loop through all elements in cohesive Elements
    for (ElementQ4Cohesive *element : myProblem->cohesiveElements) {
        // Calculate F within the element
        // DEBUG LINES
        cout << "Cohesive Kernel is: " << myProblem->_cohesiveKernel << "\n";
        element->JF(Pmat, myProblem->localJFCohesive, myProblem->localJFCohesiveSize, myProblem->_cohesiveKernel, s_tshift, t);
    }

    // myProblem->clocks[1] = clock();

    // Assemble the global Jacobian matrix
    ierr = MatAssemblyBegin(Pmat, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(Pmat, MAT_FINAL_ASSEMBLY);
    

    if (Amat != Pmat) {
        ierr = MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);
    }

    // myProblem->clocks[2] = clock();
    
    // DEBUG LINES    
    cout << "View Pmat Norm: ";
    PetscReal norm;
    MatNorm(Pmat, NORM_FROBENIUS, &norm);
    cout << norm << "\n";
    
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);

    // Output the global JF
    // cout << "Pmat: \n";
    // MatView(Pmat, PETSC_VIEWER_STDOUT_SELF);

    // myProblem->timeConsumed[0] += (double) (myProblem->clocks[1] - myProblem->clocks[0]) / CLOCKS_PER_SEC;
    // myProblem->timeConsumed[1] += (double) (myProblem->clocks[2] - myProblem->clocks[1]) / CLOCKS_PER_SEC;
    return ierr;
};


/** Prescribe slip for every node */
void Problem::prescribedSlip(double t) {
    // Currently apply constant slip everywhere and see what happens.
    for (CohesiveNode* node : cohesiveNodes) {
        slipFunction(node->getXYZ(), t);
        node->setSlip(&slip);
    }
};

/** Prescribed slip function. */
void Problem::slipFunction(const vector<double> & XYZ, double t) {
    // Prescribe slip at a given point, now hyperbolic
    // slip[0] = - (XYZ[0] - myGeometry->xRange) * XYZ[0] * 0.e-2 * t;
    slip[0] = 0.;
    slip[1] = 0.;
};


// ==================== Writing results into vtk/vtu files ============================
/** Write VTK files
 * Write into vtk nodal and cell connection values
 * At each timestep write a different file
 */
void Problem::writeVTK() const {
    string path = "./output/" + outputPrefix + to_string(stepNumber) + ".vtk";
    ofstream myFile(path);

    // Head lines
    myFile << "# vtk DataFile Version 2.0\n";        // version info
    myFile << "Time step " << stepNumber << "\n";    // title
    myFile << "ASCII"<< "\n";                        // File Format
    myFile << "DATASET UNSTRUCTURED_GRID" << "\n";   // structure type

    // Output node XYZs
    myFile << "POINTS " << _totalNofNodes << " double" << "\n";
    myFile << scientific;
    for (Node* node : upperNodes) {
        myFile << node->getXYZ()[0] << " " << node->getXYZ()[1] << " " << "0.0" << "\n";
    }
    
    
    // Output Elements
    myFile << "CELLS " << upperElements.size() << " " << 5 * upperElements.size() << "\n";
    for (ElementQ4 *element : upperElements) {        
        myFile << "4 " << element->getNID()[0]->getID() << " " 
                       << element->getNID()[1]->getID() << " " 
                       << element->getNID()[2]->getID() << " " 
                       << element->getNID()[3]->getID() << "\n";        
    }

    int I_u = 0;
    int I_v = I_u + _spaceDim;
    int I_p = I_v + _spaceDim;
    int I_e = I_p + 1;

    // Output cell type
    myFile << "CELL_TYPES " << upperElements.size() << "\n" ;
    for (int i = 0; i < upperElements.size(); i++) myFile << "9\n";


    // Output displacement
    myFile << "POINT_DATA " << _totalNofNodes << "\n";
    myFile << "VECTORS displacement double\n";
    for (Node *node : upperNodes) {
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s[I_u + i] << " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }

    // Output velocity
    myFile << "VECTORS velocity double\n";
    for (Node *node : upperNodes) {
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s[I_v + i] << " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }

    // Output pore fluid pressure
    myFile << "SCALARS pressure double 1" << "\n";
    myFile << "LOOKUP_TABLE default" << "\n";
    for (Node *node : upperNodes) {
        myFile << node->s[I_p] << "\n";
    }

    // Output volumetric strain
    myFile << "SCALARS trace_strain double 1" << "\n";
    myFile << "LOOKUP_TABLE default" << "\n";
    for (Node *node : upperNodes) {
        myFile << node->s[I_e] << "\n";
    }
};

/** Write VTU files
 * Write into vtu nodal and cell connection values
 * At each timestep write a different file
 */
void Problem::writeVTU() const {
    writeVTU_bulk();
    writeVTU_fault();
};

/** Write VTU files for the bulk
 */
void Problem::writeVTU_bulk() const {
    string path = "./output/" + outputPrefix + "Bulk" + to_string(stepNumber) + ".vtu";
    ofstream myFile(path);

    // Head lines
    myFile << "<?xml version=\"1.0\"?>" << "\n";        // version info
    myFile << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"BigEndian\">" << "\n";
    myFile << "  <UnstructuredGrid>" << "\n";    // title
    myFile << "    <Piece NumberOfPoints=\"" << upperNodes.size() + lowerNodes.size() << "\" NumberOfCells=\"" << upperElements.size() + lowerElements.size() << "\">" << "\n";                        // File Format
    
    // Output node XYZs
    myFile << "      <Points>" << "\n";
    myFile << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">" << "\n";
    myFile << scientific;
    for (Node* node : upperNodes) {
        myFile << "        " << node->getXYZ()[0] << " " << node->getXYZ()[1] << " " << "0.0" << "\n";
    }
    for (Node* node : lowerNodes) {
        myFile << "        " << node->getXYZ()[0] << " " << node->getXYZ()[1] << " " << "0.0" << "\n";
    }
    myFile << "        </DataArray>" << "\n";
    myFile << "      </Points>" << "\n";    
    
    // Output Elements
    myFile << "      <Cells>" << "\n";
    myFile << "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">" << "\n";
    for (ElementQ4 *element : upperElements) {
        myFile << "        "
               << element->getNID()[0]->getID() << " "
               << element->getNID()[1]->getID() << " "
               << element->getNID()[2]->getID() << " "
               << element->getNID()[3]->getID() << "\n";
    }
    for (ElementQ4 *element : lowerElements) {
        myFile << "        "
               << element->getNID()[0]->getID() << " "
               << element->getNID()[1]->getID() << " "
               << element->getNID()[2]->getID() << " "
               << element->getNID()[3]->getID() << "\n";
    }
    myFile << "        </DataArray>" << "\n";

    // Output offsets
    myFile << "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">" << "\n";
    myFile << "        ";
    for (int i = 0; i < upperElements.size() + lowerElements.size(); i++) myFile << 4 * (i + 1) << " ";
    myFile << "\n" << "        </DataArray>" << "\n";

    // Output cell type
    myFile << "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">" << "\n";
    myFile << "        ";
    for (int i = 0; i < upperElements.size() + lowerElements.size(); i++) myFile << "9 ";
    myFile << "\n" << "        </DataArray>" << "\n";
    myFile << "      </Cells>" << "\n";

    // ============================== Output Node results =================================
    int I_u = 0;
    int I_v = I_u + _spaceDim;
    int I_p = I_v + _spaceDim;
    int I_e = I_p + 1;
    myFile << "      <PointData>" << "\n";
    
    // Global node ID
    myFile << "        <DataArray type=\"Int32\" Name=\"GlobalNodeId\" format=\"ascii\">" << "\n";
    myFile << "        ";
    for (Node *node : upperNodes) {
        myFile << node->getID() << " ";
    }
    for (Node *node : lowerNodes) {
        myFile << node->getID() << " ";
    }
    myFile << "\n" << "        </DataArray>" << "\n";

    // Nodal displacement
    myFile << "        <DataArray type=\"Float64\" Name=\"Displacements\" NumberOfComponents=\"3\" ComponentName0=\"Ux\" ComponentName1=\"Uy\" ComponentName2=\"Uz\" Format=\"ascii\">" << "\n";
    for (Node *node : upperNodes) {
        myFile << "        ";
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s[I_u + i] << " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }
    for (Node *node : lowerNodes) {
        myFile << "        ";
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s[I_u + i] << " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }
    myFile << "        </DataArray>" << "\n";

    // Nodal displacement_t
    myFile << "        <DataArray type=\"Float64\" Name=\"Displacements_t\" NumberOfComponents=\"3\" ComponentName0=\"U_tx\" ComponentName1=\"U_ty\" ComponentName2=\"U_tz\" Format=\"ascii\">" << "\n";
    for (Node *node : upperNodes) {
        myFile << "        ";
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s_t[I_u + i] << " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }
    for (Node *node : lowerNodes) {
        myFile << "        ";
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s_t[I_u + i] << " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }
    myFile << "        </DataArray>" << "\n";

    // Output velocity
    myFile << "        <DataArray type=\"Float64\" Name=\"Velocities\" NumberOfComponents=\"3\" ComponentName0=\"Vx\" ComponentName1=\"Vy\" ComponentName2=\"Vz\" Format=\"ascii\">" << "\n";
    for (Node *node : upperNodes) {
        myFile << "        ";
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s[I_v + i] << " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }
    for (Node *node : lowerNodes) {
        myFile << "        ";
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s[I_v + i] << " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }
    myFile << "        </DataArray>" << "\n";

    // Output pore fluid pressure
    myFile << "        <DataArray type=\"Float64\" Name=\"Pressure\" NumberOfComponents=\"1\" Format=\"ascii\">" << "\n";  
    for (Node *node : upperNodes) {
        myFile << "        " << node->s[I_p] << "\n";
    }
    for (Node *node : lowerNodes) {
        myFile << "        " << node->s[I_p] << "\n";
    }
    myFile << "        </DataArray>" << "\n";

    // Output volumetric strain
    myFile << "        <DataArray type=\"Float64\" Name=\"Trace strain\" NumberOfComponents=\"1\" Format=\"ascii\">" << "\n";  
    for (Node *node : upperNodes) {
        myFile << "        " << node->s[I_e] << "\n";
    }
    for (Node *node : lowerNodes) {
        myFile << "        " << node->s[I_e] << "\n";
    }
    myFile << "        </DataArray>" << "\n";
    myFile << "      </PointData>" << "\n";
    myFile << "    </Piece>" << "\n";
    myFile << "  </UnstructuredGrid>" << "\n";
    myFile << "</VTKFile>" << "\n";
};

/** Write VTU files for the fault
 */
void Problem::writeVTU_fault() const{
    string path = "./output/" + outputPrefix + "Fault" + to_string(stepNumber) + ".vtu";
    ofstream myFile(path);
    
    // Needs to compensate for the offset because of bulk nodes
    int nodeIDOffset = upperNodes.size() + lowerNodes.size();

    // Head lines
    myFile << "<?xml version=\"1.0\"?>" << "\n";        // version info
    myFile << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"BigEndian\">" << "\n";
    myFile << "  <UnstructuredGrid>" << "\n";    // title
    myFile << "    <Piece NumberOfPoints=\"" << cohesiveNodes.size() << "\" NumberOfCells=\"" << cohesiveElements.size() << "\">" << "\n";                        // File Format
    
    // Output node XYZs
    myFile << "      <Points>" << "\n";
    myFile << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">" << "\n";
    myFile << scientific;
    for (CohesiveNode* node : cohesiveNodes) {
        myFile << "        " << node->getXYZ()[0] << " " << node->getXYZ()[1] << " " << "0.0" << "\n";
    }
    myFile << "        </DataArray>" << "\n";
    myFile << "      </Points>" << "\n";    
    
    // Output Elements
    myFile << "      <Cells>" << "\n";
    myFile << "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">" << "\n";
    for (ElementQ4Cohesive *element : cohesiveElements) {
        myFile << "        "
               << element->getNID()[0]->getID() - nodeIDOffset << " "
               << element->getNID()[1]->getID() - nodeIDOffset << "\n";
    }
    myFile << "        </DataArray>" << "\n";

    // Output offsets
    myFile << "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">" << "\n";
    myFile << "        ";
    for (int i = 0; i < cohesiveElements.size(); i++) myFile << 2 * (i + 1) << " ";
    myFile << "\n" << "        </DataArray>" << "\n";

    // Output cell type
    myFile << "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">" << "\n";
    myFile << "        ";
    for (int i = 0; i < cohesiveElements.size(); i++) myFile << "3 ";
    myFile << "\n" << "        </DataArray>" << "\n";
    myFile << "      </Cells>" << "\n";

    // ============================== Output Node results =================================
    int I_u = 0;
    int I_v = I_u + _spaceDim * 2;
    int I_p = I_v + _spaceDim * 2;
    int I_e = I_p + 2;
    int I_l = I_e + 2;
    int I_pf = I_l + 2;
    int I_psi = I_pf + 1;
    myFile << "      <PointData>" << "\n";
    
    // Global node ID
    myFile << "        <DataArray type=\"Int32\" Name=\"GlobalNodeId\" format=\"ascii\">" << "\n";
    myFile << "        ";
    for (CohesiveNode *node : cohesiveNodes) {
        myFile << node->getID() - nodeIDOffset << " ";
    }
    myFile << "\n" << "        </DataArray>" << "\n";

    // Nodal Slip
    myFile << "        <DataArray type=\"Float64\" Name=\"Slip\" NumberOfComponents=\"3\" ComponentName0=\"Slip_x\" ComponentName1=\"Slip_y\" ComponentName2=\"Slip_z\" Format=\"ascii\">" << "\n";
    for (CohesiveNode *node : cohesiveNodes) {
        myFile << "        ";
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s[I_u + _spaceDim + i] - node->s[I_u + i] << " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }
    myFile << "        </DataArray>" << "\n";

    // Nodal Slip_t
    myFile << "        <DataArray type=\"Float64\" Name=\"Slip_t\" NumberOfComponents=\"3\" ComponentName0=\"Slip_tx\" ComponentName1=\"Slip_ty\" ComponentName2=\"Slip_tz\" Format=\"ascii\">" << "\n";
    for (CohesiveNode *node : cohesiveNodes) {
        myFile << "        ";
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s_t[I_u + _spaceDim + i] - node->s_t[I_u + i] << " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }
    myFile << "        </DataArray>" << "\n";

    // Output velocity
    myFile << "        <DataArray type=\"Float64\" Name=\"Slip rate\" NumberOfComponents=\"3\" ComponentName0=\"SlipRate_x\" ComponentName1=\"SlipRate_y\" ComponentName2=\"SlipRate_z\" Format=\"ascii\">" << "\n";
    for (CohesiveNode *node : cohesiveNodes) {
        myFile << "        ";
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s[I_v + _spaceDim + i] - node->s[I_v + i]<< " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }
    myFile << "        </DataArray>" << "\n";

    // Output pore fluid pressure
    myFile << "        <DataArray type=\"Float64\" Name=\"Fault Pressure\" NumberOfComponents=\"1\" Format=\"ascii\">" << "\n";  
    for (CohesiveNode *node : cohesiveNodes) {
        myFile << "        " << node->s[I_pf] << "\n";
    }
    myFile << "        </DataArray>" << "\n";

    // Output psi
    myFile << "        <DataArray type=\"Float64\" Name=\"State Variable\" NumberOfComponents=\"1\" Format=\"ascii\">" << "\n";  
    for (CohesiveNode *node : cohesiveNodes) {
        myFile << "        " << node->s[I_psi] << "\n";
    }
    myFile << "        </DataArray>" << "\n";

    // Output lambda
    myFile << "        <DataArray type=\"Float64\" Name=\"Lagrangian Multplier\" NumberOfComponents=\"3\" Format=\"ascii\">" << "\n";  
    for (CohesiveNode *node : cohesiveNodes) {
        myFile << "        ";
        for (int i = 0; i < _spaceDim; i++) {
            myFile << node->s[I_l + i]<< " ";
        }
        if (_spaceDim == 2) myFile << "0.0 ";
        myFile << "\n";
    }

    myFile << "        </DataArray>" << "\n";
    myFile << "      </PointData>" << "\n";
    myFile << "    </Piece>" << "\n";
    myFile << "  </UnstructuredGrid>" << "\n";
    myFile << "</VTKFile>" << "\n";
};