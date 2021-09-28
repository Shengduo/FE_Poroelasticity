/** @file Problem.cc
 * Source file of class Problem.
 * Initialize, solve and write the results.
 */
#include "Problem.hh"
// Constructor
Problem::Problem(int spaceDim) {
    _spaceDim = spaceDim;
};

// Destructor
Problem::~Problem() {
    if (myGeometry) delete myGeometry;

    // Release Nodes
    deleteNodes();

    // Release Elements
    deleteElements();

    // Destroy global vectors and matrices
    VecDestroy(&globalF);
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
    vector<int> DOF_cohesive_default(spaceDim + 2, 0);
    for (int i = 0; i < myGeometry->xNodeNum; i++) {
        thisXYZ[0] = i * edgeSize[0];
        thisXYZ[1] = 0.;
        massDensity = thisXYZ[0] + thisXYZ[1];
        bodyForce  = {-thisXYZ[0], -thisXYZ[1]};
        cohesiveNodes[nodeID_in_set] = new CohesiveNode(nodeID, thisXYZ, DOF_cohesive_default, spaceDim);
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
            else {
                cohesiveNodes[i]->setDOF(j, -1); 
            }
        }
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
    
    // Output the global F
    cout << "TEST: Global F\n";
    VecView(globalF, PETSC_VIEWER_STDOUT_SELF);
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
        NID = {lowerNodes[i], lowerNodes[i + 1], upperNodes[i], upperNodes[i + 1]};
        cohesiveElements[i] = 
            new ElementQ4Cohesive(cohesiveElementST + i, NID_cohesive, NID); 
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
};

// Test all integrators
void:: Problem::testIntegrators() {
    // Test: compute body forces
    computeBodyForces();

    // Test: test gradient function
    testEvaluateF_x();

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
        cout <<"Element No." << setw(10) << element->getID() << " " << "\n";
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
        cout <<"Element No." << setw(10) << element->getID() << " " << "\n";
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
    double massDensity = 1.0;

    // Initialize some results
    vector<double> globalMassMatrix(_totalNofNodes * _totalNofNodes, 0.0);
    vector<vector<double>> eleMassMatrix(upperElements[0]->getNID().size() * upperElements[0]->getNID().size(), vector<double>(1, 0.));
    vector<vector<double>> nodalValues(upperElements[0]->getNID().size(), vector<double>(1, massDensity));

    // The global i, j indices
    int I, J;
    // Loop through upper Elements
    for (ElementQ4 *element : upperElements) {
        // Set nodal values to nodal density
        for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorNfN(eleMassMatrix, nodalValues);
        
        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalMassMatrix[_totalNofNodes * I + J] 
                    += eleMassMatrix[k * nodalValues.size() + l][0];
            }
        }
    }
    
   // Loop through lower Elements
    for (ElementQ4 *element : lowerElements) {
        // Set nodal values to nodal density
        for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorNfN(eleMassMatrix, nodalValues);
        
        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalMassMatrix[_totalNofNodes * I + J] 
                    += eleMassMatrix[k * nodalValues.size() + l][0];
            }
        }
    }
    
    // Loop through cohesive Elements
    eleMassMatrix.resize(cohesiveElements[0]->getNID().size() * cohesiveElements[0]->getNID().size());
    nodalValues.resize(cohesiveElements[0]->getNID().size());
    // Loop through cohesive Elements
    for (ElementQ4Cohesive *element : cohesiveElements) {
        // Set nodal values to nodal density
        for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorNfN(eleMassMatrix, nodalValues);
        
        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalMassMatrix[_totalNofNodes * I + J] 
                    += eleMassMatrix[k * nodalValues.size() + l][0];
            }
        }
    }

    // Printout the matrix
    ofstream myFile;
    myFile.open("Testlog.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== Global Mass Matrix ======================================" << "\n";
    printMatrix(myFile, globalMassMatrix, _totalNofNodes, _totalNofNodes);
    myFile.close();
};

/** Test integratorBfB */
void Problem::testIntegratorBfB() const {
    // Set Nodal stiffness matrix
    double massDensity = 1.0;
    vector<double> nodalD = {1., 0., 0., 1.};
    vector<vector<double>> nodalValues(4, nodalD);

    // Initialize some results
    vector<double> globalStiffMatrix(_totalNofNodes * _totalNofNodes, 0.0);
    vector<double> eleStiffMatrix(upperElements[0]->getNID().size() * upperElements[0]->getNID().size(), 0.);

    // The global i, j indices
    int I, J;
    // Loop through upper Elements
    for (ElementQ4 *element : upperElements) {
        // Set nodal values to nodal density
        // for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorBfB(eleStiffMatrix, nodalValues);
        
        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalStiffMatrix[_totalNofNodes * I + J] 
                    += eleStiffMatrix[k * nodalValues.size() + l];
            }
        }
    }
    
    // Loop through lower Elements
    for (ElementQ4 *element : lowerElements) {
        // Set nodal values to nodal density
        // for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorBfB(eleStiffMatrix, nodalValues);
        
        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalStiffMatrix[_totalNofNodes * I + J] 
                    += eleStiffMatrix[k * nodalValues.size() + l];                
            }
        }
    }
    
    // Loop through cohesive Elements
    eleStiffMatrix.resize(cohesiveElements[0]->getNID().size() * cohesiveElements[0]->getNID().size());
    nodalValues.resize(cohesiveElements[0]->getNID().size());
    for (int i = 0; i < nodalValues.size(); i++) nodalValues[i].resize(1, 1.);
    // Loop through cohesive Elements
    for (ElementQ4Cohesive *element : cohesiveElements) {
        // Set nodal values to nodal density
        // for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorBfB(eleStiffMatrix, nodalValues);

        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalStiffMatrix[_totalNofNodes * I + J] 
                    += eleStiffMatrix[k * nodalValues.size() + l];
            }
        }
    }

    // Printout the matrix
    ofstream myFile;
    myFile.open("Testlog.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== GlobalStiffMatrix ======================================" << "\n";
    printMatrix(myFile, globalStiffMatrix, _totalNofNodes, _totalNofNodes); 
    myFile.close();
};

/** Test integratorBfN */
void Problem::testIntegratorBfN() const {
    // Set Nodal BfN matrix
    double massDensity = 1.0;
    vector<double> nodalD = {1., 0.};
    vector<vector<double>> nodalValues(upperElements[0]->getNID().size(), nodalD);

    // Initialize some results
    vector<double> globalBfNMatrix(_totalNofNodes * _totalNofNodes, 0.0);
    vector<double> eleBfNMatrix(upperElements[0]->getNID().size() * upperElements[0]->getNID().size(), 0.);

    // The global i, j indices
    int I, J;
    // Loop through upper Elements
    for (ElementQ4 *element : upperElements) {
        // Set nodal values to nodal density
        // for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorBfN(eleBfNMatrix, nodalValues);
        
        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalBfNMatrix[_totalNofNodes * I + J] 
                    += eleBfNMatrix[k * nodalValues.size() + l];
            }
        }
    }
    
    // Loop through lower Elements
    for (ElementQ4 *element : lowerElements) {
        // Set nodal values to nodal density
        // for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorBfN(eleBfNMatrix, nodalValues);
        
        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalBfNMatrix[_totalNofNodes * I + J] 
                    += eleBfNMatrix[k * nodalValues.size() + l];                
            }
        }
    }
    
    // Loop through cohesive Elements
    eleBfNMatrix.resize(cohesiveElements[0]->getNID().size() * cohesiveElements[0]->getNID().size());
    nodalValues.resize(cohesiveElements[0]->getNID().size());
    for (int i = 0; i < nodalValues.size(); i++) nodalValues[i].resize(1, 1.);
    // Loop through cohesive Elements
    for (ElementQ4Cohesive *element : cohesiveElements) {
        // Set nodal values to nodal density
        // for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorBfN(eleBfNMatrix, nodalValues);

        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalBfNMatrix[_totalNofNodes * I + J] 
                    += eleBfNMatrix[k * nodalValues.size() + l];
            }
        }
    }

    // Printout the matrix
    ofstream myFile;
    myFile.open("Testlog.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== GlobalBfNMatrix ======================================" << "\n";
    printMatrix(myFile, globalBfNMatrix, _totalNofNodes, _totalNofNodes); 
    myFile.close();
};

/** Test integratorNfB */
void Problem::testIntegratorNfB() const {
    // Set Nodal NfB matrix
    double massDensity = 1.0;
    vector<double> nodalD = {1., 0.};
    vector<vector<double>> nodalValues(upperElements[0]->getNID().size(), nodalD);

    // Initialize some results
    vector<double> globalNfBMatrix(_totalNofNodes * _totalNofNodes, 0.0);
    vector<double> eleNfBMatrix(upperElements[0]->getNID().size() * upperElements[0]->getNID().size(), 0.);

    // The global i, j indices
    int I, J;
    // Loop through upper Elements
    for (ElementQ4 *element : upperElements) {
        // Set nodal values to nodal density
        // for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorNfB(eleNfBMatrix, nodalValues);
        
        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalNfBMatrix[_totalNofNodes * I + J] 
                    += eleNfBMatrix[k * nodalValues.size() + l];
            }
        }
    }
    
    // Loop through lower Elements
    for (ElementQ4 *element : lowerElements) {
        // Set nodal values to nodal density
        // for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorNfB(eleNfBMatrix, nodalValues);
        
        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalNfBMatrix[_totalNofNodes * I + J] 
                    += eleNfBMatrix[k * nodalValues.size() + l];                
            }
        }
    }
    
    // Loop through cohesive Elements
    eleNfBMatrix.resize(cohesiveElements[0]->getNID().size() * cohesiveElements[0]->getNID().size());
    nodalValues.resize(cohesiveElements[0]->getNID().size());
    for (int i = 0; i < nodalValues.size(); i++) nodalValues[i].resize(1, 1.);
    // Loop through cohesive Elements
    for (ElementQ4Cohesive *element : cohesiveElements) {
        // Set nodal values to nodal density
        // for (int n = 0; n < nodalValues.size(); n++) nodalValues[n] = {element->getNID()[n]->getMassDensity()};
       
        // Calculate the nodal values
        element->IntegratorNfB(eleNfBMatrix, nodalValues);

        for (int k = 0; k < element->getNID().size(); k++) {
            for (int l = 0; l < element->getNID().size(); l++) {
                I = element->getNID()[k]->getID();
                J = element->getNID()[l]->getID();
                globalNfBMatrix[_totalNofNodes * I + J] 
                    += eleNfBMatrix[k * nodalValues.size() + l];
            }
        }
    }

    // Printout the matrix
    ofstream myFile;
    myFile.open("Testlog.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    myFile << "\n" << "=================== GlobalNfBMatrix ======================================" << "\n";
    printMatrix(myFile, globalNfBMatrix, _totalNofNodes, _totalNofNodes); 
    myFile.close();
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
};