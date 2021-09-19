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
};

// General initialization
void Problem::initialize(const vector<double> & xRanges, const vector<int> & edgeNums) {
    
    // Initialize geometry2D
    if (_spaceDim == 2) initializeGeometry2D(xRanges, edgeNums);
    
    // Initialize nodes
    initializeNodes();

    // Initialize elements
    initializeElements();
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

    // Upper subzone nodes
    upperNodes.resize(myGeometry->nOfNodes);
    int nodeID = 0;
    int nodeID_in_set = 0;
    vector<double> thisXYZ(spaceDim);
    
    // First y
    for (int j = 0; j < myGeometry->yNodeNum; j++) {
        // Then x
        for (int i = 0; i < myGeometry->xNodeNum; i++) {
            // Reset the coordinates
            thisXYZ[0] = i * edgeSize[0];
            thisXYZ[1] = j * edgeSize[1];
            // Initialize a node
            upperNodes[nodeID_in_set] = new Node(nodeID, thisXYZ, DOF_default, spaceDim);
            if (j == myGeometry->yNodeNum - 1) // Upper surface
                upperNodes[nodeID_in_set]->setDOF(DOF_vp_fixed);
            nodeID += 1;
            nodeID_in_set += 1;
        }
    }

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
            // Initialize a node
            lowerNodes[nodeID_in_set] = new Node(nodeID, thisXYZ, DOF_default, spaceDim);
            if (j == myGeometry->yNodeNum - 1) // Lower surface
                lowerNodes[nodeID_in_set]->setDOF(DOF_vp_fixed);
            nodeID += 1;
            nodeID_in_set += 1;
        }
    }
    nodeID_in_set = 0;
    
    // Cohesive Nodes;
    cohesiveNodes.resize(myGeometry->xNodeNum);
    vector<int> DOF_cohesive_default(spaceDim + 1, 0);
    for (int i = 0; i < myGeometry->xNodeNum; i++) {
        thisXYZ[0] = i * edgeSize[0];
        thisXYZ[1] = 0.;
        cohesiveNodes[nodeID_in_set] = new CohesiveNode(nodeID, thisXYZ, DOF_cohesive_default, spaceDim);
        nodeID += 1;
        nodeID_in_set += 1;
    };

    ofstream myFile;
    myFile.open("NodeInfo.txt");

    for (Node* node : upperNodes) node->outputInfo(myFile, true);
    for (Node* node : lowerNodes) node->outputInfo(myFile, true);
    for (CohesiveNode* node : cohesiveNodes) node->outputInfo(myFile, true);
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
                   lowerNodes[(j + 1) * myGeometry->xNodeNum + i], 
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
    myFile.open("ElementInfo.txt");

    for (ElementQ4* thisElement : upperElements) thisElement->outputInfo(myFile);
    for (ElementQ4* thisElement : lowerElements) thisElement->outputInfo(myFile);
    for (ElementQ4Cohesive* thisElement : cohesiveElements) thisElement->outputInfo(myFile);
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
