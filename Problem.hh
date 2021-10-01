/** @file Problem.hh
 * Head file of class Problem
 * Initialize, solve and write the results.
 */
#include "Geometry2D.hh"
#include "ElementQ4Cohesive.hh"

/** Class Problem
 * Define a problem and solves it
 */
class Problem {
// PRIVATE MEMBERS
private:
    // Space dimension
    int _spaceDim;

    // Geometry2D
    Geometry2D* myGeometry;

    // Upper half subzone nodes
    vector<Node*> upperNodes;

    // Lower half subzone nodes
    vector<Node*> lowerNodes;

    // Cohesive nodes
    vector<CohesiveNode*> cohesiveNodes;

    // Upper half subzone elements
    vector<ElementQ4*> upperElements;

    // Lower half subzone elements
    vector<ElementQ4*> lowerElements;

    // Cohesive elements
    vector<ElementQ4Cohesive*> cohesiveElements;

    // Total degrees of freedom
    int _totalDOF;

    // Total number of nodes
    int _totalNofNodes;
    
    // Global residual vector
    Vec globalF;

    // Global solution vector
    Vec globalS;

    // Global JF
    Mat globalJF;

// PUBLIC MEMBERS
public:
    // Constructor
    Problem(int spaceDim = 2);

    // Destructor
    ~Problem();

    // Initialization
    void initialize(const vector<double> & xRanges, const vector<int> & edgeNums);

// PRIVATE MEMBERS
private:
    // Initialization of Geometry2D
    void initializeGeometry2D(const vector<double> & xRanges, const vector<int> & edgeNums);

    // Initialization of Nodes
    void initializeNodes();

    // Assign global ID for each DOF
    void assignNodalDOF();

    // TEST Try push to globalF
    void testPushGlobalF();

    // TEST Try getting from globalF
    void testFetchGlobalF();
    
    // TEST assemble globalJF
    void testAssembleGlobalJF();

    // Initialization of Elements
    void initializeElements();

    // Delete Nodes;
    void deleteNodes();

    // Delete Elements
    void deleteElements();

    // Test integrators
    void testIntegrators();

    // Compute Bodyforces
    void computeBodyForces();

    // Test integratorNfN
    void testIntegratorNfN() const;
    
    // Test integratorBfB
    void testIntegratorBfB() const;

    // Test integrator BfN
    void testIntegratorBfN() const;

    // Test integrator NfB
    void testIntegratorNfB() const;

    // TestGradient 
    void testEvaluateF_x() const;

    // Print matrix
    void printMatrix(ofstream & myFile, const vector<double>& Matrix, int nRows, int nCols) const;

// ============= Test Elastic Solution ============================================================
/** Only has 1 block upperNodes and upperElements
 * Test Petsc, Mat, Vec, KSP solver, integratorBfB
 */
// PUBLIC METHODS
public:
    // Initialize elastic problem
    void initializeElastic(const vector<double> & xRanges, const vector<int> & edgeNums);

// PRIVATE METHODS
private:
    // Initialization of Nodes
    void initializeNodesElastic();

    // Assign global ID for each DOF
    void assignNodalDOFElastic();

    // Initialization of Elements
    void initializeElementsElastic();

    // TEST Try push to globalF, push the upper surface force
    void testPushGlobalFElastic();

    // TEST Try push to globalJF, load only pushes the upper surface force
    void testPushGlobalJFElastic();

    // Linear solver
    void solveElastic();

    // Get back result from globalS
    void testFetchGlobalSElastic();
// NOT IMPLEMENTED
private:
    // Copy Constructor
    Problem(const Problem & );

    // Move-assign operator
    const Problem & operator=(const Problem & );
};
