/** @file Problem.hh
 * Head file of class Problem
 * Initialize, solve and write the results.
 */
#pragma once
#include "Geometry2D.hh"
#include "ElementQ4Cohesive.hh"
#include <time.h>

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

    // Row numbers
    PetscInt *globalRows;
    
    // Node info time
    double nodeTime;

    // Step number
    int stepNumber;
    
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

    // Global time-derivative 
    Vec globalS_t;

    // Global JF
    Mat globalJF;

    // Total non-zeros
    int _totalNonZeros;
    
    // Local residual
    double *localF;
    int localFSize;

    // Local Jacobian
    double *localJF;
    int localJFSize;

// PUBLIC MEMBERS
public:
    // Time consumed by Jacobian and the clocks
    vector<double> timeConsumed;
    vector<clock_t> clocks;
    
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

    // Non-zeros per row in global matrices
    void getNNZPerRow(PetscInt *nnz);

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

// ============= Test Poroelastic Solution ============================================================
/** Only has 1 block upperNodes and upperElements
 * Test Petsc, Mat, Vec, SNES TS solver, integratorBfB
 */
// PUBLIC METHODS
public:
    // Initialize elastic problem
    void initializePoroElastic(const vector<double> & xRanges, const vector<int> & edgeNums, double endingTime = 1.0, double dt = 0.01);

// PRIVATE METHODS
private:
    // Initialization of Nodes
    void initializeNodesPoroElastic();

    // Assign global ID for each DOF
    void assignNodalDOFPoroElastic();

    // Initialization of Elements
    void initializeElementsPoroElastic();

    // Initialize the Vecs and Mats and TS
    void initializePetsc();

    // TS solver
    void solvePoroElastic(double endingTime = 1.0, double dt = 0.01);

    // Residual function
    static PetscErrorCode IFunction(TS ts, PetscReal t, Vec s, Vec s_t, Vec F, void *ctx = NULL);
    
    // Jacobian function
    static PetscErrorCode IJacobian(TS ts, PetscReal t, Vec s, Vec s_t, PetscReal s_tshift, Mat Amat, Mat Pmat, void *ctx = NULL);

    // Write VTK files
    void writeVTK(string prefix);

    // Write VTU files
    void writeVTU(string prefix);
    
// NOT IMPLEMENTED
private:
    // Copy Constructor
    Problem(const Problem & );

    // Move-assign operator
    const Problem & operator=(const Problem & );
};
