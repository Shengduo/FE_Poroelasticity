/** @file main.cc
 * Initialize and solve the problem
 */
#include "Problem.hh"
/** Main program, now can:
 * Initialize nodes
 * ...
 */
static char help[] = "VecXY interface routines. \n\n";
int main(int argc, char **argv) {
    // Initialize problem
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, NULL, help);  if (ierr) return ierr;

    // Start the log
    PetscLogDefaultBegin();
   
    vector<double> spaceDomain = {1, 1}; 
    vector<int> nOfEdges = {1, 1};
    double endingTime = 1.0;
    double dt = 0.1;
    int BulkKernel = 1;
    int CohesiveKernel = 2; // 1 - prescribed fault, 2 - friction fault

    string outputPrefix = "Fric_Trial_Kernel2_0303";
    // Get a handle of the problem
    Problem* myProblem = new Problem();
    try {
        myProblem->initializePoroElastic(spaceDomain, nOfEdges, BulkKernel, CohesiveKernel, endingTime, dt, outputPrefix);
    }
    catch (const char* msg) {
        cerr << msg << endl;
    }

    ierr  = PetscLogView(PETSC_VIEWER_STDOUT_SELF); if (ierr) return ierr;

    /** Try out some Petsc things,
     * Start from a vector.

    
    if (ierr) return ierr;
    // Petsc vec
    Vec x, y;

    // Size of vector
    PetscInt n = 20;
    
    // Create the vector
    ierr = VecCreate(PETSC_COMM_WORLD, &x); if (ierr) return ierr;
    ierr = VecSetSizes(x, PETSC_DECIDE, n); if (ierr) return ierr;
    ierr = VecSetFromOptions(x); if(ierr) return ierr;
    ierr = VecDuplicate(x, &y); if (ierr) return ierr;
    ierr = VecSet(x, PetscScalar(1.0)); if (ierr) return ierr; 
    ierr = VecSet(y, PetscScalar(1.0)); if (ierr) return ierr;
    //ierr = VecView(x, PETSC_VIEWER_STDOUT_SELF); if (ierr) return ierr;
    //ierr = VecView(y, PETSC_VIEWER_STDOUT_SELF); if (ierr) return ierr;
    ierr = PetscFinalize(); if (ierr) return ierr;
    */
    return 0;
};