/** @file main.cc
 * Initialize and solve the problem
 */
#include "petscmat.h"
#include "Problem.hh"
/** Main program, now can:
 * Initialize nodes
 * ...
 */
static char help[] = "VecXY interface routines. \n\n";
int main(int argc, char **argv) {
    
    // Initialize problem
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, NULL, help);  if (ierr) return ierr;
    vector<double> spaceDomain = {4, 6}; 
    vector<int> nOfEdges = {4, 6};
    double endingTime = 5.0;
    double dt = 0.1;
    int nMatPartition = 8;

    
    // Get a handle of the problem
    PetscMPIInt rank, size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    try {
        Problem* myProblem = new Problem();
        myProblem->initializePoroElastic(spaceDomain, nOfEdges, endingTime, dt, nMatPartition);
        ierr = myProblem->initializePetsc(argc, argv);  if (ierr) return ierr; 
        myProblem->solvePoroElastic();    
    }
    catch (const char* msg) {
        cerr << msg << endl;
    }

    // ierr  = PetscLogView(PETSC_VIEWER_STDOUT_SELF); if (ierr) return ierr;
    ierr = PetscFinalize(); 
    return ierr;    


    // Try out some Petsc things,
    // * Start from a vector.

    /**
    if (ierr) return ierr;
    // Petsc Mat
    Mat myMat;

    // Size of vector
    PetscInt n = 8;
    
    // Local row start and end
    PetscInt rowStart, rowEnd;

    // Create the vector
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 8,  8,  2, NULL, 4, NULL, &myMat); if (ierr) return ierr;
    ierr = MatSetType(myMat, MATMPIAIJ); if (ierr) return ierr;
    // ierr = MatSetSizes(myMat, PETSC_DECIDE, PETSC_DECIDE, n, n); if(ierr) return ierr;
    ierr = MatSetFromOptions(myMat); if(ierr) return ierr;
    ierr = MatMPIAIJSetPreallocation(myMat, 2, NULL, 4, NULL);

    // Get start and end row number
    ierr = MatGetOwnershipRange(myMat, &rowStart, &rowEnd); if (ierr) return ierr;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    cout << "Rank: " << rank << "\n";
    cout << "Row Start: " << rowStart << "\n";
    cout << "Row End: " << rowEnd << "\n";

    ierr = MatAssemblyBegin(myMat, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(myMat, MAT_FINAL_ASSEMBLY); if (ierr) return ierr;
    // ierr = MatView(myMat, PETSC_VIEWER_STDOUT_SELF); if (ierr) return ierr;
    //ierr = VecView(x, PETSC_VIEWER_STDOUT_SELF); if (ierr) return ierr;
    //ierr = VecView(y, PETSC_VIEWER_STDOUT_SELF); if (ierr) return ierr;
    ierr = PetscFinalize(); if (ierr) return ierr;
    return 0;
    */
};