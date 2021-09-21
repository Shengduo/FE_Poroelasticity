/** @file main.cc
 * Initialize and solve the problem
 */
#include "Problem.hh"

/** Main program, now can:
 * Initialize nodes
 * ...
 */
int main() {
    // Initialize problem
    vector<double> spaceDomain = {1.0, 1.0};
    vector<int> nOfEdges = {1, 1};
    // Get a handle of the problem
    Problem* myProblem = new Problem();
    myProblem->initialize(spaceDomain, nOfEdges);
    return 0;
};