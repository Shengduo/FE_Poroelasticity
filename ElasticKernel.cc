/** @file ElasticKernel.cc
 * Implement 2D quasi-static linear elastic kernels
 */
#include "ElasticKernel.hh"

//======= The Residuals =================================================================================   
/** Left hand side residual 
 * F0(t, \dot{s}, s)
 */
void ElasticKernel::F0(vector<double> &F0,         // stores the result
                       int spaceDim,               // stores the dim of space
                       const vector<double> &s,    // solution vector s
                       const vector<double> &s_x,  // gradient of s
                       const vector<double> &s_t,  // time derivative of s
                       double s_tshift,            // sigma of tshift due to the time-derivative
                       const vector<double> &sOff, // offset of each solution field
                       const vector<double> &a,    // auxiliary fields
                       const vector<double> &aOff  // auxiliary fields offset
) const {

};

/** Left hand side Jacobian
 * Jf3(t, s)
 */
void ElasticKernel::Jf3(vector<double> &Jf3,        // stores the result
                        int spaceDim,               // stores the dim of space
                        const vector<double> &s,    // solution vector s
                        const vector<double> &s_x,  // gradient of s
                        const vector<double> &s_t,  // time derivative of s
                        double s_tshift,            // sigma of tshift due to the time-derivative
                        const vector<double> &sOff, // offset of each solution field
                        const vector<double> &a,    // auxiliary fields
                        const vector<double> &aOff  // auxiliary fields offset
) const {
    // Check size of Jf3
    if (Jf3.size() != spaceDim * spaceDim * (2 * spaceDim + 2) * (2 * spaceDim + 2)) throw "Jf3 size not compatible for ElasticKernel!";
    // First clear Jf3uu
    int nCols = (2 * spaceDim + 2) * spaceDim;
    for (int i = 0; i < spaceDim * spaceDim; i++) {
        for (int j = 0; j < spaceDim * spaceDim; j++)
            Jf3[i * nCols + j] = 0.;
    }
    

    // Elastic Cijkl
    // Only needs lambda and G
    // Does not depend on s, s_x, s_t, s_tshift, aOff
    int spaceDim = 2;
    int i_lambda = 3;
    int i_shearModulus = 4;
    double lambda = a[i_lambda];
    double shearModulus = a[i_shearModulus];

    // Write down Cijkl
    Jf3[0 * nCols + 0] = lambda + 2. * shearModulus;
    Jf3[3 * nCols + 3] = Jf3[0 * nCols + 0];

    Jf3[0 * nCols + 3] = lambda;  // 1122
    Jf3[3 * nCols + 0] = lambda;  // 2211

    Jf3[1 * nCols + 1] = shearModulus; // 1212
    Jf3[2 * nCols + 2] = shearModulus; // 2121
    Jf3[1 * nCols + 2] = shearModulus; // 1221
    Jf3[2 * nCols + 1] = shearModulus; // 2112
}