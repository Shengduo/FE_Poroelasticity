/** @file PoroelasticKernel.cc
 * Implement 2D quasi-static linear poroelastic kernels
 */
#include "PoroelasticKernel.hh"
/** The solution fields are
 * 1, 2 - u1, u2;
 * 3, 4 - v1, v2; (not used here)
 * 5 - pressure
 * 6 - trace-strain
 */
//======= The Residuals =================================================================================   
/** Left hand side residual 
 * F0(t, \dot{s}, s)
 */
void PoroelasticKernel::F0(vector<double> &F0,         // stores the result
                           int spaceDim,               // stores the dim of space
                           const vector<double> &s,    // solution vector s
                           const vector<double> &s_x,  // gradient of s
                           const vector<double> &s_t,  // time derivative of s
                           double s_tshift,            // sigma of tshift due to the time-derivative
                           const vector<double> &sOff, // offset of each solution field
                           const vector<double> &a,    // auxiliary fields
                           const vector<double> &aOff  // auxiliary fields offset
) {
    // Check size of F0
    if (F0.size() != 2 * spaceDim + 2) throw "F0 size not compatible for PoroelasticKernel!";

    // Clear F0
    for (int i = 0; i < F0.size(); i++) F0[i] = 0.;

    // ==================== constants ==============================
    /** Nodal properties values (
     * 0 - mass density; 
     * 1,2 - body force (force per unit volume); 
     * 3 - \lambda (drained);
     * 4 - shear modulus G;
     * 5 - Biot coefficient \alpha;
     * 6 - Biot modulus M_p;
     * 7 - Fluid mobility \kappa;
     * 8 - Fluid viscosity \mu;
     * 9 - fluid density
     * 10 - reference porosity
     * 11, 12 - fluid body force (force per unit volume)
     * 13 - fluid source density (/s) 
     * ...)
     */
    int i_bulkDensity = 0;
    int i_bodyForce = 1;
    // int i_lambda = 3;
    // int i_shearModulus = 4;
    int i_alpha = 5;
    int i_Mp = 6;
    // int i_kappa = 7;
    // int i_viscosity = 8;
    int i_source = 13;
    
    double bulkDensity = a[i_bulkDensity];
    // double lambda = a[i_lambda];
    //  double shearModulus = a[i_shearModulus];
    double alpha = a[i_alpha];
    double Mp = a[i_Mp];
    // double kappa = a[i_kappa];
    // double viscosity = a[i_viscosity];
    double source = a[i_source];
    
    // ================ f0u ======================================
    // f0u = bodyForce
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    int I_u = 0;
    for (int i = 0; i < spaceDim; i++) {
        F0[I_u + i] = bulkDensity * a[i_bodyForce + i];
    }

    // ================ f0p ======================================
    // f0p = \partial \zeta(u, p) \partial t - \gamma 
    // ** TODO ** implement \gamma
    // \zeta(u, p) = \alpha \epsilon_v + p / M_p
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    int I_p = 4;
    int I_e = 5;
    F0[I_p] = alpha * s_t[I_e] + s_t[I_p] / Mp - source;
    // ================ f0e ======================================
    // f0e = \nabla \cdot \vec{u} - \epsilon_v
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    F0[I_e] = s_x[spaceDim * I_u] + s_x[(I_u + 1) * spaceDim + 1] - s[I_e];
};

/** Left hand side residual 
 * F1(t, \dot{s}, s)
 */
void PoroelasticKernel::F1(vector<double> &F1,         // stores the result
                           int spaceDim,               // stores the dim of space
                           const vector<double> &s,    // solution vector s
                           const vector<double> &s_x,  // gradient of s
                           const vector<double> &s_t,  // time derivative of s
                           double s_tshift,            // sigma of tshift due to the time-derivative
                           const vector<double> &sOff, // offset of each solution field
                           const vector<double> &a,    // auxiliary fields
                           const vector<double> &aOff  // auxiliary fields offset
) {
    // Check size of F1
    if (F1.size() != (2 * spaceDim + 2) * spaceDim) 
        throw "F1 size not compatible for PoroelasticKernel!";

    // Clear F1
    for (int i = 0; i < F1.size(); i++) F1[i] = 0.;

    // ==================== constants ==============================
    /** Nodal properties values (
     * 0 - mass density; 
     * 1,2 - body force (force per unit volume); 
     * 3 - \lambda (drained);
     * 4 - shear modulus G;
     * 5 - Biot coefficient \alpha;
     * 6 - Biot modulus M_p;
     * 7 - Fluid mobility \kappa;
     * 8 - Fluid viscosity \mu;
     * 9 - fluid density
     * 10 - reference porosity
     * 11, 12 - fluid body force (force per unit volume)
     * 13 - fluid source density (/s) 
     * ...)
     */
    // int i_bodyForce = 1;
    int i_lambda = 3;
    int i_shearModulus = 4;
    int i_alpha = 5;
    // int i_Mp = 6;
    int i_kappa = 7;
    int i_viscosity = 8;
    int i_fluidDensity = 9;
    int i_fluidBodyForce = 11;
    
    double lambda = a[i_lambda];
    double shearModulus = a[i_shearModulus];
    double alpha = a[i_alpha];
    // double Mp = a[i_Mp];
    double kappa = a[i_kappa];
    double fluidDensity = a[i_fluidDensity];
    double viscosity = a[i_viscosity];
    
    // ================ f1u ======================================
    // f1u = - \lambda * \epsilon_v - G * (\nabla u + u \nabla)
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    int I_u = 0;
    int I_p = 4;
    int I_e = 5;
    F1[spaceDim * I_u] = -lambda * s[I_e] - 2 * shearModulus * (s_x[spaceDim * I_u]) + alpha * s[I_p];             // - \sigma_{11} + \alpha p
    F1[spaceDim * I_u + 1] = -shearModulus * (s_x[spaceDim * I_u + 1] + s_x[spaceDim * I_u + 2]); // - \sigma_{12}
    F1[spaceDim * I_u + 2] = F1[spaceDim * I_u + 1];                                              // - \sigma_{21}
    F1[spaceDim * I_u + 3] = -lambda * s[I_e] - 2 * shearModulus * (s_x[spaceDim * I_u + 3]) + alpha * s[I_p];     // - \sigma_{22} + \alpha p

    // ================ f1p ======================================
    // f1p = - q(p) = kappa / mu \cdot (\nabla p - f_fluid)
    // ** TODO ** implement f_fluid
    // \zeta(u, p) = \alpha \epsilon_v + p / M_p
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    F1[spaceDim * I_p] = kappa / viscosity * 
        (s_x[spaceDim * I_p] - fluidDensity * a[i_fluidBodyForce]);         // \partial p \partial x_1
    F1[spaceDim * I_p + 1] = kappa / viscosity * 
        (s_x[spaceDim * I_p + 1] - fluidDensity * a[i_fluidBodyForce + 1]); // \partial p \partial x_2
    
    // DEBUG LINES
    // cout << "F1[Ip] = " << s_x[spaceDim * I_p] << " " << s_x[spaceDim * I_p + 1] << "\n";
};

//======= The Jacobians ================================================================================= 
/** Left hand side Jacobian
 * Jf0(t, s)
 */
void PoroelasticKernel::Jf0(vector<double> &Jf0,        // stores the result
                            int spaceDim,               // stores the dim of space
                            const vector<double> &s,    // solution vector s
                            const vector<double> &s_x,  // gradient of s
                            const vector<double> &s_t,  // time derivative of s
                            double s_tshift,            // sigma of tshift due to the time-derivative
                            const vector<double> &sOff, // offset of each solution field
                            const vector<double> &a,    // auxiliary fields
                            const vector<double> &aOff  // auxiliary fields offset
) {
    // Check size of Jf0
    if (Jf0.size() != (2 * spaceDim + 2) * (2 * spaceDim + 2)) 
        Jf0.resize((2 * spaceDim + 2) * (2 * spaceDim + 2));

    // First clear every entry
    int nCols = (2 * spaceDim + 2);
    for (int i = 0; i < nCols; i++) {
        for (int j = 0; j < nCols; j++)
            Jf0[i * nCols + j] = 0.;
    }
    // ==================== constants ==============================
    /** Nodal properties values (
     * 0 - mass density; 
     * 1,2 - body force (force per unit volume); 
     * 3 - \lambda (drained);
     * 4 - shear modulus G;
     * 5 - Biot coefficient \alpha;
     * 6 - Biot modulus M_p;
     * 7 - Fluid mobility \kappa;
     * 8 - Fluid viscosity \mu;
     * 9 - fluid density
     * 10 - reference porosity
     * 11, 12 - fluid body force (force per unit volume)
     * 13 - fluid source density (/s) 
     * ...)
     */
    //int i_lambda = 3;
    //int i_shearModulus = 4;
    int i_alpha = 5;
    int i_Mp = 6;
    //int i_kappa = 7;
    //int i_viscosity = 8;
    //double lambda = a[i_lambda];
    //double shearModulus = a[i_shearModulus];
    double alpha = a[i_alpha];
    double Mp = a[i_Mp];
    //double kappa = a[i_kappa];
    //double viscosity = a[i_viscosity];

    // ================ Jf0pp ======================================
    // Jf0pp = s_tshift / M_p
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    int I_p = 4;
    Jf0[I_p * nCols + I_p] = s_tshift / Mp;

    // DEBUG LINES
    // cout << "Jf0pp = " << Jf0[I_p * nCols + I_p] << "\n";

    // ================ Jf0pe ======================================
    // Jf0pe = s_tshift * alpha
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    int I_e = 5;
    Jf0[I_p * nCols + I_e] = s_tshift * alpha;

    // ================ Jf0ee ======================================
    // Jf0pp = -1
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    Jf0[I_e * nCols + I_e] = -1.;
};

/** Left hand side Jacobian
 * Jf1(t, s), uses integrator NfB
 */
void PoroelasticKernel::Jf1(vector<double> &Jf1,        // stores the result
                            int spaceDim,               // stores the dim of space
                            const vector<double> &s,    // solution vector s
                            const vector<double> &s_x,  // gradient of s
                            const vector<double> &s_t,  // time derivative of s
                            double s_tshift,            // sigma of tshift due to the time-derivative
                            const vector<double> &sOff, // offset of each solution field
                            const vector<double> &a,    // auxiliary fields
                            const vector<double> &aOff  // auxiliary fields offset
) {
    // Check size of Jf1
    
    int nCols = (2 * spaceDim + 2) * spaceDim;
    int nRows = (2 * spaceDim + 2);
    if (Jf1.size() != nRows * nCols) 
        Jf1.resize(nCols * nRows);

    // First clear every entry
    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++)
            Jf1[i * nCols + j] = 0.;
    }
    // ==================== constants ==============================
    /** Nodal properties values (
     * 0 - mass density; 
     * 1,2 - body force (force per unit volume); 
     * 3 - \lambda (drained);
     * 4 - shear modulus G;
     * 5 - Biot coefficient \alpha;
     * 6 - Biot modulus M_p;
     * 7 - Fluid mobility \kappa;
     * 8 - Fluid viscosity \mu;
     * 9 - fluid density
     * 10 - reference porosity
     * 11, 12 - fluid body force (force per unit volume)
     * 13 - fluid source density (/s) 
     * ...)
     */
    // int i_lambda = 3;
    // int i_shearModulus = 4;
    // int i_alpha = 5;
    // int i_Mp = 6;
    // int i_kappa = 7;
    // int i_viscosity = 8;
    // double lambda = a[i_lambda];
    // double shearModulus = a[i_shearModulus];
    // double alpha = a[i_alpha];
    // double Mp = a[i_Mp];
    // double kappa = a[i_kappa];
    // double viscosity = a[i_viscosity];

    // ================ Jf1eu ======================================
    // Jf0pp = s_ts
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    int I_e = 5;
    int I_u = 0;
    Jf1[I_e * nCols + I_u * spaceDim] = 1.;
    Jf1[I_e * nCols + (I_u + 1) * spaceDim + 1] = 1.;
};

/** Left hand side Jacobian
 * Jf2(t, s), uses integrator BfN
 */
void PoroelasticKernel::Jf2(vector<double> &Jf2,        // stores the result
                            int spaceDim,               // stores the dim of space
                            const vector<double> &s,    // solution vector s
                            const vector<double> &s_x,  // gradient of s
                            const vector<double> &s_t,  // time derivative of s
                            double s_tshift,            // sigma of tshift due to the time-derivative
                            const vector<double> &sOff, // offset of each solution field
                            const vector<double> &a,    // auxiliary fields
                            const vector<double> &aOff  // auxiliary fields offset
) {
    // Check size of Jf2
    
    int nRows = (2 * spaceDim + 2) * spaceDim;
    int nCols = (2 * spaceDim + 2);
    if (Jf2.size() != nRows * nCols) 
        Jf2.resize(nRows * nCols);

    // First clear every entry
    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++)
            Jf2[i * nCols + j] = 0.;
    }
    // ==================== constants ==============================
    /** Nodal properties values (
     * 0 - mass density; 
     * 1,2 - body force (force per unit volume); 
     * 3 - \lambda (drained);
     * 4 - shear modulus G;
     * 5 - Biot coefficient \alpha;
     * 6 - Biot modulus M_p;
     * 7 - Fluid mobility \kappa;
     * 8 - Fluid viscosity \mu;
     * ...)
     */
    int i_lambda = 3;
    // int i_shearModulus = 4;
    int i_alpha = 5;
    // int i_Mp = 6;
    // int i_kappa = 7;
    // int i_viscosity = 8;
    double lambda = a[i_lambda];
    // double shearModulus = a[i_shearModulus];
    double alpha = a[i_alpha];
    // double Mp = a[i_Mp];
    // double kappa = a[i_kappa];
    // double viscosity = a[i_viscosity];

    // ================ Jf2up ======================================
    // Jf2up = \alpha \delta_{ij}
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    int I_p = 4;
    int I_u = 0;
    Jf2[(I_u * spaceDim) * nCols + I_p] = alpha;
    Jf2[((I_u + 1) * spaceDim + 1) * nCols + I_p] = alpha;

    // ================ Jf2ue ======================================
    // Jf2up = \alpha \delta_{ij}
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    int I_e = 5;
    Jf2[(I_u * spaceDim) * nCols + I_e] = - lambda;
    Jf2[((I_u + 1) * spaceDim + 1) * nCols + I_e] = - lambda;
};

/** Left hand side Jacobian
 * Jf3(t, s)
 */
void PoroelasticKernel::Jf3(vector<double> &Jf3,        // stores the result
                            int spaceDim,               // stores the dim of space
                            const vector<double> &s,    // solution vector s
                            const vector<double> &s_x,  // gradient of s
                            const vector<double> &s_t,  // time derivative of s
                            double s_tshift,            // sigma of tshift due to the time-derivative
                            const vector<double> &sOff, // offset of each solution field
                            const vector<double> &a,    // auxiliary fields
                            const vector<double> &aOff  // auxiliary fields offset
) {
    // Check size of Jf3
    if (Jf3.size() != spaceDim * spaceDim * (2 * spaceDim + 2) * (2 * spaceDim + 2)) 
        Jf3.resize(spaceDim * spaceDim * (2 * spaceDim + 2) * (2 * spaceDim + 2));
    // First clear Jf3uu
    int nCols = (2 * spaceDim + 2) * spaceDim;
    for (int i = 0; i < spaceDim * spaceDim; i++) {
        for (int j = 0; j < spaceDim * spaceDim; j++)
            Jf3[i * nCols + j] = 0.;
    }
    // ==================== constants ==============================
    /** Nodal properties values (
     * 0 - mass density; 
     * 1,2 - body force (force per unit volume); 
     * 3 - \lambda (drained);
     * 4 - shear modulus G;
     * 5 - Biot coefficient \alpha;
     * 6 - Biot modulus M_p;
     * 7 - Fluid mobility \kappa;
     * 8 - Fluid viscosity \mu;
     * 9 - fluid density
     * 10 - reference porosity
     * 11, 12 - fluid body force (force per unit volume)
     * 13 - fluid source density (/s) 
     * ...)
     */
    // int i_lambda = 3;
    int i_shearModulus = 4;
    // int i_alpha = 5;
    // int i_Mp = 6;
    int i_kappa = 7;
    int i_viscosity = 8;
    // double lambda = a[i_lambda];
    double shearModulus = a[i_shearModulus];
    // double alpha = a[i_alpha];
    // double Mp = a[i_Mp];
    double kappa = a[i_kappa];
    double viscosity = a[i_viscosity];

    // ================ Jf3uu ======================================
    // Poroelastic Cijkl, without lambda
    // Only needs lambda and G
    // Does not depend on s, s_x, s_t, s_tshift, aOff
    

    // Write down Cijkl
    Jf3[0 * nCols + 0] = - 2. * shearModulus;
    Jf3[3 * nCols + 3] = Jf3[0 * nCols + 0];

    // Jf3[0 * nCols + 3] = lambda;  // 1122
    // Jf3[3 * nCols + 0] = lambda;  // 2211

    Jf3[1 * nCols + 1] = - shearModulus; // 1212
    Jf3[2 * nCols + 2] = - shearModulus; // 2121
    Jf3[1 * nCols + 2] = - shearModulus; // 1221
    Jf3[2 * nCols + 1] = - shearModulus; // 2112

    // ================ Jf3pp ======================================
    // Jf3pp = \kappa / \mu * \delta_{kl} 
    /** Solution vector
     * 0, 1 - displacement
     * 2, 3 - velocity
     * 4 - pressure
     * 5 - trace-strain
     */
    int I_p = 4 * spaceDim;
    for (int i = 0; i < spaceDim; i++) {
        Jf3[(I_p + i) * nCols + (I_p + i)] = kappa / viscosity;
    }

    // DEBUG LINES
    // cout << "Jf3pp = " << Jf3[(I_p) * nCols + (I_p)] << " " << Jf3[(I_p + 1) * nCols + (I_p + 1)] << "\n";
}