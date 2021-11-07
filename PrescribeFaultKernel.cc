/** @file PrescribeFaultKernel.cc
 * Declares the quasi-static prescribed poroelastic fault with pressure diffusion
 */
#include "PrescribeFaultKernel.hh"
/** The solution fields in solution vector s are
 * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
 * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
 * 8, 9 - pressure-, pressure+;
 * 10, 11 - trace-strain-, tracestrain+
 * 12, 13 - lambda1, lambda2
 * 14 - pressure_fault
 * 15 - theta (R & S state, not used here)
 */
//======= The Residuals =================================================================================   
/** Left hand side residual 
 * F0(t, \dot{s}, s)
 */
void PrescribeFaultKernel::F0(vector<double> &F0,         // stores the result
                              int spaceDim,               // stores the dim of space
                              double t,                   // time of the simulation
                              const vector<double> &s,    // solution vector s
                              const vector<double> &s_x,  // gradient of s
                              const vector<double> &s_t,  // time derivative of s
                              double s_tshift,            // sigma of tshift due to the time-derivative
                              const vector<double> &sOff, // offset of each solution field
                              const vector<double> &a,    // auxiliary fields
                              const vector<double> &a_x, // auxiliary fields gradient
                              const vector<double> &n,    // unit normal vector
                              const vector<double> &d     // prescribed slip
) {
    // Check size of F0
    if (F0.size() != 16) F0.resize(16); 

    // Clear F0
    for (int i = 0; i < F0.size(); i++) F0[i] = 0.;

    // ==================== constants ==============================
    /** Nodal properties values (
     * 0 - mass density; 
     * 1,2 - body force (force per unit volume); 
     * 3 - fluidMobility_x
     * 4 - fluidMobility_z
     * 5 - fluidViscosity \mu
     * 6 - porosity \phi_f
     * 7 - thickness h
     * 8 - beta_p 
     * 9 - beta_sigma
     * 10 - rateStateA
     * 11 - rateStateB
     * 12 - DRateState
     * 13, 14 - fluidBodyForce
     * 15 - source 
     * ...)
     */
    int i_fluidMobilityX = 3;
    int i_fluidMobilityZ = 4;
    int i_fluidViscosity = 5;
    int i_porosity = 6;
    int i_thickness = 7;
    int i_betaP = 8;
    int i_betaSigma = 9;
    int i_fluidBodyForce = 13;
    int i_source = 15;

    // ================ f0u ======================================
    /** The solution fields in solution vector s are
     * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
     * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
     * 8, 9 - pressure-, pressure+;
     * 10, 11 - trace-strain-, tracestrain+
     * 12, 13 - lambda1, lambda2
     * 14 - pressure_fault
     * 15 - theta (R & S state, not used here)
     */
    // F0u = [-\lambda, \lambda]
    int I_u = 0;
    int I_l = 12;
    F0[I_u] = - s[I_l];
    F0[I_u + 1] = - s[I_l + 1];
    F0[I_u + 2] = s[I_l];
    F0[I_u + 3] = s[I_l + 1];

    // ================ f0p ======================================
    /** The solution fields in solution vector s are
     * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
     * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
     * 8, 9 - pressure-, pressure+;
     * 10, 11 - trace-strain-, tracestrain+
     * 12, 13 - lambda1, lambda2
     * 14 - pressure_fault
     * 15 - theta (R & S state, not used here)
     */
    // F0p = [\kappa_z / \mu (p- - pf + n \cdot f_fluid), \kappa_z / \mu (p+ - pf + n \cdot f_fluid)]
    int I_p = 8;
    int I_pf = 14; 
    F0[I_p] = a[i_fluidMobilityZ] / a[i_fluidViscosity] * 
              ((s[I_p] - s[I_pf]) / a[i_thickness] + n[0] * a[i_fluidBodyForce] + n[1] * a[i_fluidBodyForce + 1]);
    F0[I_p + 1] = a[i_fluidMobilityZ] / a[i_fluidViscosity] * 
                  ((s[I_p + 1] - s[I_pf]) / a[i_thickness] - n[0] * a[i_fluidBodyForce] - n[1] * a[i_fluidBodyForce + 1]);
    
    // ================ f0l ======================================
    /** The solution fields in solution vector s are
     * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
     * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
     * 8, 9 - pressure-, pressure+;
     * 10, 11 - trace-strain-, tracestrain+
     * 12, 13 - lambda1, lambda2
     * 14 - pressure_fault
     * 15 - theta (R & S state, not used here)
     */
    // F0l = (u+) - (u-) - d
    F0[I_l] = s[I_u + 2] - s[I_u] - d[0];
    F0[I_l + 1] = s[I_u + 3] - s[I_u + 1] - d[1];

    // ================ f0pf ======================================
    /** The solution fields in solution vector s are
     * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
     * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
     * 8, 9 - pressure-, pressure+;
     * 10, 11 - trace-strain-, tracestrain+
     * 12, 13 - lambda1, lambda2
     * 14 - pressure_fault
     * 15 - theta (R & S state, not used here)
     */
    // F0pf = ...
    F0[I_pf] = a[i_porosity] * (a[i_betaP] * (s_t[I_p] + 2 * s_t[I_pf]+ s_t[I_p + 1])
                                - a[i_betaSigma] * (n[0] * s_t[I_l] + n[1] * s_t[I_l + 1]))
               + a[i_fluidMobilityX] / a[i_fluidViscosity] * (a_x[2 * i_fluidBodyForce] + a_x[2 * i_fluidBodyForce + 3])
               - a[i_fluidMobilityZ] / a[i_fluidViscosity] * (s[I_p] - 2 * s[I_pf] + s[I_p + 1]) / pow(a[i_thickness], 2) - a[i_source];
                 
};

/** Left hand side residual 
 * F1(t, \dot{s}, s)
 */
void PrescribeFaultKernel::F1(vector<double> &F1,         // stores the result
                              int spaceDim,               // stores the dim of space
                              double t,                   // time of the simulation
                              const vector<double> &s,    // solution vector s
                              const vector<double> &s_x,  // gradient of s
                              const vector<double> &s_t,  // time derivative of s
                              double s_tshift,            // sigma of tshift due to the time-derivative
                              const vector<double> &sOff, // offset of each solution field
                              const vector<double> &a,    // auxiliary fields
                              const vector<double> &a_x, // auxiliary fields gradient
                              const vector<double> &n,    // unit normal vector
                              const vector<double> &d     // prescribed slip
) {
    // Check size of F1
    if (F1.size() != (16) * spaceDim) 
        F1.resize((16) * spaceDim);

    // Clear F1
    for (int i = 0; i < F1.size(); i++) F1[i] = 0.;

    // ==================== constants ==============================
    /** Nodal properties values (
     * 0 - mass density; 
     * 1,2 - body force (force per unit volume); 
     * 3 - fluidMobility_x
     * 4 - fluidMobility_z
     * 5 - fluidViscosity \mu
     * 6 - porosity \phi_f
     * 7 - thickness h
     * 8 - beta_p 
     * 9 - beta_sigma
     * 10 - rateStateA
     * 11 - rateStateB
     * 12 - DRateState
     * 13, 14 - fluidBodyForce
     * 15 - source 
     * ...)
     */    
    int i_fluidMobilityX = 3;
    int i_fluidViscosity = 6;
    
    // ================ f1pf ====================================== 
    // F1pf = kappa_{fx} / 4\mu * \nabla (p+ + 2 pf + p-) 
    int I_p = 8;
    int I_pf = 14;
    F1[spaceDim * I_pf] = a[i_fluidMobilityX] / a[i_fluidMobilityX] 
                          * (s_x[2 * I_p] + 2 * s_x[2 * I_pf] + s_x[2 * (I_p + 1)]);
    F1[spaceDim * I_pf + 1] = a[i_fluidMobilityX] / a[i_fluidMobilityX] 
                              * (s_x[2 * I_p + 1] + 2 * s_x[2 * I_pf + 1] + s_x[2 * (I_p + 1) + 1]);
    
    // DEBUG LINES
    // cout << "F1[Ip] = " << s_x[spaceDim * I_p] << " " << s_x[spaceDim * I_p + 1] << "\n";
};

//======= The Jacobians ================================================================================= 
/** Left hand side Jacobian
 * Jf0(t, s)
 */
void PrescribeFaultKernel::Jf0(vector<double> &Jf0,        // stores the result
                            int spaceDim,               // stores the dim of space
                            double t,                   // time of the simulation
                            const vector<double> &s,    // solution vector s
                            const vector<double> &s_x,  // gradient of s
                            const vector<double> &s_t,  // time derivative of s
                            double s_tshift,            // sigma of tshift due to the time-derivative
                            const vector<double> &sOff, // offset of each solution field
                            const vector<double> &a,    // auxiliary fields
                            const vector<double> &a_x, // auxiliary fields gradient
                            PetscBool isAssembled,      // if assembled, only calculate the time-dependent parts
                            const vector<double> &n,    // unit normal vector
                            const vector<double> &d     // prescribed slip
) {
    if (!isAssembled) {
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
        
        // Use Jf[0] to store current value of Jf0pp
        // Since only (28, 29, 35) of Jf0 are used
        Jf0[0] = s_tshift / Mp;
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

        // Use Jf[1] to store current value of Jf0pe
        // Since only (28, 29, 35) of Jf0 are used
        Jf0[1] = s_tshift * alpha;

        // ================ Jf0ee ======================================
        // Jf0pp = -1
        /** Solution vector
         * 0, 1 - displacement
         * 2, 3 - velocity
         * 4 - pressure
         * 5 - trace-strain
         */
        Jf0[I_e * nCols + I_e] = -1.;
    }
    else {
        // Check size of Jf0
        if (Jf0.size() != (2 * spaceDim + 2) * (2 * spaceDim + 2)) 
            Jf0.resize((2 * spaceDim + 2) * (2 * spaceDim + 2));
        
        int nCols = (2 * spaceDim + 2);
        /**
        // First clear every entry
        for (int i = 0; i < nCols; i++) {
            for (int j = 0; j < nCols; j++)
                Jf0[i * nCols + j] = 0.;
        }
        */
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
        Jf0[I_p * nCols + I_p] = s_tshift / Mp - Jf0[0];
        Jf0[0] = s_tshift / Mp;

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
        Jf0[I_p * nCols + I_e] = s_tshift * alpha - Jf0[1];
        Jf0[1] = s_tshift * alpha;
    }
};

/** The elements of Jf0 that requires re-assemble after the first iteration
 * non-zeros are Jf0pe, Jf0ee and Jf0pp
 * only Jf0pe, Jf0pp are time dependent
 */
const vector<int> PrescribeFaultKernel::Jf0_entries = {28, 29, 35};
const vector<int> PrescribeFaultKernel::Jf0_is = {4, 4, 5};
const vector<int> PrescribeFaultKernel::Jf0_js = {4, 5, 5};
const vector<int> PrescribeFaultKernel::Jf0_timedependent = {28, 29}; 
const vector<int> PrescribeFaultKernel::Jf0_timedependent_is = {4, 4};
const vector<int> PrescribeFaultKernel::Jf0_timedependent_js = {4, 5};

/** Left hand side Jacobian
 * Jf1(t, s), uses integrator NfB
 */
void PrescribeFaultKernel::Jf1(vector<double> &Jf1,        // stores the result
                               int spaceDim,               // stores the dim of space
                               double t,                   // time of the simulation
                               const vector<double> &s,    // solution vector s
                               const vector<double> &s_x,  // gradient of s
                               const vector<double> &s_t,  // time derivative of s
                               double s_tshift,            // sigma of tshift due to the time-derivative
                               const vector<double> &sOff, // offset of each solution field
                               const vector<double> &a,    // auxiliary fields
                               const vector<double> &a_x, // auxiliary fields gradient
                               PetscBool isAssembled,      // if assembled, only calculate the time-dependent parts
                               const vector<double> &n,    // unit normal vector
                               const vector<double> &d     // prescribed slip
) {
    if (!isAssembled) {
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
    }
};

/** The elements of Jf0 that requires re-assemble after the first iteration
 * non-zeros are Jf1eu (2 terms)
 * nothing in JF1 is time dependent
 */
const vector<int> PrescribeFaultKernel::Jf1_entries = {60, 63};
const vector<int> PrescribeFaultKernel::Jf1_is = {5, 5};
const vector<int> PrescribeFaultKernel::Jf1_js = {0, 3};
const vector<int> PrescribeFaultKernel::Jf1_timedependent = {}; 
const vector<int> PrescribeFaultKernel::Jf1_timedependent_is = {};
const vector<int> PrescribeFaultKernel::Jf1_timedependent_js = {};

/** Left hand side Jacobian
 * Jf2(t, s), uses integrator BfN
 */
void PrescribeFaultKernel::Jf2(vector<double> &Jf2,        // stores the result
                               int spaceDim,               // stores the dim of space
                               double t,                   // time of the simulation
                               const vector<double> &s,    // solution vector s
                               const vector<double> &s_x,  // gradient of s
                               const vector<double> &s_t,  // time derivative of s
                               double s_tshift,            // sigma of tshift due to the time-derivative
                               const vector<double> &sOff, // offset of each solution field
                               const vector<double> &a,    // auxiliary fields
                               const vector<double> &a_x, // auxiliary fields gradient
                               PetscBool isAssembled,      // if assembled, only calculate the time-dependent parts
                               const vector<double> &n,    // unit normal vector
                               const vector<double> &d     // prescribed slip
) {
    // If first assemble
    if (!isAssembled) {
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
    }
};

/** The elements of Jf2 that requires re-assemble after the first iteration
 * the nonzeros are Jf2up (2 terms) and Jf2ue (2 terms)
 * nothing is dependent on time
 */
const vector<int> PrescribeFaultKernel::Jf2_entries = {4, 22, 5, 23};
const vector<int> PrescribeFaultKernel::Jf2_is = {0, 3, 0, 3};
const vector<int> PrescribeFaultKernel::Jf2_js = {4, 4, 5, 5};
const vector<int> PrescribeFaultKernel::Jf2_timedependent = {}; 
const vector<int> PrescribeFaultKernel::Jf2_timedependent_is = {};
const vector<int> PrescribeFaultKernel::Jf2_timedependent_js = {};

/** Left hand side Jacobian
 * Jf3(t, s)
 */
void PrescribeFaultKernel::Jf3(vector<double> &Jf3,        // stores the result
                               int spaceDim,               // stores the dim of space
                               double t,                   // time of the simulation
                               const vector<double> &s,    // solution vector s
                               const vector<double> &s_x,  // gradient of s
                               const vector<double> &s_t,  // time derivative of s
                               double s_tshift,            // sigma of tshift due to the time-derivative
                               const vector<double> &sOff, // offset of each solution field
                               const vector<double> &a,    // auxiliary fields
                               const vector<double> &a_x, // auxiliary fields gradient
                               PetscBool isAssembled,      // if assembled, only calculate the time-dependent parts
                               const vector<double> &n,    // unit normal vector
                               const vector<double> &d     // prescribed slip
) {
    if (!isAssembled) {
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
        // Does not depend on s, s_x, s_t, s_tshift, a_x
        

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
    }
    // DEBUG LINES
    // cout << "Jf3pp = " << Jf3[(I_p) * nCols + (I_p)] << " " << Jf3[(I_p + 1) * nCols + (I_p + 1)] << "\n";
};

/** The elements of Jf3 that requires re-assemble after the first iteration
 * the non-zeros are only Jf3uu (4 terms), Jf3pp
 * nothing is time dependent
 */
const vector<int> PrescribeFaultKernel::Jf3_entries = {0, 39, 13, 14, 25, 26, 104, 117};
const vector<int> PrescribeFaultKernel::Jf3_is = {0, 3, 1, 1, 2, 2, 8, 9};
const vector<int> PrescribeFaultKernel::Jf3_js = {0, 3, 1, 2, 1, 2, 8, 9};
const vector<int> PrescribeFaultKernel::Jf3_timedependent = {}; 
const vector<int> PrescribeFaultKernel::Jf3_timedependent_is = {};
const vector<int> PrescribeFaultKernel::Jf3_timedependent_js = {};