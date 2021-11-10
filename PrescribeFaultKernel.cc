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
                              const vector<double> &n    // unit normal vector
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
     * 16, 17 - slip
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
    int i_d = 16;

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
    F0[I_l] = s[I_u + 2] - s[I_u] - a[i_d];
    F0[I_l + 1] = s[I_u + 3] - s[I_u + 1] - a[i_d + 1];

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
                              const vector<double> &n    // unit normal vector
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
    F1[spaceDim * I_pf] = a[i_fluidMobilityX] / a[i_fluidViscosity] 
                          * (s_x[2 * I_p] + 2 * s_x[2 * I_pf] + s_x[2 * (I_p + 1)]);
    F1[spaceDim * I_pf + 1] = a[i_fluidMobilityX] / a[i_fluidViscosity] 
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
                               const vector<double> &a_x,  // auxiliary fields gradient
                               PetscBool isAssembled,      // if assembled, only calculate the time-dependent parts
                               const vector<double> &n    // unit normal vector
) {
    int nCols = 16;
    // Check if the system jacobian has been assembled
    if (!isAssembled) {
        // Check size of Jf0
        if (Jf0.size() != 16 * 16) 
            Jf0.resize(16 * 16);

        // First clear every entry
        for (int i = 0; i < nCols; i++) {
            for (int j = 0; j < nCols; j++)
                Jf0[i * nCols + j] = 0.;
        }

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
        // int i_fluidMobilityX = 3;
        int i_fluidMobilityZ = 4;
        int i_fluidViscosity = 5;
        int i_porosity = 6;
        int i_thickness = 7;
        int i_betaP = 8;
        int i_betaSigma = 9;
        // int i_fluidBodyForce = 13;
        // int i_source = 15;        

        // ================ Jf0ul ======================================
        // Jf0ul = [-I; I];
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */
        int I_u = 0;
        int I_l = 12;
        for (int i = 0; i < spaceDim; i++) {
            Jf0[(I_u + i) * nCols + (I_l + i)] = -1;
            Jf0[(I_u + spaceDim + i) * nCols + (I_l + i)] = 1;
        };
        
        // DEBUG LINES
        // cout << "Jf0pp = " << Jf0[I_p * nCols + I_p] << "\n";

        // ================ Jf0pp ======================================
        // Jf0pp = [\kappa_cz / (\mu h), 0; 0, \kappa_cz / (\mu h)]
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */
        int I_p = 8;
        Jf0[I_p * nCols + I_p] = a[i_fluidMobilityZ] / a[i_fluidViscosity] / a[i_thickness];
        Jf0[(I_p + 1) * nCols + (I_p + 1)] = Jf0[I_p * nCols + I_p];
        
        // ================ Jf0ppf ======================================
        // Jf0ppf = [- \kappa_cz / (\mu h); -\kappa_cz / (\mu h)]
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */

        int I_pf = 14;
        Jf0[I_p * nCols + I_pf] = -a[i_fluidMobilityZ] / a[i_fluidViscosity] / a[i_thickness];
        Jf0[(I_p + 1) * nCols + I_pf] = Jf0[I_p * nCols + I_pf];
        
        // ================ Jf0lu ======================================
        // Jf0lu = [-I, I];
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */

        for (int i = 0; i < spaceDim; i++) {
            Jf0[(I_l + i) * nCols + (I_u + i)] = -1;
            Jf0[(I_l + i) * nCols + (I_u + spaceDim + i)] = 1;
        };
        
        // ================ Jf0pfpf ======================================
        // Jf0pfpf = 2 \kappa_z / \mu / h^2 + 2 \phi_f * \beta_P * s_tshift;
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */

        Jf0[I_pf * nCols + I_pf] = 2 * a[i_fluidMobilityZ] / a[i_fluidViscosity] / pow(a[i_thickness], 2)
                                   + 2 * a[i_porosity] * a[i_betaP] * s_tshift;
        
        // Jf0[0][0] is not used, use to store this
        Jf0[0 * nCols + 0] = Jf0[I_pf * nCols + I_pf];

        // ================ Jf0pfl ======================================
        // Jf0pfl = - \phi_f * \beta_Sigma * s_tshift * [n[0], n[1]];
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */

        Jf0[I_pf * nCols + I_l] = - a[i_porosity] * a[i_betaSigma] * s_tshift * n[0];
        Jf0[I_pf * nCols + I_l + 1] = - a[i_porosity] * a[i_betaSigma] * s_tshift * n[1];
        // Jf0[0][1] and Jf0[0][2] are not used, use to store this
        Jf0[0 * nCols + 1] = Jf0[I_pf * nCols + I_l];
        Jf0[0 * nCols + 2] = Jf0[I_pf * nCols + I_l + 1];

        // ================ Jf0pfp ======================================
        // Jf0pfp = \phi_f * \beta_P / 4 * s_tshift - \kappa_z / (\mu h^2);
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */

        Jf0[I_pf * nCols + I_p] = a[i_porosity] * a[i_betaP] * s_tshift / 4. - a[i_fluidMobilityZ] / a[i_fluidViscosity] / pow(a[i_thickness], 2);
        Jf0[I_pf * nCols + I_p + 1] = a[i_porosity] * a[i_betaP] * s_tshift / 4. - a[i_fluidMobilityZ] / a[i_fluidViscosity] / pow(a[i_thickness], 2);
        // Jf0[0][3] and Jf0[0][4] are not used, use to store this
        Jf0[0 * nCols + 3] = Jf0[I_pf * nCols + I_p];
        Jf0[0 * nCols + 4] = Jf0[I_pf * nCols + I_p + 1];
    }

    else {
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
        // int i_fluidMobilityX = 3;
        int i_fluidMobilityZ = 4;
        int i_fluidViscosity = 5;
        int i_porosity = 6;
        int i_thickness = 7;
        int i_betaP = 8;
        int i_betaSigma = 9;
        // int i_fluidBodyForce = 13;
        // int i_source = 15;
        
        // ================ Jf0pfpf ======================================
        // Jf0pfpf = 2 \kappa_z / \mu / h^2 + 2 \phi_f * \beta_Sigma * s_tshift;
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */
        int I_pf = 14;
        Jf0[I_pf * nCols + I_pf] = 2 * a[i_fluidMobilityZ] / a[i_fluidViscosity] / pow(a[i_thickness], 2)
                                   + 2 * a[i_porosity] * a[i_betaP] * s_tshift 
                                   - Jf0[0 * nCols + 0];
        
        // Jf0[0][0] is not used, use to store this
        Jf0[0 * nCols + 0] = 2 * a[i_fluidMobilityZ] / a[i_fluidViscosity] / pow(a[i_thickness], 2)
                             + 2 * a[i_porosity] * a[i_betaP] * s_tshift;
        
        // ================ Jf0pfl ======================================
        // Jf0pfl = - \phi_f * \beta_Sigma * s_tshift * [n[0], n[1]];
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */

        int I_l = 12;
        Jf0[I_pf * nCols + I_l] = - a[i_porosity] * a[i_betaSigma] * s_tshift * n[0] - Jf0[0 * nCols + 1];
        Jf0[I_pf * nCols + I_l + 1] = - a[i_porosity] * a[i_betaSigma] * s_tshift * n[1] - Jf0[0 * nCols + 2];
        // Jf0[0][1] and Jf0[0][2] are not used, use to store this
        Jf0[0 * nCols + 1] = - a[i_porosity] * a[i_betaSigma] * s_tshift * n[0];
        Jf0[0 * nCols + 2] = - a[i_porosity] * a[i_betaSigma] * s_tshift * n[1];

        // ================ Jf0pfp ======================================
        // Jf0pfp = \phi_f * \beta_P / 4 * s_tshift - \kappa_z / (\mu h^2);
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */
        
        int I_p = 8;
        Jf0[I_pf * nCols + I_p] = a[i_porosity] * a[i_betaP] * s_tshift / 4. - a[i_fluidMobilityZ] / a[i_fluidViscosity] / pow(a[i_thickness], 2) - Jf0[0 * nCols + 3];
        Jf0[I_pf * nCols + I_p + 1] = a[i_porosity] * a[i_betaP] * s_tshift / 4. - a[i_fluidMobilityZ] / a[i_fluidViscosity] / pow(a[i_thickness], 2) - Jf0[0 * nCols + 4];
        // Jf0[0][3] and Jf0[0][4] are not used, use to store this
        Jf0[0 * nCols + 3] = a[i_porosity] * a[i_betaP] * s_tshift / 4. - a[i_fluidMobilityZ] / a[i_fluidViscosity] / pow(a[i_thickness], 2);
        Jf0[0 * nCols + 4] = a[i_porosity] * a[i_betaP] * s_tshift / 4. - a[i_fluidMobilityZ] / a[i_fluidViscosity] / pow(a[i_thickness], 2);
    }
};

/** The elements of Jf0 that requires re-assemble after the first iteration
 * non-zeros are Jf0ul, Jf0pp and Jf0ppf and Jf0lu (12 non-zeros in total)
 * Jf0pfpf,  is time dependent
 */
// const vector<int> PrescribeFaultKernel::Jf0_entries = {28, 29, 35};
const vector<int> PrescribeFaultKernel::Jf0_is = {0, 1, 2, 3, 8, 9, 8, 9, 12, 12, 13, 13, 14, 14, 14, 14, 14};
const vector<int> PrescribeFaultKernel::Jf0_js = {12, 13, 12, 13, 8, 9, 14, 14, 0, 2, 1, 3, 14, 12, 13, 8, 9};
// const vector<int> PrescribeFaultKernel::Jf0_timedependent = {}; 
const vector<int> PrescribeFaultKernel::Jf0_timedependent_is = {14, 14, 14, 14, 14};
const vector<int> PrescribeFaultKernel::Jf0_timedependent_js = {14, 12, 13, 8, 9};

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
                               const vector<double> &n    // unit normal vector
) {
    if (!isAssembled) {
        // Check size of Jf1    
        int nCols = 16 * spaceDim;
        int nRows = 16;
        if (Jf1.size() != nRows * nCols) 
            Jf1.resize(nCols * nRows);

        // First clear every entry
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++)
                Jf1[i * nCols + j] = 0.;
        }
    }
};

/** The elements of Jf1 that requires re-assemble after the first iteration
 * non-zeros are NULL (0 terms)
 * nothing in JF1 is time dependent
 */
// const vector<int> PrescribeFaultKernel::Jf1_entries = {};
const vector<int> PrescribeFaultKernel::Jf1_is = {};
const vector<int> PrescribeFaultKernel::Jf1_js = {};
// const vector<int> PrescribeFaultKernel::Jf1_timedependent = {}; 
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
                               const vector<double> &a_x,  // auxiliary fields gradient
                               PetscBool isAssembled,      // if assembled, only calculate the time-dependent parts
                               const vector<double> &n     // unit normal vector
) {
    // If first assemble
    if (!isAssembled) {
        // Check size of Jf2    
        int nRows = 16 * spaceDim;
        int nCols = 16;
        if (Jf2.size() != nRows * nCols) 
            Jf2.resize(nRows * nCols);

        // First clear every entry
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++)
                Jf2[i * nCols + j] = 0.;
        }
    }
};

/** The elements of Jf2 that requires re-assemble after the first iteration
 * the nonzeros are NULL
 * nothing is dependent on time
 */
// const vector<int> PrescribeFaultKernel::Jf2_entries = {4, 22, 5, 23};
const vector<int> PrescribeFaultKernel::Jf2_is = {};
const vector<int> PrescribeFaultKernel::Jf2_js = {};
// const vector<int> PrescribeFaultKernel::Jf2_timedependent = {}; 
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
                               const vector<double> &n    // unit normal vector
) {
    if (!isAssembled) {
        // Check size of Jf3
        if (Jf3.size() != spaceDim * spaceDim * 16 * 16) 
            Jf3.resize(spaceDim * spaceDim * 16 * 16);
        // First clear Jf3uu
        int nCols = 16 * spaceDim;
        for (int i = 0; i < nCols; i++) {
            for (int j = 0; j < nCols; j++)
                Jf3[i * nCols + j] = 0.;
        }
        
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
        // int i_fluidMobilityZ = 4;
        int i_fluidViscosity = 5;
        // int i_porosity = 6;
        // int i_thickness = 7;
        // int i_betaP = 8;
        // int i_betaSigma = 9;
        // int i_fluidBodyForce = 13;
        // int i_source = 15;  

        // ================ Jf3pfpf ======================================
        // Jf3pfpf = \kappa_fx / 2\mu I
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */

        int I_pf = 14;
        Jf3[(spaceDim * I_pf) * nCols + (spaceDim * I_pf)] = a[i_fluidMobilityX] / a[i_fluidViscosity] / 2.;
        Jf3[(spaceDim * I_pf + 1) * nCols + (spaceDim * I_pf + 1)] = a[i_fluidMobilityX] / a[i_fluidViscosity] / 2.;

        // ================ Jf3pfp ======================================
        // Jf3pfp = \kappa_x / \mu / 4 I 
        /** The solution fields in solution vector s are
         * 0, 1, 2, 3 - u1-, u2-, u1+, u2+;
         * 4, 5, 6, 7 - v1-, v2-, v1+, v2+; (not used here)
         * 8, 9 - pressure-, pressure+;
         * 10, 11 - trace-strain-, tracestrain+
         * 12, 13 - lambda1, lambda2
         * 14 - pressure_fault
         * 15 - theta (R & S state, not used here)
         */
        int I_p = 8;
        Jf3[(spaceDim * I_pf) * nCols + (spaceDim * I_p)] = a[i_fluidMobilityX] / 4. / a[i_fluidViscosity];
        Jf3[(spaceDim * I_pf + 1) * nCols + (spaceDim * I_p + 1)] = a[i_fluidMobilityX] / 4. / a[i_fluidViscosity];
        Jf3[(spaceDim * I_pf) * nCols + (spaceDim * (I_p + 1))] = a[i_fluidMobilityX] / 4. / a[i_fluidViscosity];
        Jf3[(spaceDim * I_pf) * nCols + (spaceDim * (I_p + 1) + 1)] = a[i_fluidMobilityX] / 4. / a[i_fluidViscosity];
    }
    // DEBUG LINES
    // cout << "Jf3pp = " << Jf3[(I_p) * nCols + (I_p)] << " " << Jf3[(I_p + 1) * nCols + (I_p + 1)] << "\n";
};

/** The elements of Jf3 that requires re-assemble after the first iteration
 * the non-zeros are only Jf3pfpf (2 terms), Jf3pfp (4 terms)
 * nothing is time dependent
 */
// const vector<int> PrescribeFaultKernel::Jf3_entries = {0, 39, 13, 14, 25, 26, 104, 117};
const vector<int> PrescribeFaultKernel::Jf3_is = {28, 29, 28, 29, 28, 29};
const vector<int> PrescribeFaultKernel::Jf3_js = {28, 29, 16, 17, 18, 19};
// const vector<int> PrescribeFaultKernel::Jf3_timedependent = {}; 
const vector<int> PrescribeFaultKernel::Jf3_timedependent_is = {};
const vector<int> PrescribeFaultKernel::Jf3_timedependent_js = {};