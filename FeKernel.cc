/** @file FeKernel.cc,
 * Default implementation of base class FeKernel
 */
#include "FeKernel.hh"

/** Left hand side residual 
 * F0(t, \dot{s}, s)
 */
void FeKernel::F0(vector<double> & F0,                 // stores the result                   
                const vector<double> & s,              // solution vector s
                const vector<double> & s_x,            // gradient of s
                const vector<double> & s_t,            // time derivative of s
                double s_tshift,                       // sigma of tshift due to the time-derivative
                const vector<double> & sOff,           // offset of each solution field
                const vector<double> & a,              // auxiliary fields
                const vector<double> & aOff            // auxiliary fields offset
                ) const {
    
    // Default residuals are all 0.
    for (int i = 0; i < s.size(); i++) F0[i] = 0.;
};

/** Left hand side residual 
 * F1(t, \dot{s}, s)
 */
void FeKernel::F1(vector<double> & F1,                 // stores the result                   
                const vector<double> & s,              // solution vector s
                const vector<double> & s_x,            // gradient of s
                const vector<double> & s_t,            // time derivative of s
                double s_tshift,                       // sigma of tshift due to the time-derivative
                const vector<double> & sOff,           // offset of each solution field
                const vector<double> & a,              // auxiliary fields
                const vector<double> & aOff            // auxiliary fields offset
                ) const {
    
    // Default residuals are all 0.
    for (int i = 0; i < s.size(); i++) F1[i] = 0.;
};

/** Left hand side residual 
 * G0(t, s)
 */
void FeKernel::G0(vector<double> & G0,                  // stores the result                   
                const vector<double> & s,              // solution vector s
                const vector<double> & s_x,            // gradient of s
                const vector<double> & s_t,            // time derivative of s
                double s_tshift,                       // sigma of tshift due to the time-derivative
                const vector<double> & sOff,           // offset of each solution field
                const vector<double> & a,              // auxiliary fields
                const vector<double> & aOff            // auxiliary fields offset
                ) const {
    
    // Default residuals are all 0.
    for (int i = 0; i < s.size(); i++) G0[i] = 0.;
};

/** Left hand side residual 
 * G1(t, s)
 */
void FeKernel::G1(vector<double> & G1,                 // stores the result                   
                const vector<double> & s,              // solution vector s
                const vector<double> & s_x,            // gradient of s
                const vector<double> & s_t,            // time derivative of s'
                double s_tshift,                       // sigma of tshift due to the time-derivative
                const vector<double> & sOff,           // offset of each solution field
                const vector<double> & a,              // auxiliary fields
                const vector<double> & aOff            // auxiliary fields offset
                ) const {
    
    // Default residuals are all 0.
    for (int i = 0; i < s.size(); i++) G1[i] = 0.;
};
 
//======= The Jacobians =================================================================================   
/** Left hand side Jacobian
 * Jf0(t, s)
 */
void FeKernel::Jf0(vector<double> & Jf0,               // stores the result                   
                const vector<double> & s,              // solution vector s
                const vector<double> & s_x,            // gradient of s
                const vector<double> & s_t,            // time derivative of s
                double s_tshift,                       // sigma of tshift due to the time-derivative
                const vector<double> & sOff,           // offset of each solution field
                const vector<double> & a,              // auxiliary fields
                const vector<double> & aOff            // auxiliary fields offset
                ) const {

    // Default Jacobians are all 0.
    for (int i = 0; i < s.size(); i++) Jf0[i] = 0.;
};

/** Left hand side Jacobian
 * Jf1(t, s)
 */
void FeKernel::Jf1(vector<double> & Jf1,               // stores the result                   
                const vector<double> & s,              // solution vector s
                const vector<double> & s_x,            // gradient of s
                const vector<double> & s_t,            // time derivative of s
                double s_tshift,                       // sigma of tshift due to the time-derivative
                const vector<double> & sOff,           // offset of each solution field
                const vector<double> & a,              // auxiliary fields
                const vector<double> & aOff            // auxiliary fields offset
                ) const {

    // Default Jacobians are all 0.
    for (int i = 0; i < s.size(); i++) Jf1[i] = 0.;
};

/** Left hand side Jacobian
 * Jf2(t, s)
 */
void FeKernel::Jf2(vector<double> & Jf2,               // stores the result                   
                const vector<double> & s,              // solution vector s
                const vector<double> & s_x,            // gradient of s
                const vector<double> & s_t,            // time derivative of s
                double s_tshift,                       // sigma of tshift due to the time-derivative
                const vector<double> & sOff,           // offset of each solution field
                const vector<double> & a,              // auxiliary fields
                const vector<double> & aOff            // auxiliary fields offset
                ) const {

    // Default Jacobians are all 0.
    for (int i = 0; i < s.size(); i++) Jf2[i] = 0.;
};

/** Left hand side Jacobian
 * Jf3(t, s)
 */
void FeKernel::Jf3(vector<double> & Jf3,               // stores the result                   
                const vector<double> & s,              // solution vector s
                const vector<double> & s_x,            // gradient of s
                const vector<double> & s_t,            // time derivative of s
                double s_tshift,                       // sigma of tshift due to the time-derivative
                const vector<double> & sOff,           // offset of each solution field
                const vector<double> & a,              // auxiliary fields
                const vector<double> & aOff            // auxiliary fields offset
                ) const {

    // Default Jacobians are all 0.
    for (int i = 0; i < s.size(); i++) Jf3[i] = 0.;
};