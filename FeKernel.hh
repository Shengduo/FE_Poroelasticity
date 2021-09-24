/** @file FeKernel.hh
 * Define a virtual base class FeKernel
 * Allow detailed implementation of jacobian, residual functions in derived classes
 */
using namespace std;

/** Class FeKernel,
 * virtual parent class for the integration kernels of bulk elements.
 * Formulation of equations: F(t, \dot{s}, s) = G(t, s)
 * !! Only usable for implicit and IMEX scheme FOR NOW
 */
class FeKernel {
// PUBLIC METHODS
public:
//======= The Jacobians =================================================================================   
    /** Left hand side residual 
     * F0(t, \dot{s}, s)
     */
    virtual void F0(vector<double> & F0,                  // stores the result                   
                   const vector<double> & s,              // solution vector s
                   const vector<double> & s_x,            // gradient of s
                   const vector<double> & s_t,            // time derivative of s
                   const vector<double> & sOff,           // offset of each solution field
                   const vector<double> & a,              // auxiliary fields
                   const vector<double> & aOff            // auxiliary fields offset
                   ) const;
    
    /** Left hand side residual 
     * F1(t, \dot{s}, s)
     */
    virtual void F1(vector<double> & F1,                  // stores the result                   
                   const vector<double> & s,              // solution vector s
                   const vector<double> & s_x,            // gradient of s
                   const vector<double> & s_t,            // time derivative of s
                   const vector<double> & sOff,           // offset of each solution field
                   const vector<double> & a,              // auxiliary fields
                   const vector<double> & aOff            // auxiliary fields offset
                   ) const;
    
    /** Left hand side residual 
     * G0(t, s)
     */
    virtual void G0(vector<double> & G,                   // stores the result                   
                   const vector<double> & s,              // solution vector s
                   const vector<double> & s_x,            // gradient of s
                   const vector<double> & s_t,            // time derivative of s
                   const vector<double> & sOff,           // offset of each solution field
                   const vector<double> & a,              // auxiliary fields
                   const vector<double> & aOff            // auxiliary fields offset
                   ) const;

    /** Left hand side residual 
     * G1(t, s)
     */
    virtual void G1(vector<double> & G,                   // stores the result                   
                   const vector<double> & s,              // solution vector s
                   const vector<double> & s_x,            // gradient of s
                   const vector<double> & s_t,            // time derivative of s
                   const vector<double> & sOff,           // offset of each solution field
                   const vector<double> & a,              // auxiliary fields
                   const vector<double> & aOff            // auxiliary fields offset
                   ) const;
 
//======= The Jacobians =================================================================================   
    /** Left hand side Jacobian
     * Jf0(t, s)
     */
    virtual void Jf0(vector<double> & Jf0,                // stores the result                   
                   const vector<double> & s,              // solution vector s
                   const vector<double> & s_x,            // gradient of s
                   const vector<double> & s_t,            // time derivative of s
                   const vector<double> & sOff,           // offset of each solution field
                   const vector<double> & a,              // auxiliary fields
                   const vector<double> & aOff            // auxiliary fields offset
                   ) const;
    
    /** Left hand side Jacobian
     * Jf1(t, s)
     */
    virtual void Jf1(vector<double> & Jf1,                // stores the result                   
                   const vector<double> & s,              // solution vector s
                   const vector<double> & s_x,            // gradient of s
                   const vector<double> & s_t,            // time derivative of s
                   const vector<double> & sOff,           // offset of each solution field
                   const vector<double> & a,              // auxiliary fields
                   const vector<double> & aOff            // auxiliary fields offset
                   ) const;
    
    /** Left hand side Jacobian
     * Jf2(t, s)
     */
    virtual void Jf2(vector<double> & Jf2,                // stores the result                   
                   const vector<double> & s,              // solution vector s
                   const vector<double> & s_x,            // gradient of s
                   const vector<double> & s_t,            // time derivative of s
                   const vector<double> & sOff,           // offset of each solution field
                   const vector<double> & a,              // auxiliary fields
                   const vector<double> & aOff            // auxiliary fields offset
                   ) const;

    /** Left hand side Jacobian
     * Jf3(t, s)
     */
    virtual void Jf3(vector<double> & Jf3,                // stores the result                   
                   const vector<double> & s,              // solution vector s
                   const vector<double> & s_x,            // gradient of s
                   const vector<double> & s_t,            // time derivative of s
                   const vector<double> & sOff,           // offset of each solution field
                   const vector<double> & a,              // auxiliary fields
                   const vector<double> & aOff            // auxiliary fields offset
                   ) const;

// NOT IMPLEMENTED
private:
    // Copy constructor
    FeKernel(const FeKernel &); 
    
    // Move-assign operator
    const FeKernel & operator=(const FeKernel &);
};