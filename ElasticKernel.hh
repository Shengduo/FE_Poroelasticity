/** @file ElasticKernel.hh
 * Declares the quasi-static linear elastic residuals and Jacobians
 */
#include "FeKernel.hh"
using namespace std;

/** Class ElasticKernel
 * Quasi-static, linear elastic
 * only needs to redefine residual F and jacobian Jf
 */
class ElasticKernel : public FeKernel {
// PUBLIC METHODS
// PUBLIC METHODS
public:
//======= The Residuals =================================================================================   
    /** Left hand side residual 
     * F0(t, \dot{s}, s)
     */
    void F0(vector<double> &F0,         // stores the result
            int spaceDim,               // stores the dim of space
            const vector<double> &s,    // solution vector s
            const vector<double> &s_x,  // gradient of s
            const vector<double> &s_t,  // time derivative of s
            double s_tshift,            // sigma of tshift due to the time-derivative
            const vector<double> &sOff, // offset of each solution field
            const vector<double> &a,    // auxiliary fields
            const vector<double> &aOff  // auxiliary fields offset
    ) const;

    /** Left hand side residual 
     * F1(t, \dot{s}, s)
     */
    void F1(vector<double> &F1,         // stores the result
            int spaceDim,               // stores the dim of space
            const vector<double> &s,    // solution vector s
            const vector<double> &s_x,  // gradient of s
            const vector<double> &s_t,  // time derivative of s
            double s_tshift,            // sigma of tshift due to the time-derivative
            const vector<double> &sOff, // offset of each solution field
            const vector<double> &a,    // auxiliary fields
            const vector<double> &aOff  // auxiliary fields offset
    ) const;

    //======= The Jacobians =================================================================================   
    /** Left hand side Jacobian
     * Jf0(t, s)
     */
    void Jf0(vector<double> &Jf0,        // stores the result
             int spaceDim,               // stores the dim of space
             const vector<double> &s,    // solution vector s
             const vector<double> &s_x,  // gradient of s
             const vector<double> &s_t,  // time derivative of s
             double s_tshift,            // sigma of tshift due to the time-derivative
             const vector<double> &sOff, // offset of each solution field
             const vector<double> &a,    // auxiliary fields
             const vector<double> &aOff  // auxiliary fields offset
    ) const;

    /** Left hand side Jacobian
     * Jf1(t, s)
     */
    void Jf1(vector<double> &Jf1,        // stores the result
             int spaceDim,               // stores the dim of space
             const vector<double> &s,    // solution vector s
             const vector<double> &s_x,  // gradient of s
             const vector<double> &s_t,  // time derivative of s
             double s_tshift,            // sigma of tshift due to the time-derivative
             const vector<double> &sOff, // offset of each solution field
             const vector<double> &a,    // auxiliary fields
             const vector<double> &aOff  // auxiliary fields offset
    ) const;

    /** Left hand side Jacobian
     * Jf2(t, s)
     */
    void Jf2(vector<double> &Jf2,        // stores the result
             int spaceDim,               // stores the dim of space
             const vector<double> &s,    // solution vector s
             const vector<double> &s_x,  // gradient of s
             const vector<double> &s_t,  // time derivative of s
             double s_tshift,            // sigma of tshift due to the time-derivative
             const vector<double> &sOff, // offset of each solution field
             const vector<double> &a,    // auxiliary fields
             const vector<double> &aOff  // auxiliary fields offset
    ) const;

    /** Left hand side Jacobian
     * Jf3(t, s)
     */
    void Jf3(vector<double> &Jf3,        // stores the result
             int spaceDim,               // stores the dim of space
             const vector<double> &s,    // solution vector s
             const vector<double> &s_x,  // gradient of s
             const vector<double> &s_t,  // time derivative of s
             double s_tshift,            // sigma of tshift due to the time-derivative
             const vector<double> &sOff, // offset of each solution field
             const vector<double> &a,    // auxiliary fields
             const vector<double> &aOff  // auxiliary fields offset
    ) const;

    /** Default constructor */
    ElasticKernel() {};

    /** Default Virtual destructor */
    virtual ~ElasticKernel() {};

// NOT IMPLEMENTED
private:
    // Copy constructor
    ElasticKernel(const ElasticKernel &); 
    
    // Move-assign operator
    const ElasticKernel & operator=(const ElasticKernel &);
}; 