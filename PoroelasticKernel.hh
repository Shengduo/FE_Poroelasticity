/** @file PoroelasticKernel.hh
 * Declares the quasi-static linear poroelastic residuals and Jacobians
 */
#include "FeKernel.hh"
using namespace std;

/** Class ElasticKernel
 * Quasi-static, linear elastic
 * only needs to redefine residual F and jacobian Jf
 */
class PoroelasticKernel {
// PUBLIC METHODS
public:
    //======= The Residuals =================================================================================   
    /** Left hand side residual 
     * F0(t, \dot{s}, s)
     */
    static void F0(vector<double> &F0,         // stores the result
                   int spaceDim,               // stores the dim of space
                   const vector<double> &s,    // solution vector s
                   const vector<double> &s_x,  // gradient of s
                   const vector<double> &s_t,  // time derivative of s
                   double s_tshift,            // sigma of tshift due to the time-derivative
                   const vector<double> &sOff, // offset of each solution field
                   const vector<double> &a,    // auxiliary fields
                   const vector<double> &aOff  // auxiliary fields offset
    );

    /** Left hand side residual 
     * F1(t, \dot{s}, s)
     */
    static void F1(vector<double> &F1,         // stores the result
                   int spaceDim,               // stores the dim of space
                   const vector<double> &s,    // solution vector s
                   const vector<double> &s_x,  // gradient of s
                   const vector<double> &s_t,  // time derivative of s
                   double s_tshift,            // sigma of tshift due to the time-derivative
                   const vector<double> &sOff, // offset of each solution field
                   const vector<double> &a,    // auxiliary fields
                   const vector<double> &aOff  // auxiliary fields offset
    );

    //======= The Jacobians =================================================================================   
    /** Left hand side Jacobian
     * Jf0(t, s)
     */
    static void Jf0(vector<double> &Jf0,        // stores the result
                    int spaceDim,               // stores the dim of space
                    const vector<double> &s,    // solution vector s
                    const vector<double> &s_x,  // gradient of s
                    const vector<double> &s_t,  // time derivative of s
                    double s_tshift,            // sigma of tshift due to the time-derivative
                    const vector<double> &sOff, // offset of each solution field
                    const vector<double> &a,    // auxiliary fields
                    const vector<double> &aOff  // auxiliary fields offset
    );

    /** Left hand side Jacobian
     * Jf1(t, s)
     */
    static void Jf1(vector<double> &Jf1,        // stores the result
                    int spaceDim,               // stores the dim of space
                    const vector<double> &s,    // solution vector s
                    const vector<double> &s_x,  // gradient of s
                    const vector<double> &s_t,  // time derivative of s
                    double s_tshift,            // sigma of tshift due to the time-derivative
                    const vector<double> &sOff, // offset of each solution field
                    const vector<double> &a,    // auxiliary fields
                    const vector<double> &aOff  // auxiliary fields offset
    );

    /** Left hand side Jacobian
     * Jf2(t, s)
     */
    static void Jf2(vector<double> &Jf2,        // stores the result
                    int spaceDim,               // stores the dim of space
                    const vector<double> &s,    // solution vector s
                    const vector<double> &s_x,  // gradient of s
                    const vector<double> &s_t,  // time derivative of s
                    double s_tshift,            // sigma of tshift due to the time-derivative
                    const vector<double> &sOff, // offset of each solution field
                    const vector<double> &a,    // auxiliary fields
                    const vector<double> &aOff  // auxiliary fields offset
    );

    /** Left hand side Jacobian
     * Jf3(t, s)
     */
    static void Jf3(vector<double> &Jf3,        // stores the result
                    int spaceDim,               // stores the dim of space
                    const vector<double> &s,    // solution vector s
                    const vector<double> &s_x,  // gradient of s
                    const vector<double> &s_t,  // time derivative of s
                    double s_tshift,            // sigma of tshift due to the time-derivative
                    const vector<double> &sOff, // offset of each solution field
                    const vector<double> &a,    // auxiliary fields
                    const vector<double> &aOff  // auxiliary fields offset
    );

    /** Default constructor */
    PoroelasticKernel() {};

    /** Default Virtual destructor */
    virtual ~PoroelasticKernel() {};

// NOT IMPLEMENTED
private:
    // Copy constructor
    PoroelasticKernel(const PoroelasticKernel &); 
    
    // Move-assign operator
    const PoroelasticKernel & operator=(const PoroelasticKernel &);
}; 