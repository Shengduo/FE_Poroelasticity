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
public:
    //======= The Residuals =================================================================================   
 
    //======= The Jacobians =================================================================================   
    /** Left hand side Jacobian
     * Jf3(t, s)
     */
    static void Jf3(vector<double> &Jf3,    // stores the result
                    int spaceDim,           // stores the dim of space
                    const vector<double> &a // auxiliary fields
    );

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