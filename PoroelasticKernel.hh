/** @file PoroelasticKernel.hh
 * Declares the quasi-static linear poroelastic residuals and Jacobians
 */
#pragma once
#include "FeKernel.hh"
#include "Node.hh"
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
                    const vector<double> &aOff, // auxiliary fields offset
                    PetscBool isAssembled = PETSC_FALSE    // if assembled, only calculate the time-dependent parts
    );

    /** The elements of Jf0 that requires re-assemble after the first iteration
     * non-zeros are Jf0pe, Jf0ee and Jf0pp
     * only Jf0pe, Jf0pp are time dependent
     */
    const static vector<int> Jf0_entries;
    const static vector<int> Jf0_is;
    const static vector<int> Jf0_js;
    const static vector<int> Jf0_timedependent; 
    const static vector<int> Jf0_timedependent_is;
    const static vector<int> Jf0_timedependent_js;

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
                    const vector<double> &aOff, // auxiliary fields offset
                    PetscBool isAssembled = PETSC_FALSE    // if assembled, only calculate the time-dependent parts 
    );

    /** The elements of Jf0 that requires re-assemble after the first iteration
     * non-zeros are Jf1eu
     * nothing in JF1 is time dependent
     */
    const static vector<int> Jf1_entries;
    const static vector<int> Jf1_is;
    const static vector<int> Jf1_js;
    const static vector<int> Jf1_timedependent;
    const static vector<int> Jf1_timedependent_is;
    const static vector<int> Jf1_timedependent_js;

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
                    const vector<double> &aOff, // auxiliary fields offset
                    PetscBool isAssembled = PETSC_FALSE    // if assembled, only calculate the time-dependent parts
    );

    /** The elements of Jf2 that requires re-assemble after the first iteration
     * the nonzeros are Jf2up (2 terms) and Jf2ue (2 terms)
     * nothing is dependent on time
     */
    const static vector<int> Jf2_entries;
    const static vector<int> Jf2_is;
    const static vector<int> Jf2_js;
    const static vector<int> Jf2_timedependent; 
    const static vector<int> Jf2_timedependent_is;
    const static vector<int> Jf2_timedependent_js;

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
                    const vector<double> &aOff, // auxiliary fields offset
                    PetscBool isAssembled = PETSC_FALSE    // if assembled, only calculate the time-dependent parts
    );

    /** The elements of Jf3 that requires re-assemble after the first iteration
     * the non-zeros are only Jf3uu (4 terms), Jf3pp
     * nothing is time dependent
     */
    const static vector<int> Jf3_entries;
    const static vector<int> Jf3_is;
    const static vector<int> Jf3_js;
    const static vector<int> Jf3_timedependent; 
    const static vector<int> Jf3_timedependent_is;
    const static vector<int> Jf3_timedependent_js;

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