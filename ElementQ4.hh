/** @file ElementQ4.hh
 * Head file for class ElementQ4 and corresponding cohesive class ElementQ4Cohesive
 * for 2D upper/lower subzone bodies.
 */
#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include "Node.hh"
#include "ElasticKernel.hh"
#include "PoroelasticKernel.hh"

using namespace std;

/** Class ElementQ4
 * used for subzone elements
 */
class ElementQ4 {
// PRIVATE MEMBERS
private:
    /** Element ID, should start from 0 */
    int _ID;

    /** Element-node connection */
    vector<Node*> _NID;
    
    /** Faces that requires surface-load computation. */
    vector<vector<int>> _loadFaces;

    /** Shape function N, vector of 4 on 4 nodes, 
     * evaluated in base space (ksi, eta)
     */
    static vector<double> N(double ksi, double eta);

    /** Shape function N_surf, vector of 2 on 2 integral points, 
     * evaluated in base space (ksi)
     */
    static vector<double> N_surf(double ksi);

    /** Shape function N, vector of 4 on 4 nodes, 
     * evaluated in base space (ksi, eta)
     * at i, j of elemental shape matrix
     */
    double N(const vector<double> & Nvector, int i, int j) const;

    /** Gradient of shape function B, vector of 4 * 2 on 4 nodes, 2 directions, 
     * evaluated in base space (ksi, eta)
     */
    static vector<double> B(double ksi, double eta);

    /** Gradient of shape function B_x[i, j], double value, 
     * evaluated in base space (ksi, eta)
     */
    double B_x(const vector<double> & Bvector, int i, int j) const;

    /** Gradient_x of shape function B at (ksi, eta)
     * Non-static
     */
    vector<double> B_x(double ksi, double eta) const;

    /** Pre-calculate all B_x vectors at 4 integration points
     * since everything is calculated in reference config
     */
    vector<vector<double>> Bvector;

    /** Pre-calculate all N vectors at 4 integration points
     * since everything is calculated in reference config
     */
    vector<vector<double>> Nvector;

    /** Element DOF number
     * used to pre-allocate other matrices / vectors
     */
    PetscInt elementDOF;

    /** Pre-calculate the local->global indices, 
     * can call when assemblying global Jacobian and vector.
     */
    PetscInt *localGlobalIndices;

    /** Pre-calculate pointValue = w_i w_j J(i, j)
     * Used in integrators
     * Since element reference state does not change,
     * this can be pre-stored.
     */
    vector<double> pointValue;

    /** Jacobian at any given location in base space (ksi, eta), 
     * J = det(\partial (x,y) / \partial (ksi, eta))
     */
    double J(double ksi, double eta) const;

    /** calculate inverse of Jacobian at any given location in base space (ksi, eta), 
     * invJ = (\partial (x,y) / \partial (ksi, eta)) ^ (-1)
     */
    bool InvJ(vector<double> & res, double ksi, double eta) const;

// ================== Constants and structures used for JF and elementF evaluation =======
    /** Number of nodes */
    int nOfNodes;

    /** Space dimension */
    int spaceDim;

    /** Number of degree of freedom */
    int nOfDofs;

    /** Number of integration points */
    int nOfIntPts;

    /** For evaluation at each integration point, pre-allocate */
    vector<double> ss;
    vector<double> s_xs;
    vector<double> s_ts;
    vector<double> as;

    /** For collecting nodal data, pre-allocate */
    vector<vector<double> *> nodalSs;
    vector<vector<double> *> nodalS_ts;
    vector<vector<double> *> nodalAs;

    /** Pre-allocate for jacobians */
    PetscBool isJfAssembled;
    vector<vector<double>> Jf0s;
    vector<vector<double>> Jf1s;
    vector<vector<double>> Jf2s;
    vector<vector<double>> Jf3s;

    /** Pre-allocate for residuals */
    vector<vector<double>> F0s;
    vector<vector<double>> F1s;

    /** Pointer to the entire problem timer so that can call computing time */
    vector<clock_t> *clocks;
    vector<double> *timeConsumed;

// SHARED WITH COHESIVE CLASS
protected:
    /** Constants for 2-point gaussian integral */
    const static vector<double> IntPos;
    const static vector<double> IntWs;
    const static vector<double> N_surfVector;
    const static vector<int> scatter_pattern_traction;
// PUBLIC MEMBERS
public:
    /** Default Constructor */
    ElementQ4();

    /** Constructor */
    ElementQ4(int ID, 
              const vector<Node*> & NID, 
              const vector<vector<int>> & loadFaces = vector<vector<int>>(), 
              vector<clock_t> *clocks = NULL, 
              vector<double> *timeConsumed = NULL);

    /** Destructor */
    ~ElementQ4();
    
    
    /** Set element ID */
    void setID(int ID);

    /** Get element ID */
    int getID() const;
    
    /** Get elementDOF */
    PetscInt getElementDOF() const;

    /** Set element NID */
    void setNID(const vector<Node*> & NID);

    /** Get element NID */
    const vector<Node*> & getNID() const;

    /** Get face traction */
    vector<vector<double>> getFaceTraction(const vector<int>& face) const;

    /** Evaluate a function F at (i, j)th integration point */
    void evaluateF(vector<double> & res, int i, int j, 
                   const vector<vector<double> *> & NodeValues) const;

    /** Evaluate vector at (ksi, eta) in logical space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF(vector<double> & res, double ksi, double eta, 
                   const vector<vector<double> *> & NodeValues) const;
    
    /** Evaluate vector at (ksi, eta) in logical space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF(vector<double> & res, double ksi, double eta, 
                   const vector<vector<double>> & NodeValues) const;
    
    /** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at integration
     * point i, j, in LOGICAL space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF_x(vector<double> & res, int i, int j, 
                     const vector<vector<double> *> & NodeValues) const;

    /** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at (ksi, eta) 
     * in LOGICAL space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF_x(vector<double> & res, double ksi, double eta, 
                     const vector<vector<double> *> & NodeValues) const;   

    /** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at (ksi, eta) 
     * in LOGICAL space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF_x(vector<double> & res, double ksi, double eta, 
                     const vector<vector<double>> & NodeValues) const;  
//================ Integrators for ElementQ4 ==================================
    /** IntegratorNf, integrates a vector input inside an element, 
     * left side using shape function
     * RES: nOfNodes * nOfDofs
     * NODEVALUES: dim 1, nOfNodes; dim 2, nOfDofs
     * FLAG: 0 - nodevalues are given at nodes, 
     *       1 - nodevalues are given at integration points
     */
    void IntegratorNf(double *res,
                      int resSize, 
                      const vector<vector<double>> & NodeValues, 
                      int flag = 0) const;
    
    /** IntegratorNfFace, integrates a vector input over a face of the element, 
     * left side using shape function
     * RES: nOfNodes * nOfDofs
     * NODEVALUES: dim 1, nOfNodes; dim 2, nOfDofs
     * FACE: vector of the 2 nodes that define a face (order does not matter in 2D)
     * FLAG: 0 - nodevalues are given at nodes, 
     *       1 - nodevalues are given at integration points
     */
    void IntegratorNfFace(double *res,
                      int resSize, 
                      const vector<vector<double>> & NodeValues, 
                      const vector<int> & scatter_pattern, 
                      const vector<int> & Face, 
                      int flag = 0) const;

    /** IntegratorBf, integrates a vector input inside an element, 
     * left side using shape function
     * RES: nOfNodes * nOfDofs
     * NODEVALUES: dim 1, nOfNodes; dim 2, nOfDofs * spaceDim
     * FLAG: 0 - nodevalues are given at nodes, 
     *       1 - nodevalues are given at integration points
     */
    void IntegratorBf(double *res,
                      int resSize,
                      const vector<vector<double>> & NodeValues, 
                      int flag = 0) const;

    /** IntegratorNfN, integrates a vector input inside an element, 
     * both sides using shape function (nOfDofs, nOfDofs * nOfNodes)
     * first-dim: vector of nodes, second-dim: values (vector, (nOfDofs, nOfDofs))
     * FLAG: 0 - nodevalues are given at nodes, 
     *       1 - nodevalues are given at integration points
     */
    void IntegratorNfN(double *res,
                       int resSize,
                       const vector<vector<double>> & NodeValues,                         
                       const vector<int> & f_is, 
                       const vector<int> & f_js, 
                       int flag = 0) const;

    /** IntegratorBfB, integrates a vector input inside an element, 
     * both sides using gradient of shape function
     * RES:
     * first-dim: vector of nodes^2, 
     * NODEVALUES:
     * first-dim vector of nodes, 
     * second-dim (spaceDim * nDof) ^ 2 matrix.
     * FLAG: 0 - nodevalues are given at nodes
     *       1 - nodevalues are given at integration points
     */
    void IntegratorBfB(double *res,
                       int resSize, 
                       const vector<vector<double>> & NodeValues, 
                       const vector<int> & f_is, 
                       const vector<int> & f_js, 
                       int flag = 0) const;
    
    /** IntegratorBfN, integrates a vector input inside an element, 
     * left gradient of shape function, 
     * right shape function
     * RES: first-dim: vector of nOfDof^2
     * NODEVALUES: first dim: spaceDim * nOfDof, 
     * second dim: nOfDof
     * FLAG: 0 - nodevalues are given at nodes
     *       1 - nodevalues are given at integration points
     */
    void IntegratorBfN(double *res,
                       int resSize, 
                       const vector<vector<double>> & NodeValues, 
                       const vector<int> & f_is, 
                       const vector<int> & f_js, 
                       int flag = 0) const;

    /** IntegratorNfB, integrates a vector input inside an element, 
     * left side shape function, right side gradient of shape function, 
     * RES: first-dim: vector of nodes^2, second-dim: values (vector)
     * NODEVALUES: first dim: nodes, second dim: 1 by spaceDim matrix, stored as a vector)
     * FLAG: 0 - nodevalues are given at nodes, 
     *       1 - nodevalues are given at integration points
     */
    void IntegratorNfB(double *res,
                       int resSize, 
                       const vector<vector<double>> & NodeValues,                         
                       const vector<int> & f_is, 
                       const vector<int> & f_js, 
                       int flag = 0) const;

    /** Output element info */
    void outputInfo(ofstream & myFile) const;

//============ Element Jacobians and residuals =====================================================
// PUBLIC MEMBERS
public: 
    /** Calculate element jacobian JF */
    void JF(Mat & globalJF, double *localJF, int localJFSize, int Kernel, double s_tshift);

    /** Calculate element residual F */
    void elementF(Vec & globalF, double *localF, int localFSize, int Kernel, double s_tshift);

// PRIVATE MEMBERS
private:
    /** Push local Jf to global Jf */
    void JFPush(Mat & globalJF, double *elementJF, int elementJFSize) const;

    /** Push elementF to globalF */
    void elementFPush(Vec & globalF, double *elementF, int elementFSize) const;

// NOT IMPLEMENTED
private:
    /** Copy constructor */
    ElementQ4(const ElementQ4 &);

    /** Move-assign operator */
    const ElementQ4 & operator=(const ElementQ4 &);
};