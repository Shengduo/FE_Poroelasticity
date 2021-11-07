/** @file ElementQ4Cohesive.hh
 * Headfile for class ElementQ4Cohesive : ElementQ4
 */
#include "ElementQ4.hh"
//------------------------------------------------------------------------------------------------------------------------------------------------------
/** ElementQ4Cohesive
 * Cohesive version of Q4, only 2 nodes
 */

class ElementQ4Cohesive : public ElementQ4 {
// PRIVATE MEMBERS
private:
    /** Node connection */
    vector<CohesiveNode*> _NID;

    /** Shape function N, vector of 2 on 2 nodes, 
     * evaluated in base space (ksi)
     */
    static vector<double> N(double ksi);

    /** Shape function N, vector of 4 on 4 nodes, 
     * evaluated in base space (ksi)
     * at i, j of elemental shape matrix
     */
    double N(const vector<double> & Nvector, int i, int j) const;

    /** Gradient of shape function B, vector of 2 * 1 on 2 nodes, 
     * evaluated in base space (ksi)
     */
    static vector<double> B(double ksi);
    
    /** Gradient of shape function B, vector of 2 * 2 on 2 nodes, 2 spaceDims, 
     * evaluated in physical space 
     */
    vector<double> B_x(double ksi) const;

    /** Gradient of shape function B_x[i, j] at (ksi) */
    double B_x(const vector<double> & Bvector, int i, int j) const;

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

    /** Jacobian at any given location in base space (ksi, eta), J = dx / d ksi */
    double J(double ksi) const;

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

// PUBLIC MEMBERS
public:
    /** Constructor */
    ElementQ4Cohesive(int ID, const vector<CohesiveNode*> & NID, vector<clock_t> *clocks = NULL, vector<double> *timeConsumed = NULL);
    
    /** Destructor */
    ~ElementQ4Cohesive();

    /** Set element NID */
    void setNID(const vector<CohesiveNode*> & NID);

    /** Get element NID */
    const vector<CohesiveNode*> & getNID() const;

    /** Evaluate a function F at i th integration point */
    void evaluateF(vector<double> & res, int i,
                   const vector<vector<double> *> & NodeValues) const;

    /** Evaluate vector at (ksi, eta) in logical space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF(vector<double> & res, double ksi,
                   const vector<vector<double> *> & NodeValues) const;

    /** Evaluate a function at ksi */
    void evaluateF(vector<double> & res, double ksi, 
                   const vector<vector<double>> & NodeValues) const;
    
    /** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at integration
     * point i, j, in LOGICAL space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF_x(vector<double> & res, int i, 
                     const vector<vector<double> *> & NodeValues) const;

    /** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at (ksi, eta) 
     * in LOGICAL space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF_x(vector<double> & res, double ksi, 
                     const vector<vector<double> *> & NodeValues) const;   

    /** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at (ksi) 
     * in LOGICAL space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF_x(vector<double> & res, double ksi,
                   const vector<vector<double>> & NodeValues) const;

//================ Integrators for ElementQ4Cohesive ==================================
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
     */
    void IntegratorNfN(vector<double> & res,
                       const vector<vector<double>> & NodeValues) const;

    /*** IntegratorBfB, integrates a vector input inside an element, 
     * both sides using gradient of shape function
     * RES:
     * first-dim: vector of nodes^2, 
     * NODEVALUES:
     * first-dim vector of nodes, 
     * second-dim (spaceDim * nDof) ^ 2 matrix.
     */
    void IntegratorBfB(vector<double> & res,
                       const vector<vector<double>> & NodeValues) const;

    /** IntegratorBfN, integrates a vector input inside an element, 
     * left gradient of shape function, 
     * right shape function
     * RES: first-dim: vector of nOfDof^2
     * NODEVALUES: first dim: spaceDim * nOfDof, 
     * second dim: nOfDof
     */
    void IntegratorBfN(vector<double> & res,
                       const vector<vector<double>> & NodeValues) const;

    /** IntegratorNfB, integrates a vector input inside an element, 
     * left gradient of shape function, 
     * right shape function
     * RES: first-dim: vector of nOfDof^2
     * NODEVALUES: first dim: nOfDof, 
     * second dim: spaceDim * nOfDof
     */
    void IntegratorNfB(vector<double> & res,
                       const vector<vector<double>> & NodeValues) const;

    /** Output element info */
    void outputInfo(ofstream & myFile) const;

// NOT IMPLEMENTED
private:
    /** Copy constructor */
    ElementQ4Cohesive(const ElementQ4Cohesive & );

    /** Move-assign operator */
    const ElementQ4Cohesive & operator=(const ElementQ4Cohesive & ); 
};