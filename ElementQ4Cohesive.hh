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

    /** Minus and plus side node connection */
    vector<Node*> _NIDMinusPlus;

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

    /** Jacobian at any given location in base space (ksi, eta), J = dx / d ksi */
    double J(double ksi) const;



// PUBLIC MEMBERS
public:
    /** Constructor */
    ElementQ4Cohesive(int ID, const vector<CohesiveNode*> & NID, const vector<Node*> & NIDMinusPlus);
    
    /** Destructor */
    ~ElementQ4Cohesive();

    /** Set element NID */
    void setNID(const vector<CohesiveNode*> & NID);

    /** Get element NID */
    const vector<CohesiveNode*> & getNID() const;

    /** Set element NIDMinusPlus */
    void setNIDMinusPlus(const vector<Node*> & NIDMinusPlus);

    /** Get element NIDMinusPlus handle */
    const vector<Node*> & getNIDMinusPlus() const;

    /** Evaluate a function at ksi */
    void evaluateF(vector<double> & res, double ksi, 
                   const vector<vector<double>> & NodeValues) const;
    
    /** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at (ksi) 
     * in LOGICAL space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF_x(vector<double> & res, double ksi,
                   const vector<vector<double>> & NodeValues) const;

    /** IntegratorNf, integrates a vector input inside an element, 
     * left side using shape function
     * RES: nOfNodes * nOfDofs
     * NODEVALUES: dim 1, nOfNodes; dim 2, nOfDofs
     */
    void IntegratorNf(vector<double> & res, 
                      const vector<vector<double>> & NodeValues) const;
    
    /** IntegratorBf, integrates a vector input inside an element, 
     * left side using shape function
     * RES: nOfNodes * nOfDofs
     * NODEVALUES: dim 1, nOfNodes; dim 2, nOfDofs * spaceDim
     */
    void IntegratorBf(vector<double> & res, 
                      const vector<vector<double>> & NodeValues) const;

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