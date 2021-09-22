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

    /** Gradient of shape function B, vector of 2 * 1 on 2 nodes, 
     * evaluated in base space (ksi)
     */
    static vector<double> B(double ksi);

    /** Jacobian at any given location in base space (ksi, eta), J = dx / d ksi */
    double J(double ksi) const;



// PUBLIC MEMBERS
public:
    /** Constructor */
    ElementQ4Cohesive(int ID, const vector<CohesiveNode*> & NID);
    
    /** Destructor */
    ~ElementQ4Cohesive();

    /** Set element NID */
    void setNID(const vector<CohesiveNode*> & NID);

    /** Get element NID */
    const vector<CohesiveNode*> & getNID() const;

    /** Evaluate a function at ksi */
    void evaluateF(vector<double> & res, double ksi, 
                   const vector<vector<double>> & NodeValues) const;
    
    /** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at (ksi) 
     * in LOGICAL space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF_x(vector<double> & res, double ksi,
                   const vector<vector<double>> & NodeValues) const;

    /** Integrator, integrates a vector input inside an element, 
     * first-dim: vector of nodes, second-dim: values (vector)
     */
    void IntegratorNf(vector<vector<double>> & res, 
                      const vector<vector<double>> & NodeValues) const;

    /** IntegratorNfN, integrates a vector input inside an element, 
     * both sides using shape function
     * first-dim: vector of nodes^2, second-dim: values (vector)
     */
    void IntegratorNfN(vector<vector<double>> & res,
                       const vector<vector<double>> & NodeValues) const;

    /** IntegratorBfB, integrates a vector input inside an element, 
     * both sides using gradient of shape function
     * first-dim: vector of nodes^2, second-dim: values (vector)
     */
    void IntegratorBfB(vector<double> & res,
                       const vector<vector<double>> & NodeValues) const;

    /** IntegratorBfN, integrates a vector input inside an element, 
     * left side gradient of shape function, right side shape function
     * RES, first-dim: vector of nodes^2, second-dim: values (vector)
     * NODEVALUES, first dim: nodes, second dim: spaceDim by 1 matrix, stored as a vector
     */
    void IntegratorBfN(vector<double> & res,
                       const vector<vector<double>> & NodeValues) const;

    /** IntegratorNfB, integrates a vector input inside an element, 
     * left side shape function, right side gradient of shape function, 
     * RES, first-dim: vector of nodes^2, second-dim: values (vector)
     * NODEVALUES, first dim: nodes, second dim: 1 by spaceDim matrix, stored as a vector)
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