/** @file ElementQ4.hh
 * Head file for class ElementQ4 and corresponding cohesive class ElementQ4Cohesive
 * for 2D upper/lower subzone bodies.
 */
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include "Node.hh"
#include "ElasticKernel.hh"
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
    
    /** Shape function N, vector of 4 on 4 nodes, 
     * evaluated in base space (ksi, eta)
     */
    static vector<double> N(double ksi, double eta);

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

    /** Jacobian at any given location in base space (ksi, eta), 
     * J = det(\partial (x,y) / \partial (ksi, eta))
     */
    double J(double ksi, double eta) const;

    /** calculate inverse of Jacobian at any given location in base space (ksi, eta), 
     * invJ = (\partial (x,y) / \partial (ksi, eta)) ^ (-1)
     */
    bool InvJ(vector<double> & res, double ksi, double eta) const;

// SHARED WITH COHESIVE CLASS
protected:
    /** Constants for 2-point gaussian integral */
    const static vector<double> IntPos;
    const static vector<double> IntWs;

// PUBLIC MEMBERS
public:
    /** Default Constructor */
    ElementQ4();

    /** Constructor */
    ElementQ4(int ID, const vector<Node*> & NID);

    /** Destructor */
    ~ElementQ4();
    
    /** Set element ID */
    void setID(int ID);

    /** Get element ID */
    int getID() const;
    
    /** Set element NID */
    void setNID(const vector<Node*> & NID);

    /** Get element NID */
    const vector<Node*> & getNID() const;

    /** Evaluate vector at (ksi, eta) in logical space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF(vector<double> & res, double ksi, double eta, 
                   const vector<vector<double>> & NodeValues) const;
    
    /** Evaluate PHYSICAL gradient (\partial x, \partial y) of vector at (ksi, eta) 
     * in LOGICAL space with given nodal values.
     * Calculated by using shape function to map
     */
    void evaluateF_x(vector<double> & res, double ksi, double eta, 
                   const vector<vector<double>> & NodeValues) const;   

    /** IntegratorNf, integrates a vector input inside an element, 
     * both sides using shape function
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
     * both sides using shape function
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

//============ Element Jacobians and residuals =====================================================
// PUBLIC MEMBERS
public: 
    /** Calculate element jacobian JF */
    void JF(Mat & globalJF, int Kernel) const;

// PRIVATE MEMBERS
private:
    /** Push local Jf to global Jf */
    void JFPush(Mat & globalJF, const vector<double> & JF) const;

// NOT IMPLEMENTED
private:
    /** Copy constructor */
    ElementQ4(const ElementQ4 &);

    /** Move-assign operator */
    const ElementQ4 & operator=(const ElementQ4 &);
};