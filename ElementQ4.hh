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
using namespace std;

/** Class ElementQ4
 * used for subzone elements
 */
class ElementQ4 {
// PRIVATE MEMBERS
private:
    // Element ID, should start from 0
    int _ID;

    // Element-node connection
    vector<Node*> _NID;
    
    // Shape function N, vector of 4 on 4 nodes, evaluated in base space (ksi, eta)
    static vector<double> N(double ksi, double eta);

    // Gradient of shape function B, vector of 4 * 2 on 4 nodes, 2 directions, evaluated in base space (ksi, eta)
    static vector<double> B(double ksi, double eta);

    // Jacobian at any given location in base space (ksi, eta), J = det(\partial (x,y) / \partial (ksi, eta))
    double J(double ksi, double eta) const;

// PUBLIC MEMBERS
public:
    // Default Constructor
    ElementQ4();

    // Constructor
    ElementQ4(int ID, const vector<Node*> & NID);

    // Destructor
    ~ElementQ4();
    
    // Set element ID
    void setID(int ID);

    // Get element ID
    int getID() const;
    
    // Set element NID
    void setNID(const vector<Node*> & NID);

    // Get element NID
    const vector<Node*> & getNID() const;

    // Output element info
    void outputInfo(ofstream & myFile) const;

// NOT IMPLEMENTED
private:
    // Copy constructor
    ElementQ4(const ElementQ4 &);

    // Move-assign operator
    const ElementQ4 & operator=(const ElementQ4 &);
};


//------------------------------------------------------------------------------------------------------------------------------------------------------
/** ElementQ4Cohesive
 * Cohesive version of Q4, only 2 nodes
 */
class ElementQ4Cohesive : public ElementQ4 {
// PRIVATE MEMBERS
private:
    vector<CohesiveNode*> _NID;

// PUBLIC MEMBERS
public:
    // Constructor
    ElementQ4Cohesive(int ID, const vector<CohesiveNode*> & NID);
    
    // Destructor
    ~ElementQ4Cohesive();

    // Set element NID
    void setNID(const vector<CohesiveNode*> & NID);

    // Get element NID
    const vector<CohesiveNode*> & getNID() const;

    // Output element info
    void outputInfo(ofstream & myFile) const;

// NOT IMPLEMENTED
private:
    // Copy constructor
    ElementQ4Cohesive(const ElementQ4Cohesive & );

    // Move-assign operator
    const ElementQ4Cohesive & operator=(const ElementQ4Cohesive & ); 
};