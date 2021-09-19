/** @file Problem.hh
 * Head file of class Problem
 * Initialize, solve and write the results.
 */
#include "Geometry2D.hh"
#include "ElementQ4.hh"

/** Class Problem
 * Define a problem and solves it
 */
class Problem {
// PRIVATE MEMBERS
private:
    // Space dimension
    int _spaceDim;

    // Geometry2D
    Geometry2D* myGeometry;

    // Upper half subzone nodes
    vector<Node*> upperNodes;

    // Lower half subzone nodes
    vector<Node*> lowerNodes;

    // Cohesive nodes
    vector<CohesiveNode*> cohesiveNodes;

    // Upper half subzone elements
    vector<ElementQ4*> upperElements;

    // Lower half subzone elements
    vector<ElementQ4*> lowerElements;

    // Cohesive elements
    vector<ElementQ4Cohesive*> cohesiveElements;

// PUBLIC MEMBERS
public:
    // Constructor
    Problem(int spaceDim = 2);

    // Destructor
    ~Problem();

    // Initialization
    void initialize(const vector<double> & xRanges, const vector<int> & edgeNums);

// PRIVATE MEMBERS
private:
    // Initialization of Geometry2D
    void initializeGeometry2D(const vector<double> & xRanges, const vector<int> & edgeNums);

    // Initialization of Nodes
    void initializeNodes();

    // Initialization of Elements
    void initializeElements();

    // Delete Nodes;
    void deleteNodes();

    // Delete Elements
    void deleteElements();

// NOT IMPLEMENTED
private:
    // Copy Constructor
    Problem(const Problem & );

    // Move-assign operator
    const Problem & operator=(const Problem & );
};
