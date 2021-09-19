/** @file Geometry2D.hh
 * Class Geometry2D
 * 
 * Give the length in x and y direction for a subzone
 * x in [0, xRange]
 * y in [0, yRange] -- upper subzone
 *      [0, - yRange] -- lower subzone 
 * 
 * Give the number of element edges in x and y
 * Uniformly divide x and y ranges such that all elements are rectangles.
 */

/** Class Geometry2D
 * Generate 2D mesh based on 4 values above
 */
class Geometry2D {
    friend class Problem;
// PROTECTED MEMBERS
protected:
    // Ranges of x and y domain for a subdomain
    double xRange, yRange;
    
    // Number of element edges for a subdomain
    int xEdgeNum, yEdgeNum;

    // Number of element edges for a subdomain
    int xNodeNum, yNodeNum;

    // Element numbers
    int nOfElements;

    // Node numbers
    int nOfNodes;

// PUBLIC MEMBERS
public:
    // Constructor
    Geometry2D(double xRange, double yRange, int xEdgeNum, int yEdgeNum) :
        xRange(xRange), 
        yRange(yRange), 
        xEdgeNum(xEdgeNum), 
        yEdgeNum(yEdgeNum), 
        xNodeNum(xEdgeNum + 1), 
        yNodeNum(yEdgeNum + 1), 
        nOfElements(xEdgeNum * yEdgeNum), 
        nOfNodes(xNodeNum * yNodeNum) {};

    // Destructor
    ~Geometry2D(){};

// NOT IMPLEMENTED
private:
    // Copy constructor
    Geometry2D(const Geometry2D &);

    // Move-assign operator
    const Geometry2D & operator=(const Geometry2D &);
};