/** @file Node.hh
 * Declaration of class Node and class CohesiveNode
 */
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
using namespace std;

//------------------------------------------------------------------------------------------------------------------------------------------------------
/** Class Node, declares node
 * contains vector XYZ (2 or 3 dim)
 * node ID
 * nodal degree of freedom
 * and its match to the global one
 */
class Node {
// PRIVATE MEMBERS
protected:
    // Space dim
    int _spaceDim;

    // ID number of the node
    int _ID;

    // Nodal coordinates
    vector<double> _nodalXYZ;

    // Nodal variables existence, false OR true, 
    // [displacement, velocity, pressure, trace_strain]
    // vector<string> *_varNames;

    // Nodal degree of freedom, 0 free, 1 fixed
    // displacement(spaceDim), velocity(spaceDim), pressure(1) 
    vector<int> _nodalDOF;

    /** Nodal properties values (0 - mass density; 
     * 1,2 - body force (force per unit volume); 
     * 3 - \lambda (drained);
     * 4 - shear modulus G;
     * 5 - Biot coefficient \alpha;
     * 6 - Biot modulus M_p;
     * 7 - Fluid mobility \kappa;
     * 8 - Fluid viscosity \mu;
     * ...)
     */
    vector<double> _nodalProperties;

    /** Current vector [displacement, velocity, pressure, tracestrain] */
    vector<double> s;
    
    /** Current vector time derivative d/dt [displacement, velocity, pressure, tracestrain] */
    vector<double> s_t;
    
// PUBLIC MEMBERS
public:
    /** TO DO 
     * Change this nodalBodyForce to private, currently for debugging purpose..
     */
    vector<double> _nodalBodyForce;
    
    // Default Constructor
    Node();

    // Constructor 1
    Node(int ID, int spaceDim = 2);

    // Constructor 2
    Node(int ID, 
         const vector<double> & XYZ, 
         const vector<int> & DOF, 
         int spaceDim = 2, 
         double density = 1., 
         const vector<double> *bodyForce = NULL, 
         double lambda = 0., 
         double shearModulus = 0., 
         double biotAlpha = 0., 
         double biotMp = 0., 
         double fluidMobility = 0., 
         double fluidViscosity = 0.);

    // Destructor
    ~Node();

    // Set spaceDim
    void setSpaceDim(int spaceDim);

    // Set ID
    void setID(int ID);

    // Set nodal coordinates XYZ
    void setXYZ(const vector<double> & XYZ);

    // Set DOF 1
    void setDOF(const vector<int> & DOF);

    // Set DOF 2
    void setDOF(int index, int DOF);

    // Set mass density
    void setMassDensity(double density);

    // Set body force
    void setBodyForce(const vector<double> *bodyForce);

    // Set lambda
    void setLambda(double lambda) {
        if (_nodalProperties.size() < 4) _nodalProperties.resize(4);
        _nodalProperties[3] = lambda;
    };

    // Set shear modulus
    void setShearModulus(double shearModulus) {
        if (_nodalProperties.size() < 5) _nodalProperties.resize(5);
        _nodalProperties[4] = shearModulus;
    };

    // Set Biot coefficient alpha
    void setBiotAlpha(double biotAlpha) {
        if (_nodalProperties.size() < 6) _nodalProperties.resize(6);
        _nodalProperties[5] = biotAlpha;
    };

    // Set Biot modulus M_p
    void setBiotMp(double biotMp) {
        if (_nodalProperties.size() < 7) _nodalProperties.resize(7);
        _nodalProperties[6] = biotMp;
    };

    // Set fluid mobility kappa
    void setFluidMobility(double fluidMobility) {
        if (_nodalProperties.size() < 8) _nodalProperties.resize(8);
        _nodalProperties[7] = fluidMobility;
    };

    // Set fluid viscosity mu
    void setFluidViscosity(double fluidViscosity) {
        if (_nodalProperties.size() < 9) _nodalProperties.resize(9);
        _nodalProperties[8] = fluidViscosity;
    };

    // Get spaceDim
    int getSpaceDim() const;

    // Get ID
    int getID() const;

    // Get XYZ
    const vector<double> & getXYZ() const;

    // Get DOF 1
    const vector<int> & getDOF() const;

    // Get DOF 2
    int getDOF(int index) const;

    // Get mass density
    double getMassDensity() const;

    // Get body force
    vector<double> getBodyForce() const;

    // Get lambda
    double getLambda() const {
        if (_nodalProperties.size() < 4) throw("Lambda not initialized!");
        return _nodalProperties[3];
    };

    // Get shear modulus
    double getShearModulus() const {
        if (_nodalProperties.size() < 5) throw("Shear modulus not initialized!");
        return _nodalProperties[4];
    };

    // Get Biot coefficient alpha
    double getBiotAlpha() const {
        if (_nodalProperties.size() < 6) throw("Biot Alpha not initialized!");
        return _nodalProperties[5];
    };

    // Set Biot modulus M_p
    double getBiotMp() const {
        if (_nodalProperties.size() < 7) throw("Biot Mp not initialized!");
        return _nodalProperties[6];
    };

    // Set fluid mobility kappa
    double getFluidMobility() const {
        if (_nodalProperties.size() < 8) throw("Fluid mobility not initialized!");
        return _nodalProperties[7];
    };

    // Set fluid viscosity mu
    double getFluidViscosity() const {
        if (_nodalProperties.size() < 9) throw("Fluid viscosity not initialized!");
        return _nodalProperties[8];
    };

    /** Initialize s = [displacement, velocity, pressure, trace_strain] */
    void initializeS(const vector<double> & initialS);

    /** Push initial s to the global vector s */
    void pushS(vector<double> & globalS) const;

    /** Get current s from the global vector s */
    void fetchS(const vector<double> & globalS);

    /** Get current s_t from the global vector s_t */
    void fetchS_t(const vector<double> & globalS_t);

    // Output nodal information to a file
    void outputInfo(ofstream & myFile, bool outputElse = false) const;

// NOT IMPLEMENTED
private:
    Node(const Node &); ///< Copy constructor, not implemented
    
    const Node& operator=(const Node &); ///< Move-assign operator, not implemented
};

//------------------------------------------------------------------------------------------------------------------------------------------------------
/** Class CohesiveNode, 
 * deals with the cohesive nodes initialization and calculation.
 * has DOF: {lambda(spaceDim), pressure_fault, theta}
 */
class CohesiveNode : public Node {
// PUBLIC MEMBERS
public:
    // Constructor 1
    CohesiveNode(int ID, int spaceDim = 2);

    // Constructor 2
    CohesiveNode(int ID, const vector<double> & XYZ, const vector<int> & DOF, int spaceDim = 2);

    /** Initialize s */
    void initializeS(const vector<double> & initialS);

// NOT IMPLEMENTED
private:
    // Copy constructor
    CohesiveNode(const CohesiveNode &);

    // Move-assign operator
    const CohesiveNode & operator=(const CohesiveNode &);
};