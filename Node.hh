/** @file Node.hh
 * Declaration of class Node and class CohesiveNode
 */
#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include "petsc.h"
#include "petscksp.h"
#include "petscts.h"
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
    
// PUBLIC MEMBERS
public:
    /** TO DO 
     * Change this nodalBodyForce to private, currently for debugging purpose..
     */
    vector<double> _nodalBodyForce;
    
    /** Current vector [displacement, velocity, pressure, tracestrain] */
    vector<double> s;
    
    /** Current vector time derivative d/dt [displacement, velocity, pressure, tracestrain] */
    vector<double> s_t;
    
    /** Nodal properties values (0 - mass density; 
     * 1,2 - body force (force per unit volume); 
     * 3 - \lambda (drained);
     * 4 - shear modulus G;
     * 5 - Biot coefficient \alpha;
     * 6 - Biot modulus M_p;
     * 7 - Fluid mobility \kappa;
     * 8 - Fluid viscosity \mu;
     * 9 - fluid density
     * 10 - reference porosity
     * 11, 12 - fluid body force (force per unit volume)
     * 13 - fluid source density (/s) 
     * ...)
     */
    vector<double> _nodalProperties;

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
         double fluidViscosity = 0.,
         double fluidDensity = 1., 
         double porosity = 0., 
         const vector<double> *fluidBodyForce = NULL, 
         double source = 0.);

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

    // Set fluid density
    void setFluidDensity(double fluidDensity) {
        if (_nodalProperties.size() < 10) _nodalProperties.resize(10);
        _nodalProperties[9] = fluidDensity;
    };

    // Set porosity
    void setPorosity(double porosity) {
        if (_nodalProperties.size() < 10) _nodalProperties.resize(11);
        _nodalProperties[10] = porosity;
    };

    // Set body force
    void setFluidBodyForce(const vector<double> *fluidBodyForce) {
        if (_nodalProperties.size() < 13) _nodalProperties.resize(13);
        if (fluidBodyForce) {
            for (int i = 0; i < _spaceDim; i++) {
                _nodalProperties[11 + i] = (*fluidBodyForce)[i];
            }
        }
        else {
            for (int i = 0; i < _spaceDim; i++) {
                _nodalProperties[11 + i] = 0.;
            }
        }
    };

    // Set fluid source
    void setSource(double source) {
        if (_nodalProperties.size() < 14) _nodalProperties.resize(14);
        _nodalProperties[13] = source;
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
        if (_nodalProperties.size() < 4) throw "Lambda not initialized!";
        return _nodalProperties[3];
    };

    // Get shear modulus
    double getShearModulus() const {
        if (_nodalProperties.size() < 5) throw "Shear modulus not initialized!";
        return _nodalProperties[4];
    };

    // Get Biot coefficient alpha
    double getBiotAlpha() const {
        if (_nodalProperties.size() < 6) throw "Biot Alpha not initialized!";
        return _nodalProperties[5];
    };

    // Get Biot modulus M_p
    double getBiotMp() const {
        if (_nodalProperties.size() < 7) throw "Biot Mp not initialized!";
        return _nodalProperties[6];
    };

    // Get fluid mobility kappa
    double getFluidMobility() const {
        if (_nodalProperties.size() < 8) throw "Fluid mobility not initialized!";
        return _nodalProperties[7];
    };

    // Get fluid viscosity mu
    double getFluidViscosity() const {
        if (_nodalProperties.size() < 9) throw "Fluid viscosity not initialized!";
        return _nodalProperties[8];
    };

    // Get fluid density
    double getFluidDensity() const {
        if (_nodalProperties.size() < 10) throw "Fluid density not initialized!";
        return _nodalProperties[9];
    };

    // Get porosity
    double getPorosity() const {
        if (_nodalProperties.size() < 10) throw "Porosity not initialized!";
        return _nodalProperties[10];
    };

    // Get fluid body force
    vector<double> getFluidBodyForce() const {
        if (_nodalProperties.size() < 13) 
            throw "Fluid bodyforce not initialized!";
        vector<double> res = {_nodalProperties[11], _nodalProperties[12]};
        return res;
    };

    // Get fluid source
    double getSource(double source) {
        if (_nodalProperties.size() < 14) throw "Fluid source not initialized!";
        return _nodalProperties[13];
    };

    /** Initialize s = [displacement, velocity, pressure, trace_strain] */
    void initializeS(const vector<double> & initialS);

    /** Push initial s to the global vector s */
    void pushS(Vec & globalS) const;

    /** Get current s from the global vector s */
    void fetchS(const Vec & globalS);

    /** Get current s_t from the global vector s_t */
    void fetchS_t(const Vec & globalS_t);

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
// PRIVATE MEMBERS
    // Vector of points to lower and upper nodes.
    vector<Node*> _lowerUpperNodes;

// PUBLIC MEMBERS
public:
    // Constructor 1
    CohesiveNode(int ID, int spaceDim = 2);

    /** Nodal properties values (0 - mass density; 
     * 1,2 - body force (force per unit volume); 
     * 3 - fluidMobility_x
     * 4 - fluidMobility_z
     * 5 - fluidViscosity \mu
     * 6 - porosity \phi_f
     * 7 - thickness h
     * 8 - beta_p 
     * 9 - beta_sigma
     * 10 - rateStateA
     * 11 - rateStateB
     * 12 - DRateState
     * 13, 14 - fluidBodyForce
     * 15 - source 
     * ...)
     */

    // Constructor 2
    CohesiveNode(int ID,
                 const vector<double> &XYZ,
                 const vector<int> &DOF,
                 const vector<Node*> &upperLowerNodes, 
                 int spaceDim = 2,
                 int density = 0.0,
                 const vector<double> *bodyForce = NULL,
                 double fluidMobility_x = 0.,
                 double fluidMobility_z = 0., 
                 double fluidViscosity = 1.,
                 double porosity = 1.0,
                 double thickness = 0.01, 
                 double beta_p = 1.0,
                 double beta_sigma = 1.0,
                 double rateStateA = 1.0, 
                 double rateStateB = 1.0, 
                 double DRateState = 1.0,
                 const vector<double> *fluidBodyForce = NULL,
                 double source = 0.);

    // Set ID
    void setID(int ID);

    // Get ID
    int getID() const {
        return _ID;
    }

    // Set nodal coordinates XYZ
    void setXYZ(const vector<double> & XYZ);

    // Get XYZ
    const vector<double> & getXYZ() const {
        return _nodalXYZ;
    };

    // Set DOF 1 [lambda, pressure, theta]
    void setDOF(const vector<int> & DOF) {
        if (_nodalDOF.size() < 2 * (2 * _spaceDim + 2) + _spaceDim + 2) 
            _nodalDOF.resize(2 * (2 * _spaceDim + 2) + _spaceDim + 2);
        if (DOF.size() != _nodalDOF.size()) {
            throw "Input DOF incompatible with nodal DOF!";
        }
        for(int i = 0; i < _nodalDOF.size(); i++) _nodalDOF[i] = DOF[i];
    };    

    // Get DOF 1
    const vector<int> & getDOF() const {
        return _nodalDOF;
    };

    // Set DOF 2
    void setDOF(int index, int DOF) {
        if (_nodalDOF.size() < index) throw "Nodal DOF does not contain this index!";
        _nodalDOF[index] = DOF;
    }

    // Get DOF2
    int getDOF(int index) const {
        if (_nodalDOF.size() < index) throw "Nodal DOF does not contain this index!";
        return _nodalDOF[index];
    };

    // UpdateDOF
    void updateDOF() {
        // u(-+)
        _nodalDOF[0] = _lowerUpperNodes[0]->getDOF(0);
        _nodalDOF[1] = _lowerUpperNodes[0]->getDOF(1);
        _nodalDOF[2] = _lowerUpperNodes[1]->getDOF(0);
        _nodalDOF[3] = _lowerUpperNodes[1]->getDOF(1);
        
        // v(-+) 
        _nodalDOF[4] = _lowerUpperNodes[0]->getDOF(2);
        _nodalDOF[5] = _lowerUpperNodes[0]->getDOF(3);
        _nodalDOF[6] = _lowerUpperNodes[1]->getDOF(2);
        _nodalDOF[7] = _lowerUpperNodes[1]->getDOF(3);
        
        // p(-+) 
        _nodalDOF[8] = _lowerUpperNodes[0]->getDOF(4);
        _nodalDOF[9] = _lowerUpperNodes[1]->getDOF(4);

        // trace_strain(-+)
        _nodalDOF[10] = _lowerUpperNodes[0]->getDOF(5);
        _nodalDOF[11] = _lowerUpperNodes[1]->getDOF(5);
    };

    // Set lowerUpperNodes
    void setLowerUpperNodes(const vector<Node*> & lowerUpperNodes) {
        if (lowerUpperNodes.size() < 2) throw "lowerUpperNodes size not compatible with this cohesive node!";
        _lowerUpperNodes.resize(2);
        _lowerUpperNodes[0] = lowerUpperNodes[0];
        _lowerUpperNodes[1] = lowerUpperNodes[1];
    };

    // Get upperLowerNodes
    const vector<Node*> & getLowerUpperNodes() const {
        return _lowerUpperNodes;
    };

    // Set spaceDim
    void setSpaceDim(int spaceDim);

    // Get spaceDim
    int getSpaceDim() const {
        return _spaceDim;
    };

    // Set mass density
    void setMassDensity(double density);

    // Get massDensity
    double getMassDensity() const {
        if (_nodalProperties.size() < 1) throw "Mass density not initialized!";
        return _nodalProperties[0];
    };

    // Set body force
    void setBodyForce(const vector<double> *bodyForce);

    // Get body force
    vector<double> getBodyForce() const {
        if (_nodalProperties.size() < 3) throw "Body force not initialized!";
        vector<double> res(_spaceDim, 0.);
        for (int i = 0; i < _spaceDim; i++) {
            res[i] = _nodalProperties[1 + i];
        }
        return res;
    };

    // Set fluid mobility in x-dir
    void setFluidMobility_x(double fluidMobility_x) {
        if (_nodalProperties.size() < 4) _nodalProperties.resize(4);
        _nodalProperties[3] = fluidMobility_x;
    };

    // Get fluid mobility in x-dir
    double getFluidMobility_x() const {
        if (_nodalProperties.size() < 4) throw "FluidMobility_x not initialized!";
        return _nodalProperties[3];
    }

    // Set fluid mobility in z-dir, 
    void setFluidMobility_z(double fluidMobility_z) {
        if (_nodalProperties.size() < 5) _nodalProperties.resize(5);
        _nodalProperties[4] = fluidMobility_z;
    };

    // Get fluid mobility in z-dir
    double getFluidMobility_z() const {
        if (_nodalProperties.size() < 5) throw "FluidMobility_z not initialized!";
        return _nodalProperties[4];
    }

    // Set fluid viscosity \mu
    void setFluidViscosity(double fluidViscosity) {
        if (_nodalProperties.size() < 6) _nodalProperties.resize(6);
        _nodalProperties[5] = fluidViscosity;
    };

    // Get fluid viscosity \mu
    double getFluidViscosity() const {
        if (_nodalProperties.size() < 6) throw "FluidViscosity not initialized!";
        return _nodalProperties[5];
    }

    // Set porosity \phi_f
    void setPorosity(double porosity) {
        if (_nodalProperties.size() < 7) _nodalProperties.resize(7);
        _nodalProperties[6] = porosity;
    };

    // Get porosity \phi_f
    double getPorosity() const {
        if (_nodalProperties.size() < 7) throw "Porosity not initialized!";
        return _nodalProperties[6];
    };

    // Set thickness of the shear layer
    void setThickness(double thickness) {
        if (_nodalProperties.size() < 8) _nodalProperties.resize(8);
        _nodalProperties[7] = thickness;
    };

    // Get thickness of the shear layer
    double getThickness() const {
        if (_nodalProperties.size() < 8) throw "Thickness not initialized!";
        return _nodalProperties[7];
    }

    // Set beta_p
    void setBeta_p(double beta_p) {
        if (_nodalProperties.size() < 9) _nodalProperties.resize(9);
        _nodalProperties[8] = beta_p;
    };

    // Get beta_p
    double getBeta_p() const {
        if (_nodalProperties.size() < 9) throw "Beta_p not initialized!";
        return _nodalProperties[8];
    };

    // Set beta_sigma
    void setBeta_sigma(double beta_sigma) {
        if (_nodalProperties.size() < 10) _nodalProperties.resize(10);
        _nodalProperties[9] = beta_sigma;
    };

    // Get beta_sigma
    double getBeta_sigma() const {
        if (_nodalProperties.size() < 10) throw "Beta_sigma not initialized!";
        return _nodalProperties[9];
    };

    // Set rate and state constant a
    void setRateStateA(double rateStateA) {
        if (_nodalProperties.size() < 11) _nodalProperties.resize(11);
        _nodalProperties[10] = rateStateA;
    };

    // Get rate and state constant a
    double getRateStateA() const {
        if (_nodalProperties.size() < 11) throw "Rate and State A not initialized!";
        return _nodalProperties[10];
    };

    // Set rate and state constant b
    void setRateStateB(double rateStateB) {
        if (_nodalProperties.size() < 12) _nodalProperties.resize(12);
        _nodalProperties[11] = rateStateB;
    };

    // Get rate and state constant b
    double getRateStateB() const {
        if (_nodalProperties.size() < 12) throw "Rate and State B not initialized!";
        return _nodalProperties[11];
    };

    // Set rate and state D_rs
    void setDRateState(double DRateState) {
        if (_nodalProperties.size() < 13) _nodalProperties.resize(13);
        _nodalProperties[12] = DRateState;
    };

    // Get rate and state D_rs
    double getDRateState() const {
        if (_nodalProperties.size() < 13) throw "D Rate and State not initialized!";
        return _nodalProperties[12];
    };

    // Set fluid body force
    void setFluidBodyForce(const vector<double> *fluidBodyForce) {
        if (_nodalProperties.size() < 15) _nodalProperties.resize(15);
        if (fluidBodyForce) {
            for (int i = 0; i < _spaceDim; i++) {
                _nodalProperties[13 + i] = (*fluidBodyForce)[i];
            }
        }
        else {
            for (int i = 0; i < _spaceDim; i++) {
                _nodalProperties[13 + i] = 0.;
            }
        }
    };
    
    // Get fluid body force
    vector<double> getFluidBodyForce() const {
        if (_nodalProperties.size() < 15) throw "Fluid body force not initialized!";
        vector<double> res = {_nodalProperties[13], _nodalProperties[14]};
        return res;
    }

    // Set fluid source
    void setSource(double source) {
        if (_nodalProperties.size() < 16) _nodalProperties.resize(16);
        _nodalProperties[15] = source;
    };

    // Get fluid source
    double getSource() const {
        if (_nodalProperties.size() < 16) throw "Source not initialized!";
        return _nodalProperties[16];
    };

    /** Initialize s */
    void initializeS(const vector<double> & initialS);

    /** Push initial s to the global vector s */
    void pushS(Vec & globalS) const;

// NOT IMPLEMENTED
private:
    // Copy constructor
    CohesiveNode(const CohesiveNode &);

    // Move-assign operator
    const CohesiveNode & operator=(const CohesiveNode &);
};