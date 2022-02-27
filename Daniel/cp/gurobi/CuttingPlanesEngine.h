/*----------------------------------------------------------------------------------------------+
|               Cutting Planes Engine to solve an LP adding constraints one by one              |
------------------------------------------------------------------------------------------------+
See accompanying .cpp file for licence info     */

#ifndef CUTPLANESENGINE_H_INCLUDED
#define CUTPLANESENGINE_H_INCLUDED


#include "gurobi_c++.h"
#include <climits>



//given current solution primals (dimension nrVars), find a new cut newCut and the
//associated rightHandVal. Writing a=newCut et b=rightHandVal, the added cut is:
//  maximize=1   =>              a^T x <= b   => Return value  b - a^T x
//  maximize=0   =>              a^T x >= b   => Return value: a^T x - b
//The returned violation is a negative value if x violates newCut, positive otherwise
//If the cut separator realizes the gap is closed, it returns INT_MAX and cutPlanes stops
//Attention: one can modify primals, but it is not recommended
typedef double (*cutSeprtSimple_t) (const int nrVars, double*primals, 
                                  double * newCut, double&newRightHand);
class CuttingPlanesEngine{
    public:
    //constructs object using nrVars unbounded variables
    CuttingPlanesEngine(int nrVars, cutSeprtSimple_t pSimple);
    //classical destructor
    ~CuttingPlanesEngine();

    void setVarBounds(double*varLB, double*varUB);
    //the lower and upper bound on each variable is defined by doubles varLB and resp varUb
    void setVarBounds(int start_idx, int end_idx, double varLB, double varUB);
    void setVarBounds(double varLB, double varUB);
    void turnAllVarsInteger();
    //the lower bound on each variable is defined varLB , impose no upper bounding
    //Set the objective function for minimization
    //All cuts have form ...>=rightHand
    void setObjCoefsMinimize(double*    coefs);
    void setObjCoefsMinimize(int*    coefs);
    //Set the objective function for maximization
    //All cuts have form ...<=rightHand
    void setObjCoefsMaximize(double*    coefs);
    void setObjCoefsMaximize(int*    coefs);
    double getObjVal();
    //Add to the model a new cut 
    void modelAddCut(int * coefs, double rightHand);
    void modelAddCut(double * coefs, double rightHand);
    //Add to the model an equality, return id of the new constraint
    void modelAddEquality(double * coefs, double rightHand);
    //solve over current columns and return the value of the objective function (restricted model)
    double solve();
    //puts the dual variables in *yStar (if non null)
    void getPrimals(double*ystar);
    //puts xstart in the primal variables that are a starting point

    //Solve column generator using a cutSeprt routine called 'cutSeprt' 
    //Last 2 args (filled by the routine) indicate the nr. of iterations and time required to finish
    //Returns EXIT_SUCCESS (0) in case of success or -1 otherwise
    int runCutPlanes(int& it, double&tm);
    int runCutPlanes(int, double, int&, double&);
    private:
    int n;                          //the total number of decision variables
    double * primals;
    double currObj;                 //the optimum of the restricted program 
    int maximize ;
    cutSeprtSimple_t internalCutSeprtSimple;

    GRBEnv env ;
    GRBModel model ;
    GRBVar* vars ;

    //below is called by modelAddEquality and modelAddCut above
    void modelAddCut(int * coefs, double rightHand, int sense);
    void modelAddCut(double * coefs, double rightHand, int sense);
};
#endif
