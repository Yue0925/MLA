/*----------------------------------------------------------------------------------------------+
|               Cutting Planes Engine to solve an LP adding constraints one by one              |
|                   - using the cheshire cat technique, no need to include cplex                | 
|                     .h headers in all files that include CuttingPlanesEngine.h,               |
|                     significantly speeding-up compilation                                     |    
-----------------+---------------------------------------------------------------+--------------+
                 | Author: Daniel Porumbel 2018       daniel.porumbel@cnam.fr    |
                 |License: Creative Commons Attribution 3.0 Unported License     |
                 |         http://creativecommons.org/licenses/by/3.0/           |             
                 +--------------------------------------------------------------*/

#include "CuttingPlanesEngine.h"
#include <cstdlib>
#include <vector>
#include <cmath>
using namespace std;

#ifdef USE_CPU_TIME_FUNCTION
#include "general.h"
#else
#include <time.h>
double getCPUTime(){
    return clock()/CLOCKS_PER_SEC; //attention to the wrap around issue, in man 3 clock
}
#endif


//void CuttingPlanesEngine::turnVarInteger(int i)
//{
//}
//void CuttingPlanesEngine::turnAllVarsInteger()
//{
//}
//void CuttingPlanesEngine::freeData(double*&a,double*&b, double**&c,int nr)
//{
//}

CuttingPlanesEngine::CuttingPlanesEngine(int nrVars, cutSeprtSimple_t cutSeprtSimple): 
                           n(nrVars), 
                           internalCutSeprtSimple(cutSeprtSimple),
                           env(),
                           model(env)
{
    vars = new GRBVar[n];
    for(int i=0;i<n;i++){
        char valName[10];
        sprintf(valName,"x%d",i+1);
        vars[i] = model.addVar(0,GRB_INFINITY, 0, GRB_CONTINUOUS, valName);
    }
    primals = new double[n];
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);
}

CuttingPlanesEngine::~CuttingPlanesEngine()
{ 
    //for(int i=noRows-1;i>=0;i--) modelDelCut(i);
    delete[] primals;
}
void CuttingPlanesEngine::setVarBounds(int start, int end, double varLb, double  varUb)
{
    for (int i=start;i<=end;i++){
        model.addConstr(vars[i]>=varLb);
        model.addConstr(vars[i]<=varUb);
    }
}
void CuttingPlanesEngine::turnAllVarsInteger(){
    for(int i=0;i<n;i++)
        vars[i].set(GRB_CharAttr_VType, GRB_INTEGER);
}

void CuttingPlanesEngine::setVarBounds(double varLb, double  varUb)
{
    for (int i=0;i<n;i++){
        model.addConstr(vars[i]>=varLb);
        model.addConstr(vars[i]<=varUb);
    }
}
#define ADD_OBJ(MinOrMax)                                    \
    GRBLinExpr expr = 0;                                     \
    for(int i=0;i<n;i++)                                     \
            expr+=coefs[i]* vars[i];                         \
    model.setObjective(expr,MinOrMax);                       
    //expr.end();
void CuttingPlanesEngine::setObjCoefsMaximize(int * coefs)
{
    ADD_OBJ(GRB_MAXIMIZE);
    maximize = 1;
    currObj  = INT_MAX;
}
void CuttingPlanesEngine::setObjCoefsMinimize(int * coefs)
{
    ADD_OBJ(GRB_MINIMIZE);
    maximize = 0;
}
void CuttingPlanesEngine::setObjCoefsMaximize(double * coefs)
{
    ADD_OBJ(GRB_MAXIMIZE);
    maximize = 1;
    currObj  = INT_MAX;
}
void CuttingPlanesEngine::setObjCoefsMinimize(double * coefs)
{
    ADD_OBJ(GRB_MINIMIZE);
    maximize = 0;
}

void CuttingPlanesEngine::modelAddCut(int * coefs, double rightHand, int sense)
{
    int i;
    GRBLinExpr expr = 0;
    for(i=0;i<n;i++){
            expr+=coefs[i]* vars[i];
    }
    if(sense==1)
        model.addConstr(expr<=rightHand);
    if(sense==0)
        model.addConstr(expr==rightHand);
    if(sense==-1)
        model.addConstr(expr>=rightHand);
    //expr.end();    
}
void CuttingPlanesEngine::modelAddCut(double * coefs, double rightHand, int sense)
{
    int i;
    GRBLinExpr expr = 0;
    for(i=0;i<n;i++){
            expr+=coefs[i]* vars[i];
    }
    if(sense==1)
        model.addConstr(expr<=rightHand);
    if(sense==0)
        model.addConstr(expr==rightHand);
    if(sense==-1)
        model.addConstr(expr>=rightHand);
    //expr.end();    
}
void CuttingPlanesEngine::modelAddCut(int * coefs, double rightHand)
{
    if(maximize)
        modelAddCut(coefs, rightHand,1);
    else
        modelAddCut(coefs, rightHand,-1);
}
void CuttingPlanesEngine::modelAddCut(double * coefs, double rightHand)
{
    if(maximize)
        modelAddCut(coefs, rightHand,1);
    else
        modelAddCut(coefs, rightHand,-1);
}
void CuttingPlanesEngine::modelAddEquality(double * coefs, double rightHand)
{
    modelAddCut(coefs, rightHand,0);
}

double CuttingPlanesEngine::getObjVal()
{
    return currObj;
}
double CuttingPlanesEngine::solve()
{
    model.optimize();
    currObj = model.get(GRB_DoubleAttr_ObjVal);
    return currObj;
}//end solve
void CuttingPlanesEngine::getPrimals(double *ystar)
{
    for(int i=0;i<n;i++)
        ystar[i] = vars[i].get(GRB_DoubleAttr_X);
}


int CuttingPlanesEngine::runCutPlanes(const int itMax, const double tmMax, int& it, double&tm)
{
    double rHand;
    it = 0;
    solve();
    while(true){
        it++;
        double* newCut = new double[n];
        getPrimals(primals);
        if(internalCutSeprtSimple(n,primals,newCut,rHand)>=0)
            break;
        modelAddCut(newCut,rHand);
        delete[] newCut;
        solve();
        if(it==10)
            model.write("test.lp");
    }
    return 0;
    

    return 0;
}

int CuttingPlanesEngine::runCutPlanes(int& it, double&tm)
{
    return runCutPlanes(INT_MAX, INT_MAX, it, tm);
}
