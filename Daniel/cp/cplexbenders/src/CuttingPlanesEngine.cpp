/*----------------------------------------------------------------------------------------------+
|               Cutting Planes Engine to solve an LP adding constraints one by one              |
|                   - using the cheshire cat technique, no need to include cplex                | 
|                     .h headers in all files that include CuttingPlanesEngine.h,               |
|                     significantly speeding-up compilation                                     |    
|                   - @Daniel Porumbel 2018, Creative Commons Attribution 3.0 Unported License  |                                                   |
-----------------+---------------------------------------------------------------+-------------*/

#include "CuttingPlanesEngine.h"
#include <cstdlib>
#include <vector>
#include <cmath>
#include <ilcplex/ilocplex.h>                    
using namespace std;



double lowerBound;                             //a global visible anywhere
double upperBound;                             //a global visible anywhere
int    switchToIntVarsNow;                     //a global usable from cutSeprt
class cplexCheshireData{                      //From other files, one can refer to this struct
    public:
    IloEnv         env;                        //without needing to include cplex headers
    IloModel       model;
    IloNumVarArray vars;
    IloCplex       cplex ;
    IloNumArray    lb ;
    IloNumArray    ub ;
    IloRangeArray cuts;                     
    #ifdef TMP_CUT_ERASER
    IloRange rrecConstr;
    #endif
    cplexCheshireData(int n):
            env(), model(env), vars(env,n,0,10000000), cplex (env),
            lb(env,n), ub (env,n), cuts(env){
        //it can assign -EPS to some variables >0
        //additional constraints below help avoiding this
        for(int i=0;i<n;i++){
            model.add(vars[i]>=0);
            lb[i]=0;
        }
    };

};

void CuttingPlanesEngine::removeVar(int outIdx)
{
    d.vars[outIdx].end();
    d.vars.remove(outIdx);
    d.lb.remove(outIdx);
    d.ub.remove(outIdx);
    n--;
    #ifdef TMP_CUT_ERASER
    if(outIdx==recordedConstr){
        CPLOG("Erasing"<<d.rrecConstr<<endl);
        d.rrecConstr.end();
        //d.model.remove(d.rrecConstr);
    }
    #endif

}
void CuttingPlanesEngine::addVar(double objCoef, double min, double max, int swapOK)
{
    IloObjective obj = d.cplex.getObjective();
    d.vars.add(obj(objCoef));
    n++;
    d.lb.add(min);
    d.ub.add(max);

    d.model.add(d.vars[n-1]>=d.lb[n-1]);
    d.model.add(d.vars[n-1]<=d.ub[n-1]);
    if(swapOK){
        IloNumVar tmp = d.vars[n-1];
        d.vars[n-1] = d.vars[n-2];
        d.vars[n-2] = tmp; 
        double  tmpd= d.lb[n-1];
        d.lb[n-1]   = d.lb[n-2];
        d.lb[n-2]   = tmpd;
        tmpd= d.ub[n-1];
        d.ub[n-1]   = d.ub[n-2];
        d.ub[n-2]   = tmpd;
    }
    //delete[] primals;
    if(primals!=NULL)
        delete[] primals;
    primals = new double[n];
}
void CuttingPlanesEngine::addVar(double objCoef, double min, double max)
{
    return addVar(objCoef,min,max,0);
}
void CuttingPlanesEngine::init()
{
    //by default there is no log output
    cplog.setstate(std::ios_base::failbit);
    ::upperBound = INT_MAX;
    ::lowerBound = INT_MIN;
    ::switchToIntVarsNow = 0;
    timeoutSet = -1;

    d.cplex.extract(d.model);
    d.cplex.setParam(IloCplex::Threads, THREADS);
    //For numerical precision in solving LPs
    d.cplex.setParam(IloCplex::NumericalEmphasis, true);

    #ifdef NO_CPLEX_OUTPUT
    d.cplex.setOut(d.env.getNullStream());
    d.cplex.setError(d.env.getNullStream());
    //This says when it is MIP
    d.cplex.setWarning(d.env.getNullStream());
    #endif 
    primals  = new double[n];
    noRows = 0;
    tmOnlySolve = 0;
    intVars = 0;
    intStatus = new int[n];
    totalNrCoefs      = 0;
    maximize = 0;
}
int  CuttingPlanesEngine::nbIntVars()
{
    return intVars;
}
void CuttingPlanesEngine::setToleranceParamToEpsilon()
{
     //the first below seems an old version (<Cplex 2.6) of second below
     #ifndef CPLEXVER
        d.cplex.setParam(IloCplex::EpInt, EPS);
     #else
        //CPLEXVER passed in command line via -D. It should contain 3 digits
        #if (CPLEXVER>=126)
            d.cplex.setParam(IloCplex::Param::MIP::Tolerances::Integrality,EPS);//you said cplex >=12.6
        #endif
     #endif
}
void CuttingPlanesEngine::turnVarInteger(int i)
{
     if(intStatus[i]==1){
        CPLOG("Variable " << i <<" is already integer, can not re-turn it integer.");
        return;
     }
     CPLOG("Turning variable " << i <<" integer...");
     d.model.add(IloConversion(d.env, d.vars[i], ILOINT));
     setToleranceParamToEpsilon();
     intVars ++;
     intStatus[i] = 1;
     CPLOG("Done\n");
}
void CuttingPlanesEngine::turnAllVarsInteger()
{

    if(maximize) 
         upperBound=INT_MAX;
    else
         lowerBound=INT_MIN;

    if(intVars==0){
        d.model.add(IloConversion(d.env, d.vars, ILOINT));
        setToleranceParamToEpsilon();
        intVars =n ;
        for(int i=0;i<n;i++)
             intStatus[i] = 1;
        return ;
    }
    if(intVars<n){
        CPLOG("\n\nMaking all variables DISCRETE******************************************************** \n\n\n");
        //d.model.add(IloConversion(d.env, d.vars, ILOINT));
        for(int i=0;i<n;i++)
            turnVarInteger(i);
    }
}

void CuttingPlanesEngine::activateLog()
{
    cplog.clear();                   //it clear fail bit, if previously set
}
void CuttingPlanesEngine::printProgressMsg(ostream& cstream)
{
    cplog.rdbuf(cstream.rdbuf());
}
CuttingPlanesEngine::CuttingPlanesEngine(int nrVars, cutSeprtSimple_t cutSeprtSimple): 
                           cplog(clog.rdbuf()),//disabled by default via failbit
                           n(nrVars), currObj(0), d(*new cplexCheshireData(n)),
                           internalCutSeprtSimple(cutSeprtSimple),
                           internalCutSeprtSolver(NULL),
                           internalCutSeprtExtended(NULL),
                           maxMoreConstr(0),
                           turnIntegerEnd(0)
{
    init();
}
CuttingPlanesEngine::CuttingPlanesEngine(int nrVars, cutSeprtSolver_t cutSeprtSolver): 
                           cplog(clog.rdbuf()),//disabled by default via failbit
                           n(nrVars), currObj(0), d(*new cplexCheshireData(n)),
                           internalCutSeprtSimple(NULL),
                           internalCutSeprtSolver(cutSeprtSolver),
                           internalCutSeprtExtended(NULL),
                           maxMoreConstr(0),
                           turnIntegerEnd(0)
{
    init();
}
CuttingPlanesEngine::CuttingPlanesEngine(int nrVars, cutSeprtExtended_t cutSeprtExtended, int maxMoreConstrToReturn): 
                           cplog(clog.rdbuf()),//disabled by default via failbit
                           n(nrVars), currObj(0), d(*new cplexCheshireData(n)),
                           internalCutSeprtSimple(NULL),
                           internalCutSeprtSolver(NULL),
                           internalCutSeprtExtended(cutSeprtExtended),
                           maxMoreConstr(maxMoreConstrToReturn),
                           turnIntegerEnd(0)
{
    init();
}

CuttingPlanesEngine::~CuttingPlanesEngine()
{ 
    d.lb.end();
    d.ub.end();
    for(int i=noRows-1;i>=0;i--)
        modelDelCut(i);
    d.cuts.end();
    d.vars.end();
    d.model.end();
    d.env.end();
    delete[] primals;
    delete[] intStatus;
}
void CuttingPlanesEngine::setVarBounds(double * varLb, double * varUb)
{
    for (int i=0;i<n;i++){
        d.model.add(d.vars[i]>=varLb[i]);
        d.model.add(d.vars[i]<=varUb[i]);
        d.lb[i] = varLb[i];
        d.ub[i] = varUb[i];
    }
    d.vars.setBounds(d.lb,d.ub);
}
void CuttingPlanesEngine::setVarBounds(int start, int end, double varLb, double varUb){
    for (int i=start;i<=end;i++){
        d.lb[i] = varLb;
        d.ub[i] = varUb;
        d.model.add(d.vars[i]>=d.lb[i]);
        d.model.add(d.vars[i]<=d.ub[i]);
    }
}
void CuttingPlanesEngine::setVarBounds(double varLb, double  varUb)
{
    for (int i=0;i<n;i++){
        d.lb[i] = varLb;
        d.ub[i] = varUb;
        d.model.add(d.vars[i]>=d.lb[i]);
        d.model.add(d.vars[i]<=d.ub[i]);
    }
    d.vars.setBounds(d.lb,d.ub);
}
void CuttingPlanesEngine::setVarLowerBoundsOnly(double varLb)
{
    for (int i=0;i<n;i++){
        d.lb[i] = varLb;
        d.ub[i] = IloInfinity;
        d.model.add(d.vars[i]>=d.lb[i]);
    }
    d.vars.setBounds(d.lb,d.ub);
}
#define ADD_OBJ(MinOrMax)                                    \
    IloExpr expr(d.env);                                     \
    for(int i=0;i<n;i++)                                     \
        if(abs(coefs[i])>EPS)                            \
            expr+=coefs[i]* d.vars[i];                       \
    IloObjective obj(d.env, expr, IloObjective::MinOrMax);   \
    d.model.add(obj);                                        \
    expr.end();
void CuttingPlanesEngine::setObjCoefsMaximize(int * coefs)
{
    ADD_OBJ(Maximize);
    maximize = 1;
    currObj  = INT_MAX;
}
void CuttingPlanesEngine::setObjCoefsMaximize(double * coefs)
{
    ADD_OBJ(Maximize);
    maximize = 1;
    currObj  = INT_MAX;
}
void CuttingPlanesEngine::setObjCoefsMinimize(int * coefs)
{
    ADD_OBJ(Minimize);
    maximize = 0;
}
void CuttingPlanesEngine::setObjCoefsMinimize(double * coefs)
{
    ADD_OBJ(Minimize);
    maximize = 0;
}
void CuttingPlanesEngine::alwaysTurnIntegerInTheEnd()
{
    turnIntegerEnd = 1;
}
int CuttingPlanesEngine::getNbCuts()
{
    return noRows;
}
double CuttingPlanesEngine::getCutDualVal(int i)
{
    return d.cplex.getDual(d.cuts[i]);
}
double CuttingPlanesEngine::getCutRightHand(int i)
{
    return d.cuts[i].getLB();
}
void CuttingPlanesEngine::modelDelCut(int i)
{
    d.cuts[i].end();
}
int CuttingPlanesEngine::modelAddCut(double * coefs, double rightHand)
{
    int i;
    IloExpr expr(d.env);
    for(i=0;i<n;i++)
        if(abs(coefs[i])>EPS_CUT_COEFS){
            totalNrCoefs+=1;
            expr+=coefs[i]* d.vars[i];
        }
    #ifdef TMP_CUT_ERASER
    if((!recordedConstr)&&(coefs[n-3]>0)&&(coefs[n-3]<1)){
        recordedConstr = n-3;
        d.rrecConstr = IloRange(expr<=rightHand);
        d.model.add(d.rrecConstr);
        expr.end();
        assert(1==0);
        return -1;
    }
    #endif
    CPLOG("Cut Generated:");
    #ifndef NDEBUG
    for(i=0;i<n;i++)
        if(abs(coefs[i])>EPS_CUT_COEFS){
                    CPLOG(coefs[i]<<"*y["<<i<<"]+");
        }
    if(maximize){
        CPLOG("<=");
    } else{
        CPLOG(">=");
    }
    CPLOG(rightHand<<endl);
    #endif
    if(maximize)
        d.cuts.add(expr<=rightHand);
    else
        d.cuts.add(expr>=rightHand);
    noRows++;
    d.model.add(d.cuts[noRows-1]);
    expr.end();    
    return noRows-1;
}

int CuttingPlanesEngine::modelAddCut(int * coefs, double rightHand, int sense)
{
    int i;
    IloExpr expr(d.env);
    for(i=0;i<n;i++){
            expr+=coefs[i]* d.vars[i];
            totalNrCoefs++;
    }
    if(sense==1)
        d.cuts.add(expr<=rightHand);
    if(sense==0)
        d.cuts.add(expr==rightHand);
    if(sense==-1)
        d.cuts.add(expr>=rightHand);

    noRows++;
    d.model.add(d.cuts[noRows-1]);
    expr.end();    
    return noRows-1;
}
int CuttingPlanesEngine::modelAddCut(int * coefs, double rightHand)
{
    if(maximize)
        return modelAddCut(coefs, rightHand,1);
    else
        return modelAddCut(coefs, rightHand,-1);
}
int CuttingPlanesEngine::modelAddEquality(int * coefs, double rightHand)
{
    return modelAddCut(coefs, rightHand,0);
}

double CuttingPlanesEngine::getObjVal()
{
    return currObj;
}
//set ::lowerBound = INT_MAX if error (not enough time as set by //setTimeoutSolve())
double CuttingPlanesEngine::solve()
{
    double objVal=INT_MIN;
    try{
        double startTmSolve = getCPUTime();
        cout<<"Start solve master..."<<flush;
        d.cplex.solve();
        cout<<"...done master"<<endl;
        tmOnlySolve+=(getCPUTime()-startTmSolve);
        objVal = d.cplex.getObjValue();   //Extract solution
        //CPLOG("Cut Planes objVal="<<objVal<<endl);
        //CPLOG("Status:"<<d.cplex.getCplexStatus()<<endl);
        //if(d.cplex.getCplexStatus() ==IloCplex::AbortTimeLim) CPLOG("Time limit\n");
        if(d.cplex.getCplexStatus() ==IloCplex::AbortTimeLim){
                 cerr<<"\n\n\nATTENTION: time limit "<<timeoutSet<<" secs exceeded in ConstrGener. I will report this and\n\
                       and probably stop afterwords. Did you use setTimeoutSolve()?\n\n";
                 if(maximize){
                    ::upperBound = INT_MIN;
                    currObj      = INT_MIN;
                    return  INT_MIN;
                 }else{
                    ::lowerBound = INT_MAX;
                    currObj      = INT_MAX;
                    return  INT_MAX;
                 }
        }
    }catch (IloCplex::Exception e){
        if(d.cplex.getStatus()==IloAlgorithm::InfeasibleOrUnbounded){
            cerr<<"Attention: InfeasibleOrUnbounded. Last time I had this error, it was because\n \
                   the program you served me was infeasible, it was not possible to find the \n \
                   configurations that do satisfy all set-covering constraints\n\
                   You can think of temporarily adding the 1 1 1 ... 1 columns with cost INT_MAX\n\
                   Later, you can erase it with modelDelCut(i) where i is the return value of \n\
                   the call modelAddCut for the above temporarily column. You do this when \n\
                   this column reaches a dual value of 0, via getCutCoefSol(...)\n\
                   Model saved to trycatch.lp";
            cerr<<"I do not call exit but set currObj=INT_MAX, because this situation might be normal.  You should check this.\n"; 
            d.cplex.exportModel("trycatch.lp");
            if(maximize){
                objVal = INT_MIN;
                ::upperBound = INT_MIN;
                currObj = INT_MIN;
            }else{
                objVal = INT_MAX;
                ::lowerBound = INT_MAX;
                currObj = INT_MAX;
            }
        }
    }
    if(d.cplex.getStatus()==IloAlgorithm::Infeasible){
        d.cplex.exportModel("trycatch.lp");
        cerr<<"ATTENTION on Cutting-Planes: finding infeasible solution.  Model saved to trycatch.lp. \n";
        cerr<<"I do not call exit but set currObj=INT_MAX, because this situation might be normal.  You should check this.\n"; 
        if(maximize){
            objVal = INT_MIN;
            ::upperBound = INT_MIN;
            currObj = INT_MIN;
        }else{
            objVal = INT_MAX;
            ::lowerBound = INT_MAX;
            currObj = INT_MAX;
        }
    }

    IloNumArray solution(d.env);
    d.cplex.getValues (d.vars,solution);
    for(int i=0;i<n;i++){
            primals[i] = solution[i];
            //CPLOG(primals[i]<<",");
            //assert(primals[i]>=d.lb[i]);
            //assert(d.lb[i]==0);
            //use below if you don't set NumericalEmphasis true
            assert(primals[i]>=d.lb[i]-EPS);
            //rounding the primals to the correct values
            //if(primals[i]<d.lb[i]+EPS)primals[i]=d.lb[i];
    }
    //CPLOG("\n");
    
    if(maximize){
        if(objVal<currObj){                    //Update bound
            currObj      = objVal;
            ::upperBound = objVal;
        }
    }else{
        if(objVal>currObj){                    //Update bound
            currObj = objVal;
            ::lowerBound = objVal;
        }
    }
    return currObj;
}//end solve
void CuttingPlanesEngine::getPrimals(double *ystar)
{
    int i;
    if(ystar!=NULL)
        for(i=0;i<n;i++)
            ystar[i] = primals[i];
}

void CuttingPlanesEngine::setPrimals(double *startPrimals, double objValInit)
{
     IloNumArray startValues(d.env);
     IloNumVarArray varCpy(d.env);
     if(maximize)
        upperBound = objValInit;
     else
        lowerBound = objValInit;
     for (int i = 0; i < n; ++i){
             startValues.add(startPrimals[i]);
             varCpy.add(d.vars[i]);
             primals[i] = startPrimals[i];
     }
     try{
        //d.cplex.addMIPStart(varCpy, startValues);
        d.cplex.setStart(startValues,NULL,varCpy,NULL,NULL,NULL);
     }catch (IloCplex::Exception e){
        cerr<<"\n\nConcert Exception in when setting initial values for variables:"<<e.getMessage()<<endl;
        exit(EXIT_FAILURE);
     }
     varCpy.end();
     startValues.end();
}

double CuttingPlanesEngine::coefCutContrib(int coef, int idx)
{
    if(abs(coef)<EPS) //Useless when coef is integer
        return 0;
   return (double)coef * primals[idx]; 
}

double CuttingPlanesEngine::leftHandSum(int * coefs)
{
    double totalSum = 0;
    for (int i=0;i<n;i++)
        totalSum += coefCutContrib(coefs[i],i);
    return totalSum;
}

int CuttingPlanesEngine::violatedCut(int * coefs, double rightHand)
{
    double redCost = leftHandSum(coefs);
    if(maximize)
        return (redCost > rightHand + EPS);
    else
        return (redCost < rightHand - EPS);
}
void CuttingPlanesEngine::exportModel(char* filename)
{
    d.cplex.exportModel(filename);
}
void CuttingPlanesEngine::freeData(double*&a,double*&b, double**&c,int nr)
{
    delete[] a;
    if(internalCutSeprtExtended!=NULL){
        delete[] b;
        int i;
        for(i=0;i<nr;i++)
            delete[] c[i];
        delete[] c;
    }
}

double CuttingPlanesEngine::runSelectedCutSeprt(const int n, double*primals, double * newCut,double&newRightHand,
                                     int it, double tm, double ** newCutMore,double*newRightHandMore,
                                     int&newMore, int maxMoreConstr)
{
      if(internalCutSeprtExtended!=NULL)
           return internalCutSeprtExtended(n,primals,newCut, newRightHand,it,tm,
                                   newCutMore,newRightHandMore,newMore, maxMoreConstr);
      if(internalCutSeprtSolver!=NULL)
           return internalCutSeprtSolver(n,primals,newCut, newRightHand,it,tm);
      return internalCutSeprtSimple(n,primals,newCut, newRightHand);
} 
void CuttingPlanesEngine::setTimeoutSolve(double timeOut)
{
        d.cplex.setParam(IloCplex::ClockType,1);//timeout expressed in CPU time
        d.cplex.setParam(IloCplex::TiLim,timeOut);
        timeoutSet = timeOut;
}

int CuttingPlanesEngine::runCutPlanes(const int itMax, const double tmMax, int& it, double&tm)
{
  try{
       it = 0;
       tm = 0;
       double newRightHand, newViolation,startTm = getCPUTime();
       double* newCut = new double[n];
   
       double** newCutMore = NULL;        //not necessarily used 
       double* newRightHandMore=NULL;     //only for cutSeprtExtended
       int newMore;
       if(internalCutSeprtExtended!=NULL){
           newCutMore = new double*[maxMoreConstr];
           newRightHandMore = new double[maxMoreConstr];
           int iii;
           for(iii=0;iii<maxMoreConstr;iii++)
                newCutMore[iii] = new double[n];
       }

       //if(!d.cplex.isDualFeasible())                        //if unbounded
       //    setVarBounds(-IloInfinity,INT_MAX/((double)n));  //Add some bounds on the variables

       if((maximize) &&(upperBound==INT_MAX))//unless set via setPrimals
          solve();
       if((!maximize)&&(lowerBound==INT_MIN))//unless set via setPrimals
          solve();
       do{
           try{
           newViolation = runSelectedCutSeprt(n,primals,newCut, newRightHand,it,tm,
                                newCutMore,newRightHandMore, newMore, maxMoreConstr);
           if(newViolation==INT_MAX){//gap closed
                 freeData(newCut,newRightHandMore,newCutMore,maxMoreConstr);
                 tm = getCPUTime() - startTm;
                 return EXIT_SUCCESS;
           }
           }catch (IloCplex::Exception e){
                cerr<<"\n\nConcert Exception in your cutSeprt called by the Cutting-Planes:"<<e.getMessage()<<endl;
                cerr<<"I do not call exit, because this might be normal.  You should know better, maybe catch the exception in your cutSeprt.\n\n\n"; 
                if(maximize){
                    ::upperBound = INT_MIN;
                    currObj = INT_MIN;
                }else{
                    ::lowerBound = INT_MAX;
                    currObj = INT_MAX;
                }
                freeData(newCut,newRightHandMore,newCutMore,maxMoreConstr);
                tm = getCPUTime() - startTm;
                return EXIT_FAILURE;
           }
           //max: newViolation = rightHand - a^T x
           //min: newViolation = a^T x   - rightHand
           if(newViolation >= -EPS){        //no violation
                 if(turnIntegerEnd){
                     if(intVars==n){           //all integer => exit
                          freeData(newCut,newRightHandMore,newCutMore,maxMoreConstr);
                          tm = getCPUTime() - startTm;
                          return EXIT_SUCCESS;
                     }
                     turnAllVarsInteger();
                 }else{
                      freeData(newCut,newRightHandMore,newCutMore,maxMoreConstr);
                      tm = getCPUTime() - startTm;
                      return EXIT_SUCCESS;
                 }
           }
           modelAddCut(newCut,newRightHand);
           if(internalCutSeprtExtended!=NULL){
                int ii;
                for(ii=0;ii<newMore;ii++)
                    modelAddCut(newCutMore[ii], newRightHandMore[ii]);
           }
           if(::switchToIntVarsNow)  //the cutSeprt can set this global variable,
                turnAllVarsInteger();//but it can't call turnAllVarsInteger() directly
           solve();
           if(maximize)
           if(::upperBound==INT_MAX){//solve sets ::upperBound=INT_MAX if infeasible
                cerr<<"\n\n\nATTENTION: THE CONSTRAINT/COL GENERATOR LP BECOME INFEASIBLE. STOP HERE\n \
                       without calling exit because this situation might be normal. Check it. \n\n\n";
                break;
           }
           if(!maximize)
           if(::lowerBound==INT_MAX){//solve sets ::lowerBound=INT_MAX if infeasible
                cerr<<"\n\n\nATTENTION: THE CONSTRAINT/COL GENERATOR LP BECOME INFEASIBLE. STOP HERE\n \
                       without calling exit because this situation might be normal. Check it. \n\n\n";
                break;
           }
           //CPLOG("  -> "<<upBound<<endl);
           it ++;
           tm = getCPUTime() - startTm;
           CPLOG("ITER="<<it<<" TIME="<<tm<<" LBND="<<currObj<<"        //message from Cutting Planes engine"<<endl);
       }while((tm<=tmMax)&& it<=itMax);
       freeData(newCut,newRightHandMore,newCutMore,maxMoreConstr);
       return EXIT_FAILURE;                                 //not success
  }
  catch (IloCplex::Exception e){
       d.cplex.exportModel("trycatch.lp");
       cerr<<"\n\nConcert Exception Cutting-Planes:"<<e.getMessage()<<endl;
       cerr<<"Model saved to trycatch.lp. I. I think the problem is infeasible. \n";
       cerr<<"I do not call exit, because this infeasibility/error might be normal.  You should check this.\n"; 
       return EXIT_FAILURE;
  }
}

int CuttingPlanesEngine::runCutPlanes(int& it, double&tm)
{
    return runCutPlanes(INT_MAX, INT_MAX, it, tm);
}
double CuttingPlanesEngine::getTmOnlySolve()
{
    return tmOnlySolve;
}
long CuttingPlanesEngine::getTotalNrCoefs()
{
    return totalNrCoefs;
}




extern "C"{
  /*
   * Author:  David Robert Nadeau
   * Site:    http://NadeauSoftware.com/
   * License: Creative Commons Attribution 3.0 Unported License
   *          http://creativecommons.org/licenses/by/3.0/deed.en_US
   */
  
  #if defined(_WIN32)
  #include <Windows.h>
  
  #elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
  #include <unistd.h>
  #include <sys/resource.h>
  #include <sys/times.h>
  #include <time.h>
  
  #else
  #error "Unable to define getCPUTime( ) for an unknown OS."
  #endif
  
  /**
   * Returns the amount of CPU time used by the current process,
   * in seconds, or -1.0 if an error occurred.
   */
  double getCPUTime( )
  {
  #if defined(_WIN32)
  	/* Windows -------------------------------------------------- */
  	FILETIME createTime;
  	FILETIME exitTime;
  	FILETIME kernelTime;
  	FILETIME userTime;
  	if ( GetProcessTimes( GetCurrentProcess( ),
  		&createTime, &exitTime, &kernelTime, &userTime ) != -1 )
  	{
  		SYSTEMTIME userSystemTime;
  		if ( FileTimeToSystemTime( &userTime, &userSystemTime ) != -1 )
  			return (double)userSystemTime.wHour * 3600.0 +
  				(double)userSystemTime.wMinute * 60.0 +
  				(double)userSystemTime.wSecond +
  				(double)userSystemTime.wMilliseconds / 1000.0;
  	}
  
  #elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
  	/* AIX, BSD, Cygwin, HP-UX, Linux, OSX, and Solaris --------- */
  
  #if _POSIX_TIMERS > 0
  	/* Prefer high-res POSIX timers, when available. */
  	{
  		clockid_t id;
  		struct timespec ts;
  #if _POSIX_CPUTIME > 0
  		/* Clock ids vary by OS.  Query the id, if possible. */
  		if ( clock_getcpuclockid( 0, &id ) == -1 )
  #endif
  #if defined(CLOCK_PROCESS_CPUTIME_ID)
  			/* Use known clock id for AIX, Linux, or Solaris. */
  			id = CLOCK_PROCESS_CPUTIME_ID;
  #elif defined(CLOCK_VIRTUAL)
  			/* Use known clock id for BSD or HP-UX. */
  			id = CLOCK_VIRTUAL;
  #else
  			id = (clockid_t)-1;
  #endif
  		if ( id != (clockid_t)-1 && clock_gettime( id, &ts ) != -1 )
  			return (double)ts.tv_sec +
  				(double)ts.tv_nsec / 1000000000.0;
  	}
  #endif
  
  #if defined(RUSAGE_SELF)
  	{
  		struct rusage rusage;
  		if ( getrusage( RUSAGE_SELF, &rusage ) != -1 )
  			return (double)rusage.ru_utime.tv_sec +
  				(double)rusage.ru_utime.tv_usec / 1000000.0;
  	}
  #endif
  
  #if defined(_SC_CLK_TCK)
  	{
  		const double ticks = (double)sysconf( _SC_CLK_TCK );
  		struct tms tms;
  		if ( times( &tms ) != (clock_t)-1 )
  			return (double)tms.tms_utime / ticks;
  	}
  #endif
  
  #if defined(CLOCKS_PER_SEC)
  	{
  		clock_t cl = clock( );
  		if ( cl != (clock_t)-1 )
  			return (double)cl / (double)CLOCKS_PER_SEC;
  	}
  #endif
  
  #endif
  
  	return -1.0;		/* Failed. */
  }
}
