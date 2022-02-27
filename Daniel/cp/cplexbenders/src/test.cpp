#include "CuttingPlanesEngine.h"                               
#include <iostream>
#include <cstdlib>
#include <ilcplex/ilocplex.h>                    
#define EPSILON 1.0e-6
#define COUT(x) {}
#define COUT2(x) {cout<<x;}
using namespace std;

////int n=6;                     //number of vars including w
////int f[] = {7, 2, 2, 7, 7, 1};//last entry is the coef of w
////int c[] = {66, 5, 4, 3, 2};
////int d   = 3;
int n,d;
int *f,*c;

void loadData(int nn, int benders){
    n=nn;
    d=n/2;
    f = new int[n+1];
    c = new int[n];
    f[0] = 7;
    c[0] = 8;
    for(int i=1;i<n;i++){
        f[i] = (f[i-1]*f[0])%159;
        c[i] = (c[i-1]*c[0])%61;
    }
    if(benders){
        f[n] = 1;
        n++;
    }
}
//étant donnée la solution courante y, trouver une contrainte violée
//a^T y >= d et renvoyer a^T y - d. Remplir vecteur a et valeur d!!
double separateFunc (const int n, double*y, double * a, double&rHand)
{ 
    cout<<"Starting sep function 1...";
    cout.flush();
    COUT("Cutting y=");
    for(int i=0;i<n;i++)
        COUT(y[i]<<" ");
    COUT(endl);
    IloEnv env;               //environnement toujours nécessaire
    IloModel model(env);      //le modèle contient les contraintes+l'objectif
    IloCplex cplex (model);   //L'objet Cplex s'occupe de l'algorithmique

    IloNumVarArray v(env, n, 0, INT_MAX);    
    IloExpr obj(env);
    obj =  d*v[n-1];            //v[n-1] is d in the slides
    for(int i=0;i<n-1;i++)
        obj+= -y[i] * v[i];
    IloObjective ilo_obj(env, obj, IloObjective::Maximize);
    model.add(ilo_obj);

    //ajouter le code ci-dessous pour demander à cplex de minimiser
    //les valeurs des v lorsqu'il y a plusieurs solutions v optimales
    //// IloNumArray newvals(env,n);
    //// for(int i=0;i<n-1;i++){
    ////     if(y[i]<0.1)
    ////         newvals[i] = -1;
    ////     else
    ////         newvals[i] = 0;
    //// }
    //// newvals[n-1]=0;
    //// ilo_obj.setLinearCoefs(v,newvals);
    //// obj.end();
    //// model.add(ilo_obj);
    //// cout<<"solving again.."<<flush;
    //// cplex.solve();
    //// cout<<"...done"<<endl;
    //// ilo_obj.end();
    //// newvals.end();
    
    //renvoyer a^T y - d (négatif si a^T y >= d est violée)
}
int main(int argc, char**  argv)
{
    if(argc==1){
        cerr<<"Give me a size"<<endl;
        return 0;
    }
    int nn = atoi(argv[1]);
    loadData(nn,1);
    int iters;               
    double finalObj, time; 

    //construire l'objet de plans coupants avec fonc de séparation
    CuttingPlanesEngine cutPlanes(n,separateFunc);      
    cutPlanes.setObjCoefsMinimize(f);
    cutPlanes.turnAllVarsInteger();         
    //cutPlanes.activateLog();              //voir messages de log

    int * y0_greater_y1 = new int[n];
    int * y0_greater_y2 = new int[n];
    int *        sum_y  = new int[n];
    cout<<"n="<<n<<endl;
    for(int i=0;i<n;i++){
        y0_greater_y1[i] = 0;
        y0_greater_y2[i] = 0;
               sum_y [i] = 1;
    }
    y0_greater_y1[0]   = 1;
    y0_greater_y2[0]   = 1;
    y0_greater_y1[1]   = -1;
    y0_greater_y2[2]   = -1;
           sum_y [n-1] = 0;
    
    cutPlanes.modelAddCut(y0_greater_y1, 0);
    cutPlanes.modelAddCut(y0_greater_y2, 0);
    cutPlanes.modelAddCut(       sum_y , d);
    //fix =d for(int i=0;i<n-1;i++) sum_y [i] = - 1; cutPlanes.modelAddCut(sum_y , -d);
    cout<<"new cuts..."<<flush;
    cutPlanes.setVarBounds(0,n-2,0,1);
    cout<<"done";

    ////cout<<"Affectation costs c=";
    ////for(int i=0;i<n-1;i++)
    ////    cout<<c[i]<<" ";
    ////cout<<endl;
    //Faire tourner les plans coupants
    cutPlanes.runCutPlanes(iters, time);    //renvoie 0 si échec



    finalObj = cutPlanes.getObjVal();

    //Imprimer la valeur objectif finale 
    cout<<"\nObj="<<finalObj<<", after "<<iters<<" iters.\n";
    double*y = new double[n]; 
    cutPlanes.getPrimals(y); 

    cout<<"Open facilities   :";

    for(int i=0;i<((n-1<60)?n-1:60);i++)
        if(y[i]>1.0e-6)
                cout<<1;
            else
                cout<<0;
    cout<<endl; 
    exit(1);
    //cout<<y[n-1];
    //cout<<"d="<<d<<endl;
    //for(int i=0;i<n-1;i++)
    //    if(y[i]>EPSILON)
    //        cout<<f[i]<<"/"<<c[i]<<endl;
}
