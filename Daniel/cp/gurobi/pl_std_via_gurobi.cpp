#include "gurobi_c++.h"
using namespace std;

int main()
{
  try {
    int n = 3;
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    GRBVar* vars = new GRBVar[n];
    for(int i=0;i<n;i++)
        vars[i] = model.addVar(0,GRB_INFINITY, 0, GRB_CONTINUOUS, "x");

    GRBLinExpr expr = 0;
    for(int i=0;i<n;i++)
        expr += vars[i];
    model.setObjective(expr, GRB_MINIMIZE);


    GRBLinExpr constr = 0;
    for(int i=0;i<n;i++)
        constr += vars[i];
    model.addConstr(constr >= 4, "c0");

    model.optimize();

    model.addConstr(constr >= 5, "c0");

    model.optimize();

    if (model.get(GRB_IntAttr_SolCount) > 0) {
        for (int i = 0; i < n; i++){
            double s;
            s = vars[i].get(GRB_DoubleAttr_X);
            cout<<"vars["<<i<<"]="<<s<<endl;
        }
    }
    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}
