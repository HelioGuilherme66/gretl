#include <oxstd.h>
#import <packages/dpd/dpd>

main()
{
    decl dpd = new DPD();
    decl vcv, uhat;

    dpd.Load("unbal3.csv");   
    dpd.SetIndex("id"); 
    dpd.SetYear("year"); 
    dpd.Info();

    dpd.Select(Y_VAR, {"y", 0, 1});
    // dpd.Select(X_VAR, {"x", 0, 0});
    // dpd.Select(I_VAR, {"x", 0, 0});

    dpd.Gmm("y", 2, 99); 
    dpd.SetDummies(D_CONSTANT);
    //    dpd.SetOptions(-1,-1,-1,TRUE);
    dpd.Estimate();            // 1-step estimation
    //vcv = dpd.GetCovar();
    //print(vcv);

    dpd.SetDummies(D_NONE);
    dpd.Estimate();

    delete dpd;                    
}
