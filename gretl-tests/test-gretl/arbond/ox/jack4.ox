#include <oxstd.h>
#import <packages/dpd/dpd>

main()
{
    decl dpd = new DPD();
    decl vcv, uhat;

    dpd.Load("jack4.csv");   
    dpd.SetIndex("id"); 
    dpd.SetYear("year"); 
    dpd.Info();

    dpd.Select(Y_VAR, {"y", 0, 1});

    dpd.Gmm("y", 2, 99); 
    //dpd.SetDummies(D_CONSTANT);
    //    dpd.SetOptions(-1,-1,-1,TRUE);
    //dpd.Estimate();            // 1-step estimation
    //vcv = dpd.GetCovar();
    //print(vcv);

    dpd.SetDummies(D_NONE);
    dpd.Estimate();
    uhat = dpd.GetResiduals();
    print("resids", uhat);

    delete dpd;                    
}
