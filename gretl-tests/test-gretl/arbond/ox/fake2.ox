#include <oxstd.h>
#import <packages/dpd/dpd>

main()
{
    decl dpd = new DPD();
    decl vcv, uhat;

    dpd.Load("fakedata2.csv");   // load data
    dpd.SetYear("year");      // specify columns with years

    dpd.Select(Y_VAR, {"y", 0, 1});     // formulate model
    dpd.Select(X_VAR, {"x", 0, 0});
    dpd.Select(I_VAR, {"x", 0, 0});

    dpd.Gmm("y", 2, 99);           // GMM-type instrument
    dpd.SetDummies(D_NONE);        // specify dummies
    dpd.SetTest(1,2);
    dpd.Estimate();                // 1-step estimation
    vcv = dpd.GetCovar();
    print(vcv);
    uhat = dpd.GetResiduals();
    print(uhat);

    delete dpd;                     // finished with object
}
