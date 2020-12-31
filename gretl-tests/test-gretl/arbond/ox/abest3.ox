#include <oxstd.h>
#import <packages/dpd/dpd>

main()
{
    decl dpd = new DPD();

    dpd.Load("abdata.in7");                   // load data
    dpd.SetYear("YEAR");     // specify columns with years

    dpd.Select(Y_VAR, {"n", 0, 2});     // formulate model
    dpd.Select(X_VAR, {"w", 0, 1, "k", 0, 0, "ys", 0, 1});
    dpd.Select(I_VAR, {"w", 0, 1, "k", 0, 0, "ys", 0, 1});

    dpd.Gmm("n", 2, 3);            // GMM-type instrument
    dpd.SetDummies(D_CONSTANT + D_TIME);// specify dummies
    dpd.SetTest(1, 2);// specification,Sargan,AR 1-2 tests
    dpd.Estimate();                   // 1-step estimation

    delete dpd;                     // finished with object
}
