#include <oxstd.h>
#import <packages/dpd/dpd>

main()
{
    decl dpd = new DPD();

    dpd.Load("ardatamiss.csv");   // load data
    dpd.SetYear("year");      // specify columns with years

    dpd.Select(Y_VAR, {"y", 0, 1});     // formulate model

    dpd.Gmm("y", 2, 99);           // GMM-type instrument
    dpd.SetDummies(D_NONE);        // specify dummies
    dpd.SetTest(1,2);
    dpd.Estimate();                // 1-step estimation

    delete dpd;                     // finished with object
}
