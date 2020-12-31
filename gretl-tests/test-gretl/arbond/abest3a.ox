#include <oxstd.h>
#import <packages/dpd/dpd>

main()
{
    decl dpd = new DPD();

    dpd.Load("abdata.in7");
    dpd.SetYear("YEAR");
    dpd.SetGroup("IND");
	
    //----------------- column (c) ------------------------
    dpd.SetYear("YEAR");      // specify columns with years
    dpd.Select(Y_VAR, {"n", 0, 2});      // formulate model
    dpd.Select(X_VAR, {"w", 0, 1, "k", 0, 0, "ys", 0, 1});
    dpd.Select(I_VAR, {"ys", 0, 1});
    dpd.SetDummies(D_CONSTANT + D_TIME); // specify dummies

    dpd.Gmm("n", 2, 99);             // GMM-type instrument
    dpd.Gmm("w", 2, 3);              // GMM-type instrument
    dpd.Gmm("k", 2, 3);              // GMM-type instrument

    print("\n***** Arellano & Bond (1991), Table 4 (c)\n");
    print("        (but using different instruments!!)\n");
    dpd.SetMethod(M_2STEP); // changed here for now
    dpd.Estimate();                    // 2-step estimation

    delete dpd;                      // finished with object
}
