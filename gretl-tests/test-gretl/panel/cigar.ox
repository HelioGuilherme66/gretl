#include <oxstd.h>
#import <packages/dpd/dpd>

main()
{
    decl dpd = new DPD();

    dpd.Load("cigar.csv");

    dpd.SetYear("year");
    dpd.SetGroup("state");
    dpd.SetTest(1, 2); // Sargan, AR 1-2 tests

    dpd.Select(Y_VAR, {"lnc", 0, 1});
    dpd.Select(X_VAR, {"lnp", 0, 0, "lnpn", 0, 0, "lny", 0, 0});
    dpd.Select(I_VAR, {"lnp", 0, 0, "lnpn", 0, 0, "lny", 0, 0});
    dpd.Gmm("lnc", 2, 99); 
    dpd.SetDummies(D_CONSTANT + D_TIME);     
    dpd.SetMethod(M_2STEP);

    dpd.Estimate();
    
    delete dpd; 
}
