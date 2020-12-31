#include <oxstd.h>
#import <packages/dpd/dpd>

main()
{
  decl dpd = new DPD();
  dpd.Load("abdata.in7");
  dpd.SetYear("YEAR");

  // model-specific code here
  dpd.Select(Y_VAR, {"n", 0, 1});
  dpd.Select(X_VAR, {"w", 0, 1, "k", 0, 1});
  dpd.SetDummies(D_CONSTANT + D_TIME);

  print("\n\n***** Blundell & Bond (1998), Table 4: 1976-86 GMM-DIF");
  dpd.Gmm("n", 2, 99);
  dpd.Gmm("w", 2, 99);
  dpd.Gmm("k", 2, 99);
  dpd.SetMethod(M_2STEP);
  dpd.Estimate();

  print("\n\n***** Blundell & Bond (1998), Table 4: 1976-86 GMM-SYS");
  dpd.GmmLevel("n", 1, 1);
  dpd.GmmLevel("w", 1, 1);
  dpd.GmmLevel("k", 1, 1);
  dpd.SetMethod(M_2STEP);
  dpd.Estimate();

  delete dpd;
}
