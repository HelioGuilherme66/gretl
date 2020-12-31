#include <oxstd.h>
#import <packages/dpd/dpd>

main()
{
    decl dpd = new DPD();

    dpd.Load("alessia.csv");
    dpd.SetYear("year");
    dpd.SetGroup("panelid");
	
    dpd.Select(Y_VAR, {"lprod", 0, 1});      
    dpd.Select(X_VAR, {"lemp", 0, 0, "lva", 0, 0});
    dpd.Select(I_VAR, {"lemp", 0, 0, "lva", 0, 0});
    dpd.SetDummies(D_CONSTANT);
    dpd.Gmm("lprod", 2, 10);
    dpd.SetMethod(M_2STEP);
	dpd.SetOptions(FALSE); // asymptotic std errors
    dpd.Estimate(); 

    delete dpd; 
}
