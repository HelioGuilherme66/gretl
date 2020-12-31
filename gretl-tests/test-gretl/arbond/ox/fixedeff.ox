#include <oxstd.h>
#import <packages/dpd/dpd>

main()
{
    decl dpd = new DPD();

    dpd.Load("SeatBelts.csv");
    dpd.SetYear("year");
    dpd.SetOptions(FALSE);
    dpd.Select(Y_VAR, {"fatalityrate", 0, 0}); 
    dpd.Select(X_VAR, {"sb_useage", 0, 0,
                       "speed65", 0, 0, 
                       "speed70", 0, 0,
                       "drinkage21", 0, 0,
                       "ba08", 0, 0});

    dpd.SetDummies(D_CONSTANT+D_TIME);
    dpd.SetTransform(T_WITHIN);
    dpd.SetTest(1,2);
    dpd.Estimate(); 

    delete dpd;
}
