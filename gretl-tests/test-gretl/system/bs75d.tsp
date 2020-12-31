options double crt nwidth=14,signif=7;
name bs75d 'Berndt and Savin (1975) - 3 equation system';
? Berndt and Savin 1975, 3 equation singular model.
? Reproduces results for 3 equations with cross-equation restrictions.
? Also demonstrates invariance to which equation is dropped.
? by Clint Cummins, 10/1995
? revised to add comparison table 9/2004
load;
smpl 30,71;   ? drop first observation
?
title 'unrestricted equations';
frml u1 y1 = a1 + g11*lx1 + g12*lx2 + g13*lx3;
frml u2 y2 = a2 + g21*lx1 + g22*lx2 + g23*lx3;
frml u3 y3 = (1-a1-a2) + (-g11-g21)*lx1 + (-g12-g22)*lx2 + (-g13-g23)*lx3;
param a1-a3 g11-g13 g21-g23 g31-g33;
lsq u1 u2;
lsq(wname=own) u1 u3;  ? demonstrate invariance to dropped eqn
lsq(wname=own) u2 u3;
copy @logl loglu;   ? save LogL value for LR test
frml mu1 a3 = (1-a1-a2);
frml mu2 g31 = -g11-g21;
frml mu3 g32 = -g12-g22;
frml mu4 g33 = -g13-g23;
analyz mu1-mu4;
?
? Benchmark
? Berndt & Savin "Estimation and Hypothesis Testing in Singular Equation
?  Systems with Autoregressive Disturbances", Econometrica 43, 1975
?  Table I, Model 8, p.952  (model with no autoregressive disturbances)
?  Table II, Model 8, p.954  (LR test)
?
?         Berndt &    TSP 4.5
?         Savin
? G11      .220      .21994420 
?(SE)     (.018)    (.01770392)
? G12     -.148     -.14768561
?         (.017)    (.01663795)
? G13     -.072     -.07225858
?         (.012)    (.01160467)
? G22      .050      .04968651
?         (.025)    (.02444804)
? G23      .098      .09799911
?         (.013)    (.01247807)
? G33     -.026     -.02574052
?         (.011)    (.01069897)
? A1                -.01881157
?                   (.06733621)
? A2                 .89265980
?                   (.05485714)
? A3                 .12615177
?                   (.05793081)
? LogL  289.5129  211.02127
? LR     47.36     47.3594284
? Notes:
? 1. All coefficients and SEs are reproduced to the given 2-3 digits,
?    except the SE for G22, which is good to 1.9 digits.
? 2. The TSP coefficients and SEs should be good for the 7-8 significant
?    digits shown.
? 3. Results are invariant which equation is dropped.
? 4. The coefficients and SEs for G13 G23 G33 are obtained by the delta
?    method.
? 5. The LogL reported by Berndt & Savin is simply off by an additive
?    constant, since their LR test statistic is reproduced.
? 6. Berndt designed the algorithm which TSP's LSQ command uses for
?    estimation, and he used an early TSP version (probably 2.3, which
?    did calculations in single precision) on the IBM 360/70 (p.949).
?    The results for nonlinear estimation in 1974 are holding up very well.
?
title 'restricted equations - symmetry';
frml r1 y1 = a1 + g11*lx1 + g12*lx2 + (-g11-g12)*lx3;
frml r2 y2 = a2 + g12*lx1 + g22*lx2 + (-g12-g22)*lx3;
frml r3 y3 = (1-a1-a2) + (-g11-g12)*lx1 + (-g12-g22)*lx2 + (g11+2*g12+g22)*lx3;
lsq(tol=1e-6,wname=own,maxit=50) r1 r2;
lsq(tol=1e-6,wname=own) r1 r3;  ? demonstrate invariance to dropped eqn
lsq(tol=1e-6,wname=own) r2 r3;
frml mr2 g31 = -g11-g12;
frml mr4 g33 = g11+2*g12+g22;
analyz mu1 mr2 mu3 mr4;  ? estimates and SEs from delta method
set lrstat = 2*(loglu-@logl);
print lrstat;

title 'check results using FIML';
? The coefficients and LogL are the same, and the SEs differ somewhat
? (as expected).
fiml(endog=(y1,y2),tol=1e-6,hcov=c) r1 r2;
end;
?--------------------------------------------------------------
freq a;
smpl 29,71;  ? data from Berndt + Savin, Econometrica 1975, p.956
noprint;
read y1 y2 lx1 lx2 lx3;
.554779 .266051 9.87415 8.60371 3.75673
.516336 .297473 9.68621 8.50492 3.80692
.508201 .324199 9.47010 8.27596 3.80437
.500482 .339580 9.24692 7.99221 3.76251
.501695 .276974 9.35356 7.98968 3.69071
.591426 .287008 9.42084 8.03620 3.63564
.565047 .244096 9.53972 8.15803 3.59285
.605954 .239187 9.69140 8.32009 3.56678
.620161 .218232 9.76817 8.42001 3.57381
.592091 .243161 9.51295 8.19771 3.60240
.613115 .217042 9.68047 8.30987 3.58147
.632455 .215269 9.78417 8.49624 3.57744
.663139 .184409 10.0606 8.69868 3.60950
.698179 .164348 10.2892 8.84523 3.66664
.704590 .146865 10.4731 8.93024 3.65388
.694067 .161722 10.4465 8.96044 3.62434
.674668 .197231 10.2790 8.82522 3.61489
.635916 .204232 10.1536 8.77547 3.67562
.642855 .187224 10.2053 8.77481 3.82632
.641063 .186566 10.2227 8.83821 3.96038
.646317 .203646 10.1127 8.82364 4.04470
.645476 .187497 10.2067 8.84161 4.08128
.638030 .197361 10.2773 8.94010 4.10951
.634626 .209992 10.2830 9.01603 4.16930
.631144 .219287 10.3217 9.06317 4.21727
.593088 .235335 10.2101 9.05664 4.25670
.607360 .227035 10.2720 9.07566 4.29193
.607204 .246631 10.2743 9.12407 4.32252
.586994 .256784 10.2396 9.15880 4.37792
.548281 .271022 10.1248 9.14025 4.42641
.553401 .261815 10.2012 9.15980 4.43460
.552105 .275137 10.1846 9.19297 4.43173
.544133 .280783 10.1479 9.19533 4.44407
.553820 .281286 10.1970 9.21544 4.45074
.549951 .283030 10.2036 9.22841 4.46403
.547204 .291287 10.2271 9.23954 4.48447
.555110 .281313 10.2882 9.26531 4.52057
.558182 .280151 10.3530 9.31675 4.58156
.545735 .294385 10.3351 9.35382 4.65983
.538964 .294593 10.3525 9.38361 4.71804
.542764 .299927 10.3676 9.40725 4.76329
.534595 .315319 10.2968 9.39139 4.81136
.545591 .315828 10.2592 9.34121 4.84082
;
y3 = 1 - y1 - y2;
