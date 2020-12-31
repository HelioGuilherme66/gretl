options double crt nwidth=17,signif=9;
name penngrow 'Penn World Table data - growth model';

? Dynamic Panel Model, balanced data:
? OLS, Fixed Effects, Random Effects (FGLS and Conditional ML), Between,
? and Anderson-Hsiao 2SLS (AH-L and AH-D).
? Follows:
? Nerlove, Marc, "Properties of Alternative Estimators of
?   Dynamic Panel Models:  An Empirical Analysis of Cross-Country
?   Data for the Study of Economic Growth", in C. Hsiao, K. Lahiri,
?   L.-F. Lee, and M. H. Pesaran, eds, "Analysis of Panels and
?   Limited Dependent Variable Models," Cambridge University Press,
?   August 1999.  (348p., 0-521-63169-6)
? See also:
? Nerlove, Marc, "Likelihood Inference for Dynamic Panel Models,"
?   Annales d'Economie et Statistiques, vol. 54, 2000.
?   http://www.ensae.fr/annalesEnglish/
? Nerlove, Marc, "Likelihood Inference in Econometrics,"
?   Academic Press, forthcoming 2000.
? Data used:  Penn World Table, 1960-85:
?   y = log(real per capita GDP)
?   x = log(savings) - log(population growth)  (5 year averages)
? by Clint Cummins 10/1999

freq(panel,n=94,t=6,id=country,time=year) n;
set nt = 94*6;
smpl 1,nt;
?read(file='penngrow.txt') country year y x;
read(file='penngrow.wks');

title '94 country sample';
select year>1960;

? Benchmarks:
?---------
? 1. Fixed Effects OLS / Within
?           Nerlove(1999)      TSP
? Y(-1)        .7204       .7203687510
? (SEasy)     (.0211)     (.0211205374)    asymptotic SE (comparable to ML SEs)
? (SEdf)                  (.0236765322)    SE with df correction
? X            .1656       .1655889301
? (SEasy)     (.0172)     (.0172453860)
? (SEdf)                  (.0193324123)
? SSR/NT       .0113       .0112885550
? S2                       .0141861520
?---------
? 2. Feasible GLS / Random Effects / Variance Components (one-way)
?           Nerlove(1999)      TSP
? RHO used     .2675       .268312999      (checking method used in TSP...)
? Y(-1)        .9130       .9129384943
? (SEasy)     (.0119)     (.0113539782)
? X            .1520       .1520979768
? (SEasy)     (.0135)     (.0128574848)
? SSR/NT       .0213
? S2                       .0200371810
?---------
? 3. Pooled OLS / Total
?           Nerlove(1999)      TSP
? Y(-1)        .9487       .9487096052
? (SEasy)     (.0090)
? (SEdf)                  (.009071038586)
? X            .1244       .1243758682
? (SEasy)     (.0108)
? (SEdf)                  (.0107937856)
? SSR/NT       .0193
? S2                       .0193882794
?---------
? 4. Country Means OLS / Between
?           Nerlove(1999)      TSP
? Y(-1)        .9817       .9817017898
? (SEasy)     (.0112)
? (SEdf)                  (.0113862874)
? X            .0919       .0918577087
? (SEasy)     (.0140)
? (SEdf)                  (.0140213718)
? SSR/NT       .0047
? S2                       .004872144122
?---------
? 5. Conditional ML - Random Effects, y0 fixed
?           Nerlove(1999)      TSP
? RHO          .1133       .1133500346
? (SEhess)    (.0497)     (.0498165772)   SE from analytic second derivatives
? Y(-1)        .9339       .9338582881
? (SEhess)    (.0122)     (.0122508692)
? X            .1370       .1369797209
? (SEhess)    (.0131)     (.0131629846)
? SIGMA2       .0194       .0193751264
? (SEhess)    (.0013)     (.001308721996)
? LogLike   264.927     264.9274773

panel y y(-1) x;

title 'fixed effects SEs without DF correction';
mat vcov = @vcovw*(@nob-94-2)/@nob;
tstats(names=@rnmsw) @coefw vcov;
set resvar = @ssrw/@nob;
print resvar;

title 'intraclass correlation coefficient used';
set vbet = @s2t - @s2w;
set rho = vbet/(vbet + @s2w);
print rho;

title 'conditional ML random effects - y0 fixed';
?--------------------------------
? starting values from pooled OLS
unmake @coeft gamma beta alpha;
set sigma2 = @s2t;
param gamma beta alpha sigma2 rho,.5;

? replicate individual means across all time periods
select year=1960;
unmake @mean yibar yLibar xibar;
select year>1960;
yL = y(-1);
dot y yL x;
   .ibar = .ibar(-1);
   frml .star (. - .ibar)/sqrt(eta) + .ibar/sqrt(xi);
enddot;

const t,5;
frml condml logl = -log(sigma2)/2 + lnorm(e/sqrt(sigma2))
                   -log(xi)/(2*T) - log(eta)*((T-1)/(2*T));
frml xi  1 + rho*T - rho;
frml eta 1 - rho;
frml e ystar - gamma*yLstar - beta*xstar - alpha/sqrt(xi);
eqsub condml e ystar yLstar xstar xi eta;
ML(hiter=n,hcov=nbw,tol=1e-11) condml;
?----------------------------------

title '22 country sample';
? Reproduces Nerlove(1999), same as the 94 country sample
set nt22 = 22*6;
smpl 1,nt22;
select year>1960;
panel y y(-1) x;

smpl 1,nt;

title 'Anderson-Hsiao estimates';
? Proposed Benchmarks (not in Nerlove(1999)):
?---------
? 6. AH-L / using Y(-2) as instrument for DY(-1)
?                TSP
? Y(-1)      .6869656513
? (SEdf)    (.0644896228)
? X          .1484237304
? (SEdf)    (.0211451444)
? Num.obs.       376
?---------
? 7. AH-D / using DY(-2) as instrument for DY(-1)
?                TSP
? Y(-1)      .5843594276
? (SEdf)    (.0717878873)
? X          .1347469433
? (SEdf)    (.0218479837)
? Num.obs.       282

? Create first differenced variables
select year > 1960;
dot y x;
  d. = . - .(-1);
enddot;

select year > 1965;
title 'AH-L';
2sls(inst=( y(-2) dx c)) dy dy(-1) dx;

select year > 1970;
title 'AH-D';
2sls(inst=(dy(-2) dx c)) dy dy(-1) dx;
