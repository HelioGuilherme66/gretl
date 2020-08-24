author = Riccardo "Jack" Lucchetti and Sven Schreiber
email = r.lucchetti@univpm.it
tags = C32
version = 1.92
date = 2020-08-24
description = Structural VARs
public = SVAR_setup SVAR_restrict SVAR_ident SVAR_estimate \
    SVAR_cumulate SVAR_boot SVAR_hd SVAR_coint \
    GetShock IRFplot IRFsave FEVD \
    GUI_SVAR GUI_plot FEVDplot FEVDsave HDplot HDsave \
    SVAR_bundle_print \
    SVAR_SRplain IRF_plotdata \
    SVAR_SRexotic SVAR_spagplot SVAR_SRfull \
    SVAR_SRirf SVAR_SRdraw SRgetbest \
    SVAR_getshock SVAR_HD
   # SVAR_setidIRF
gui-main = GUI_SVAR

bundle-plot = GUI_plot
bundle-print = SVAR_bundle_print
label = Structural VARs

help = SVAR.pdf

sample-script = examples/simple_C.inp
min-version = 2017a
data-requirement = needs-time-series-data
data-files = examples
