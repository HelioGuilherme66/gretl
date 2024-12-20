author = Riccardo "Jack" Lucchetti and Sven Schreiber
email = r.lucchetti@univpm.it
tags = C32
version = @VERSION@
date = 2024-10-20
description = Structural VARs
public = SVAR_setup SVAR_restrict SVAR_namedrestrict \
    SVAR_ident SVAR_estimate \
    SVAR_cumulate SVAR_boot SVAR_hd SVAR_coint \
    GetShock IRFplot IRFsave FEVD \
    GUI_SVAR GUI_plot FEVDplot FEVDsave HDplot HDsave \
    SVAR_bundle_print \
    SVAR_SRplain \
    SVAR_SRexotic SVAR_spagplot SVAR_SRfull \
    SVAR_SRirf SVAR_SRdraw SVAR_SRgetbest SRgetbest \
    SVAR_getshock SVAR_HD \
    SVAR_SRresetalpha
   # SVAR_setidIRF IRF_plotdata 
gui-main = GUI_SVAR

bundle-plot = GUI_plot
bundle-print = SVAR_bundle_print
label = Structural VARs

help = SVAR.pdf

sample-script = examples/Traditional/simple_C.inp
data-files = examples

min-version = @VERSION@
depends = extra
data-requirement = needs-time-series-data

