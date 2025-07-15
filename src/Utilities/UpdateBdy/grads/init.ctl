dset ^init.dat
options sequential big_endian
title model output from grapes
undef -9.99E+33
xdef  500 linear  70.0000    0.1500
ydef  330 linear  15.0000    0.1500
zdef   33 linear 1 1
tdef   81 linear 02z04MAR2007   60mn
vars 5
pi   33 99 PI
th   32 99 potential tempertature
u   31 99 U wind component
v   31 99 V wind component
q   32 99 mixing ratio
endvars
