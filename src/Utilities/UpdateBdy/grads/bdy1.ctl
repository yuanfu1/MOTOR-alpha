dset ^bdy1.dat
options sequential big_endian
title model output from grapes
undef -9.99E+33
xdef  1201 linear    70.0000    0.1500
ydef  6 linear  15.0000    0.1500
zdef   67 linear 1 1
tdef   1 linear 02z04MAR2007   60mn
vars 10
pi   67 99 pressure
pt   67 99 pressure
th   67 99 perturbed PI
tht   67 99 perturbed potential tempertature
u   67 99 U wind component
ut   67 99 U wind component
v   67 99 V wind component
vt   67 99 V wind component
q   67 99 terrain following coordinate W wind component 
qt   67 99 terrain following coordinate W wind component 
endvars
