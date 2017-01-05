from sympy import *
init_printing(use_unicode=True)
from sympy.physics.quantum import TensorProduct

Ein_s,Ein_1,Ein_2,Ein_d = symbols('psi^in_s psi^in_1 psi^in_2 psi^in_d')
Eout_s,Eout_1,Eout_2,Eout_d = symbols('psi^out_s psi^out_1 psi^out_2 psi^out_d')
Uss,Us1,Us2,Usd,U1s,U11,U12,U1d,U2s,U21,U22,U2d,Uds,Ud1,Ud2,Udd=symbols('U_ss U_s1 U_s2 U_sd U_1s U_11 U_12 U_1d U_2s U_21 U_22 U_2d U_ds U_d1 U_d2 U_dd')

Ein = Matrix([[Ein_s],[Ein_1],[Ein_2],[Ein_d]])
Eout = Matrix([[Eout_s],[Eout_1],[Eout_2],[Eout_d]])
U = Matrix([[Uss,Us1,Us2,Usd],[U1s,U11,U12,U1d],[U2s,U21,U22,U2d],[Uds,Ud1,Ud2,Udd]])

Ein_ideal = Matrix([[Ein_s],[0],[0],[0]])



