#universal constants
h = 6.62607015*(10**-34)
hbar = 1.054571816*(10**-34)
eV = 1.602176634*(10**-19) #eV to J
eps0 = 8.854187817*(10**-12)
m0 = 9.10938356*(10**-31)
Ry = 13.605693009 #ev
coulomb = 1.4399 #evnm

#WSe2 constants
class WSe2:
    a = .3286 #nm
    me = 0.33 #eV .15-.25
    mh = 0.45
    M = me+mh #eV, 
    Eb = 0.45 #eV, .3-.6 depending on env
    r0 = 4.5 #nm, 4-5
    dielec = 4.5 #1-4.5 based on env, tune further
    Lm = 10 #nm, mismatch, test 5-30
    V0 = .008 #eV
    psi = 120
    eoff = 0 #V
    hoff = 0 #V

class WS2: #UPDATE ALL PARAMETERS
    a = .3157 #nm
    me = 0.27 #eV .15-.25
    mh = .35
    M = me+mh #eV, 
    Eb = 0.45 #eV, .3-.6 depending on env
    r0 = 4.5 #nm, 4-5
    dielec = 4.5 #1-4.5 based on env, tune further
    Lm = 10 #nm, mismatch, test 5-30
    V0 = .008 #eV
    psi = 0
    eoff = -.3 #V
    hoff = -.4 #V
'''
class MoS2:
    a= 3.16
    me= .35*m0
    M= .44*m0
    Eb= 0.44
    r0= 4.2
    dielec= 2.5
    Lm= 10
    V0= 30

'''
