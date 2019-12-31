import math
import numpy as np
                       
k = 8.98755 * (10**9) #Which is written in Nm^2c^(-2)
KConstantHCC = 5.540 #mdyne / A
KConstant = 5.562
KConstantCC = 4.135
KConstantCCC = 4.14
#KTORSION = LL
#n = # Where this is the amplitude
#P = # Angle that we want for torsion
#D = # Phase shift of the torsion

#Electrostatic is in
#Bond length would be in mdyne * A
# Bond Angle would be the angle multiplying with mdyne*A
# 
#

#Also have to convert everything to the same units

def BondLengthCH(x):
    A1 =  (20 * KConstant * ((x - 1.09))) * (1*10**(-8)) * (1*10**(-10)) * .00024 * 6.022*10**(-23)
    return A1

def DERDERBondLengthCH(x):
    AA1 = 20 * KConstant 
    return AA1

def OGBondLengthCH(x):
    A22 = 10 *  KConstant * ((x - 1.09)**2) 
    return A22

def BondLengthCC1(x):
    CC1 = 6 * KConstantCC * ((x - 1.54)) * (1*10**(-8)) * (1*10**(-10)) * .00024 * 6.022*10**(-23)
    return CC1

def DERDERBondLengthCC1(x):
    ACC = 6 * KConstantCC 
    return ACC
    
def OGBondLengthCC(x):
    F22 = 3 * KConstant * ((x - 1.54)**2)
    return F22

def AngleCCC2(x):
    CCC2 = 4 * KConstantCCC * (x - 110)  #Need to do angles
    return CCC2

def DERDERAngleCCC2(x):
    ACCC2 = 4 * KConstantCCC
    return ACCC2
    
def OGBondAngleCCC(x):
    D22 = 2 * KConstantCCC * ((x - 110)**2)
    return D22

def AngleHCC1(x):
    HCC1 = 16 *  KConstantHCC * (x - 110)
    return HCC1

def DERDERAngleHCC1(x):
    AHCC1 = 16 * KConstantHCC
    return AHCC1
    
def OGBondAngleHCC(x):
    L22 = 8 * KConstant * ((x - 110)**2)
    return L22

def ES143(x):
    QQ143 = (-1*k * (((1.441958814*10**(-20))*(-5.767835256*10**(-20))) / ((x - 1.54)**2) )) * (1*10**(10)) * .00024 * (6.022*10**(-23))
    return QQ143
##All of the Electro. should be in J hence why multiplied by 10^10^-1
def DERDERES143(x):
    AQQ143 = -1 * k * (1.6634 * 10**(-39) / ((x-1.54)**3) ) 
    return AQQ143
    
def OGBondESHC(x):
    MM22 =  k * (((1.441958814*10**(-20))*(-5.767835256*10**(-20))) / (x-1.54) )
    return MM22
###2.52 is equil. distance between the H-H ones assuming 110 bond between them
def ES56(x):
    QQ56 = (-1*k * ((((1.441958814*10**(-20)))*((1.441958814*10**(-20)))) / ((x-2.52)**2) ) ) * (1*10**(10)) * .00024 * (6.022*10**(-23))
    return QQ56

def DERDERES56(x):
    AQQ56 = k * ( (4.15849 * 10 ** (-40)) / ((x-2.52)**3) )  
    return AQQ56
    
def OGBondESHH(x):
    L22 =  k * (((1.441958814*10**(-20))*((1.441958814*10**(-20)))) / (x-2.52) )
    return L22

def VDW510(x):
    VW510 = ((481.14*((x-2.67269)**6) - 350751) / ((x-2.67269)**13))  
    return VW510

def DERDERVDW510(x):
    AVDW510 = ((4559763 - (3367.98*((x-2.67269)**6))) / ((x-2.67269)**14))  
    return AVDW510
    
def OGVDW(x):
    VDW11 = ((29229.255) / (x-2.67269)**12) - ((80.19) / (x - 2.67269)**6)
    return VDW11
    


    
L1 = []
L2 = []
Step = .001
tol = .00001
K = 1000

def SteepestDescent(C, Initial, K, tol): # B is not needed 
    counter = 1
    for i in range(1, K):
        New = Initial - Step*C(Initial)
        if abs(C(Initial)) < tol:
            #print(i, New)
            L1.append(New)
            L2.append(C(New))
            break
        else:
            Initial = New
            L1.append(Initial)
            L2.append(C(New))
            counter += 1
            #print(i, New)

    return Initial
#Change Initial here to whatever fits best for out approximation
BondCH = SteepestDescent(BondLengthCH, 1.00, K, tol)
BondCC = SteepestDescent(BondLengthCC1, 1.45, K, tol)
AngleCCC = SteepestDescent(AngleCCC2, 100, K, tol)
AngleHCC = SteepestDescent(AngleHCC1,  100, K, tol)
ElectroCH = SteepestDescent(ES143,  1.45, K, tol)
ElectroHH = SteepestDescent(ES56,  2.4, K, tol)
VanDerWaal = SteepestDescent(VDW510, 2.4, K, tol)
Torsion = 0 # Because of the fact that we are calculating for Anti... the optimal Torsion Energy is 0
# Where Torsion is zero becuase in optimal optimized postion of 180 degrees it has 0 PE
        
PotentialEnergy = (OGBondLengthCH(BondCH)) + (OGBondLengthCC(BondCC)) + (OGBondAngleCCC(AngleCCC)) + (OGBondAngleHCC(AngleHCC)) + (14*OGBondESHC(ElectroCH)) + (24*OGBondESHH(ElectroHH)) + (6 * OGVDW(VanDerWaal)) + Torsion

CHLengthEn = OGBondLengthCH(BondCH)
HHLengthEn = OGBondLengthCC(BondCC)
CCCAngleEn = OGBondAngleCCC(AngleCCC)
HCCAngleEn = OGBondAngleHCC(AngleHCC)
CHElectrostaticEn = 14*OGBondESHC(ElectroCH)
HHElectrostaticEn = 24*OGBondESHH(ElectroHH)
VDWEn = 6 * OGVDW(VanDerWaal)

print('-------------------------------------------------')

F1 = []
F2 = []
tolZ = .00001
KZ = 400



def NewtonMethod(A, B, InitialZ, KZ, tolZ):
    counter = 1
    for i in range(1, KZ):
        if B(InitialZ) == 0:
            return InitialZ
        NewZ = InitialZ - A(InitialZ) / B(InitialZ) #Where the A is the first derivative and
        if abs(NewZ-InitialZ) < tol: #               The B is teh second Derivative
            #print(i, NewZ)
            F1.append(NewZ)
            F2.append(InitialZ - (A(InitialZ) / B(InitialZ)))
            break
        else:
            InitialZ = NewZ
            F1.append(InitialZ)
            F2.append(A(NewZ))
            counter += 1
            #print(i, InitialZ)

    return NewZ

NEWBondCH = NewtonMethod(BondLengthCH,DERDERBondLengthCH, 1.01, KZ, tolZ)
NEWBondCC = NewtonMethod(BondLengthCC1,DERDERBondLengthCC1, 1.45, KZ, tolZ)
NEWAngleCCC = NewtonMethod(AngleCCC2,DERDERAngleCCC2, 100, KZ, tolZ)
NEWAngleHCC = NewtonMethod(AngleHCC1,DERDERAngleHCC1, 100, KZ, tolZ)
NEWElectroCH = NewtonMethod(ES143,DERDERES143, 1.45, KZ, tolZ)
NEWElectroHH = NewtonMethod(ES56,DERDERES56, 2.4, KZ, tolZ)
NEWVanDerWaal = NewtonMethod(VDW510,DERDERVDW510, 2.4, KZ, tolZ)
NEWTorsion = 0 #STILL NEED TO PUT IN OWN INITIALZ VALUE HERE

NEWTPotentialEnergy = (OGBondLengthCH(NEWBondCH)) + (OGBondLengthCC(NEWBondCC)) + (OGBondAngleCCC(NEWAngleCCC)) + (OGBondAngleHCC(NEWAngleHCC)) + (14*OGBondESHC(NEWElectroCH)) + (24*OGBondESHH(NEWElectroHH)) + (6* OGVDW(NEWVanDerWaal)) + NEWTorsion

print('The minimized Energy with Steepesty Descent will be: ', PotentialEnergy)
print('The minimized Energy with Newtons Method will be: ', NEWTPotentialEnergy)

print('Now for Steepest Descent')

print('The Bond Length CH Energy is:', CHLengthEn) 
print('The Bond Length CC Energy is:', HHLengthEn) 
print('The Bond Angle CCC Energy is:', CCCAngleEn) 
print('The Bond Length HCC Energy is:', HCCAngleEn ) 
print('The Electrostatic CH Energy is:,' , CHElectrostaticEn)  
print('The Electrostatic HH Energy is: ',  HHElectrostaticEn )
print('The  VDW Energy is:',  VDWEn) 

print('--------------------------------------------------------------')
print('Now for Newtons Descent')

print('The Bond Length CH Energy is:', OGBondLengthCH(NEWBondCH)) 
print('The Bond Length CC Energy is:', OGBondLengthCC(NEWBondCC) )
print('The Bond Angle CCC Energy is:', OGBondAngleCCC(NEWAngleCCC))
print('The Bond Length HCC Energy is:', OGBondAngleHCC(NEWAngleHCC) ) 
print('The Electrostatic CH Energy is:,' , 14*OGBondESHC(NEWElectroCH))  
print('The Electrostatic HH Energy is: ', 24*OGBondESHH(NEWElectroHH) )
print('The  VDW Energy is:',  6* OGVDW(NEWVanDerWaal)) 

