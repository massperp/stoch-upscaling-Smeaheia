import numpy as np
import sys
import pickle
import numpy as np

def krg(So, Sor, bc_exponent):
  So_eff = (So - Sor) / (1.0 - Sor)
  return np.power( 1 - So_eff, 2.0) * (1.0 - np.power( So_eff, ((2.0 + bc_exponent) / bc_exponent)))
  
def krog(So, Sor, bc_exponent):
  if (So < Sor):
    return 0.0
  elif (So == Sor):
    return 0.0001
  
  So_eff = (So - Sor) / (1.0 - Sor)
  return np.power(So_eff, (2.0 + 3.0 * bc_exponent) / bc_exponent)

def pcog(So, pc_entry, Sor, bc_exponent):
  if (So <= Sor):
    return 0.0
  
  So_eff = (So - Sor) / (1.0 - Sor)
  
  return pc_entry * np.power( So_eff, -1.0/bc_exponent)
    

def writeSOF2(f, sor, exponent):
  smin = 0.0
  smax = 0.95 #0.9757 #1.0 - data.get("Sor", 1.0) 
  n = 10
  for s in np.linspace(smin, smax, n):
    kr = krog(s, sor, exponent)
    f.write( ("{:f} {:f} \n").format(s, kr))
    #print(("{:f} {:f} \n").format(s, kr))
  
  f.write("/\n")
  
def writeSGFN(f, sor, maxkrg, pc_entry, exponent):
  smin = 0.0
  smax = 1.0-sor-0.01 #0.9757 #1.0 - data.get("Sor", 1.0) 
  n = 10
  for s in np.linspace(smin, smax, n):
    kr = maxkrg * krg(1.0 - s, sor, exponent)
    pc = pcog(1.0 - s, pc_entry, sor, exponent)
    f.write( ("{:f} {:f} {:f} \n").format(s, kr, pc))
    #print(("{:f} {:f} {:f} \n").format(s, kr, pc))

  f.write("/\n")

def writeSOF2_F(f, S, Kw, swi):
  if S[0] > 0:
    f.write( ("{:f} {:f} \n").format(0.0, 0.0))
  for s, knw in zip(S,Kw):
    seff = s*(1 - swi)+swi
    f.write( ("{:f} {:f} \n").format(s, knw))
  if S[-1] < 1:
    f.write( ("{:f} {:f} \n").format(1.0, 1.0))
  f.write("/\n")


def writeSGFN_F(f, S, Kw, Pc, swi, maxkrg):
  #if S[0] > 0:
  #  f.write( ("{:f} {:f} {:f}\n").format(0.0, 0.0, Pc[len(Pc)-1]/1e5))
  for s, kw, pc in zip(reversed(S),reversed(Kw),reversed(Pc)):
    seff = s*(1 - swi)+swi
    f.write( ("{:f} {:f} {:f} \n").format(1 - seff, maxkrg * kw, pc/1e5))
 
  f.write("/\n")
