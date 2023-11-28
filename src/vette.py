import numpy as np
from opm.io.ecl import EGrid
from analytical_effective_perm import damage_zone_permeability
from analytical_effective_perm import damage_zone_width_full
import sys

def isvette(xyz, ijk):
  return xyz[1][0]>6728000 and xyz[1][0]<6750000 and xyz[2][0] < 1350 #and ijk[2] == 0

def istrollold(xyz):
  return xyz[1][0]>6710000 and xyz[1][0]<6722000

def istroll(xyz):
  return xyz[1][0]>6710000 and xyz[1][0]<6734000 and xyz[2][0] > 1500 #1330 #from https://www.norskpetroleum.no/en/facts/field/troll/


def throw(y):
  maxPoint = 6739819.56
  minPoint = 6667544.54
  return 500 * (1.0 - (abs(maxPoint - y) / (maxPoint - minPoint)) )
 
fn = "input/vette.aqucon"
fn2 = "input/vette.multx"
if(len(sys.argv) > 1):
    ind = str(sys.argv[1])
    fn = "input/vette_"+ ind + ".aqucon"
    fn2 = "input/vette_"+ ind + ".multx" 

kfakm = 0.1
if(len(sys.argv) > 2):
    kfakm = float(sys.argv[2])

level = "0"
if(len(sys.argv) > 3):
    level = sys.argv[3]
    
casename ="SMEAHEIA_BOX_"+ level + ".EGRID"
grid = EGrid(casename)

aquiferfile = "aquifer_" + level + ".inc"
boundarytxt = np.genfromtxt(aquiferfile, skip_header=1, skip_footer=1, dtype=(float, float, float, float))
aquifers = []
for txt in boundarytxt:
  aquifers.append(txt[0])
  aquifers.append(txt[1])
  aquifers.append(txt[2])
  aquifers.append(txt[3])

h = open(fn, "w")
h.write("AQUCON\n")
f = open(fn2, "w")
#h.write("-- Index  I1  I2  J1  J2  K1  K2  Face (I/J/K) \n")
#h.write("-- Index  I1  I2  J1  J2  K1  K2  Face (I/J/K) TransMult Trans.opt\n")
# 1   1   1     11   11     1   27      I+     / 
aqunum = 1
connum = 0
#print(len(aquifers))
#volumes = grid.cellvolumes()
#vettevolumes = []
global_index = -1
for aquifer in aquifers:
  global_index +=1

  if (aquifer == 0):
    continue
  
  ijk = grid.ijk_from_global_index(global_index)
  #if (ijk[2] > 0):
  #  continue
  
  xyz = grid.xyz_from_ijk(ijk[0], ijk[1], ijk[2])
  
  transMult = 1.0
  if (kfakm > 0):
    t = throw( xyz[1][0])
    perms = damage_zone_permeability(t, None,kfakm)
    dl = damage_zone_width_full(t)
    X = dl / 200;
    transMult = perms[0] / (perms[0] * (1.0 - X) + X)
 
  if (isvette(xyz, ijk)):
    h.write(str(1) + " ")
  elif (istroll(xyz)):
    h.write(str(2) + " ")
  else:
    continue
    
  h.write(str(ijk[0] +1 ) + " " + str(ijk[0] +1 ) + " " + str(ijk[1] +1 ) + " " + str(ijk[1] +1 ))
  
  #h.write(" 1 27 ")
  h.write(" " + str(ijk[2] +1 ) + " " + str(ijk[2] +1 ) + " ")
  h.write("'I-' " + str(1.0) + " ")
  
  if (isvette(xyz, ijk)):
    h.write(str(0) + " /\n")
  elif (istroll(xyz)):
    h.write(str(1) + " /\n")
  else:
    continue

  if (isvette(xyz, ijk)):
    f.write("BOX\n")
    f.write(str(ijk[0] +1 ) + " " + str(ijk[0] +1 ) + " " + str(ijk[1] +1 ) + " " + str(ijk[1] +1 ))
    #only pick top
    #f.write(" 1 27 ")
    f.write(" " + str(ijk[2] +1 ) + " " + str(ijk[2] +1 ) + " ")
    f.write("/\n\nMULTX- \n "+str(1)+"*" + str(transMult) + "/\n\n")
    f.write("ENDBOX\n\n")

h.write("/\n")
h.close()
