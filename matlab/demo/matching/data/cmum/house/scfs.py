from os import system
import sys

for q in range(1,112):
  q2 = (3-len(str(q)))*'0' + str(q)
  system("./scfeatures none house/houses/house" + str(q) + " house/house" + q2 + ".scf")
