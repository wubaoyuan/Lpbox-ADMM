from os import system

for i in range(1,112):
  q = (3-len(str(i)))*'0' + str(i)
  system("python delaunay.py houses/house" + str(i) + " > house" + q + ".adj")
