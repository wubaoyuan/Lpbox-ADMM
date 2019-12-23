import sys
from os import system


for baseline in [0,10,20,30,40,50,60,70,80,90]:
  outfile = "pairs/houses.txt"
  outtrain = "pairs/houses" + str(baseline) + "_train.txt"
  outvalidate = "pairs/houses" + str(baseline) + "_validate.txt"
  outtest = "pairs/houses" + str(baseline) + "_test.txt"

  f = open(outfile, 'w')

  fs = [open(outtrain, 'w'),open(outvalidate, 'w'),open(outtest, 'w')]

  def tofile(f, i):
    st = "0"*(3-len(str(i))) + str(i)
    corner = "../Data/house/houses/house" + str(i)
    feature = "../Data/house/house" + st + ".scf"
    image = "../Data/house/images/house.seq" + str(i-1) + ".png"
    adj = "../Data/house/house" + st + ".adj"
    f.write(corner + " " + feature + " " + image + " " + adj + "\n")


  for i in range(1,112):
    #system("wget http://www.vasc.ri.cmu.edu//idb/images/motion/house/house.seq"+str(i-1)+ ".png")
    tofile(f,i)
    for j in range(i,112):
      if (j-i == baseline):
        tofile(fs[(i-1) % 3], i)
        tofile(fs[(i-1) % 3], j)

  f.close()
  [f1.close() for f1 in fs]
