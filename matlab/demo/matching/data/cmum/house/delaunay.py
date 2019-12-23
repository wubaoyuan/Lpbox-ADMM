from __future__ import division
from math import *
from os import system
from sets import Set
import sys
import numpy
import random

N = 30

mat = []
for i in range(N):
  mat.append([0]*N)

fname = sys.argv[1]

a = open(fname, 'r').readlines()
verts = [tuple([float(c) for c in b.split()]) for b in a]
verts = verts[:N]

triangles = Set()
edges = dict()

def dist((x1,y1),(x2,y2)):
  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))

def det(m):
  return numpy.linalg.det(m)

def dictadd(d, a):
  if (d.has_key(a)): d[a] += 1
  else: d[a] = 1


def centrerad((Ax, Ay), (Bx, By), (Cx, Cy)):
  A = (Ax, Ay)
  B = (Bx, By)
  C = (Cy, Cy)
  modA = Ax*Ax + Ay*Ay
  modB = Bx*Bx + By*By
  modC = Cx*Cx + Cy*Cy
  Sx = 0.5*det(numpy.array([[modA, Ay, 1], [modB, By, 1], [modC, Cy, 1]]))
  Sy = 0.5*det(numpy.array([[Ax, modA, 1], [Bx, modB, 1], [Cx, modC, 1]]))
  a = det(numpy.array([[Ax, Ay, 1], [Bx, By, 1], [Cx, Cy, 1]]))
  b = det(numpy.array([[Ax, Ay, modA], [Bx, By, modB], [Cx, Cy, modC]]))

  return ((Sx/a, Sy/a), sqrt(b/a + (Sx*Sx + Sy*Sy)/(a*a)))

verts.append((0, 0))
verts.append((0, 10000))
verts.append((10000, 0))

triangles.add((N, N+1, N+2))

order = range(N)
random.shuffle(order)

for i in order:
  for (A,B,C) in list(triangles):
    ((Sx, Sy), rad) = centrerad(verts[A], verts[B], verts[C])
    if (dist((Sx, Sy), verts[i]) < rad):
      e1 = [A,B]
      e2 = [B,C]
      e3 = [C,A]
      e1.sort()
      e2.sort()
      e3.sort()
      dictadd(edges, tuple(e1))
      dictadd(edges, tuple(e2))
      dictadd(edges, tuple(e3))
      triangles.remove((A,B,C))
  for (v1,v2) in edges.keys():
    if (edges[(v1,v2)] > 1): continue
    tri = [i,v1,v2]
    tri.sort()
    triangles.add(tuple(tri))
  edges = dict()

for (v1,v2,v3) in triangles:
  if (v1 >= N or v2 >= N or v3 >= N): continue
  mat[v1][v2] = 1
  mat[v2][v1] = 1
  mat[v2][v3] = 1
  mat[v3][v2] = 1
  mat[v3][v1] = 1
  mat[v1][v3] = 1

for m in mat:
  for m1 in m:
    print m1,
  print


