// g++ -o scfeatures scfeatures.cpp -lm

//#include "CImg.h"

#include <vector>
#include <utility>
#include <iostream>

#include "math.h"
//using namespace cimg_library;
using namespace std;

#define PI 3.14159265

double dist(pair<double, double> p1, pair<double, double> p2)
{
  double x1 = p1.first;
  double y1 = p1.second;
  double x2 = p2.first;
  double y2 = p2.second;

  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

int main(int argc, char** argv)
{
  //char* imname = argv[1];
  char* cornername = argv[2];
  char* outfile = argv[3];

  FILE* f = fopen(cornername, "r");
  double x, y;
  vector<pair<double, double> > corners;

  while (fscanf(f, "%lf %lf", &x, &y) == 2)
  {
    pair<double, double> p(x,y);
    corners.push_back(p);
  }

  fclose(f);

  f = fopen(outfile, "w");

  // Shape context...
  vector <int*> scf;

  double alpha = 0;
  for (int i = 0; i < (int) corners.size(); i ++)
    for (int j = 0; j < (int) corners.size(); j ++)
      alpha += dist(corners[i], corners[j]);
  alpha /= (corners.size()*corners.size());

  for (int i = 0; i < (int) corners.size(); i++) {
    int* feat = new int [60];
    for (int j = 0; j < 60; j ++) {
      feat[j] = 0;
    }

    for (int j = 0; j < (int) corners.size(); j++) {
      if (j == i) {
	continue;
      }

      double xco = corners[j].first - corners[i].first;
      double yco = corners[j].second - corners[i].second;
      // convert to polar coordinates...

      double theta = atan2(yco, xco) + PI;
      double r = dist(corners[i], corners[j])/alpha;

      int bin1 = 0;
      int bin2 = 0;

      if      (theta < 2*PI/12)  bin1 = 0;
      else if (theta < 4*PI/12)  bin1 = 1;
      else if (theta < 6*PI/12)  bin1 = 2;
      else if (theta < 8*PI/12)  bin1 = 3;
      else if (theta < 10*PI/12) bin1 = 4;
      else if (theta < 12*PI/12) bin1 = 5;
      else if (theta < 14*PI/12) bin1 = 6;
      else if (theta < 16*PI/12) bin1 = 7;
      else if (theta < 18*PI/12) bin1 = 8;
      else if (theta < 20*PI/12) bin1 = 9;
      else if (theta < 22*PI/12) bin1 = 10;
      else                       bin1 = 11;

      if      (r < 0.125) bin2 = 0;
      else if (r < 0.25)  bin2 = 1;
      else if (r < 0.5)   bin2 = 2;
      else if (r < 1)     bin2 = 3;
      else if (r < 2)     bin2 = 4;
      else continue;

      feat[bin2*12 + bin1] ++;
    }

    scf.push_back(feat);
    for (int j = 0; j < 60; j ++)
      fprintf(f, "%d ", feat[j]);
    fprintf(f, "\n");
  }

  // There are no features based on the image for now.
  //CImg<unsigned char> image(imname);

  fclose(f);

  return 0;
}
