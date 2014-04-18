#include <iostream>
#include "../src/data3d.h"
#include "../src/Point.h"

using namespace std;

int main(){
   data3d<double> ez;
   Point p(0,5,2);
   ez.create3DArray(10,10,10,0.0);
   cout<<ez[p]<<'\t';
   cout<<ez.p[p.x][p.y][p.z]<<endl;
   return 0;
}

