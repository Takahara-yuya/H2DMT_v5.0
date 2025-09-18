#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int main() {
    double max_depth = 50000.0;  // 50 km以内
    double depth_intervals[] = {0.3, 0.5, 3.0};  // 不同深度范围的最大间隔

    double current_depth = 20;

    double mesh=5;


    ofstream fout("depth.dat");
    double in1=1.2;
    double in2=1.5;
    double in3=2.0;

	fout<<current_depth/1000<<" "<<mesh/1000<<endl;

    while(current_depth<25000)
    {
	mesh=mesh+pow(10,in1);
	current_depth=current_depth+mesh;
	fout<<current_depth/1000<<" "<< log(mesh) <<endl;
    }

      while(current_depth<45000)
    {
        mesh=mesh+pow(10,in2);
        current_depth=current_depth+mesh;
        fout<<current_depth/1000<<" "<<log(mesh)<<endl;
    }

      while(current_depth<60000)
    {
        mesh=mesh+pow(10,in3);
        current_depth=current_depth+mesh;
        fout<<current_depth/1000<<" "<< log(mesh) <<endl;
    }


}

