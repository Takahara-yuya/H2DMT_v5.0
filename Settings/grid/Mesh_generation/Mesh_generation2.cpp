#include <algorithm>
#include <vector>
#include <string>
#include <iostream> // 
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
using namespace std;
int main()
{
	ifstream fin("90z.dat");
	vector<double> thick, node;
	node.resize(90);
	thick.resize(89);
	for (int i = 0; i < 90; i++)
	{
		fin >> node[i];
	}
	for (int i = 0; i < 89; i++)
	{
		thick[i] = node[i + 1] - node[i];
	}
	for (int i = 0; i < 88; i++)
	{
		cerr << log10(thick[i + 1] - thick[i]) << endl;
	}
	/*int layer = 5;
	ofstream fout("zmesh.dat");
	double y1 = 0.01;
	double ys;
	ys = y1 * layer;
	double y2;
	y2 = y1;
	vector<double> yNode, ythick;
	while (1)
	{
		if (ys > 300.)break;
		y1 *= 1.1;
		ys += y1;
		ythick.push_back(y1);
	}
	yNode.resize(ythick.size() + 1 + layer);
	yNode[0] = 0;
	for (int i = 1; i < layer + 1; i++)
	{
		yNode[i] = yNode[i - 1] + y2;
	}
	for (int i = layer + 1; i < ythick.size() + 1 + layer; i++)
	{
		yNode[i] = yNode[i - 1] + ythick[i - layer - 1];
	}
	for (int i = 0; i < yNode.size(); i++)
	{
		fout << yNode[i] * 1000. << endl;
	}*/
}