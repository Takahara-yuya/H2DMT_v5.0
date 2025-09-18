#pragma once
#ifndef MTMESH2D_H_
#define MTMESH2D_H_

#include "CommonHeaders.h"
#include"MTdata2D.h"

class MTMesh2D
{
public:
	//Basic Mesh
	const double rho_air = 10.;
	const double rho_water = -0.52;

	bool marine;

	int Num_xCen, Num_yCen;

	int Num_xNode, Num_yNode, Num_Node, Num_Cen;

	double sinValue, xmin, xmax, ymin, ymax;

	vector<double> xNode, xCen, yNode, yCen, xTopo;
	vector<vector<int>> fix;
	vector<int> site_index, upbound;
	VectorXd sgm1D;

	void Init_MTMesh2D_File(string Filename, MTData2D data_obs);

	void Init_MTMesh2D_Empty(int a, int b, vector<double> a1, vector<double> b1, double c);

	void ReadMeshFile(string Filename, MTData2D data_obs);

	void OutputMeshFile(string Filename);

	void AutoMesh(MTData2D data_obs, string yFilename, string topo, double value, int expand, int grid_expand, int refinement, bool site_Ref, bool Marine);

	void Siteindex(MTData2D data_cal);

	int findCloseElements(vector<double> a, double b);

protected:
	const int sitref = 2;//left and right
	const double disref = 100.;//50m

	double _linearinter1D(vector<double> x_data, vector<double> y_data, double x);
};

#endif