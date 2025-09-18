/*
***********************************************************************

MTmesh2D.cpp	(Constructing MT model)
This file is part of H2DMT.

***********************************************************************

Sep 01, 2023
Copyright 2023

Zuwei Huang
hzw1498218560@tongji.edu.cn
School of Ocean and Earth Science, Tongji University
Integrated Geophysics Group

Luolei Zhang
zhangluolei@hotmail.com
School of Ocean and Earth Science, Tongji University
Integrated Geophysics Group

version 5.6.0
Dec 19, 2024


***********************************************************************
*/

#include "../include/MTmesh2D.h"

void MTMesh2D::Init_MTMesh2D_File(string Filename, MTData2D data_obs)
{
	MTMesh2D::ReadMeshFile(Filename, data_obs);
	Num_Cen = Num_xCen * Num_yCen;
	Num_Node = Num_xNode * Num_yNode;
	//check the observation point , if not located at the center of the grid, exit!
}

void MTMesh2D::Init_MTMesh2D_Empty(int a, int b, std::vector<double> a1, std::vector<double> b1, double c)
{
	Num_xCen = a;
	Num_yCen = b;
	sinValue = c;
	Num_xNode = Num_xCen + 1;
	Num_yNode = Num_yCen + 1;

	Num_Cen = Num_xCen * Num_yCen;
	Num_Node = Num_xNode * Num_yNode;

	if (a1.size() != a + 1)
	{
		std::cerr << "Demension Misses!" << std::endl;
		exit(1);
	}
	if (b1.size() != b + 1)
	{
		std::cerr << "Demension Misses!" << std::endl;
		exit(1);
	}

	xNode.resize(Num_xNode);
	xNode.resize(Num_yNode);
	xCen.resize(Num_xCen);
	yCen.resize(Num_yCen);

	xNode.assign(a1.begin(), a1.end());
	yNode.assign(b1.begin(), b1.end());

	xmin = xNode.front(); xmax = xNode.back();
	ymin = yNode.front(); ymax = yNode.back();

	for (int i = 0; i < Num_xCen; i++)
	{
		xCen[i] = 0.5 * (xNode[i] + xNode[i + 1]);
	}

	for (int i = 0; i < Num_yCen; i++)
	{
		yCen[i] = 0.5 * (yNode[i] + yNode[i + 1]);
	}

	sgm1D.resize(Num_Cen);
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		for (int iy = 0; iy < Num_yCen; iy++)
		{
			sgm1D(iy + ix * Num_yCen) = -c;
		}
	}

}

void MTMesh2D::ReadMeshFile(string Filename, MTData2D data_obs)
{
	double tmp; string tmpp;
	std::ifstream fin(Filename);
	if (!fin) {
		std::cerr << "Can't open MeshFile!\n";
		return;
	}
	fin >> tmpp;
	if (tmpp == "MARINE")marine = true;
	else marine = false;
	fin >> Num_xCen >> Num_yCen;

	Num_xNode = Num_xCen + 1;
	Num_yNode = Num_yCen + 1;
	Num_Cen = Num_xCen * Num_yCen;
	Num_Node = Num_xNode * Num_yNode;

	xNode.resize(Num_xNode);
	yNode.resize(Num_yNode);
	xCen.resize(Num_xCen);
	xTopo.resize(Num_xCen);
	yCen.resize(Num_yCen);
	upbound.resize(Num_xCen);

	for (int i = 0; i < Num_xCen; i++)
	{
		fin >> xCen[i];
		xCen[i] *= 1000.;
	}

	for (int i = 1; i < Num_xNode - 1; i++)xNode[i] = 0.5 * (xCen[i - 1] + xCen[i]);
	xNode[0] = xCen[0] - (xNode[1] - xCen[0]);
	xNode[Num_xNode - 1] = xCen[Num_xCen - 1] + (xCen[Num_xCen - 1] - xNode[Num_xNode - 2]);

	for (int i = 0; i < Num_yCen; i++)
	{
		fin >> yCen[i];
		yCen[i] *= 1000.;
	}

	for (int i = 1; i < Num_yNode - 1; i++)yNode[i] = 0.5 * (yCen[i - 1] + yCen[i]);
	yNode[0] = yCen[0] - abs(yCen[0] - yNode[1]);
	yNode[Num_yNode - 1] = yCen[Num_yCen - 1] + (yCen[Num_yCen - 1] - yNode[Num_yNode - 2]);

	for (int i = 0; i < Num_xCen; i++)
	{
		fin >> xTopo[i];
		xTopo[i] *= 1000.;
	}
	//
	fix.resize(Num_xCen);
	for (int i = 0; i < Num_xCen; ++i) {
		fix[i].resize(Num_yCen);
	}
	//
	sgm1D.resize(Num_Cen);
	for (int i = 0; i < Num_xCen; i++)
	{
		for (int j = 0; j < Num_yCen; j++)
		{
			fin >> tmp;
			if (tmp > 9.9)
				fix[i][j] = 1;
			else if (abs(tmp - rho_water) < 0.000001 && (yCen[j] <= xTopo[i]))
				fix[i][j] = 2;
			else
				fix[i][j] = 0;
			sgm1D(j + i * Num_yCen) = -tmp;
		}
		for (int j = 0; j < Num_yCen - 1; j++)
		{
			if (fix[i][0] == 0)
			{
				upbound[i] = 0;
				break;
			}
			if (fix[i][j] != 0 && fix[i][j + 1] == 0)
			{
				upbound[i] = j + 1;
				break;
			}

		}
	}
	fin.close();

	/*Siteindex(data_obs);*/
	site_index.resize(data_obs.nSite);
	for (int i = 0; i < data_obs.nSite; i++)
	{
		site_index[i] = findCloseElements(xCen, data_obs.xSite[i]);
	}

	xmin = xNode[0]; xmax = xNode[Num_xNode - 1];
	ymin = yNode[0]; ymax = yNode[Num_yNode - 1];
}

void MTMesh2D::OutputMeshFile(std::string Filename)
{
	ofstream fout(Filename + ".dat");
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		for (int iy = 0; iy < Num_yCen; iy++)
		{
			fout << xCen[ix] / 1000. << " " << yCen[iy] / 1000. << " " << -sgm1D(iy + ix * Num_yCen) << endl;
		}
	}
	fout.close();
	//.mesh
	string f2 = Filename + ".mesh";
	ofstream fout2(f2);
	if (marine)
		fout2 << "MARINE" << endl;
	else
		fout2 << "LAND" << endl;

	fout2 << Num_xCen << " " << Num_yCen << endl;
	for (int i = 0; i < Num_xCen; i++)
	{
		fout2 << xCen[i] * 0.001 << " ";
	}
	fout2 << endl;
	for (int i = 0; i < Num_yCen; i++)
	{
		fout2 << yCen[i] * 0.001 << " ";
	}
	fout2 << endl;
	for (int i = 0; i < Num_xCen; i++)
	{
		fout2 << xTopo[i] * 0.001 << " ";
	}
	fout2 << endl;
	for (int i = 0; i < Num_xCen; i++)
	{
		for (int j = 0; j < Num_yCen; j++)
		{
			fout2 << -sgm1D(j + i * Num_yCen) << " ";
		}
		fout2 << endl;
	}
	//gmt
	string f3 = Filename + ".gmt";
	ofstream fout3(f3);
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		for (int iy = 0; iy < Num_yCen; iy++)
		{
			if (fix[ix][iy] == 0)
			{
				fout3 << "> -Z" << -sgm1D(iy + ix * Num_yCen) << endl;
				fout3 << xNode[ix] / 1000. << " " << yNode[iy] / 1000. << endl;
				fout3 << xNode[ix] / 1000. << " " << yNode[iy + 1] / 1000. << endl;
				fout3 << xNode[ix + 1] / 1000. << " " << yNode[iy + 1] / 1000. << endl;
				fout3 << xNode[ix + 1] / 1000. << " " << yNode[iy] / 1000. << endl;
			}
		}
	}
	fout3.close();
}

void MTMesh2D::AutoMesh(MTData2D data_obs, string yFilename, string topo,
	double value, int expand, int grid_expand, int refinement, bool site_Ref, bool& local_Ref, double& lr, bool Marine)
{
	lr *= 1000.;
	if (Marine)
		marine = true;
	else
		marine = false;
	ifstream fin(yFilename);
	if (!fin) {
		std::cerr << "Can't open Depth_meshFile!" << endl;
		exit(0);
	}
	while (!fin.eof())
	{
		double yy;
		fin >> yy;
		if (fin.fail())break;
		yNode.push_back(yy);
	}
	double expand2 = expand;
	expand = expand + grid_expand;
	Num_yCen = yNode.size() - 1;
	Num_yNode = Num_yCen + 1;
	yCen.resize(Num_yCen); site_index.resize(data_obs.nSite); yNode.resize(Num_yNode);
	for (int i = 0; i < Num_yCen; i++)yCen[i] = 0.5 * (yNode[i] + yNode[i + 1]);
	//find min
	double sm, sm_last, x11, x22;
	sm_last = 0.5 * abs(data_obs.xSite[1] - data_obs.xSite[0]);
	xNode.push_back(data_obs.xSite[0] + sm_last);
	xNode.push_back(data_obs.xSite[0] - sm_last);
	for (int i = 1; i < data_obs.nSite - 1; i++)
	{
		x11 = 0.5 * abs(data_obs.xSite[i + 1] - data_obs.xSite[i]);
		x22 = abs(data_obs.xSite[i] - data_obs.xSite[i - 1]) - sm_last;
		if (x22 >= 2 * x11)sm = x11;
		else sm = x22;
		xNode.push_back(data_obs.xSite[i] + sm);
		xNode.push_back(data_obs.xSite[i] - sm);
		sm_last = sm;
	}
	xNode.push_back(data_obs.xSite[data_obs.nSite - 1] + abs(data_obs.xSite[data_obs.nSite - 1] - data_obs.xSite[data_obs.nSite - 2]) - sm);
	sort(xNode.begin(), xNode.end());
	auto last = unique(xNode.begin(), xNode.end());
	xNode.erase(last, xNode.end());
	vector<double> xNode_new, xNode_new2;
	xNode_new = xNode;
	xCen.resize(xNode.size() - 1);
	for (int i = 0; i < xNode.size() - 1; i++)
	{
		xCen[i] = 0.5 * (xNode[i] + xNode[i + 1]);
	}
	if (local_Ref)
	{
		xNode_new.insert(xNode_new.end(), xCen.begin(), xCen.end());
		sort(xNode_new.begin(), xNode_new.end());
		xNode_new2 = xNode_new;
		for (int i = 0; i < xNode_new.size() - 1; i++)
		{
			int n2 = int((abs(xNode_new[i + 1] - xNode_new[i])-0.1) / lr);
			for (int j = 0; j < n2; j++)
			{
				xNode_new2.push_back(xNode_new[i] + (j + 1) * abs(xNode_new[i + 1] - xNode_new[i]) / (n2 + 1.0));
			}
		}
		std::unordered_set<double> toRemove(xCen.begin(), xCen.end());
		// 使用 std::remove_if 删除 vec1 中在 vec2 中出现的元素
		xNode_new2.erase(
			std::remove_if(xNode_new2.begin(), xNode_new2.end(), [&](double x) {
				return toRemove.find(x) != toRemove.end();
				}),
			xNode_new2.end()
		);
		sort(xNode_new2.begin(), xNode_new2.end());
		last = unique(xNode_new2.begin(), xNode_new2.end());
		xNode_new2.erase(last, xNode_new2.end());
		xNode = xNode_new2;
	}
	//expand
	double xl = abs(xNode[1] - xNode[0]);
	double xr = abs(xNode[xNode.size() - 1] - xNode[xNode.size() - 2]);
	vector<double> xWidth_l, xWidth_r;
	xWidth_l.resize(expand); xWidth_r.resize(expand);
	if (expand > 0)xWidth_l[0] = xl;
	if (expand > 0)xWidth_r[0] = xr;
	double xs = xNode[0]; double xe = xNode[xNode.size() - 1];
	for (int i = 1; i < xWidth_l.size(); i++)
	{
		if (i > 1)
		{
			xWidth_l[i] += xWidth_l[i - 1] + 1.35 * (xWidth_l[i - 1] - xWidth_l[i - 2]);
			xWidth_r[i] += xWidth_r[i - 1] + 1.35 * (xWidth_r[i - 1] - xWidth_r[i - 2]);
		}
		else
		{
			xWidth_l[i] += xWidth_l[i - 1] + 1.35 * (xWidth_l[i - 1] - 0);
			xWidth_r[i] += xWidth_r[i - 1] + 1.35 * (xWidth_r[i - 1] - 0);
		}
	}
	for (int i = 0; i < xWidth_l.size(); i++)
	{
		xNode.push_back(xs - xWidth_l[i]);
		xNode.push_back(xe + xWidth_r[i]);
	}
	sort(xNode.begin(), xNode.end());
	//
	Num_xNode = xNode.size(); Num_xCen = Num_xNode - 1;
	xCen.resize(Num_xCen);
	for (int i = 0; i < Num_xCen; i++)
	{
		xCen[i] = 0.5 * (xNode[i] + xNode[i + 1]);
	}
	Num_Cen = Num_xCen * Num_yCen;
	Num_Node = Num_xNode * Num_yNode;
	sgm1D.resize(Num_xCen * Num_yCen); 	xTopo.resize(Num_xCen); upbound.resize(Num_xCen);
	for (int i = 0; i < data_obs.nSite; i++)
	{
		site_index[i] = findCloseElements(xCen, data_obs.xSite[i]);
	}
	unordered_set<double> seen;
	bool hasDuplicates = false;
	for (double value : site_index) {
		if (seen.find(value) != seen.end()) {
			hasDuplicates = true;
		}
		else {
			seen.insert(value);
		}
	}
	if (hasDuplicates) {
		cerr<< "Siteindex Error!" << endl;
		exit(0);
	}
	///add topo
	fix.resize(Num_xCen);
	for (int i = 0; i < Num_xCen; ++i) {
		fix[i].resize(Num_yCen);
	}
	//fileread
	vector<double> xtopo_tmp, ytopo_tmp;
	ifstream fint(topo);
	if (!fint) {
		std::cerr << "Can't open topography file!\n";
		exit(0);
	}
	while (!fint.eof())
	{
		double xx1, yy1;
		fint >> xx1 >> yy1;
		xx1 *= 1000.; yy1 *= 1000.;
		if (Marine && yy1 < 0)
		{
			cerr << "Marine topography must >= 0 !" << endl;
			exit(0);
		}
		xtopo_tmp.push_back(xx1);
		ytopo_tmp.push_back(yy1);
	}
	fint.close();
	xtopo_tmp.push_back(999999999.); xtopo_tmp.insert(xtopo_tmp.begin(), -999999999.);
	double tt1 = ytopo_tmp.front(); double tt2 = ytopo_tmp.back();
	ytopo_tmp.push_back(tt2); ytopo_tmp.insert(ytopo_tmp.begin(), tt1);
	//interpolate
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		xTopo[ix] = _linearinter1D(xtopo_tmp, ytopo_tmp, xCen[ix]);
	}
	//
	if (!Marine)
	{
		auto minElement = std::min_element(xTopo.begin(), xTopo.end());
		ymin = *minElement;
		for (int iy = 0; iy < Num_yNode; iy++)
		{
			yNode[iy] += ymin;
			if (iy < Num_yCen)
				yCen[iy] += ymin;
		}
	}
	else
	{
		ymin = 0;
	}
	//
	for (int i = 0; i < Num_xCen; i++)
	{
		for (int j = 0; j < Num_yCen; j++)
		{
			if (yCen[j] <= xTopo[i])
			{
				if (!Marine)
					fix[i][j] = 1;
				else
					fix[i][j] = 2;
			}
			else fix[i][j] = 0;
		}
		for (int j = 0; j < Num_yCen - 1; j++)
		{
			if (fix[i][0] == 0)
			{
				upbound[i] = 0;
				break;
			}
			if (fix[i][j] != 0 && fix[i][j + 1] == 0)
			{
				upbound[i] = j + 1;
				break;
			}
		}
	}
	///
	sgm1D.resize(Num_Cen);
	for (int i = 0; i < Num_xCen; i++)
	{
		for (int j = 0; j < Num_yCen; j++)
		{
			if (fix[i][j] == 0)
				sgm1D(j + i * Num_yCen) = -value;
			else
			{
				if (fix[i][j] == 1)
					sgm1D(j + i * Num_yCen) = -rho_air;
				else
					sgm1D(j + i * Num_yCen) = -rho_water;
			}
		}
	}
	xmin = xNode[0]; xmax = xNode[Num_xNode - 1];
	ymin = yNode[0]; ymax = yNode[Num_yNode - 1];
	//0113

}

void MTMesh2D::Siteindex(MTData2D data_cal)
{
	site_index.resize(data_cal.nSite);
	for (int i = 0; i < data_cal.nSite; i++)
	{
		site_index[i] = findCloseElements(xCen, data_cal.xSite[i]);
	}
}

int MTMesh2D::findCloseElements(vector<double> a, double b)
{
	for (int i = 0; i < a.size(); ++i)
	{
		//cerr << a[i]<<" "<<b<<" "<<abs(a[i] - b) << endl;
		if (std::abs(a[i] - b) < 5)
		{
			return i;
		}
	}
	return 0;
}

double MTMesh2D::_linearinter1D(vector<double> x_data, vector<double> y_data, double x)
{
	if (x < x_data.front()) {
		cerr << x << " " << x_data.front() << endl;
		cerr << "Interpolation value is outside the range of known data" << endl;
		exit(1);
	}


	int i = 0;
	while (x > x_data[i + 1])
	{
		++i;
	}

	// 1D
	if (x <= x_data.back())
	{
		double x0 = x_data[i];
		double x1 = x_data[i + 1];
		double y0 = y_data[i];
		double y1 = y_data[i + 1];

		return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
	}
	else
	{
		double x0 = x_data[i - 1];
		double x1 = x_data[i];
		double y0 = y_data[i - 1];
		double y1 = y_data[i];

		return y1 + (y1 - y0) * (x - x1) / (x1 - x0);
	}

}