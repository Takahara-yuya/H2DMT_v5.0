#pragma once
#ifndef MTFWD2D_H_
#define MTFWD2D_H_
#include "MTdata2D.h"
#include "MTmesh2D.h"
const double PI = 3.1415926535897932;
const double EPS = 8.854187817 * (1E-12);
const double MIU = 4 * PI * (1E-7);
const double RHO_AIR = 1E+10;

class MTFwd2D
{
public:

	//! Constructor
	//!
	MTFwd2D() {

	};
	//! Destructor
	~MTFwd2D() {

	};
	
	//mpi parameters//
	int np, myid;

	//Basic Parameter
	string OutputFile_Root, OutputTmp_Root;
	int Num_AirL_TE, Numx_expand;

	MatrixXd Jacobi;// 

	vector<double> Air_Thick;
	vector<double> X_expand, X_expandR;
	vector<double> Y_expand;
	vector<double> Rhosite, RhoL_Xexpand;

	MatrixXd k2e;

	void Init_Module(string Filename, MTData2D data_fwd, MTMesh2D model);

	void Fwd(MTMesh2D& model, MTData2D& data_fwd, int TETM, int inv, int print);

protected:

	void _Fwd2DMT(MTMesh2D& model, MTData2D& data_fwd, int idx, int inv, int TETM, int print);

	void _XYI8(int NX, int NY, vector<double> X, vector<double> Y,
		vector< vector<int> >& I8, vector< vector<double> >& XY, int idx);

	void _K88(double A, double B, complex<double> lamda, complex<double> tao, complex<double> k,
		MatrixXcd& KE, MatrixXcd& KE3);

	void _E2K2(int NX, int NY, double freq, vector< vector<int> > I8, vector< vector<double> > XY,
		SparseMatrix < complex< double> >& KT, vector< vector<double> > RO, int idx);

	void _Init_Kpara(complex<double>& k, complex<double>& tao, complex<double>& lamda,
		double freq, double Rho, int idx);

	void _Update_para(int nx, int ny, int nsite, vector< vector<double> > Cen_val_RhoL, vector<int> site_index, vector<int> upbound);
};
#endif