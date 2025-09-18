/*
***********************************************************************

MTfwd2D.cpp (MT Forward Modeling)
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

* **********************************************************************
*
*
 *
 *  For TE-mode,  E=(Ex,0,0),  H=(0,Hy,Hz);
 *  For TM-mode,  E=(0,Ey,Ez), H=(Hx,0,0);
 *       z (strike)
 *      /
 *    /
 * o/-------- x
 *  |
 *  |
 *  |
 *  |
 *  y
*/
#include "../include/MTfwd2D.h"

void MTFwd2D::Init_Module(string Filename, MTData2D data_fwd, MTMesh2D model)
{
	ifstream fin(Filename);
	string input_line, tmp1, tmp2, tmp3;
	fin >> tmp1 >> tmp2 >> tmp3 >> OutputFile_Root;
	fin >> tmp1 >> tmp2 >> tmp3 >> OutputTmp_Root;
	fin >> tmp1 >> tmp2 >> tmp3 >> Num_AirL_TE;
	Air_Thick.resize(Num_AirL_TE);
	fin >> tmp1 >> tmp2 >> tmp3;
	for (int i = 0; i < Num_AirL_TE; i++)
	{
		fin >> Air_Thick[i];
		Air_Thick[i] *= 1000.;
	}
	fin >> tmp1 >> tmp2 >> Numx_expand;
	X_expand.resize(Numx_expand);
	X_expandR.resize(Numx_expand);
	fin >> tmp1 >> tmp2 >> tmp3;
	for (int i = 0; i < Numx_expand; i++)
	{
		fin >> X_expand[i];
		X_expand[i] *= 1000.;
		X_expand[i] = X_expand[i] + model.xNode[0];
		//	if (i == 0)X_expand[i] = -1.5 * abs(model.xNode[1] - model.xNode[0]) + model.xNode[0];
		//	else X_expand[i] = -1.5 * abs(X_expand[i] - X_expand[i - 1]) + X_expand[i - 1];
	}
	fin >> tmp1 >> tmp2 >> tmp3;
	for (int i = 0; i < Numx_expand; i++)
	{
		fin >> X_expandR[i];
		X_expandR[i] *= 1000.;
		X_expandR[i] = X_expandR[i] + model.xNode[model.Num_xNode - 1];
		//	if (i == 0)X_expandR[i] = 1.5 * abs(model.xNode[model.Num_xNode - 1] - model.xNode[model.Num_xNode - 2]) + model.xNode[model.Num_xNode - 1];
		//	else X_expandR[i] = 1.5 * abs(X_expandR[i] - X_expandR[i - 1]) + X_expandR[i - 1];
	}
	RhoL_Xexpand.resize(model.Num_yCen);

	k2e.resize(8, 8);
	k2e.setZero();
	k2e(0, 0) = 6.;
	k2e(1, 0) = 2.; k2e(1, 1) = 6.;
	k2e(2, 0) = 3.; k2e(2, 1) = 2.; k2e(2, 2) = 6.;
	k2e(3, 0) = 2.; k2e(3, 1) = 3.; k2e(3, 2) = 2.; k2e(3, 3) = 6.;
	k2e(4, 0) = -6.; k2e(4, 1) = -6.; k2e(4, 2) = -8.; k2e(4, 3) = -8.; k2e(4, 4) = 32.;
	k2e(5, 0) = -8.; k2e(5, 1) = -6.; k2e(5, 2) = -6.; k2e(5, 3) = -8.; k2e(5, 4) = 20.; k2e(5, 5) = 32.;
	k2e(6, 0) = -8.; k2e(6, 1) = -8.; k2e(6, 2) = -6.; k2e(6, 3) = -6.; k2e(6, 4) = 16.; k2e(6, 5) = 20.; k2e(6, 6) = 32.;
	k2e(7, 0) = -6.; k2e(7, 1) = -8.; k2e(7, 2) = -8.; k2e(7, 3) = -6.; k2e(7, 4) = 20.; k2e(7, 5) = 16.; k2e(7, 6) = 20.; k2e(7, 7) = 32.;

	MatrixXd k22 = k2e.transpose();
	for (int i = 0; i < 8; i++)k22(i, i) = 0;
	k2e = k2e + k22;

	if (model.marine)Num_AirL_TE = 0;
}

void MTFwd2D::_Fwd2DMT(MTMesh2D& model, MTData2D& data_fwd, int idx, int inv, int TETM, int print)
{
	//
	if ((print == 1) && (myid == 0))
	{
		cerr << "-------------FORWARD MODELING---------------" << endl;
		cerr << "NUMBER OF FREQUENCY: " << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	auto start_time = chrono::high_resolution_clock::now();

	int NX = model.Num_xCen + 2 * Numx_expand;
	int NY = model.Num_yCen + idx * Num_AirL_TE;
	int ND = 3 * NX * NY + 2 * NX + 2 * NY + 1;
	int NE = NX * NY;
	double DL = abs(model.yNode[3] - model.yNode[0]);
	double DL1 = abs(model.yNode[1] - model.yNode[0]);
	double A, B;
	complex<double> WU, WU1;

	vector< vector<double> >  Cen_val_RhoL;
	Cen_val_RhoL.resize(model.Num_xCen);
	for (int i = 0; i < Cen_val_RhoL.size(); i++)
	{
		Cen_val_RhoL[i].resize(model.Num_yCen);
	}
	for (int ix = 0; ix < model.Num_xCen; ix++)
	{
		for (int iy = 0; iy < model.Num_yCen; iy++)
		{
			Cen_val_RhoL[ix][iy] = pow(10, -1 * model.sgm1D(iy + ix * model.Num_yCen));
		}
	}
	//upbound expand
	vector<int>upb;
	upb.resize(NX);
	for (int i = 0; i < Numx_expand; i++)
	{
		upb[i] = model.upbound[0];
	}
	for (int i = Numx_expand; i < Numx_expand + model.Num_xCen; i++)
	{
		upb[i] = model.upbound[i - Numx_expand];
	}
	for (int i = Numx_expand + model.Num_xCen; i < 2 * Numx_expand + model.Num_xCen; i++)
	{
		upb[i] = model.upbound[model.Num_xCen - 1];
	}
	///
	VectorXcd AF, PF, AFD, PFD, source;
	MatrixXcd Z;
	AF.resize(data_fwd.nSite); PF.resize(ND);
	AFD.resize(ND); PFD.resize(ND); source.resize(ND);
	Z.resize(data_fwd.nSite, data_fwd.nFre);
	vector< vector<int> > I8;  vector< vector<double> > XY;
	SparseMatrix < complex< double> > K;
	VectorXcd a, b;
	VectorXcd a1, b1, c1;
	a.resize(4); b.resize(4); b.setZero(); a.setZero();
	//a.resize(2); b.resize(2); b.setZero(); a.setZero();
	a1.resize(8); b1.resize(8); c1.resize(8);
	//starting forward
	_Update_para(model.Num_xCen, model.Num_yCen, data_fwd.nSite, Cen_val_RhoL, model.site_index, model.upbound);
	_XYI8(NX, NY, model.xNode, model.yNode, I8, XY, idx);
	SparseLU<SparseMatrix<complex< double>>, COLAMDOrdering<int>> solver;
	for (int ifre = myid; ifre < data_fwd.nFre; ifre += np)
	{
		//Finite element
		if (print == 1)
		{
			cerr << ifre + 1 << " ";
		}
		_E2K2(NX, NY, data_fwd.FreqL[ifre], I8, XY, K, Cen_val_RhoL, idx);
		source.setZero();
		//add boundary condition
		for (int ix = 0; ix < NX; ix++)
		{
			if (ix < NX - 1)
			{
				int d1, d2, L;
				if (model.marine)//Marine boundary condition
				{
					L = ix * NY;
				}
				else
				{
					if (idx == 0)L = ix * NY + upb[ix];//upper bound
					else L = ix * NY;
				}
				d1 = I8[0][L];
				d2 = I8[7][L];
				K.coeffRef(d1, d1) *= 1.0E+10;
				K.coeffRef(d2, d2) *= 1.0E+10;
				source(d1) = K.coeffRef(d1, d1);
				source(d2) = K.coeffRef(d2, d2);
			}
			else
			{
				int d1, d2, d3, L;
				if (model.marine)//Marine boundary condition
				{
					L = ix * NY;
				}
				else
				{
					if (idx == 0)L = ix * NY + upb[ix];//upper bound
					else L = ix * NY;
				}
				d1 = I8[0][L];
				d2 = I8[7][L];
				d3 = I8[3][L];
				K.coeffRef(d1, d1) *= 1.0E+10;
				K.coeffRef(d2, d2) *= 1.0E+10;
				K.coeffRef(d3, d3) *= 1.0E+10;
				source(d1) = K.coeffRef(d1, d1);
				source(d2) = K.coeffRef(d2, d2);
				source(d3) = K.coeffRef(d3, d3);
			}
		}
		//solver.compute(K);
		if (ifre == myid)solver.analyzePattern(K);
		solver.factorize(K);
		PF = solver.solve(source);
		//Calculate Aux Field
		for (int isite = 0; isite < data_fwd.nSite; isite++)
		{
			int site_cen = NY * (model.site_index[isite] + Numx_expand) + idx * Num_AirL_TE + model.upbound[model.site_index[isite]];
			int site_node = I8[7][site_cen];
			if (idx == 0)//TM mode
			{
				WU = 1. / Rhosite[isite];
				a.coeffRef(0) = 1. / WU * (-11.) / (2. * DL);
				a.coeffRef(1) = 1. / WU * (18.) / (2. * DL);
				a.coeffRef(2) = 1. / WU * (-9.) / (2. * DL);
				a.coeffRef(3) = 1. / WU * (2.) / (2. * DL);
				b(0) = 1.0;
			}
			else//TE mode
			{
				WU = complex<double>(0, 2 * PI * data_fwd.FreqL[ifre] * MIU);
				b.coeffRef(0) = 1. / WU * (-11.) / (2. * DL);
				b.coeffRef(1) = 1. / WU * (18.) / (2. * DL);
				b.coeffRef(2) = 1. / WU * (-9.) / (2. * DL);
				b.coeffRef(3) = 1. / WU * (2.) / (2. * DL);
				a(0) = 1.0;
			}
			complex<double>av = a.transpose() * PF.segment(site_node, 4);
			complex<double>bv = b.transpose() * PF.segment(site_node, 4);
			Z(isite, ifre) = av / bv;
			double Rho = 1. / (2 * PI * data_fwd.FreqL[ifre] * MIU) * abs(pow(Z(isite, ifre), 2));
			double Phase = atan(Z(isite, ifre).imag() / Z(isite, ifre).real());
			if (idx == 0)//TM mode
			{
				AF.coeffRef(isite) = av;//AF=EX PF=HZ
				data_fwd.Rho_TM[isite][ifre] = log10(Rho);
				data_fwd.Pha_TM[isite][ifre] = 180 + Phase * 180. / PI; //atan-->atan2
			}
			else//TE mode
			{
				AF.coeffRef(isite) = bv;//AF=HX PF=EZ
				data_fwd.Rho_TE[isite][ifre] = log10(Rho);
				data_fwd.Pha_TE[isite][ifre] = Phase * 180. / PI;
			}
		}
		//Jacobian calculate
		if (inv == 1)
		{
			//add upbound condition
			int d1, d2, d3, L;
			//calculate PFD and AFD
			for (int isite = 0; isite < data_fwd.nSite; isite++)
			{
				if (idx == 0)
				{
					WU = 1. / Rhosite[isite];
					WU1 = complex<double>(0, -2 * PI * data_fwd.FreqL[ifre] * MIU);
				}
				else
				{
					WU = complex<double>(0, 2 * PI * data_fwd.FreqL[ifre] * MIU);
					WU1 = 1. / Rhosite[isite];
				}
				int site_cen = NY * (model.site_index[isite] + Numx_expand) + idx * Num_AirL_TE + model.upbound[model.site_index[isite]];
				int site_node = I8[7][site_cen];
				//PFD
				if (idx == 0)
				{
					if (!model.marine)
						PFD.setZero();
					else
					{
						source.setZero();
						source.coeffRef(site_node) = complex<double>(1., 0);
						PFD = solver.solve(source);
					}
				}
				else
				{
					source.setZero();
					source.coeffRef(site_node) = complex<double>(1., 0);
					PFD = solver.solve(source);
				}
				//AFD
				source.setZero();
				if (!model.marine)
					source.coeffRef(site_node) = 1. / WU * (11.) / (1. * DL) * double(idx);
				else
					source.coeffRef(site_node) = 1. / WU * (11.) / (1. * DL);
				source.coeffRef(site_node + 1) = 1. / WU * (-18.) / (1. * DL);
				source.coeffRef(site_node + 2) = 1. / WU * (9.) / (1. * DL);
				source.coeffRef(site_node + 3) = 1. / WU * (-2.) / (1. * DL);
				AFD = solver.solve(source);
				//genarate Jacobian matrix
				for (int ix = Numx_expand; ix < model.Num_xCen + Numx_expand; ix++)
				{
					for (int iy = Num_AirL_TE * idx; iy < NY; iy++)
					{
						double rhoc = Cen_val_RhoL[ix - Numx_expand][iy - Num_AirL_TE * idx];
						int KK = iy + ix * (NY);
						int KK2 = iy - Num_AirL_TE * idx + (ix - Numx_expand) * model.Num_yCen;
						A = abs(XY[0][int(I8[3][KK])] - XY[0][int(I8[0][KK])]);
						B = abs(XY[1][int(I8[0][KK])] - XY[1][int(I8[1][KK])]);
						for (int i = 0; i < 8; i++)
						{
							a1.coeffRef(i) = PF.coeff(I8[i][KK]);
							if (idx == 0)a1.coeffRef(i) *= WU1 * rhoc;//TM*iwu/σ
						}
						for (int i = 0; i < 8; i++)b1.coeffRef(i) = PFD.coeff(I8[i][KK]);
						for (int i = 0; i < 8; i++)c1.coeffRef(i) = AFD.coeff(I8[i][KK]);
						complex<double> AFDC, PFDC, ZDC;
						double RhoDC, PhaDC;
						PFDC = A * B / 180. * b1.transpose() * k2e * a1;
						AFDC = -A * B / 180. * c1.transpose() * k2e * a1;
						if (idx == 0)
						{
							if (ix == (model.site_index[isite] + Numx_expand) && iy == 0)
							{
								AFDC = AFDC + AF(isite) * WU;
							}
						}
						if (idx == 0)ZDC = 1. / pow((PF.coeff(site_node)), 2) * (AFDC * PF.coeff(site_node) - PFDC * AF.coeff(isite));
						else ZDC = 1. / pow((AF.coeff(isite)), 2) * (PFDC * AF.coeff(isite) - AFDC * PF.coeff(site_node));
						RhoDC = 1. / (2 * PI * data_fwd.FreqL[ifre] * MIU) * ((ZDC * conj(Z(isite, ifre)) + conj(ZDC) * Z(isite, ifre))).real();
						PhaDC = (ZDC.imag() * Z(isite, ifre).real() - ZDC.real() * Z(isite, ifre).imag()) * 1. / (pow(Z(isite, ifre).real(), 2))
							* 1. / (1. + pow(Z(isite, ifre).imag() / Z(isite, ifre).real(), 2));
						if (model.fix[ix - Numx_expand][iy - Num_AirL_TE * idx] != 0)
						{
							RhoDC = 0; PhaDC = 0;
						}
						//combine Jacobian
						int idata = ifre + isite * data_fwd.nFre;
						int ndata = data_fwd.nSite * data_fwd.nFre;
						if (idx == 0)//TM mode
						{
							int tmp = int(TETM / 3.0);
							Jacobi(idata + 2 * ndata * tmp, KK2) = RhoDC * 1. / rhoc / (pow(10, data_fwd.Rho_TM[isite][ifre]));
							Jacobi(idata + ndata + 2 * ndata * tmp, KK2) = 180. / 3.1415926535 * PhaDC * 1. / rhoc * 2.3025850929940;
						}
						else
						{
							Jacobi(idata, KK2) = RhoDC * 1. / rhoc / (pow(10, data_fwd.Rho_TE[isite][ifre]));
							Jacobi(idata + ndata, KK2) = 180. / 3.1415926535 * PhaDC * 1. / rhoc * 2.3025850929940;
						}
					}
				}
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if ((print == 1) && (myid == 0))
	{
		cerr << endl;
		auto end_time = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
		cerr << "Forward modeling time: " << duration.count() << " s" << endl;
	}
	Cen_val_RhoL.clear();
	///Ruduce_and_Bcast
	//data and Jacobi
	MPI_Barrier(MPI_COMM_WORLD);
}

void MTFwd2D::_Init_Kpara(complex<double>& k, complex<double>& tao, complex<double>& lamda, double freq, double Rho, int idx)
{
	if (idx == 0)//TM
	{
		double SGM = 1. / Rho;
		double W = (2. * PI * freq);
		tao = 1. / complex<double>(SGM, -W * EPS);
		lamda = complex<double>(0., W * MIU);
		k = sqrt(complex<double>(0., -W * SGM * MIU));
	}
	else//TE
	{
		double SGM = 1. / Rho;
		double W = (2. * PI * freq);
		tao = 1. / complex<double>(0., (W * MIU));
		lamda = complex<double>(SGM, -W * EPS);
		k = sqrt(complex<double>(0., -W * SGM * MIU));
	}
}

void MTFwd2D::_XYI8(int NX, int NY, vector<double> X, vector<double> Y,
	vector< vector<int> >& I8, vector< vector<double> >& XY, int idx)
{
	VectorXd y1, x1;
	int ND = 3 * NX * NY + 2 * NX + 2 * NY + 1;
	int NE = NX * NY;
	int N, N1;
	XY.resize(2);
	I8.resize(8);

	x1.resize(NX + 1); y1.resize(NY + 1);
	for (int i = 0; i < Numx_expand; i++)x1(i) = X_expand[i];
	for (int i = Numx_expand; i < NX + 1 - Numx_expand; i++)x1(i) = X[i - Numx_expand];
	for (int i = NX + 1 - Numx_expand; i < NX + 1; i++)x1(i) = X_expandR[i - (NX + 1 - Numx_expand)];

	for (int i = 0; i < Num_AirL_TE * idx; i++)y1(i) = Air_Thick[i];
	for (int i = Num_AirL_TE * idx; i < NY + 1; i++)y1(i) = Y[i - Num_AirL_TE * idx];

	for (int i = 0; i < 8; i++)
	{
		if (i < 2)XY[i].resize(ND);
		I8[i].resize(NE);
	}

	for (int ix = 0; ix < NX; ix++)
	{
		for (int iy = 0; iy < NY; iy++)
		{
			N = ix * NY + iy;
			N1 = ix * (3 * NY + 2) + 2 * iy;
			I8[0][N] = N1; I8[1][N] = N1 + 2; I8[2][N] = I8[1][N] + 3 * NY + 2; I8[3][N] = N1 + 3 * NY + 2;
			I8[4][N] = N1 + 1; I8[5][N] = I8[1][N] + 2 * NY - iy; I8[6][N] = I8[4][N] + 3 * NY + 2; I8[7][N] = I8[5][N] - 1;
		}
	}

	for (int ix = 0; ix < NX + 1; ix++)
	{
		for (int iy = 0; iy < NY; iy++)
		{
			N = (ix) * (3 * NY + 2) + 2 * iy;
			N1 = N + 1;
			XY[0][N] = x1[ix];
			XY[1][N] = y1[iy];
			XY[0][N1] = x1[ix];
			XY[1][N1] = y1[iy] + (y1[iy + 1] - y1[iy]) * 0.5;
		}
		XY[0][N + 2] = x1(ix);
		XY[1][N + 2] = y1(NY);
	}

	for (int ix = 0; ix < NX; ix++)
	{
		for (int iy = 0; iy < NY + 1; iy++)
		{
			N = ix * (3 * NY + 2) + 2 * NY + iy + 1;
			XY[0][N] = x1[ix] + 0.5 * (x1[ix + 1] - x1[ix]);
			XY[1][N] = y1[iy];
		}
	}

}

void MTFwd2D::_K88(double a, double b, complex<double> lamda, complex<double> tao, complex<double> k,
	MatrixXcd& KE, MatrixXcd& KE3)
{
	complex<double> alpha1 = tao * b / (a * 90.);
	complex<double> alpha2 = a * b * lamda / 180.;
	complex<double> beta = tao * a / (b * 90.);
	complex<double> beta2 = tao * k * b / 30.;

	KE.resize(8, 8);
	KE3.resize(8, 8); KE3.setZero();
	MatrixXcd k1; k1.resize(8, 8); k1.setZero();
	MatrixXcd k2; k2.resize(8, 8); k2.setZero();

	k1(0, 0) = 52. * alpha1 + 52. * beta; k1(1, 0) = 17. * alpha1 + 28. * beta;
	k1(2, 0) = 23. * alpha1 + 23. * beta; k1(3, 0) = 28. * alpha1 + 17. * beta;
	k1(4, 0) = 6. * alpha1 - 80. * beta; k1(5, 0) = -40. * alpha1 - 6. * beta;
	k1(6, 0) = -6. * alpha1 - 40. * beta; k1(7, 0) = -80. * alpha1 + 6. * beta;
	k1(1, 1) = k1(0, 0); k1(2, 1) = k1(3, 0);
	k1(3, 1) = k1(2, 0); k1(4, 1) = k1(4, 0);
	k1(5, 1) = k1(7, 0); k1(6, 1) = k1(6, 0);
	k1(7, 1) = k1(5, 0); k1(2, 2) = k1(0, 0);
	k1(3, 2) = k1(1, 0); k1(4, 2) = k1(6, 0);
	k1(5, 2) = k1(7, 0); k1(6, 2) = k1(4, 0);
	k1(7, 2) = k1(5, 0); k1(3, 3) = k1(0, 0);
	k1(4, 3) = k1(6, 0); k1(5, 3) = k1(5, 0);
	k1(6, 3) = k1(4, 0); k1(7, 3) = k1(7, 0);
	k1(4, 4) = 48. * alpha1 + 160. * beta; k1(5, 4) = 0;
	k1(6, 4) = -48. * alpha1 + 80. * beta; k1(7, 4) = 0;
	k1(5, 5) = 160. * alpha1 + 48. * beta; k1(6, 5) = 0;
	k1(7, 5) = 80. * alpha1 - 48. * beta; k1(6, 6) = k1(4, 4);
	k1(7, 6) = 0; k1(7, 7) = k1(5, 5);

	MatrixXcd k11 = k1.transpose();
	for (int i = 0; i < 8; i++)k11(i, i) = 0;
	k1 = k1 + k11;

	k2(0, 0) = 6.;
	k2(1, 0) = 2.; k2(1, 1) = 6.;
	k2(2, 0) = 3.; k2(2, 1) = 2.; k2(2, 2) = 6.;
	k2(3, 0) = 2.; k2(3, 1) = 3.; k2(3, 2) = 2.; k2(3, 3) = 6.;
	k2(4, 0) = -6.; k2(4, 1) = -6.; k2(4, 2) = -8.; k2(4, 3) = -8.; k2(4, 4) = 32.;
	k2(5, 0) = -8.; k2(5, 1) = -6.; k2(5, 2) = -6.; k2(5, 3) = -8.; k2(5, 4) = 20.; k2(5, 5) = 32.;
	k2(6, 0) = -8.; k2(6, 1) = -8.; k2(6, 2) = -6.; k2(6, 3) = -6.; k2(6, 4) = 16.; k2(6, 5) = 20.; k2(6, 6) = 32.;
	k2(7, 0) = -6.; k2(7, 1) = -8.; k2(7, 2) = -8.; k2(7, 3) = -6.; k2(7, 4) = 20.; k2(7, 5) = 16.; k2(7, 6) = 20.; k2(7, 7) = 32.;

	MatrixXcd k22 = k2.transpose();
	for (int i = 0; i < 8; i++)k22(i, i) = 0;
	k2 = k2 + k22;

	KE = k1 - k2 * alpha2;

	KE3(0, 0) = 4.;
	KE3(1, 0) = -1.; KE3(1, 1) = 4.;
	KE3(4, 0) = 2.; KE3(4, 1) = 2.; KE3(4, 4) = 16.;

	MatrixXcd k33 = KE3.transpose();
	for (int i = 0; i < 8; i++)k33(i, i) = 0;

	KE3 = KE3 + k33;
	KE3 *= beta2;

}

void MTFwd2D::_E2K2(int NX, int NY, double freq, vector< vector<int> > I8, vector< vector<double> > XY,
	SparseMatrix < complex< double> >& KT, vector< vector<double> > RO, int idx)//idx=0 TM-mode
{
	int ND = 3 * NX * NY + 2 * NX + 2 * NY + 1;
	int NE = NX * NY;

	int L, NJ, NK;
	double a, b;
	complex<double>	lamda, k, tao;

	MatrixXcd KE3, KE;
	vector<Eigen::Triplet<complex<double>>> triplets;

	for (int ix = 0; ix < NX; ix++)
	{
		for (int iy = 0; iy < NY; iy++)
		{
			L = ix * NY + iy;
			a = abs(XY[0][int(I8[3][L])] - XY[0][int(I8[0][L])]);
			b = abs(XY[1][int(I8[0][L])] - XY[1][int(I8[1][L])]);
			if (iy < Num_AirL_TE * idx)
			{
				_Init_Kpara(k, tao, lamda, freq, RHO_AIR, idx);
				_K88(a, b, lamda, tao, k, KE, KE3);
			}
			else
			{
				if (ix >= Numx_expand && ix < NX - Numx_expand)
				{
					_Init_Kpara(k, tao, lamda, freq, RO[ix - Numx_expand][iy - Num_AirL_TE * idx], idx);
					_K88(a, b, lamda, tao, k, KE, KE3);
				}
				else
				{
					_Init_Kpara(k, tao, lamda, freq, RhoL_Xexpand[iy - Num_AirL_TE * idx], idx);
					_K88(a, b, lamda, tao, k, KE, KE3);
				}
			}
			for (int j = 0; j < 8; j++)
			{
				NJ = int(I8[j][L]);
				for (int k = 0; k < 8; k++)
				{
					NK = int(I8[k][L]);
					if (iy == NY - 1)
					{
						triplets.push_back(Triplet<complex<double>>(NJ, NK, KE(j, k) + KE3(j, k)));
					}
					else
					{
						triplets.push_back(Triplet<complex<double>>(NJ, NK, KE(j, k)));
					}
				}
			}
		}
	}
	KT.resize(ND, ND);
	KT.setFromTriplets(triplets.begin(), triplets.end());
	KT.makeCompressed();
}

void MTFwd2D::_Update_para(int nx, int ny, int nsite, vector< vector<double> > Cen_val_RhoL, vector<int> site_index, vector<int> upbound)
{
	//update rosite roexdpand
	Rhosite.resize(nsite);
	for (int i = 0; i < ny; i++)
	{
		double Rhoave = 0;
		/*	for (int j = 0; j < nx; j++)
			{
				Rhoave += Cen_val_RhoL[j][i];
			}
			Rhoave /= nx;*/
		RhoL_Xexpand[i] = Cen_val_RhoL[0][i];
	}
	for (int i = 0; i < nsite; i++)
	{
		Rhosite[i] = Cen_val_RhoL[site_index[i]][upbound[site_index[i]]];
	}

}

void MTFwd2D::Fwd(MTMesh2D& model, MTData2D& data_fwd, int TETM, int inv, int print)
{
	for (int ifre = 0; ifre < data_fwd.nFre; ifre++)
	{
		for (int isite = 0; isite < data_fwd.nSite; isite++)
		{
			data_fwd.Rho_TE[isite][ifre] = 0;
			data_fwd.Rho_TM[isite][ifre] = 0;
			data_fwd.Pha_TE[isite][ifre] = 0;
			data_fwd.Pha_TM[isite][ifre] = 0;
		}
	}
	if (inv == 1)
	{
		if (TETM == 3) Jacobi.resize(4 * data_fwd.nSite * data_fwd.nFre, model.Num_Cen);
		else Jacobi.resize(2 * data_fwd.nSite * data_fwd.nFre, model.Num_Cen);
		Jacobi.setZero();
	}

	if (TETM == 1)//TE
	{
		if ((print == 1) && (myid == 0))cerr << endl;
		_Fwd2DMT(model, data_fwd, 1, inv, TETM, print);
	}
	if (TETM == 2)//TM
	{
		if ((print == 1) && (myid == 0))cerr << endl;
		_Fwd2DMT(model, data_fwd, 0, inv, TETM, print);
	}
	if (TETM == 3)//TE+TM
	{
		if ((print == 1) && (myid == 0))cerr << endl;
		_Fwd2DMT(model, data_fwd, 1, inv, TETM, print);
		_Fwd2DMT(model, data_fwd, 0, inv, TETM, print);
	}
	//
	vector<double> Rho_TE_buf(data_fwd.nSite);
	vector<double> Rho_TM_buf(data_fwd.nSite);
	vector<double> Pha_TE_buf(data_fwd.nSite);
	vector<double> Pha_TM_buf(data_fwd.nSite);
	vector<double> Rho_TE_result(data_fwd.nSite);
	vector<double> Rho_TM_result(data_fwd.nSite);
	vector<double> Pha_TE_result(data_fwd.nSite);
	vector<double> Pha_TM_result(data_fwd.nSite);
	for (int ifre = 0; ifre < data_fwd.nFre; ifre++)
	{
		vector<double> Rho_TE_buf(data_fwd.nSite);
		vector<double> Rho_TM_buf(data_fwd.nSite);
		vector<double> Pha_TE_buf(data_fwd.nSite);
		vector<double> Pha_TM_buf(data_fwd.nSite);
		for (int isite = 0; isite < data_fwd.nSite; isite++)
		{
			Rho_TE_buf[isite] = data_fwd.Rho_TE[isite][ifre];
			Rho_TM_buf[isite] = data_fwd.Rho_TM[isite][ifre];
			Pha_TE_buf[isite] = data_fwd.Pha_TE[isite][ifre];
			Pha_TM_buf[isite] = data_fwd.Pha_TM[isite][ifre];
		}
		MPI_Allreduce(Rho_TE_buf.data(), Rho_TE_result.data(), data_fwd.nSite, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(Rho_TM_buf.data(), Rho_TM_result.data(), data_fwd.nSite, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(Pha_TE_buf.data(), Pha_TE_result.data(), data_fwd.nSite, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(Pha_TM_buf.data(), Pha_TM_result.data(), data_fwd.nSite, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		for (int isite = 0; isite < data_fwd.nSite; isite++)
		{
			data_fwd.Rho_TE[isite][ifre] = Rho_TE_result[isite];
			data_fwd.Rho_TM[isite][ifre] = Rho_TM_result[isite];
			data_fwd.Pha_TE[isite][ifre] = Pha_TE_result[isite];
			data_fwd.Pha_TM[isite][ifre] = Pha_TM_result[isite];
		}
	}
	// 
	if (inv)//jacobi
	{
		vector<double> row_sendbuf(Jacobi.cols());
		vector<double> row_recvbuf(Jacobi.cols());
		for (int i = 0; i < Jacobi.rows(); i++)
		{
			for (int j = 0; j < Jacobi.cols(); j++)
			{
				row_sendbuf[j] = Jacobi.coeffRef(i, j);
			}
			MPI_Allreduce(row_sendbuf.data(), row_recvbuf.data(), Jacobi.cols(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			for (int j = 0; j < Jacobi.cols(); j++)
			{
				Jacobi.coeffRef(i, j) = row_recvbuf[j];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
