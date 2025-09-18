/*
***********************************************************************

MTinverse2D.cpp	(MT inversion)
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

* ***********************************************************************
*/

#include "../include/MTinverse2D.h"
#include "../include/mathfunction.h"

void MTInv2D::Init_Module(string Filename)
{
	ifstream fin(Filename);
	if (!fin) {
		if (myid == 0)cerr << "Error Opening " << Filename << " !" << endl;
		exit(1);
	}
	///
	string tmp1, tmp2, tmp3;
	fin >> tmp1 >> tmp2 >> tmp3 >> OutputFile_Root;
	fin >> tmp1 >> tmp2 >> tmp3 >> OutputTmp_Root;
	fin >> tmp1 >> tmp2 >> tmp3 >> LogFile_Root;
	fin >> tmp1 >> tmp2 >> static_shift;
	fin >> tmp1 >> tmp2 >> alpha;
	fin >> tmp1 >> tmp2 >> damp;
	fin >> tmp1 >> tmp2 >> dt;
	fin >> tmp1 >> tmp2 >> stabilizer;
	fin >> tmp1 >> tmp2 >> ee;
	fin >> tmp1 >> tmp2 >> v2h;
	fin >> tmp1 >> tmp2 >> smdw;
	fin >> tmp1 >> tmp2 >> gamma;
	fin >> tmp1 >> tmp2 >> err_floor_RTE >> err_floor_PTE >> err_floor_RTM >> err_floor_PTM;
	fin >> tmp1 >> tmp2 >> target_misfit;
	fin >> tmp1 >> tmp2 >> max_nls;
	fin >> tmp1 >> tmp2 >> c1 >> c2;
	fin >> tmp1 >> tmp2 >> pre_type;
	fin >> tmp1 >> tmp2 >> cg_type;
	fin >> tmp1 >> tmp2 >> cg_break;
	fin >> tmp1 >> tmp2 >> Max_iteration;
	fin.close();
	if (myid == 0)
	{
		cerr << "--------------------------------Inversion parameters---------------------------------" << endl;
		cerr << "Max iteration:" << Max_iteration << endl;
		cerr << "Starting regularization trad-off: " << alpha << endl;
		cerr << "Data mistfit reduction threshold: " << dt << "%" << endl;
		if (stabilizer == 1)cerr << "Stabilizer: First-order smooth matrix" << endl;
		if (stabilizer == 2)cerr << "Stabilizer: Second-order smooth matrix" << endl;
		if (stabilizer == 3)cerr << "Stabilizer: Minimum support (MS) stabilizer" << endl;
		if (stabilizer == 4)cerr << "Stabilizer: Minimum gradient support (MGS) stabilizer" << endl;
		if (stabilizer > 2)cerr << "Focusing parameter: " << ee << endl;
		cerr << "Regularization trad-off attenuation rate: " << damp << endl;
		cerr << "Precondition gamma: " << gamma << endl;
		cerr << "Wolfe linesearch parameter (c1, c2): " << c1 << " " << c2 << endl;
		if (cg_type == 1)cerr << "Conjugate type: Mix Fletcher-Reeves and Polak-Ribiere method" << endl;
		if (cg_type == 2)cerr << "Conjugate type: Fletcher-Reeves method" << endl;
		if (cg_type == 3)cerr << "Conjugate type: Polak-Ribiere method" << endl;
		if (cg_type == 4)cerr << "Conjugate type: Dai-Kou method" << endl;
		if (pre_type == 1)cerr << "Precondition type: Nonlinear-conjugate damping matrix (identity)" << endl;
		if (pre_type == 2)cerr << "Precondition type: Iterative Gauss-Newton incomplete Hessian matrix" << endl;
		if (pre_type == 3)cerr << "Precondition type: Diagnoal Hessian matrix" << endl;
		if (pre_type == 4)cerr << "Precondition type: Gauss-Newton Hessian matrix" << endl;
		if (myid == 0)cerr << "-------------------------------------------------------------------------------------" << endl;
	}
}

void MTInv2D::NonLinearConjugateGradient(MTFwd2D& fwd, MTMesh2D& model, MTMesh2D& mapr,
	MTData2D& data_obs, MTData2D& data_cal, int TETM)
{
	//mapr
//	if (stabilizer == 3)
//	{
	mapr.sgm1D = model.sgm1D;
	//	}
	//	else
	//	{
	//		mapr.sgm1D.setZero();
	//	}
		//flog
	ofstream flog;
	if (myid == 0)
	{
		flog.open(LogFile_Root);
		flog << left << setw(16) << "niter" << " " << setw(16) << "data_misfit" << " " << setw(16) << "fcost" << " " << setw(16) << "data_misfit(TE)"
			<< " " << setw(16) << "data_misfit(TM)" << " " << setw(16) << "alpha" << endl;
	}
	//err-floor
	for (int i = 0; i < data_obs.nSite; i++)
	{
		for (int j = 0; j < data_obs.nFre; j++)
		{
			if (data_obs.Err_RhoTE[i][j] < err_floor_RTE * 0.00434)data_obs.Err_RhoTE[i][j] = err_floor_RTE * 0.00434;
			if (data_obs.Err_PhaTE[i][j] < err_floor_PTE * 0.285)data_obs.Err_PhaTE[i][j] = err_floor_PTE * 0.285;
			if (data_obs.Err_RhoTM[i][j] < err_floor_RTM * 0.00434)data_obs.Err_RhoTM[i][j] = err_floor_RTM * 0.00434;
			if (data_obs.Err_PhaTM[i][j] < err_floor_PTM * 0.285)data_obs.Err_PhaTM[i][j] = err_floor_PTM * 0.285;
		}
	}
	//
	_Data_Weighting_Matrix_Cons(data_obs, TETM);//Wd
	_Stabilizer_Function(model, mapr);//We
	int n = model.Num_Cen;                             // Dimension
	//
	MTData2D data_cal_err = data_obs;
	MTMesh2D model_ref = model;
	MTData2D data_cal_test = data_cal;
	//
	VectorXd  h(n), r, g(n), g_prev(n), h_prev(n), d(n), d_prev(n), delta_g(n), ddg(n);
	double Data_Misfit, TE_Misfit, TM_Misfit, DF, MF, Data_Misfit_last, DF_last, MF_last, f_ref, eps, alpha0;
	MatrixXd H, Jacobi;
	bool lsh = true;
	//
	double step, beta, beta_PR, beta_FR, beta_DK, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, v;
	step = 0; beta = 0; Bad_reduction = 0;
	step_gn = 1;
	d_prev.setZero();
	double vb, tb, wb;
	//
	//init
	int nit = 0;
	fwd.Fwd(model, data_cal, TETM, 1, 1);
	_Misfit_Calculate(r, Data_Misfit, TM_Misfit, TE_Misfit, data_obs, data_cal, data_cal_err,
		model, mapr, TETM);//rn=d-F(m)

	while (true)
	{
		nit++;
		if (_Iteration_Termination_Condition(nit, Data_Misfit))break;
		//change Sensitivity
	/*	if (TETM == 3)
		{
			for (int ifre = 0; ifre < data_obs.nFre; ifre++)
			{
				for (int isite = 0; isite < data_obs.nSite; isite++)
				{
					int idata = ifre + isite * data_cal.nFre;
					int ndata = data_cal.nSite * data_cal.nFre;
					for (int ii = 0; ii < model.Num_Cen; ii++)
					{
						fwd.Jacobi(idata + 2 * ndata, ii) *= TM_Misfit / (TE_Misfit + 0.001);
						fwd.Jacobi(idata + 3 * ndata, ii) *= TM_Misfit / (TE_Misfit + 0.001);
					}
				}
			}
		}*/
		//aplha
		//nit decrease
		if (nit == 1)alpha0 = alpha;
		else alpha = alpha0 / pow(nit, 1);
		//Lellivier 2012
		/*if (nit > 1)
		{
			wb = Data_Misfit / target_misfit;
			vb = 1 + 1. * abs(wb - 1);
			alpha /= vb;
		}	*/
		//
		if (nit == 1)
		{
			_Quadratic_Approximation(g, H, f_ref, DF, MF, r, fwd.Jacobi, model.sgm1D, mapr.sgm1D, 1);
			eps = DF / alpha - MF;
		}
		else _Quadratic_Approximation(g, H, f_ref, DF, MF, r, fwd.Jacobi, model.sgm1D, mapr.sgm1D, 0);
		_Precondition(h, g, H, fwd.Jacobi, model.Num_xCen, model.Num_yCen);
		///
		DF_last = DF;
		MF_last = MF;
		if (myid == 0)
		{
			_Print_Information(nit, Data_Misfit, f_ref, DF, MF, TE_Misfit, TM_Misfit, model, data_cal, data_cal_err, fwd);
			flog << left << setw(16) << nit << " " << setw(16) << Data_Misfit << " " << setw(16) << f_ref << " " << setw(16) << TE_Misfit
				<< " " << setw(16) << TM_Misfit << " " << setw(16) << alpha << endl;
		}
		if (nit > 1)
		{
			if (lsh)
			{
				Bad_reduction = 0;
				//restart
				v = abs(double(g_prev.transpose() * g)) / double(g.transpose() * g);
				//if (myid == 0)cerr << "Orthogonal condition number: " << v << endl;
				if ((v >= 0.2 && cg_break == 1) || (pre_type == 4))
				{
					if (myid == 0 && cg_break == 1)cerr << "Orthogonal condition failed, restart conjugate at: " << v << endl;
					beta = 0;
				}
				else
				{
					//beta_PR
					delta_g = g - g_prev;
					_Precondition(ddg, delta_g, H, fwd.Jacobi, model.Num_xCen, model.Num_yCen);
					tmp1 = g.transpose() * ddg;
					tmp2 = g_prev.transpose() * h_prev;
					beta_PR = tmp1 / tmp2;
					//beta_FR
					tmp1 = g.transpose() * h;
					tmp2 = g_prev.transpose() * h_prev;
					beta_FR = tmp1 / tmp2;
					//beta_DK
					tmp1 = g.transpose() * ddg;
					tmp2 = delta_g.transpose() * d_prev;
					tmp3 = delta_g.transpose() * ddg;
					tmp4 = d_prev.transpose() * g;
					beta_DK = tmp1 / tmp2 - (tmp3 * tmp4) / (tmp2 * tmp2);
					//test which to use
					if (cg_type == 1)
					{
						if (beta_PR < -beta_FR)
						{
							if (myid == 0)cerr << "Conjugate condition 1, Fletcher-Reeves conjugate method......" << endl;
							beta = -beta_FR;
						}
						else if (abs(beta_PR) <= beta_FR)
						{
							if (myid == 0)cerr << "Conjugate condition 2, Polak-Ribiere conjugate method......" << endl;
							beta = beta_PR;
						}
						else if (beta_PR > beta_FR)
						{
							if (myid == 0)cerr << "Conjugate condition 3, Fletcher-Reeves conjugate method......" << endl;
							beta = beta_FR;
						}
					}
					else if (cg_type == 2)
					{
						if (myid == 0)cerr << "Fletcher-Reeves conjugate method......" << endl;
						beta = beta_FR;
					}
					else if (cg_type == 3)
					{
						if (myid == 0)cerr << "Polak-Ribiere conjugate method......" << endl;
						beta = beta_PR;
					}
					else if (cg_type == 4)
					{
						if (myid == 0)cerr << "Dai-Kou(2013) conjugate method......" << endl;
						beta = beta_DK;
					}
				}
			}
			else
			{
				beta = 0;
			}
			d = -h + beta * d_prev;
		}
		else
			d = -h;
		//
		_Descent_direction_output(nit, d, g, model);
		//update
		d_prev = d; h_prev = h; g_prev = g;
		//steplength
		_Steplength_wolfe(step, lsh, n, f_ref, g, d, r, H, Jacobi, data_obs, data_cal_test, data_cal_err, model_ref, mapr, fwd, TETM, 1);
		//update
		//if (pre_type == 4 && step < 0.1)step = 0.1;
		if (stabilizer > 2)mapr.sgm1D = model.sgm1D;
		model.sgm1D = model.sgm1D + step * d;
		/*if (step <= 1 && step > 0)step_gn = step;
		else step_gn = 1;*/
		step_gn = 1;
		_add_constrain(model);
		if (stabilizer > 2)_Stabilizer_Function(model, mapr);//We
		model_ref.sgm1D = model.sgm1D;
		///
		fwd.Fwd(model, data_cal, TETM, 0, 1);
		_Misfit_Calculate(r, Data_Misfit, TM_Misfit, TE_Misfit, data_obs, data_cal, data_cal_err,
			model, mapr, TETM);//rn=d-F(m)
		fwd.Jacobi = Jacobi;
		if ((Data_Misfit_last - Data_Misfit) < dt * 0.01 * Data_Misfit_last && nit > 1)
		{
			if (myid == 0)
			{
				cerr << "The error reduction rate is less than " << dt << "%, reduce the regularization trade-off!";
				cerr << endl;
			}
			lsh = false;
			Bad_reduction += 1;
			alpha0 *= damp;
		}
		Data_Misfit_last = Data_Misfit;
	}
	//
	if (myid == 0)flog.close();
}

void MTInv2D::_Data_Weighting_Matrix_Cons(MTData2D data_obs, int TETM)
{
	if (myid == 0)cerr << "Constructing Data Weighting Matrix...... " << endl;
	int ndata = data_obs.nFre * data_obs.nSite;
	vector<Triplet<double> > triplets;
	if (TETM == 1)//TE
	{
		Wd.resize(ndata * 2, ndata * 2);
		Wd2n.resize(ndata * 2, ndata * 2);
		for (int isite = 0; isite < data_obs.nSite; isite++)
		{
			for (int ifre = 0; ifre < data_obs.nFre; ifre++)
			{
				int tmp = ifre + isite * data_obs.nFre;
				triplets.push_back(Triplet<double>(tmp, tmp, 1. / (
					data_obs.Err_RhoTE[isite][ifre])));
				triplets.push_back(Triplet<double>(tmp + ndata, tmp + ndata, 1. / (
					data_obs.Err_PhaTE[isite][ifre])));
			}
		}
		Wd.setZero(); // 
		for (const auto& triplet : triplets) {
			Wd.insert(triplet.row(), triplet.col()) = triplet.value();
		}
		Wd.makeCompressed(); // 
	}
	if (TETM == 2)//TM
	{
		Wd.resize(ndata * 2, ndata * 2);
		Wd2n.resize(ndata * 2, ndata * 2);
		for (int isite = 0; isite < data_obs.nSite; isite++)
		{
			for (int ifre = 0; ifre < data_obs.nFre; ifre++)
			{
				int tmp = ifre + isite * data_obs.nFre;
				triplets.push_back(Triplet<double>(tmp, tmp, 1. /
					(data_obs.Err_RhoTM[isite][ifre])));
				triplets.push_back(Triplet<double>(tmp + ndata, tmp + ndata, 1. /
					(data_obs.Err_PhaTM[isite][ifre])));
			}
		}
		Wd.setZero(); // 
		for (const auto& triplet : triplets) {
			Wd.insert(triplet.row(), triplet.col()) = triplet.value();
		}
		Wd.makeCompressed(); // 
	}
	if (TETM == 3)//TE+TM
	{
		Wd.resize(ndata * 4, ndata * 4);
		Wd2n.resize(ndata * 4, ndata * 4);
		for (int isite = 0; isite < data_obs.nSite; isite++)
		{
			for (int ifre = 0; ifre < data_obs.nFre; ifre++)
			{
				int tmp = ifre + isite * data_obs.nFre;
				triplets.push_back(Triplet<double>(tmp, tmp, 1. / (
					data_obs.Err_RhoTE[isite][ifre])));
				triplets.push_back(Triplet<double>(tmp + ndata, tmp + ndata, 1. / (
					data_obs.Err_PhaTE[isite][ifre])));
				triplets.push_back(Triplet<double>(tmp + 2 * ndata, tmp + 2 * ndata, 1. / (
					data_obs.Err_RhoTM[isite][ifre])));
				triplets.push_back(Triplet<double>(tmp + 3 * ndata, tmp + 3 * ndata, 1. /
					(data_obs.Err_PhaTM[isite][ifre])));
			}
		}
		Wd.setZero(); // 
		Wd.setFromTriplets(triplets.begin(), triplets.end());
		Wd.makeCompressed(); // 
	}
}

void MTInv2D::_Stabilizer_Function(MTMesh2D model, MTMesh2D model_apr)
{
	//if (myid == 0)cerr << "Constructing Model Penalty Matrix...... " << endl;
	vector<Triplet<double> > triplets, triplets2, triplets3;
	double dm, mm;
	if (stabilizer == 0)
	{
		Wey.resize(model.Num_Cen, model.Num_Cen);
		Wey.setZero();
		Wey.makeCompressed();
		//
		for (int i = 0; i < model.Num_xCen; i++)
		{
			for (int j = 0; j < model.Num_yCen; j++)
			{
				int np = j + i * model.Num_yCen;
				if (model.fix[i][j] == 0)
				{
					triplets.push_back(Triplet<double>(np, np, 1));//EVENLY WEIGHT
				}
			}
		}
		Wex.resize(model.Num_Cen, model.Num_Cen);
		Wex.setFromTriplets(triplets.begin(), triplets.end());
		Wex.makeCompressed();
	}
	else if (stabilizer == 1)
	{
		//vertical smooth
		for (int i = 0; i < model.Num_xCen; i++)
		{
			for (int j = 0; j < model.Num_yCen - 1; j++)
			{
				int np = j + i * model.Num_yCen;
				if (model.fix[i][j] == 0 && model.fix[i][j + 1] == 0)
				{
					if (static_shift == 1)
					{
						if (!findNumber(model.site_index, i) && j != model.upbound[i])
						{
							triplets.push_back(Triplet<double>(np, np, 1. * v2h));//EVENLY WEIGHT
							triplets.push_back(Triplet<double>(np, np + 1, -1. * v2h));//EVENLY WEIGHT
						}
					}
					else
					{
						triplets.push_back(Triplet<double>(np, np, 1. * v2h));//EVENLY WEIGHT
						triplets.push_back(Triplet<double>(np, np + 1, -1. * v2h));//EVENLY WEIGHT
					}
				}
			}
		}
		Wey.resize(model.Num_Cen, model.Num_Cen);
		Wey.setFromTriplets(triplets.begin(), triplets.end());
		Wey.makeCompressed();
		//horizontal smooth
		double v, h, ssm;
		for (int nl = 0; nl < model.Num_yCen; nl++)//L layer
		{
			v = abs(model.yNode[nl + 1] - model.yNode[nl]);
			for (int np = 0; np < model.Num_xCen - 1; np++)
			{
				int tmp = nl + np * model.Num_yCen;
				h = abs(model.xNode[np + 1] - model.xNode[np]);
				ssm = max(smdw * v / h, 1.0);
				if (model.fix[np][nl] == 0 && model.fix[np + 1][nl] == 0)
				{
					if (static_shift == 1)
					{
						if (!findNumber(model.site_index, np) && !findNumber(model.site_index, np + 1) && nl != model.upbound[np])
						{
							triplets2.push_back(Triplet<double>(np, np, 1. * ssm));//EVENLY WEIGHT
							triplets2.push_back(Triplet<double>(np, np + 1, -1. * ssm));//EVENLY WEIGHT
						}
					}
					else
					{
						triplets2.push_back(Triplet<double>(tmp, tmp, 1. * ssm));
						triplets2.push_back(Triplet<double>(tmp, tmp + model.Num_yCen, -1. * ssm));
					}
				}
			}
		}
		Wex.resize(model.Num_Cen, model.Num_Cen);
		Wex.setFromTriplets(triplets2.begin(), triplets2.end());
		Wex.makeCompressed();
	}
	else if (stabilizer == 2)
	{
		//vertical smooth
		for (int i = 0; i < model.Num_xCen; i++)
		{
			for (int j = 0; j < model.Num_yCen - 1; j++)
			{
				int np = j + i * model.Num_yCen;
				if (model.fix[i][j] == 0 && model.fix[i][j + 1] == 0)
				{
					if (static_shift == 1)
					{
						if (!findNumber(model.site_index, i) && j != model.upbound[i])
						{
							triplets.push_back(Triplet<double>(np, np, 1. * v2h));//EVENLY WEIGHT
							triplets.push_back(Triplet<double>(np, np + 1, -1. * v2h));//EVENLY WEIGHT
						}
					}
					else
					{
						triplets.push_back(Triplet<double>(np, np, 1. * v2h));//EVENLY WEIGHT
						triplets.push_back(Triplet<double>(np, np + 1, -1. * v2h));//EVENLY WEIGHT
					}
				}
			}
		}
		Wey.resize(model.Num_Cen, model.Num_Cen);
		Wey.setFromTriplets(triplets.begin(), triplets.end());
		Wey.makeCompressed();
		Wey = Wey.transpose() * Wey;
		//horizontal smooth
		double v, h, ssm;
		for (int nl = 0; nl < model.Num_yCen; nl++)//L layer
		{
			v = abs(model.yNode[nl + 1] - model.yNode[nl]);
			for (int np = 0; np < model.Num_xCen - 1; np++)
			{
				int tmp = nl + np * model.Num_yCen;
				h = abs(model.xNode[np + 1] - model.xNode[np]);
				ssm = max(smdw * v / h, 1.0);
				if (model.fix[np][nl] == 0 && model.fix[np + 1][nl] == 0)
				{
					if (static_shift == 1)
					{
						if (!findNumber(model.site_index, np) && !findNumber(model.site_index, np + 1) && nl != model.upbound[np])
						{
							triplets.push_back(Triplet<double>(np, np, 1. * ssm));//EVENLY WEIGHT
							triplets.push_back(Triplet<double>(np, np + 1, -1. * ssm));//EVENLY WEIGHT
						}
					}
					else
					{
						triplets2.push_back(Triplet<double>(tmp, tmp, 1. * ssm));
						triplets2.push_back(Triplet<double>(tmp, tmp + model.Num_yCen, -1. * ssm));
					}
				}
			}
		}
		Wex.resize(model.Num_Cen, model.Num_Cen);
		Wex.setFromTriplets(triplets2.begin(), triplets2.end());
		Wex = Wex.transpose() * Wex;
		Wex.makeCompressed();
	}
	else if (stabilizer == 3)
	{
		Wey.resize(model.Num_Cen, model.Num_Cen);
		Wey.setZero();
		Wey.makeCompressed();
		//
		for (int i = 0; i < model.Num_xCen; i++)
		{
			for (int j = 0; j < model.Num_yCen; j++)
			{
				int np = j + i * model.Num_yCen;
				dm = pow(model.sgm1D(np) - model_apr.sgm1D(np), 2);
				if (model.fix[i][j] == 0)
				{
					triplets.push_back(Triplet<double>(np, np, 1. / sqrt(ee * ee + dm)));//EVENLY WEIGHT
				}
			}
		}
		Wex.resize(model.Num_Cen, model.Num_Cen);
		Wex.setFromTriplets(triplets.begin(), triplets.end());
		Wex.makeCompressed();
	}
	else if (stabilizer == 4)
	{
		//vertical 
		for (int i = 0; i < model.Num_xCen; i++)
		{
			for (int j = 0; j < model.Num_yCen - 1; j++)
			{
				int np = j + i * model.Num_yCen;
				dm = model.sgm1D(np) - model.sgm1D(np + 1);
				mm = model.sgm1D(np) - model_apr.sgm1D(np);
				if (model.fix[i][j] == 0 && model.fix[i][j + 1] == 0)
				{
					//	triplets.push_back(Triplet<double>(np, np, dm / sqrt(dm * dm + ee * ee) / sqrt(mm * mm + 0.1)));//EVENLY WEIGHT
					triplets.push_back(Triplet<double>(np, np, -1. * sqrt(1. / (dm * dm + ee * ee))));//EVENLY WEIGHT
					triplets.push_back(Triplet<double>(np, np + 1, 1. * sqrt(1. / (dm * dm + ee * ee))));//EVENLY WEIGHT
				}
			}
		}
		Wey.resize(model.Num_Cen, model.Num_Cen);
		Wey.setFromTriplets(triplets.begin(), triplets.end());
		Wey.makeCompressed();
		//horizontal 
		double v, h, ssm;
		for (int nl = 0; nl < model.Num_yCen; nl++)//L layer
		{
			v = abs(model.yNode[nl + 1] - model.yNode[nl]);
			for (int np = 0; np < model.Num_xCen - 1; np++)
			{
				int tmp = nl + np * model.Num_yCen;
				h = abs(model.xNode[np + 1] - model.xNode[np]);
				ssm = max(smdw * v / h, 1.0);
				dm = model.sgm1D(tmp) - model.sgm1D(tmp + model.Num_yCen);
				mm = model.sgm1D(tmp) - model_apr.sgm1D(tmp);
				if (model.fix[np][nl] == 0 && model.fix[np + 1][nl] == 0)
				{
					//	triplets2.push_back(Triplet<double>(tmp, tmp, dm / sqrt(dm * dm + ee * ee) / sqrt(mm * mm + 0.1)));
					triplets2.push_back(Triplet<double>(tmp, tmp, -1. * sqrt(1. / (dm * dm + ee * ee))));
					triplets2.push_back(Triplet<double>(tmp, tmp + model.Num_yCen, 1. * sqrt(1. / (dm * dm + ee * ee))));
				}
			}
		}
		Wex.resize(model.Num_Cen, model.Num_Cen);
		Wex.setFromTriplets(triplets2.begin(), triplets2.end());
		Wex.makeCompressed();
	}
	//	We = Wex.transpose() * Wex + Wey.transpose() * Wey;
}

void MTInv2D::_Misfit_Calculate(VectorXd& Rn, double& Data_Misfit, double& TM_Misfit, double& TE_Misfit,
	MTData2D data_obs, MTData2D data_cal, MTData2D& data_cal_err, MTMesh2D model, MTMesh2D mapr, int TETM)
{
	FTE_rho_err.resize(data_obs.nSite); FTM_rho_err.resize(data_obs.nSite); FTE_pha_err.resize(data_obs.nSite); FTM_pha_err.resize(data_obs.nSite);
	for (int isite = 0; isite < data_obs.nSite; isite++)
	{
		FTE_rho_err[isite].resize(data_obs.nFre); FTM_rho_err[isite].resize(data_obs.nFre); FTE_pha_err[isite].resize(data_obs.nFre); FTM_pha_err[isite].resize(data_obs.nFre);
	}
	site_err.resize(data_obs.nSite); site_err_TE.resize(data_obs.nSite); site_err_TM.resize(data_obs.nSite);
	site_err.assign(site_err.size(), 0.0); site_err_TE.assign(site_err_TE.size(), 0.0); site_err_TM.assign(site_err_TM.size(), 0.0);
	Data_Misfit = 0; TM_Misfit = 0; TE_Misfit = 0;
	int exclude = 0; int exclude1 = 0; int exclude2 = 0;
	int ndata = data_obs.nFre * data_obs.nSite;
	if (TETM == 1)//TE only
	{
		Rn.resize(ndata * 2);
		for (int isite = 0; isite < data_obs.nSite; isite++)
		{
			for (int ifre = 0; ifre < data_obs.nFre; ifre++)
			{
				int tmp = ifre + isite * data_obs.nFre;
				Rn(tmp) = -data_cal.Rho_TE[isite][ifre] + data_obs.Rho_TE[isite][ifre];
				Rn(tmp + ndata) = -PI * (data_cal.Pha_TE[isite][ifre] + data_obs.Pha_TE[isite][ifre]) / 180.;
				data_cal_err.Rho_TE[isite][ifre] = -data_cal.Rho_TE[isite][ifre] + data_obs.Rho_TE[isite][ifre];
				data_cal_err.Pha_TE[isite][ifre] = -1. * (data_cal.Pha_TE[isite][ifre] + data_obs.Pha_TE[isite][ifre]);
				if (data_obs.Err_RhoTE[isite][ifre] > 9000.)exclude += 1;
				if (data_obs.Err_PhaTE[isite][ifre] > 9000.)exclude += 1;
				///
				FTE_rho_err[isite][ifre] = abs(Rn(tmp) * Wd.coeff(tmp, tmp)); FTE_pha_err[isite][ifre] = abs(Rn(tmp + ndata) * Wd.coeff(tmp + ndata, tmp + ndata));
				FTM_rho_err[isite][ifre] = 0; FTM_pha_err[isite][ifre] = 0;
				site_err_TE[isite] += pow(Rn(tmp) * Wd.coeff(tmp, tmp), 2) + pow(Rn(tmp + ndata) * Wd.coeff(tmp + ndata, tmp + ndata), 2);
			}
			site_err_TE[isite] = sqrt(site_err_TE[isite] / (2 * data_obs.nFre));
			site_err[isite] = site_err_TE[isite];
		}
		for (int i = 0; i < ndata * 2; i++)
		{
			Data_Misfit += pow(Rn(i) * Wd.coeff(i, i), 2);
		}
		Data_Misfit = sqrt(Data_Misfit / (ndata * 2 - exclude));
		TE_Misfit = Data_Misfit;
	}
	else if (TETM == 2)//TM only
	{
		Rn.resize(ndata * 2);
		for (int isite = 0; isite < data_obs.nSite; isite++)
		{
			for (int ifre = 0; ifre < data_obs.nFre; ifre++)
			{
				int tmp = ifre + isite * data_obs.nFre;
				Rn(tmp) = -data_cal.Rho_TM[isite][ifre] + data_obs.Rho_TM[isite][ifre];
				data_cal_err.Rho_TM[isite][ifre] = -data_cal.Rho_TM[isite][ifre] + data_obs.Rho_TM[isite][ifre];
				Rn(tmp + ndata) = -(data_cal.Pha_TM[isite][ifre] - (180 - data_obs.Pha_TM[isite][ifre]));// *PI / 180.;
				data_cal_err.Pha_TM[isite][ifre] = 180 + (data_cal.Pha_TM[isite][ifre] - (180 - data_obs.Pha_TM[isite][ifre])) * PI / 180.;
				if (data_obs.Err_RhoTM[isite][ifre] > 9000.)exclude += 1;
				if (data_obs.Err_PhaTM[isite][ifre] > 9000.)exclude += 1;
				///
				FTE_rho_err[isite][ifre] = 0; FTE_pha_err[isite][ifre] = 0;
				FTM_rho_err[isite][ifre] = abs(Rn(tmp) * Wd.coeff(tmp, tmp)); FTM_pha_err[isite][ifre] = abs(Rn(tmp + ndata) * Wd.coeff(tmp + ndata, tmp + ndata));
				site_err_TM[isite] += pow(Rn(tmp) * Wd.coeff(tmp, tmp), 2) + pow(Rn(tmp + ndata) * Wd.coeff(tmp + ndata, tmp + ndata), 2);
			}
			site_err_TM[isite] = sqrt(site_err_TM[isite] / (2 * data_obs.nFre));
			site_err[isite] = site_err_TM[isite];
		}
		for (int i = 0; i < ndata * 2; i++)
		{
			Data_Misfit += pow(Rn(i) * Wd.coeff(i, i), 2);
		}
		Data_Misfit = sqrt(Data_Misfit / (ndata * 2 - exclude));
		TM_Misfit = Data_Misfit;
	}
	else if (TETM == 3)//TE+TM
	{
		Rn.resize(ndata * 4);

		for (int isite = 0; isite < data_obs.nSite; isite++)
		{
			for (int ifre = 0; ifre < data_obs.nFre; ifre++)
			{
				int tmp = ifre + isite * data_obs.nFre;
				Rn(tmp) = -data_cal.Rho_TE[isite][ifre] + data_obs.Rho_TE[isite][ifre];
				Rn(tmp + ndata) = -1. * (data_cal.Pha_TE[isite][ifre] + data_obs.Pha_TE[isite][ifre]);
				Rn(tmp + ndata * 2) = -data_cal.Rho_TM[isite][ifre] + data_obs.Rho_TM[isite][ifre];
				Rn(tmp + ndata * 3) = -(data_cal.Pha_TM[isite][ifre] - (180 - data_obs.Pha_TM[isite][ifre]));
				data_cal_err.Rho_TE[isite][ifre] = -data_cal.Rho_TE[isite][ifre] + data_obs.Rho_TE[isite][ifre];
				data_cal_err.Pha_TE[isite][ifre] = -(data_cal.Pha_TE[isite][ifre] + data_obs.Pha_TE[isite][ifre]);
				data_cal_err.Rho_TM[isite][ifre] = -data_cal.Rho_TM[isite][ifre] + data_obs.Rho_TM[isite][ifre];
				data_cal_err.Pha_TM[isite][ifre] = 180 + (data_cal.Pha_TM[isite][ifre] - (180 - data_obs.Pha_TM[isite][ifre]));
				if (data_obs.Err_RhoTM[isite][ifre] > 9000.)exclude2 += 1;
				if (data_obs.Err_PhaTM[isite][ifre] > 9000.)exclude2 += 1;
				if (data_obs.Err_RhoTE[isite][ifre] > 9000.)exclude1 += 1;
				if (data_obs.Err_PhaTE[isite][ifre] > 9000.)exclude1 += 1;
				///
				FTE_rho_err[isite][ifre] = abs(Rn(tmp) * Wd.coeff(tmp, tmp)); FTE_pha_err[isite][ifre] = abs(Rn(tmp + ndata) * Wd.coeff(tmp + ndata, tmp + ndata));
				FTM_rho_err[isite][ifre] = abs(Rn(tmp + 2 * ndata) * Wd.coeff(tmp + 2 * ndata, tmp + 2 * ndata));
				FTM_pha_err[isite][ifre] = abs(Rn(tmp + 3 * ndata) * Wd.coeff(tmp + 3 * ndata, tmp + 3 * ndata));
				site_err_TE[isite] += pow(Rn(tmp) * Wd.coeff(tmp, tmp), 2) + pow(Rn(tmp + ndata) * Wd.coeff(tmp + ndata, tmp + ndata), 2);
				site_err_TM[isite] += pow(Rn(tmp + ndata * 2) * Wd.coeff(tmp + ndata * 2, tmp + ndata * 2), 2) + pow(Rn(tmp + ndata * 3) * Wd.coeff(tmp + ndata * 3, tmp + ndata * 3), 2);
				site_err[isite] = site_err_TE[isite] + site_err_TM[isite];
			}
			site_err_TM[isite] = sqrt(site_err_TM[isite] / (2 * data_obs.nFre));
			site_err_TE[isite] = sqrt(site_err_TE[isite] / (2 * data_obs.nFre));
			site_err[isite] = sqrt(site_err[isite] / (4 * data_obs.nFre));
		}
		for (int i = 0; i < ndata * 4; i++)Data_Misfit += pow(Rn(i) * Wd.coeff(i, i), 2);
		Data_Misfit = sqrt(Data_Misfit / (ndata * 4 - exclude1 - exclude2));
		for (int i = 0; i < ndata * 2; i++)TE_Misfit += pow(Rn(i) * Wd.coeff(i, i), 2);
		TE_Misfit = sqrt(TE_Misfit / (ndata * 2 - exclude1));
		for (int i = ndata * 2; i < ndata * 4; i++)TM_Misfit += pow(Rn(i) * Wd.coeff(i, i), 2);
		TM_Misfit = sqrt(TM_Misfit / (ndata * 2 - exclude2));
	}
}

void MTInv2D::_Quadratic_Approximation(VectorXd& g, MatrixXd& H, double& Fai_ref, double& DF, double& MF,
	VectorXd Rn, MatrixXd& J, VectorXd m, VectorXd mapr, int calc)
{
	double tmp1, tmp2;
	/*tmp1 = Rn.transpose() * Rn;
	tmp2 = Rn.transpose() * Wd.transpose() * Wd * Rn;*/
	Wd2n = Wd.transpose() * Wd;// *tmp1 / tmp2;
	Fai_ref = 0; DF = 0; MF = 0;
	MF = alpha * (m - mapr).transpose() * Wex.transpose() * Wex * (m - mapr);
	MF += alpha * (m - mapr).transpose() * Wey.transpose() * Wey * (m - mapr);
	DF = Rn.transpose() * Wd2n * Rn;
	Fai_ref = DF + MF;
	if (calc == 1)
	{
		//Hessian
		H = 2 * J.transpose() * Wd2n * J;
		H += 2 * alpha * (Wex.transpose() * Wex + Wey.transpose() * Wey);
		//Gradient
		g = -2 * J.transpose() * Wd2n * Rn + 2 * alpha * (Wex.transpose() * Wex + Wey.transpose() * Wey) * (m - mapr);
		//Fai
	}
}

void MTInv2D::_Print_Information(double nit, double Data_Misfit, double f_ref, double DF, double MF, double TE_Misfit, double TM_Misfit,
	MTMesh2D model, MTData2D data_cal, MTData2D data_cal_err, MTFwd2D fwd)
{
	cerr << endl;
	cerr << "********************************INVERSION********************************" << endl;
	cerr << "NUMBER OF ITERATION: " << nit << endl;
	cerr << "Data Misfit: " << Data_Misfit << " " << "Function cost: " << f_ref << endl;
	cerr << "Data funtion: " << DF << " " << "Model function: " << MF << endl;
	cerr << "TE Misfit: " << TE_Misfit << " " << "TM Misfit: " << TM_Misfit << endl;
	cerr << "Reguralization trad-off: " << alpha << " " << "Damping factor: " << damp << endl;
	cerr << "***************************************************************************" << endl;
	model.OutputMeshFile(OutputFile_Root + "Inv_model");
	model.OutputMeshFile(OutputTmp_Root + "inv_nlcg_itr" + to_string(int(nit)));
	data_cal.OutputFile(fwd.OutputFile_Root + "tem-in-cal.dat");
	data_cal.OutputFile(fwd.OutputTmp_Root + "tem-in-cal-itr" + to_string(int(nit)) + ".dat");
	data_cal_err.OutputFile(fwd.OutputFile_Root + "tem-in-cal_err.dat");
	data_cal_err.OutputFile(fwd.OutputTmp_Root + "tem-in-cal-itr" + to_string(int(nit)) + "_err.dat");
	_site_err_output(fwd.OutputTmp_Root + "itr" + to_string(int(nit)) + "_site_err.dat", data_cal);
	_site_err_output(fwd.OutputFile_Root + "site_err.dat", data_cal);
	_fre_err_output(fwd.OutputTmp_Root + "itr" + to_string(int(nit)) + "_freq_err.dat", data_cal);
	_fre_err_output(fwd.OutputFile_Root + "freq_err.dat", data_cal);
}

bool MTInv2D::_Iteration_Termination_Condition(double nit, double misfit)
{
	if (nit > Max_iteration)
	{
		if (myid == 0)cerr << "Maximum number of iterations reached, iteration terminated!" << endl;
		return true;
	}
	else if (misfit < target_misfit)
	{
		if (myid == 0)cerr << "Target misfit achieved, iteration terminated!" << endl;
		return true;
	}
	/*else if (alpha < 0.1)
	{
		if (myid == 0)cerr << "Regularization trad-off less than 0.1, iteration terminated!" << endl;
		return true;
	}*/
	else if (Bad_reduction > 300)
	{
		if (myid == 0)cerr << "Ten consecutive iterations failed to meet the error reduction criterion, iteration terminated!" << endl;
		return true;
	}
	else
	{
		return false;
	}
}

void MTInv2D::_Precondition(VectorXd& h, VectorXd& g, MatrixXd& H, MatrixXd& J, int nx, int ny)
{
	if (pre_type == 0)
	{
		h = g * 0.01;
		//constraint
		VectorXd abs_h = h.array().abs();
		double mean = abs_h.mean();
		double stddev = sqrt((abs_h.array() - mean).square().sum() / abs_h.size());
		for (int i = 0; i < h.size(); i++)
		{
			if (abs(h(i)) > (mean + 5 * stddev))h(i) = abs(h(i) + 1e-6) / (h(i) + 1e-6) * (mean + 5 * stddev);
		}
	}
	else if (pre_type == 1)
	{
		SparseMatrix<double> Cl;//preconditionier
		vector<Triplet<double> > triplets;
		for (int i = 0; i < g.size(); i++)
		{
			triplets.push_back(Triplet<double>(i, i, gamma));
		}
		Cl.resize(g.size(), g.size());
		Cl.setFromTriplets(triplets.begin(), triplets.end());
		Cl.makeCompressed();
		Cl = Cl + alpha * (Wex.transpose() * Wex + Wey.transpose() * Wey);
		SparseLU< SparseMatrix<double>>solver;
		solver.compute(Cl);
		h = solver.solve(g);
		//constraint
		VectorXd abs_h = h.array().abs();
		double mean = abs_h.mean();
		double stddev = sqrt((abs_h.array() - mean).square().sum() / abs_h.size());
		for (int i = 0; i < h.size(); i++)
		{
			if (abs(h(i)) > (mean + 5 * stddev))h(i) = abs(h(i) + 1e-6) / (h(i) + 1e-6) * (mean + 5 * stddev);
		}
	}
	else if (pre_type == 2)
	{
		BiCGSTAB< MatrixXd, DiagonalPreconditioner<double>> solver;
		solver.setTolerance(1e-1);
		solver.compute(H);
		h = solver.solve(g);
	}
	else if (pre_type == 3)
	{
		for (int ix = 1; ix < nx - 1; ix++)
		{
			for (int iy = 0; iy < ny; iy++)
			{
				int i = iy + ix * ny;
				h(i) = g(i) / (H.coeffRef(i, i));
			}
		}
		for (int iy = 0; iy < ny; iy++)
		{
			int i1 = iy + 0 * ny;
			int i2 = iy + 1 * ny;
			int i3 = iy + (nx - 1) * ny;
			int i4 = iy + (nx - 2) * ny;
			h(i1) = h(i2);
			h(i3) = h(i4);
		}
		//constraint
		VectorXd abs_h = h.array().abs();
		double mean = abs_h.mean();
		double stddev = sqrt((abs_h.array() - mean).square().sum() / abs_h.size());
		for (int i = 0; i < h.size(); i++)
		{
			if (abs(h(i)) > (mean + 5 * stddev))h(i) = abs(h(i) + 1e-6) / (h(i) + 1e-6) * (mean + 5 * stddev);
		}
	}
	else
	{
		h = H.ldlt().solve(g);
		VectorXd abs_h = h.array().abs();
		double mean = abs_h.mean();
		double stddev = sqrt((abs_h.array() - mean).square().sum() / abs_h.size());
		for (int i = 0; i < h.size(); i++)
		{
			if (abs(h(i)) > (mean + 5 * stddev))h(i) = abs(h(i) + 1e-6) / (h(i) + 1e-6) * (mean + 5 * stddev);
		}
	}
}

void MTInv2D::_Steplength_wolfe(double& steplength, bool& lsh, int n, double finit, VectorXd& g, VectorXd d, VectorXd& r, MatrixXd& H, MatrixXd& Jacobi,
	MTData2D data_obs, MTData2D& data_cal, MTData2D& data_cal_err, MTMesh2D& m, MTMesh2D mapr, MTFwd2D fwd, int TETM, int print)
{
	double amax = 1.0;
	double Data_Misfit, TE_Misfit, TM_Misfit, DF, MF;
	double f_prev, df_prev, fk, dfk, f_cub, df0, f0;
	double a_prev, a, abset;
	MTMesh2D mk = m;
	VectorXd gk, rk, sgm1D_prev; MatrixXd Hk;
	gk = g;
	Hk = H;
	int nls = 0;
	a = 0;
	a_prev = a;
	f_prev = finit; f0 = finit; f_cub = finit;
	df0 = g.transpose() * d;
	df_prev = g.transpose() * d;
	lsh = false;
	//Init
	if (myid == 0)
	{
		cerr << "Bracketing the step length in Wolfe-linesearch......" << endl;
		cerr << left << setw(16) << "nls" << " " << setw(16) << "steplength" << " " << setw(16) << "fcost" << " " << setw(16) << "d_fcost" << " " << setw(16) << "fcost_prev" << endl;
		cerr << left << setw(16) << nls << " " << setw(16) << a << " " << setw(16) << f0 << " " << setw(16) << df0 << " " << setw(16) << f_prev << endl;
	}
	if (pre_type < 4)a = a - double(gk.transpose() * d) / double(d.transpose() * Hk * d);//init_value
	else a = step_gn;
	if ((a < 0.01 || a>10) && pre_type > 1)
	{
		if (myid == 0)cerr << "Warning!! Abnormal steplength status!! Please adjust the gamma !!" << endl;
	}
	//
	while (nls < max_nls)
	{
		nls++;
		//
		sgm1D_prev = mk.sgm1D;
		mk.sgm1D = m.sgm1D + a * d;
		_add_constrain(mk);
		if (stabilizer > 2)_Stabilizer_Function(mk, mapr);//We
		//fwd
		fwd.Fwd(mk, data_cal, TETM, 1, 0);
		_Misfit_Calculate(rk, Data_Misfit, TM_Misfit, TE_Misfit, data_obs, data_cal, data_cal_err,
			mk, mapr, TETM);//rn=d-F(m)
		_Quadratic_Approximation(gk, Hk, fk, DF, MF, rk, fwd.Jacobi, mk.sgm1D, mapr.sgm1D, 1);//new g,H,f
		dfk = gk.transpose() * d;
		if (myid == 0)cerr << left << setw(16) << nls << " " << setw(16) << a << " " << setw(16) << fk << " " << setw(16) << dfk << " " << setw(16) << f_prev << endl;
		//
		if ((fk > f0 + c1 * a * df0) || ((fk >= f_prev) && nls > 1))
		{
			if (myid == 0)
			{
				cerr << "Linesearched a mis-match the Sufficient decrease condition(1)......" << endl;
				cerr << "enter cube interpolation in alo: " << a_prev << " flo: " << f_prev << " ahi: " << a << " fhi: " << fk << endl;
			}
			abset = _Zoom(f0, df0, a, a_prev, fk, f_prev, dfk, df_prev, fwd, m, mk, mapr, data_obs, data_cal, data_cal_err, d, gk, Hk, Jacobi, TETM);
			H = Hk;
			g = gk;
			Jacobi = fwd.Jacobi;
			lsh = true;
			break;
		}
		if (abs(dfk) <= -c2 * df0)
		{
			if (myid == 0)cerr << "Linesearch satisfied Wolfe standard, accept the steplength......" << endl;
			abset = a;
			H = Hk;
			g = gk;
			Jacobi = fwd.Jacobi;
			lsh = true;
			break;
		}
		if (dfk >= 0)
		{
			if (myid == 0)
			{
				cerr << "Function dirivate lager than 0, enter cube interpolation......" << endl;
				cerr << "enter cube interpolation in alo: " << a << " flo: " << fk << " ahi: " << a_prev << " fhi: " << f_prev << endl;
			}
			abset = _Zoom(f0, df0, a_prev, a, f_prev, fk, df_prev, dfk, fwd, m, mk, mapr, data_obs, data_cal, data_cal_err, d, gk, Hk, Jacobi, TETM);
			H = Hk;
			g = gk;
			Jacobi = fwd.Jacobi;
			lsh = true;
			break;
		}
		//a = a - double(gk.transpose() * d) / double(d.transpose() * Hk * d); ...Rodi iteration...
		a_prev = a;
		f_prev = fk;
		df_prev = dfk;
		a = a * 1.5;
		//	
	}
	steplength = abset;
	if (nls > max_nls - 1)
	{
		if (myid == 0)cerr << "Linesearch Failed! Reached Max Linesearch count!" << endl;
		steplength = abset;
	}
}

double MTInv2D::_Cubic_Interpolate(double x0, double x1, double f0, double f1, double fp0, double fp1)
{
	double x2;
	double d1 = fp0 + fp1 - 3 * (f0 - f1) / (x0 - x1);
	double d2 = sign(x1 - x0) * pow(d1 * d1 - fp1 * fp0, 0.5);
	//interpolation
	x2 = x1 - (x1 - x0) * (fp1 + d2 - d1) / (fp1 - fp0 + 2 * d2);
	return x2;
}

double MTInv2D::_Zoom(double f0, double df0, double ahi, double alo, double f_ahi, double f_alo, double df_ahi, double df_alo,
	MTFwd2D fwd, MTMesh2D m, MTMesh2D mk, MTMesh2D mapr, MTData2D data_obs, MTData2D data_cal, MTData2D data_cal_err,
	VectorXd d, VectorXd& g, MatrixXd& H, MatrixXd& Jacobi, int TETM)
{
	double Data_Misfit, TE_Misfit, TM_Misfit, f_aj, df_aj, DF, MF;
	VectorXd gk, rk, gk_lo;
	MatrixXd Hk, Hk_lo, Jacobi_lo;
	Hk_lo = H; gk_lo = g; Jacobi_lo = Jacobi;
	double aj;
	double abest;
	int count = 0;
	if (myid == 0)cerr << "-------------------------------------------------------------------------------------" << endl;
	while (true)
	{
		count++;
		if (count > max_nls)
		{
			aj = alo;
			if (pre_type != 1)H = Hk_lo;
			g = gk_lo;
			Jacobi = Jacobi_lo;
			break;
		}
		if (abs(ahi - alo) < 0.05 * alo || ahi < 1e-4)
		{
			if (myid == 0)cerr << "Alpha bracketing zoom is too small, accept the lower bound steplength at: " << alo << endl;
			aj = alo;
			H = Hk_lo;
			g = gk_lo;
			Jacobi = Jacobi_lo;
			break;
		}
		aj = _Cubic_Interpolate(alo, ahi, f_alo, f_ahi, df_alo, df_ahi);
		if (myid == 0)cerr << "Cubic interpolate times: " << count << " at steplength: " << aj << endl;
		mk.sgm1D = m.sgm1D + aj * d;
		_add_constrain(mk);
		if (stabilizer > 2)_Stabilizer_Function(mk, mapr);//We
		//fwd
		fwd.Fwd(mk, data_cal, TETM, 1, 0);
		_Misfit_Calculate(rk, Data_Misfit, TM_Misfit, TE_Misfit, data_obs, data_cal, data_cal_err,
			mk, mapr, TETM);//rn=d-F(m)
		_Quadratic_Approximation(gk, Hk, f_aj, DF, MF, rk, fwd.Jacobi, mk.sgm1D, mapr.sgm1D, 1);//new g,H,f
		df_aj = gk.transpose() * d;
		///
		if ((f_aj > f0 + c1 * aj * df0) || (f_aj >= f_alo))
		{
			if (myid == 0)cerr << "Alpha mis-match the Sufficient Decrease condition(1), adjust the higher boundary ahi: " << aj << " fai: " << f_aj << endl;
			ahi = aj;
			df_ahi = df_aj;
			f_ahi = f_aj;
		}
		else
		{
			if (abs(df_aj) <= -c2 * df0)
			{
				if (myid == 0)cerr << "Alpha satisfied the Wolfe condition at steplength: " << aj << " fcost: " << f_aj << endl;
				abest = aj;
				H = Hk;
				g = gk;
				Jacobi = fwd.Jacobi;
				break;
			}
			if (myid == 0)cerr << "Alpha mis-match the Curtrave condition(2), adjust the lower boundary alo: " << aj << " flo: " << f_aj << endl;
			if (df_aj * (ahi - alo) >= 0)
				ahi = alo;
			alo = aj;
			df_alo = df_aj;
			f_alo = f_aj;
			Hk_lo = Hk;
			gk_lo = gk;
			Jacobi_lo = fwd.Jacobi;
		}
	}
	if (myid == 0)cerr << "-------------------------------------------------------------------------------------" << endl;
	return aj;
}

void MTInv2D::_site_err_output(string Filename, MTData2D data_cal)
{
	ofstream fout(Filename);
	fout << "xsite err err_TE err_TM" << endl;
	for (int i = 0; i < data_cal.nSite; i++)
	{
		fout << data_cal.xSite[i] * 0.001 << " " << site_err[i] << " " << site_err_TE[i] << " " << site_err_TM[i] << endl;
	}
	fout.close();
}

void MTInv2D::_fre_err_output(string Filename, MTData2D data_cal)
{
	ofstream fout(Filename);
	for (int i = 0; i < data_cal.nSite; i++)
	{
		for (int j = 0; j < data_cal.nFre; j++)
		{
			fout << data_cal.xSite[i] * 0.001 << " " << data_cal.Freq[j] << " " << FTE_rho_err[i][j] << " " << FTE_pha_err[i][j] << " " << FTM_rho_err[i][j] << " " << FTM_pha_err[i][j] << endl;
		}
	}
	fout.close();
}

void MTInv2D::_add_constrain(MTMesh2D& model)
{
	int np;
	for (int ix = 0; ix < model.Num_xCen; ix++)
	{
		for (int iy = 0; iy < model.Num_yCen; iy++)
		{
			np = iy + ix * model.Num_yCen;
			if (model.fix[ix][iy] == 1)model.sgm1D(np) = -model.rho_air;
			else if (model.fix[ix][iy] == 2)model.sgm1D(np) = -model.rho_water;
			else
			{
				if (model.sgm1D(np) < -5)model.sgm1D(np) = -5;
				if (model.sgm1D(np) > 2)model.sgm1D(np) = 2;
			}
		}
	}
}

void MTInv2D::_lanczos(MatrixXd& A, int k, VectorXd& eigenvalues, MatrixXd& eigenvectors)
{
	int n = A.rows();
	std::vector<Eigen::VectorXd> q; // 保存正交基
	Eigen::VectorXd q0 = Eigen::VectorXd::Zero(n);
	Eigen::VectorXd q1 = Eigen::VectorXd::Random(n).normalized();
	q.push_back(q1);

	Eigen::VectorXd v(n), w(n);
	std::vector<double> alpha, beta;

	for (int i = 0; i < k; ++i) {
		v = A * q.back();
		if (i > 0) v -= beta.back() * q[q.size() - 2];
		alpha.push_back(q.back().dot(v));
		v -= alpha.back() * q.back();
		beta.push_back(v.norm());

		if (beta.back() > 1e-10) {
			q.push_back(v / beta.back());
		}
		else {
			break; // 基向量不足
		}
	}

	// 构造三对角矩阵
	Eigen::MatrixXd T = Eigen::MatrixXd::Zero(k, k);
	for (int i = 0; i < k; ++i) {
		T(i, i) = alpha[i];
		if (i < k - 1) {
			T(i, i + 1) = beta[i];
			T(i + 1, i) = beta[i];
		}
	}

	// 计算T的特征值和特征向量
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(T);
	eigenvalues = solver.eigenvalues();
	eigenvectors = solver.eigenvectors();
}

void MTInv2D::_Descent_direction_output(double nit, VectorXd  d, VectorXd g, MTMesh2D model)
{
	int np;
	ofstream fout(OutputTmp_Root + "gradient_nlcg_itr" + to_string(int(nit)) + ".dat");
	for (int ix = 0; ix < model.Num_xCen; ix++)
	{
		for (int iy = 0; iy < model.Num_yCen; iy++)
		{
			np = iy + ix * model.Num_yCen;
			fout << model.xCen[ix] / 1000. << " " << model.yCen[iy] / 1000. << " " << g(np) << " " << d(np) << endl;
		}
	}
	fout.close();
}