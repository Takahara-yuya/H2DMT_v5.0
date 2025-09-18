#pragma once
#ifndef MTINVERSE2D_H_
#define MTINVERSE2D_H_
#include "CommonHeaders.h"
#include"MTdata2D.h"
#include "MTmesh2D.h"
#include "MTfwd2D.h"

class MTInv2D
{
public:
	//! Constructor
	MTInv2D() {

	};
	//! Destructor
	~MTInv2D() {

	};
	//
	double ee = 0.1;
	//MPI
	int np, myid;

	string OutputFile_Root, LogFile_Root, OutputTmp_Root;
	int stabilizer, Max_iteration, max_nls, pre_type, cg_type, cg_break;
	double target_misfit, err_floor_RTE, err_floor_PTE, err_floor_RTM, err_floor_PTM;
	vector<double> site_err, site_err_TE, site_err_TM;
	vector<vector<double>> FTE_rho_err, FTM_rho_err, FTE_pha_err, FTM_pha_err;
	double alpha, damp, dt, v2h, gamma, smdw, t1, t2, static_shift, step_gn;
	double c1, c2;
	VectorXd Rn;//data misfit
	SparseMatrix<double> Wd, Wd2n, Wm, We, Wex, Wey;//data weighting ;//model weighting matrix
	//terminate condition
	int Bad_reduction;

	void Init_Module(string Filename);

	void NonLinearConjugateGradient(MTFwd2D& fwd, MTMesh2D& model, MTMesh2D& mapr,
		MTData2D& data_obs, MTData2D& data_cal, int TETM);

	void _Data_Weighting_Matrix_Cons(MTData2D data_obs, int TETM);

	void _Stabilizer_Function(MTMesh2D model, MTMesh2D model_apr);

	void _Misfit_Calculate(VectorXd& Rn, double& Data_Misfit, double& TM_Misfit, double& TE_Misfit
		, MTData2D data_obs, MTData2D data_cal, MTData2D& data_cal_err, MTMesh2D model, MTMesh2D mapr, int TETM);

	void _Precondition(VectorXd& h, VectorXd& g, MatrixXd& H, MatrixXd& J, int nx, int ny);

	void _Steplength_wolfe(double& steplength, bool& lsh, int n, double finit, VectorXd& g, VectorXd d, VectorXd& r, MatrixXd& H, MatrixXd& Jacobi,
		MTData2D data_obs, MTData2D& data_cal, MTData2D& data_cal_err, MTMesh2D& m, MTMesh2D mapr, MTFwd2D fwd, int TETM, int print);

	bool _Iteration_Termination_Condition(double nit, double misfit);

	void _Print_Information(double nit, double Data_Misfit, double f_ref, double DF, double MF, double TE_Misfit, double TM_Misfit,
		MTMesh2D model, MTData2D data_cal, MTData2D data_cal_err, MTFwd2D fwd);

	void _Quadratic_Approximation(VectorXd& g, MatrixXd& H, double& Fai_ref, double& DF, double& MF,
		VectorXd Rn, MatrixXd& J, VectorXd m, VectorXd mapr, int calc);

	double _Cubic_Interpolate(double x0, double x1, double f0, double f1, double fp0, double fp1);

	double _Zoom(double f0, double df0, double ahi, double alo, double f_ahi, double f_alo, double df_ahi, double df_alo,
		MTFwd2D fwd, MTMesh2D m, MTMesh2D mk, MTMesh2D mapr, MTData2D data_obs, MTData2D data_cal, MTData2D data_cal_err,
		VectorXd d, VectorXd& g, MatrixXd& H, MatrixXd& Jacobi, int TETM);

	void _site_err_output(string Filename, MTData2D data_cal);

	void _fre_err_output(string Filename, MTData2D data_cal);

	void _add_constrain(MTMesh2D& model);

	void _lanczos(MatrixXd& A, int k, VectorXd& eigenvalues, MatrixXd& eigenvectors);

	void _Descent_direction_output(double nit, VectorXd  d, VectorXd g, MTMesh2D model);

protected:

	bool findNumber(vector<int>& vec, int target)
	{
		for (int num : vec) {
			if (num == target) {
				return true; // 
			}
		}
		return false; // 
	}

};
#endif
