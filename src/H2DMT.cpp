/*
***********************************************************************

H2DMT.cpp
This file is the main program of H2DMT.

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
#include "../include/CommonHeaders.h"
#include "../include/MTmesh2D.h"
#include "../include/MTdata2D.h"
#include "../include/MTfwd2D.h"
#include "../include/MTfileread.h"
#include "../include/MTinverse2D.h"

int main(int argc, char** argv)
{
	//mpi initialize//
	int myid, np;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	// Copyright and license notice:
	if (myid == 0)
	{
		cerr << endl;
		cerr << "***************************************************************************" << endl << endl;
		cerr << "                              H2DMT  V5.6.0                   " << endl << endl;
		cerr << "                       Copyright 2024, Zuwei Huang         " << endl << endl;

		cerr << "  H2DMT is free software: you can redistribute it and/or modify it under the " << endl <<
			"  terms of the GNU Lesser General Public License as published by the Free " << endl <<
			"  Software Foundation, version 3 of the License. " << endl << endl;

		cerr << "  H2DMT is distributed in the hope that it will be useful, but WITHOUT ANY " << endl <<
			"  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS " << endl <<
			"  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for " << endl <<
			"  more details. " << endl << endl;
		cerr << "***************************************************************************" << endl << endl;
	}
	std::this_thread::sleep_for(std::chrono::seconds(1));
	// The user should specify the configuration file as input to H2DMT:
	//************************** Read Startup File *************************************
	int TETM, Auto, site_Expand, grid_Expand, refine, Fwd_only;
	string Fwd_Setting_File, Inverse_Setting_File, Mesh_File;
	string Model_File, Observed_data_File, Site_File, Zdepth_mesh_File, topo_File;
	double Automesh_value,lr;
	bool Marine, site_Ref, local_Ref ;
	string startup = "Startup";

	MTstartupRead(startup, TETM, Auto, site_Expand, grid_Expand, refine,
		Fwd_only, Automesh_value, site_Ref, local_Ref, lr, Marine,
		Fwd_Setting_File, Inverse_Setting_File, Mesh_File,
		Observed_data_File, Site_File, Zdepth_mesh_File, topo_File);
	if ((TETM == 1) && (myid == 0))cerr << "Polarization mode: TE mode" << endl;
	if ((TETM == 2) && (myid == 0))cerr << "Polarization mode: TM mode" << endl;
	if ((TETM == 3) && (myid == 0))cerr << "Polarization mode: TE+TM mode" << endl;
	//****************************************************************************
	//****************************Parameter Initialization*********************************
	MTData2D MTdata_obs, MTdata_cal;
	MTMesh2D MTmodel, MTmapr;
	MTFwd2D MTFwd;
	MTInv2D MTInv;
	//MPI
	MTFwd.myid = myid; MTFwd.np = np;
	MTInv.myid = myid; MTInv.np = np;
	///
	MTdata_obs.Init_MTData2D_File(Observed_data_File, Site_File);
	MTdata_cal.Init_MTData2D_Empty(MTdata_obs.nFre, MTdata_obs.nSite, MTdata_obs.xSite,
		MTdata_obs.ySite, MTdata_obs.Freq);
	//
	if (Auto == 1)
	{
		if (myid == 0)cerr << "Mesh settings: Automatic mesh generation..." << endl;
		MTmodel.AutoMesh(MTdata_obs, Zdepth_mesh_File, topo_File
			, Automesh_value, site_Expand, grid_Expand, refine, site_Ref, local_Ref, lr, Marine);
		MTmapr.Init_MTMesh2D_Empty(MTmodel.Num_xCen, MTmodel.Num_yCen, MTmodel.xNode, MTmodel.yNode, 0.0);
	}
	else
	{
		if (myid == 0)cerr << "Mesh settings from meshfile: " << Mesh_File << endl;
		MTmodel.Init_MTMesh2D_File(Mesh_File, MTdata_obs);
		MTmapr.Init_MTMesh2D_Empty(MTmodel.Num_xCen, MTmodel.Num_yCen, MTmodel.xNode, MTmodel.yNode, 0.0);
	}
	if (myid == 0)
	{
		cerr << "------------------------------Grid and Data parameters-------------------------------" << endl;
		cerr << "Grid nx center:" << MTmodel.Num_xCen << endl;
		cerr << "Grid ny center:" << MTmodel.Num_yCen << endl;
		cerr << "Observed data site number :" << MTdata_obs.nSite << " Frequency number :" << MTdata_obs.nFre << endl;
	}
	//
	MTFwd.Init_Module(Fwd_Setting_File, MTdata_obs, MTmodel);
	MTInv.Init_Module(Inverse_Setting_File);
	//-------------------------------------------------------------------------------------------------//
	//
	if (Fwd_only == 1)
	{
		if (myid == 0)cerr << "FORWARD MODELING ONLY! OUTPUT DATA: " << MTFwd.OutputFile_Root << "tem-in-fwd.dat" << endl;
		MTFwd.Fwd(MTmodel, MTdata_cal, TETM, 0, 1);
		MTdata_cal.OutputFile(MTFwd.OutputFile_Root + "tem-in-fwd.dat");
		MTmodel.OutputMeshFile(MTFwd.OutputFile_Root + "model");
	}
	else
	{
		if (myid == 0)cerr << "Inversion Method : WNLCG ( Jorge Nocedal <<Numerical Optimization 2nd Edition>> )" << endl;
		if (myid == 0)cerr << "-------------------------------------------------------------------------------------" << endl;
		if (myid == 0)cerr << endl;
		MTInv.NonLinearConjugateGradient(MTFwd, MTmodel, MTmapr, MTdata_obs, MTdata_cal, TETM);
	}
	MPI_Finalize();
	exit(0);
}
