#pragma once
#ifndef FILEREAD_H_
#define FILEREAD_H_

#include "CommonHeaders.h"
#include "MTdata2D.h"
#include "MTmesh2D.h"
#include "MTfwd2D.h"

void MTstartupRead(string filename, int& TETM, int& Auto, int& site_expand, int& grid_expand, int& refine,
	double& Fwd_only, double& Automesh_value, bool& site_Ref, bool &Marine,
	string& Fwd_Setting_File, string& Inverse_Setting_File, string& Mesh_File,
	string& Observed_data_File, string& Site_File, string& Zdepth_mesh_File, string& topo_File);

void MTAutomeshRead(string filename, string& Zdepth_mesh_File, string& topo_File, int& site_Expand, int& grid_Expand, int& refine, double& Automesh_value, bool& site_Ref, bool &Marine);

#endif