#include "../include/MTfileread.h"

void MTstartupRead(string filename, int& TETM, int& Auto, int& site_expand, int& grid_expand, int& refine,
	int& Fwd_only, double& Automesh_value, bool& site_Ref, bool& local_Ref, double& lr, bool& Marine,
	string& Fwd_Setting_File, string& Inverse_Setting_File, string& Mesh_File,
	string& Observed_data_File, string& Site_File, string& Zdepth_mesh_File, string& topo_File)
{
	ifstream fin(filename);

	string input_line, tmp1, tmp2, tmp3;
	fin >> tmp1 >> tmp2 >> input_line;

	for (char& c : input_line) {
		c = std::toupper(c);
	}

	if (input_line == "TE")
	{
		TETM = 1;
	}
	else if (input_line == "TM")
	{
		TETM = 2;
	}
	else if (input_line == "TE+TM" || input_line == "TM+TE" || input_line == "TMTE" || input_line == "TETM")
	{
		TETM = 3;
	}
	else
	{
		exit(1);
	}

	fin >> tmp1 >> tmp2 >> input_line;
	for (char& c : input_line) {
		c = std::toupper(c);
	}
	if (input_line == "AUTO" || input_line == "AUTOMESH")Auto = 1;
	else
	{
		Auto = 2;
	}

	fin >> tmp1 >> tmp2 >> tmp3 >> input_line;
	for (char& c : input_line) {
		c = std::toupper(c);
	}
	if (input_line == "FWD" || input_line == "FORWARD" || input_line == "FWDONLY")Fwd_only = 1;
	else
	{
		Fwd_only = 2;
	}

	fin >> tmp1 >> tmp2 >> tmp3 >> Fwd_Setting_File;
	fin >> tmp1 >> tmp2 >> tmp3 >> Inverse_Setting_File;
	fin >> tmp1 >> tmp2 >> tmp3 >> Mesh_File;
	fin >> tmp1 >> tmp2 >> tmp3 >> Observed_data_File;
	fin >> tmp1 >> tmp2 >> tmp3 >> Site_File;
	if (Auto == 1)MTAutomeshRead(Mesh_File, Zdepth_mesh_File, topo_File, site_expand, grid_expand, refine, Automesh_value, site_Ref, local_Ref, lr, Marine);

	fin.close();
}

void MTAutomeshRead(string filename, string& Zdepth_mesh_File, string& topo_File,
	int& site_Expand, int& grid_Expand, int& refine, double& Automesh_value, bool& site_Ref, bool& local_Ref, double& lr, bool& Marine)
{
	string tmp1, tmp2, tmp3, tmp4;
	ifstream fin2(filename);
	fin2 >> tmp1 >> tmp2 >> tmp3;
	for (char& c : tmp3) {
		c = std::toupper(c);
	}
	if (tmp3 == "MARINE")
		Marine = true;
	else
		Marine = false;
	fin2 >> tmp1 >> tmp2 >> tmp3 >> tmp4;
	for (char& c : tmp4) {
		c = std::toupper(c);
	}
	if (tmp4 == "YES")
		site_Ref = true;
	else
		site_Ref = false;
	fin2 >> tmp1 >> tmp2 >> tmp3 >> tmp4;
	for (char& c : tmp4) {
		c = std::toupper(c);
	}
	if (tmp4 == "YES")
		local_Ref = true;
	else
		local_Ref = false;
	fin2 >> tmp1 >> tmp2 >> tmp3 >> lr;
	fin2 >> tmp1 >> tmp2 >> Zdepth_mesh_File;
	fin2 >> tmp1 >> tmp2 >> topo_File;
	fin2 >> tmp1 >> tmp2 >> site_Expand;
	fin2 >> tmp1 >> tmp2 >> grid_Expand;
	fin2 >> tmp1 >> tmp2 >> refine;
	fin2 >> tmp1 >> tmp2 >> Automesh_value;
	fin2.close();
}
