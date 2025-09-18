/*
***********************************************************************

MTdata2D.h	(Constructing MT data)
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

version 5.0.0
Dec 04, 2024


***********************************************************************
*/

#include "../include/MTdata2D.h"

void MTData2D::Init_MTData2D_File(std::string Filename1, std::string Filename2)
{
	MTData2D::ReadDataFile(Filename1);
	MTData2D::ReadSiteFile(Filename2);
}

void MTData2D::Init_MTData2D_Empty(int nfre, int nsite, std::vector<double> xsite, std::vector<double> ysite, std::vector<double> freq)
{
	nSite = nsite; nFre = nfre;

	xSite.resize(nSite); ySite.resize(nSite); Freq.resize(nFre); FreqL.resize(nFre);
	Rho_TE.resize(nSite); Rho_TM.resize(nSite); Pha_TE.resize(nSite); Pha_TM.resize(nSite);

	for (int i = 0; i < nSite; i++)
	{
		Rho_TE[i].resize(nFre); Rho_TM[i].resize(nFre); Pha_TE[i].resize(nFre); Pha_TM[i].resize(nFre);
	}

	xSite.assign(xsite.begin(), xsite.end());
	ySite.assign(ysite.begin(), ysite.end());
	Freq.assign(freq.begin(), freq.end());

	for (int i = 0; i < nFre; i++)
	{
		FreqL[i] = pow(10., Freq[i]);
	}

	for (int i = 0; i < nSite; i++)
	{
		for (int j = 0; j < nFre; j++)
		{
			Rho_TE[i][j] = 0; Rho_TM[i][j] = 0; Pha_TE[i][j] = 0; Pha_TM[i][j] = 0;
		}
	}

}

void MTData2D::ReadDataFile(std::string Filename)
{
	std::ifstream fin(Filename);
	if (!fin) {
		std::cerr << "Can't open observed data file!\n";
		exit(1);
	}
	fin >> nSite >> nFre;

	xSite.resize(nSite); ySite.resize(nSite); Freq.resize(nFre); FreqL.resize(nFre);
	Rho_TE.resize(nSite); Rho_TM.resize(nSite); Pha_TE.resize(nSite); Pha_TM.resize(nSite);
	Err_RhoTE.resize(nSite); Err_RhoTM.resize(nSite);  Err_PhaTE.resize(nSite);  Err_PhaTM.resize(nSite);

	for (int i = 0; i < nSite; i++)
	{
		Rho_TE[i].resize(nFre); Rho_TM[i].resize(nFre); Pha_TE[i].resize(nFre); Pha_TM[i].resize(nFre);
		Err_RhoTE[i].resize(nFre); Err_RhoTM[i].resize(nFre);  Err_PhaTE[i].resize(nFre);  Err_PhaTM[i].resize(nFre);
	}

	for (int i = 0; i < nSite; i++)
	{
		for (int j = 0; j < nFre; j++)
		{
			fin >> xSite[i] >> Freq[j] >> Rho_TE[i][j] >> Pha_TE[i][j] >> Rho_TM[i][j] >> Pha_TM[i][j] >> 
				Err_RhoTE[i][j] >> Err_RhoTM[i][j] >> Err_PhaTE[i][j] >> Err_PhaTM[i][j];
			//Check Data
			//check phase
			// all change to the first quadrant
			if (Pha_TE[i][j] > 90 && Pha_TE[i][j] <= 180)Pha_TE[i][j] = 180 - Pha_TE[i][j];
			if (Pha_TE[i][j] > 180 && Pha_TE[i][j] <= 270)Pha_TE[i][j] = Pha_TE[i][j] - 180;
			if (Pha_TE[i][j] > 270 && Pha_TE[i][j] <= 360)Pha_TE[i][j] = 360 - Pha_TE[i][j];
			if (Pha_TM[i][j] > 90 && Pha_TM[i][j] <= 180)Pha_TM[i][j] = 180 - Pha_TM[i][j];
			if (Pha_TM[i][j] > 180 && Pha_TM[i][j] <= 270)Pha_TM[i][j] = Pha_TM[i][j] - 180;
			if (Pha_TM[i][j] > 270 && Pha_TM[i][j] <= 360)Pha_TM[i][j] = 360 - Pha_TM[i][j];
			//
		}
	}

	for (int i = 0; i < nFre; i++)
	{
		FreqL[i] = pow(10., Freq[i]);
	}

	fin.close();
}

void MTData2D::ReadSiteFile(std::string Filename)
{
	std::ifstream fin(Filename);
	if (!fin) {
		std::cerr << "Can't open site file!\n";
		exit(1);
	}

	for (int i = 0; i < nSite; i++)
	{
		std::string data;
		std::getline(fin, data);
		std::stringstream iss(data);
		iss >> xSite[i] >> ySite[i];
		xSite[i] *= 1000.;
	}

	fin.close();
}

void MTData2D::OutputFile(std::string Filename)
{
    ofstream fout(Filename);
	for (int i = 0; i < nSite; i++)
	{
		for (int j = 0; j < nFre; j++)
		{
			fout << xSite[i] / 1000. << " " << Freq[j] << " " << Rho_TE[i][j] << " " << abs(Pha_TE[i][j])
				<< " " << Rho_TM[i][j] << " " << 180 - Pha_TM[i][j] << endl;
		}
	}

}

