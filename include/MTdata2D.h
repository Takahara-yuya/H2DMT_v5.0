#pragma once
#ifndef MTDATA2D_H_
#define MTDATA2D_H_

#include "CommonHeaders.h"

class MTData2D
{
public:

	int nFre, nSite;

	std::vector<double> xSite, ySite, Freq, AF, FreqL;

	std::vector< std::vector<double> > Rho_TE, Rho_TM, Pha_TE, Pha_TM, Err_RhoTE, Err_RhoTM, Err_PhaTE, Err_PhaTM;

	void Init_MTData2D_File(std::string Filename1, std::string Filename2);

	void Init_MTData2D_Empty(int nfre, int nsite, std::vector<double> xsite, std::vector<double> ysite, std::vector<double> freq);

	void ReadDataFile(std::string Filename);

	void ReadSiteFile(std::string Filename);

	void OutputFile(std::string Filename);

};
#endif
