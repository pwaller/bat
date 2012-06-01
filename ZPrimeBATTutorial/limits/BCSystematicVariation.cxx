#include "BCSystematicVariation.h"

#include <TH1D.h>

// ---------------------------------------------------------
BCSystematicVariation::BCSystematicVariation(const char* channelname, const char* systematicname, int nprocesses)
{  
	fChannelName = channelname;
	fSystematicName = systematicname;
	for (int i = 0; i < nprocesses; ++i) {
		fHistogramUpContainer.push_back(0);
		fHistogramDownContainer.push_back(0);
	}
};

// ---------------------------------------------------------
BCSystematicVariation::~BCSystematicVariation()
{
};

// ---------------------------------------------------------
