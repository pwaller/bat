#include "BCTemplate.h"

#include <TH1D.h>

// ---------------------------------------------------------
BCTemplate::BCTemplate(const char* channelname, const char* processname) : fEfficiency(0)
																																				 , fHistogram(0)
																																				 , fNBins(0)
{  
	fChannelName = channelname;
	fProcessName = processname;
};

// ---------------------------------------------------------
BCTemplate::~BCTemplate()
{
	// debugKK
	//	if (fHistogram)
	//		delete fHistogram; 
};

// ---------------------------------------------------------
void BCTemplate::SetHistogram(TH1D* hist)
{
	fHistogram = hist; 
	if (hist)
		fNBins = fHistogram->GetNbinsX();
}; 

// ---------------------------------------------------------
