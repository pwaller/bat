#include "BCTemplate.h"

#include <TH1D.h>
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
void BCTemplate::SetNextHistogram(TH1D* hist)
{
	fNextHistogram = hist; 
}; 

// ---------------------------------------------------------
void BCTemplate::SetTemplateHistograms(TObjArray* obj)
{
	fTemplateHistograms = obj; 
	if ( obj->GetEntries() > 0 )
	   fHistogram = ((TH1D*)(TObjArray*)obj->At(0));
	if (fHistogram)
		fNBins = fHistogram->GetNbinsX();
}; 
