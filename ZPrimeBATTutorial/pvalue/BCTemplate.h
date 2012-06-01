#ifndef __BCTEMPLATE__H
#define __BCTEMPLATE__H

#include <string>
#include <TObjArray.h>

class TH1D;

// ---------------------------------------------------------
class BCTemplate
{
 public:
	
	// Constructors and destructor
	BCTemplate(const char* channelname, const char* processname);
	~BCTemplate();
	
	// setters
	
	// set efficiency
	void SetEfficiency(double eff)
	{ fEfficiency = eff; }; 
	
	// set histogram
	void SetHistogram(TH1D* hist); 

	// set histogram
	void SetNextHistogram(TH1D* hist); 

	// set histograms
	void SetTemplateHistograms(TObjArray* obj); 
	
	// getters
	
	// return the name of the channel
	std::string GetChannelName() {
		return fChannelName; }; 

	// return the name of the process
	std::string GetProcessName() {
		return fProcessName; }; 

	// return efficiency
	double GetEfficiency()
	{ return fEfficiency; };

	// return histogram	
	TH1D* GetHistogram()
	{ return fHistogram; }; 

	// return next histogram	
	TH1D* GetNextHistogram()
	{ return fNextHistogram; }; 

	TObjArray* GetTemplateHistograms()
	{ return fTemplateHistograms; }; 

	// return the number of bins
	int GetNBins()
	{ return fNBins; }; 

 private:

		// the efficiency of the contribution
		double fEfficiency; 

		// the template histogram
		TH1D* fHistogram;

		// the next template histogram
		TH1D* fNextHistogram;

		// the template histograms
		TObjArray* fTemplateHistograms;

		// number of bins in the histogram
		int fNBins;

		// channel name
		std::string fChannelName;

		// process name
		std::string fProcessName;

};
// ---------------------------------------------------------

#endif

