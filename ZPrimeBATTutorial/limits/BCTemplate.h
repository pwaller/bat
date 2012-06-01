#ifndef __BCTEMPLATE__H
#define __BCTEMPLATE__H

#include <string>

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

	// return the number of bins
	int GetNBins()
	{ return fNBins; }; 

 private:

		// the efficiency of the contribution
		double fEfficiency; 

		// the template histogram
		TH1D* fHistogram;

		// number of bins in the histogram
		int fNBins;

		// channel name
		std::string fChannelName;

		// process name
		std::string fProcessName;

};
// ---------------------------------------------------------

#endif

