#include <BCMultiTemplateFitter.h>
#include <BCChannel.h>
#include <BCProcess.h>
#include <BCTemplate.h>
#include <BCSystematic.h>
#include <BCSystematicVariation.h>

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

#include <TCanvas.h>
#include <THStack.h>
#include <TH1D.h>
#include <TAxis.h> 

#include <iostream>
#include <fstream>

// ---------------------------------------------------------
BCMultiTemplateFitter::BCMultiTemplateFitter() : BCModel()
																							 , fNChannels(0)
																							 , fNProcesses(0)
																							 , fNSystematics(0)
																							 , fFlagEfficiencyConstraint(false)
{  
};

// ---------------------------------------------------------
BCMultiTemplateFitter::~BCMultiTemplateFitter()
	// default destructor
{
	for (int i = 0; i < fNChannels; ++i) 
		delete fChannelContainer.at(i);
};

// ---------------------------------------------------------
int BCMultiTemplateFitter::GetChannelIndex(const char* name)
{
	// loop over all channels and compare names
	for (int i = 0; i < fNChannels; ++i) {
		// get channel
		BCChannel* channel = GetChannel(i); 

		// compare names
		if (!channel->GetName().compare(name))
			return i; 
	}
	
	// if channel does not exist, return -1
	return -1;
}

// ---------------------------------------------------------
int BCMultiTemplateFitter::GetProcessIndex(const char* name)
{
	// loop over all processs and compare names
	for (int i = 0; i < fNProcesses; ++i) {
		// get process
		BCProcess* process = GetProcess(i); 

		// compare names
		if (!process->GetName().compare(name))
			return i; 
	}

	// if process does not exist, return -1
	return -1;
}

// ---------------------------------------------------------
int BCMultiTemplateFitter::GetSystematicIndex(const char* name)
{
	// loop over all systematics and compare names
	for (int i = 0; i < fNSystematics; ++i) {
		// get systematic
		BCSystematic* systematic = GetSystematic(i); 

		// compare names
		if (!systematic->GetName().compare(name))
			return i; 
	}

	// if process does not exist, return -1
	return -1;
}

// ---------------------------------------------------------
int BCMultiTemplateFitter::SetTemplate(const char* channelname, const char* processname, TH1D hist, double efficiency)
{
	// get channel index
	int channelindex = GetChannelIndex(channelname);

	// check if channel exists
	if (channelindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Channel does not exist.");
			return -1;
	}

	// get process index
	int processindex = GetProcessIndex(processname);

	// check if process exists
	if (processindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Process does not exist.");
			return -1;
	}

	// get channel
	BCChannel* channel = GetChannel(channelindex);
	
	// get template
	BCTemplate* bctemplate = channel->GetTemplate(processindex);

	// normalize histogram
	if (hist.Integral())
		hist.Scale(1.0 / hist.Integral()); 

	// remove statistics box
	hist.SetStats(kFALSE);

	// set color and fill style
	hist.SetFillColor(2 + processindex);
	hist.SetFillStyle(1001);

	// create new histogram
	TH1D* temphist = new TH1D(hist);

	// set histogram
	bctemplate->SetHistogram(temphist);
	
	// set efficiency
	bctemplate->SetEfficiency(efficiency);

	// no error
	return 1;
}

// ---------------------------------------------------------
int BCMultiTemplateFitter::SetData(const char* channelname, TH1D hist)
{
	int channelindex = GetChannelIndex(channelname);

	// check if channel exists
	if (channelindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Channel does not exist.");
			return -1;
	}

	// get channel
	BCChannel* channel = GetChannel(channelindex);
	
	// get template
	BCTemplate* data = channel->GetData();

	// remove statistics box
	hist.SetStats(kFALSE);

	// set marker
	hist.SetMarkerStyle(20);
	hist.SetMarkerSize(1.1);

	// remove old data set if it exists
	if (data->GetHistogram()) {
		delete data->GetHistogram();
		data->SetHistogram(0);
	}

	// set histogram
	data->SetHistogram(new TH1D(hist));
	
	// no error
	return 1;	
}

// ---------------------------------------------------------
int BCMultiTemplateFitter::AddChannel(const char* name)
{  
	// check if channel exists
	for (int i = 0; i < fNChannels; ++i) { 
		// compare names
		if (GetChannelIndex(name) >= 0) {
      BCLog::OutWarning("BCMultitemplateFitter::AddChannel() : Channel with this name exists already.");
			return -1;
		}
	}

	// create new channel
	BCChannel* channel = new BCChannel(name);

	// create new data
	BCTemplate* bctemplate = new BCTemplate(channel->GetName().c_str(), "data");

	// add data
	channel->SetData(bctemplate); 

	// add process templates
	for (int i = 0; i < fNProcesses; ++i) {
		// get process
		BCProcess* process = GetProcess(i); 

		// create new template
		BCTemplate* bctemplate = new BCTemplate(name, process->GetName().c_str()); 

		// add template
		channel->AddTemplate(bctemplate);
	}

	// loop over all systematics
	for (int i = 0; i < fNSystematics; ++i) {
		// get systematic
		BCSystematic* systematic = GetSystematic(i); 

		// create new systematic variation
		BCSystematicVariation* variation = new BCSystematicVariation(name, systematic->GetName().c_str(), fNProcesses); 

		// add systematic variation
		channel->AddSystematicVariation(variation);
	}

	// add channel
	fChannelContainer.push_back(channel);

	// increase number of channels
	fNChannels++;

	// no error
	return 1;
};

// ---------------------------------------------------------
int BCMultiTemplateFitter::AddProcess(const char* name, double nmin, double nmax)
{  
	// check if process exists
	for (int i = 0; i < fNProcesses; ++i) { 
		// compare names
		if (GetProcessIndex(name) >= 0) {
      BCLog::OutWarning("BCMultitemplateFitter::AddProcess() : Process with this name exists already.");
			return -1;
		}
	}

	// create new process
	BCProcess* process = new BCProcess(name);

	// add process
	fProcessContainer.push_back(process);

	// add process templates
	for (int i = 0; i < fNChannels; ++i) {
		// get channel
		BCChannel* channel = GetChannel(i); 

		// create new template
		BCTemplate* bctemplate = new BCTemplate(channel->GetName().c_str(), name); 

		// add template
		channel->AddTemplate(bctemplate);

		// loop over all systematic 
		for (int j = 0; j < fNSystematics; ++j) {
			// get systematic variation
			BCSystematicVariation* variation = channel->GetSystematicVariation(j);
			
			// add histogram
			variation->AddHistograms(0, 0);
		}
	}

	// increase number of processes
	fNProcesses++;

	// add parameter index to container
	fProcessParIndexContainer.push_back(GetNParameters());

	// add parameter
	AddParameter(name, nmin, nmax);

	// no error
	return 1;
};

// ---------------------------------------------------------
int BCMultiTemplateFitter::AddSystematic(const char* name, double min, double max)
{
	// check if systematic exists
	for (int i = 0; i < fNSystematics; ++i) { 
		// compare names
		if (GetSystematicIndex(name) >= 0) {
      BCLog::OutWarning("BCMultitemplateFitter::AddSystematic() : Systematic with this name exists already.");
			return -1;
		}
	}

	// create new systematic
	BCSystematic* systematic = new BCSystematic(name);

	// add systematic
	fSystematicContainer.push_back(systematic);

	// add systematic variations
	for (int i = 0; i < fNChannels; ++i) {
		// get channel
		BCChannel* channel = GetChannel(i); 

		// create new systematic variation
		BCSystematicVariation* variation = new BCSystematicVariation(channel->GetName().c_str(), name, fNProcesses); 

		// add systematic variation
		channel->AddSystematicVariation(variation);
	}
	// ...

	// increase number of systematices
	fNSystematics++;

	// add parameter index to container
	fSystematicParIndexContainer.push_back(GetNParameters());

	// add parameter
	AddParameter(name, min, max);

	// no error
	return 1;

}

// ---------------------------------------------------------
int BCMultiTemplateFitter::SetSystematicVariation(const char* channelname, const char* processname,  const char* systematicname, TH1D hist_up, TH1D hist_down)
{
	// get channel index
	int channelindex = GetChannelIndex(channelname);

	// check if channel exists
	if (channelindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Channel does not exist.");
			return -1;
	}

	// get process index
	int processindex = GetProcessIndex(processname);

	// check if process exists
	if (processindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Process does not exist.");
			return -1;
	}

	// get systematic index
	int systematicindex = GetSystematicIndex(systematicname);

	// check if systematic exists
	if (systematicindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Systematic does not exist.");
			return -1;
	}

	// get channel
	BCChannel* channel = GetChannel(channelindex);
	
	// get systematic variation
	BCSystematicVariation* variation = channel->GetSystematicVariation(systematicindex);

	// set histogram
	variation->SetHistograms(processindex, new TH1D(hist_up), new TH1D(hist_down));
	
	// no error
	return 1;
}

// ---------------------------------------------------------
int BCMultiTemplateFitter::SetSystematicVariation(const char* channelname, const char* processname,  const char* systematicname, TH1D hist, TH1D hist_up, TH1D hist_down)
{
	// get number of bins
	int nbins = hist.GetNbinsX(); 

	TH1D* hist_up_rel = new TH1D(hist);
	TH1D* hist_down_rel = new TH1D(hist);

	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin) {
		hist_up_rel->SetBinContent(ibin, (hist_up.GetBinContent(ibin) - hist.GetBinContent(ibin)) / hist.GetBinContent(ibin));
		hist_down_rel->SetBinContent(ibin, (hist.GetBinContent(ibin) - hist_down.GetBinContent(ibin)) / hist.GetBinContent(ibin));
	}
	
	// set the systematic variation
	return SetSystematicVariation(channelname, processname, systematicname, *hist_up_rel, *hist_down_rel);
}

// ---------------------------------------------------------
int BCMultiTemplateFitter::PrintSummary(const char* filename)
{
	// open file
	std::ofstream ofi(filename);
	ofi.precision(3);
	
	// check if file is open
	if(!ofi.is_open()) {
		BCLog::OutWarning(Form("BCMultitemplateFitter::PrintSummary() : Could not open file %s", filename));
		return 0;
	}
	
	ofi
		<< " Multi template fitter summary " << std::endl
		<< " ----------------------------- " << std::endl
		<< std::endl
		<< " Number of channels      : " << fNChannels << std::endl
		<< " Number of processes     : " << fNProcesses << std::endl
		<< " Number of systematics   : " << fNSystematics << std::endl
		<< std::endl;
		
	ofi
		<< " Channels :" << std::endl;
		for (int i = 0; i < GetNChannels(); ++i) {
			ofi 
				<< " " << i 
				<< " : \"" << GetChannel(i)->GetName().c_str() << "\""
				<< std::endl;
		}
		ofi 
			<< std::endl;

	ofi
		<< " Processes :" << std::endl;
		for (int i = 0; i < GetNProcesses(); ++i) {
			ofi 
				<< " " << i 
				<< " : \"" << GetProcess(i)->GetName().c_str()  << "\""
				<< " (par index " << GetParIndexProcess(i) << ")" 
				<< std::endl;
		}
		ofi 
			<< std::endl;

	ofi
		<< " Systematics :" << std::endl;
		for (int i = 0; i < GetNSystematics(); ++i) {
			ofi 
				<< " " << i 
				<< " : \"" << GetSystematic(i)->GetName().c_str()  << "\""
				<< " (par index " << GetParIndexSystematic(i) << ")" 
				<< std::endl;
		}
		ofi 
			<< std::endl;
		if (GetNSystematics() == 0)
			ofi 
				<< " - none - " << std::endl;

		ofi 
			<< " Goodness-of-fit: " << std::endl;
		for (int i = 0; i < GetNChannels(); ++i) {
			ofi 
				<< " i : \"" << GetChannel(i)->GetName().c_str() << "\" : chi2 = " 
				<< CalculateChi2( i, GetBestFitParameters() )
				<< std::endl;
		}
		ofi 
			<< std::endl;
		

	// close file
	ofi.close();

	// no error
	return 1;
}

// ---------------------------------------------------------
double BCMultiTemplateFitter::Expectation(int channelindex, int binindex, std::vector<double>& parameters)
{
	double expectation = 0.;

	// loop over all processes
	for (int i = 0; i < fNProcesses; ++i) {
		// get efficiency
		double efficiency = Efficiency(channelindex, i, binindex, parameters); 

		// get probability
		double probability = Probability(channelindex, i, binindex, parameters);

		// get parameter index
		int parindex = fProcessParIndexContainer[i]; 

		// add to expectation
		expectation += parameters[parindex]
			* efficiency
			* probability; 
	}

	// check if expectation is positive
	if (expectation < 0)
		expectation = 0.;

	return expectation;
}

// ---------------------------------------------------------
double BCMultiTemplateFitter::Efficiency(int channelindex, int processindex, int binindex, std::vector<double>& parameters)
{
	// get channel
	BCChannel* channel = fChannelContainer[channelindex];

	double efficiency = channel->GetTemplate(processindex)->GetEfficiency();

	double defficiency = 1.;

	// loop over all systematics
	for (int i = 0; i < fNSystematics; ++i) {
		if (!(fSystematicContainer[i]->GetFlagSystematicActive()))
			continue;
		
		// get parameter index
		int parindex = fSystematicParIndexContainer[i]; 
		
		// get histogram
		TH1D* hist = 0; 
		
		if (parameters[parindex] > 0)
			hist = channel->GetSystematicVariation(i)->GetHistogramUp(processindex);
		else
			hist = channel->GetSystematicVariation(i)->GetHistogramDown(processindex);
		
		// check if histogram exists
		if (!hist) 
			continue;

		// multiply efficiency
		defficiency += parameters[parindex] * hist->GetBinContent(binindex);
	}

	// calculate efficiency
	efficiency *= defficiency;

	// make sure efficiency is between 0 and 1
	if (fFlagEfficiencyConstraint) {
		if (efficiency < 0.)
			efficiency = 0.;
		if (efficiency > 1.)
			efficiency = 1.;
	}

	return efficiency;
}

// ---------------------------------------------------------
double BCMultiTemplateFitter::Probability(int channelindex, int processindex, int binindex, std::vector<double>& parameters)
{
	// get histogram
	TH1D* hist = fChannelContainer[channelindex]->GetTemplate(processindex)->GetHistogram();

	// this needs to be fast
	if (!hist)
		return 0.;
	
	return hist->GetBinContent(binindex);
}

// ---------------------------------------------------------
int BCMultiTemplateFitter::PrintStack(const char* channelname, std::vector<double> parameters, const char* filename, const char* options) 
{
	int index = GetChannelIndex(channelname);

	return PrintStack(index, parameters, filename, options);
}

// ---------------------------------------------------------
int BCMultiTemplateFitter::PrintStack(int channelindex, std::vector<double> parameters, const char* filename, const char* options)
{
	// check if parameters are filled
	if (!parameters.size())
		return -1;

	// check options
	bool flag_logx   = false; // plot x-axis in log-scale
	bool flag_logy   = false; // plot y-axis in log-scale

	if (std::string(options).find("logx") < std::string(options).size())
		flag_logx = true;

	if (std::string(options).find("logy") < std::string(options).size())
		flag_logy = true;

	// get channel
	BCChannel* channel = GetChannel(channelindex);

	// create canvas
	TCanvas* c1 = new TCanvas(); 
	c1->cd();

	// set log or linear scale
	if (flag_logx)
		c1->SetLogx();
	
	if (flag_logy)
		c1->SetLogy();
	
	// create stack
	THStack* stack = new THStack("", ""); 

	// create a container of temporary histograms
	std::vector<TH1D*> histcontainer;

	// get number of templates
	unsigned int ntemplates = GetNProcesses();

	// loop over all templates
	for (unsigned int i = 0; i < ntemplates; ++i) {

		// get histogram
		TH1D* temphist = channel->GetTemplate(i)->GetHistogram(); 

		// create new histogram
		TH1D* hist;

		if (temphist)
			hist = new TH1D( *(temphist) );
		else
			continue;

		// get number of bins
		int nbins = hist->GetNbinsX(); 

		// scale histogram
		for (int ibin = 1; ibin <= nbins; ++ibin) {

			// get efficiency
			double efficiency = Efficiency(channelindex, i, ibin, parameters); 

			// get probability
			double probability = Probability(channelindex, i, ibin, parameters);

			// get parameter index
			int parindex = GetParIndexProcess(i); 

			// add to expectation
			double expectation = parameters.at(parindex)
				* efficiency
				* probability;
			
			// set bin content
			hist->SetBinContent(ibin, expectation);
		}

		// add histogram to container (for memory management)
		histcontainer.push_back(hist);

		// add histogram to stack
		stack->Add(hist);
	}

	//draw data
	channel->GetData()->GetHistogram()->Draw("PE");
	
	// set range user
	channel->GetData()->GetHistogram()->GetYaxis()->SetRangeUser(0., channel->GetData()->GetHistogram()->GetMaximum() + 2. * sqrt(channel->GetData()->GetHistogram()->GetMaximum()));

	// draw stack
	stack->Draw("SAMEHIST");

	//draw data again
	channel->GetData()->GetHistogram()->Draw("SAMEPE");

	// redraw the axes
	gPad->RedrawAxis();

	// print
	c1->Print(filename);

	// free memory
	for (unsigned int i = 0; i < histcontainer.size(); ++i) {
		TH1D* hist = histcontainer.at(i);
		delete hist;
	}
	delete stack;
	delete c1; 

	// no error
	return 1;
}

// ---------------------------------------------------------
double BCMultiTemplateFitter::CalculateChi2(int channelindex, std::vector<double> parameters)
{
	if (parameters.size() == 0)
		return -1;

	double chi2 = 0;
	
	// get channel
	BCChannel* channel = GetChannel(channelindex);
	
	// get data
	BCTemplate* data = channel->GetData();
	
	// get histogram
	TH1D* hist = data->GetHistogram(); 
	
	// check if histogram exists
	if (hist) {
		// get number of bins in data
		int nbins = hist->GetNbinsX();
		
		// loop over all bins
		for (int ibin = 1; ibin <= nbins; ++ibin) {
			// get expectation value
			double expectation = Expectation(channelindex, ibin, parameters);
			
			// get observation
			double observation = hist->GetBinContent(ibin);
			
			// add Poisson term
			chi2 += (expectation - observation) * (expectation - observation) / expectation;
		}
	}
	
	// return chi2;
	return chi2;
}

// ---------------------------------------------------------
double BCMultiTemplateFitter::CalculateChi2(std::vector<double> parameters)
{
	if (parameters.size() == 0)
		return -1;

	double chi2 = 0;

	// get number of channels
	int nchannels = GetNChannels();

	// loop over all channels
	for (int i = 0; i < nchannels; ++i) {
		chi2 += CalculateChi2(i, parameters);
	}
	
	// return chi2
	return chi2;
}

// ---------------------------------------------------------
double BCMultiTemplateFitter::CalculateCash(int channelindex, std::vector<double> parameters)
{
	if (parameters.size() == 0)
		return -1;

	double cash = 0;
	
	// get channel
	BCChannel* channel = GetChannel(channelindex);
	
	// get data
	BCTemplate* data = channel->GetData();
	
	// get histogram
	TH1D* hist = data->GetHistogram(); 
	
	// check if histogram exists
	if (hist) {
		// get number of bins in data
		int nbins = hist->GetNbinsX();
		
		// loop over all bins
		for (int ibin = 1; ibin <= nbins; ++ibin) {
			// get expectation value
			double expectation = Expectation(channelindex, ibin, parameters);
			
			// get observation
			double observation = hist->GetBinContent(ibin);
			
			// calculate Cash statistic
			cash += 2. * (expectation - observation); 

			// check negative log
			if (observation > 0)
				cash += 2. * observation * log (observation/expectation);
		}
	}
	
	// return cash;
	return cash;
	
}

// ---------------------------------------------------------
double BCMultiTemplateFitter::CalculateCash(std::vector<double> parameters)
{
	if (parameters.size() == 0)
		return -1;

	double cash = 0;

	// get number of channels
	int nchannels = GetNChannels();

	// loop over all channels
	for (int i = 0; i < nchannels; ++i) {
		cash += CalculateCash(i, parameters);
	}
	
	// return cash
	return cash;
}

// ---------------------------------------------------------
double BCMultiTemplateFitter::LogLikelihood(std::vector<double> parameters)
{
	double logprob = 0.;
	
	// loop over all channels
	for (int ichannel = 0; ichannel < fNChannels; ++ichannel) {
		
		// get channel
		BCChannel* channel = fChannelContainer[ichannel];
		
		// check if channel is active
		if (!(channel->GetFlagChannelActive()))
			continue;

		// get data
		BCTemplate* data = channel->GetData();

		// get histogram
		TH1D* hist = data->GetHistogram(); 

		// check if histogram exists
		if (!hist) 
			continue;

		// get number of bins in data
		int nbins = data->GetNBins();

		// loop over all bins
		for (int ibin = 1; ibin <= nbins; ++ibin) {
			// get expectation value
			double expectation = Expectation(ichannel, ibin, parameters);
			
			// get observation
			double observation = hist->GetBinContent(ibin);

			// add Poisson term
			logprob += BCMath::LogPoisson(observation, expectation);
		}
	}

	return logprob;
}

// ---------------------------------------------------------
