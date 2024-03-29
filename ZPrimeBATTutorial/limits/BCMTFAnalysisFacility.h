#ifndef __BCMTFANALYSISFACILITY__H
#define __BCMTFANALYSISFACILITY__H

#include <string>
#include <vector>

class BCMultiTemplateFitter;
class TTree; 
class TH1D; 
class TRandom3; 

// ---------------------------------------------------------
class BCMTFAnalysisFacility
{

	public:

		// Constructors and destructor
		BCMTFAnalysisFacility(BCMultiTemplateFitter* mtf);
		~BCMTFAnalysisFacility();

		// setters
		
		// set MTF
		void SetMTF(BCMultiTemplateFitter* mtf)
		{ fMTF = mtf; }; 

		// set flag for using MCMC (true) or not (false)
		void SetFlagMCMC(bool flag)
		{ fFlagMCMC = flag; }; 

		// getters
		
		// return MTF
		BCMultiTemplateFitter* GetMTF()
		{ return fMTF; }; 

		// misc

		// perform the full set of single channel analyses and the
		// combination
		int PerformSingleChannelAnalyses(const char* dirname, const char* options = "");

		// perform the analysis using one systematic at a time, without
		// systematic and with all systematics
		// flag_nuisance: use nuisance parameters (true) or delta method (false)
		int PerformSingleSystematicAnalyses(const char* dirname, const char* options = "");

		// perform calibration curve
		int PerformCalibrationAnalysis(const char* dirname, std::vector<double> default_parameters, int index, std::vector<double> parametervalues, int nensembles = 1000);

		// perform full ensemble test
		int PerformEnsembleTest(std::vector<double> parameters);

		// build a single ensemble based on a single set of parameters
		std::vector<TH1D> BuildEnsemble(std::vector<double> parameters);

		// build ensembles based on a single set of parameters
		TTree* BuildEnsembles(std::vector<double> parameters, int nensembles);

		// build ensembles based on a varying sets of parameters, e.g., using the prior or posterior
		TTree* BuildEnsembles(TTree* tree, int nensembles);
		
		// perform ensemble test based on one set of parameters
		TTree* PerformEnsembleTest(std::vector<double> parameters, int nensembles); 

		// perform ensemble test based on varying sets of parameters
		TTree* PerformEnsembleTest(TTree* tree, int nensembles); 

		// transform a matrix to a set of histograms
		std::vector<TH1D> MatrixToHistograms(std::vector< std::vector<double> > matrix);

 private:

		// the multi template fitter
		BCMultiTemplateFitter* fMTF;

		// a random number generator
		TRandom3* fRandom; 

		// flag: use MCMC for analysis
		bool fFlagMCMC;
};
// ---------------------------------------------------------

#endif

