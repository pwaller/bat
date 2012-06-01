#ifndef __BCCHANNEL__H
#define __BCCHANNEL__H

#include <string>
#include <vector>

class BCTemplate;
class BCSystematicVariation; 

// ---------------------------------------------------------
class BCChannel
{

	public:

		// Constructors and destructor
		BCChannel(const char* name);
		~BCChannel();

		// setters
		
		// set name
		void SetName(const char* name)
		{ fName = name; }; 

		// set data
		void SetData(BCTemplate* bctemplate)
		{ fData = bctemplate; }; 

		// set a flag for using this channel or not
		void SetFlagChannelActive(bool flag)
		{ fFlagChannelActive = flag; }; 

		// getters
		std::string GetName()
			{ return fName; }; 

		// return a template
		BCTemplate* GetTemplate(int index)
		{ return fTemplateContainer.at(index); }; 

		// return a systematicvariation
		BCSystematicVariation* GetSystematicVariation(int index)
		{ return fSystematicVariationContainer.at(index); }; 

		// return the data
		BCTemplate* GetData()
		{ return fData; }; 

		// return flag
		bool GetFlagChannelActive()
		{ return fFlagChannelActive; }; 

		// misc

		// add a template
		void AddTemplate(BCTemplate* bctemplate)
		{ fTemplateContainer.push_back(bctemplate); }; 

		// add a systematic variation
		void AddSystematicVariation(BCSystematicVariation* variation)
		{ fSystematicVariationContainer.push_back(variation); }; 

		// print templates
		void PrintTemplates(const char* filename);

		// print template with systematics
		void PrintTemplate(int index, const char* filename);

 private:

		// name of the channel
		std::string fName;

		// the data set
		BCTemplate* fData;

		// a container of templates
		std::vector<BCTemplate*> fTemplateContainer;

		// a container of systematics
		std::vector<BCSystematicVariation*> fSystematicVariationContainer;

		// flag: channel is used (true) or not (false) in fit
		bool fFlagChannelActive; 

};
// ---------------------------------------------------------

#endif

