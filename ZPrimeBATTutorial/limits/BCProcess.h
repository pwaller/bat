#ifndef __BCPROCESS__H
#define __BCPROCESS__H

#include <string>

// ---------------------------------------------------------
class BCProcess
{
	public:

		// Constructors and destructor
		BCProcess(const char* name);
		~BCProcess();

		// setters
		
		// set name
		void SetName(const char* name)
		{ fName = name; }; 

		// getters
		std::string GetName()
			{ return fName; }; 

 private:

		// name of the channel
		std::string fName;

};
// ---------------------------------------------------------

#endif

