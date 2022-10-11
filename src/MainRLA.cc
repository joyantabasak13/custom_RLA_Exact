/*
 * MainRLA.cpp
 *
 *  Created on: Aug 5, 2014
 *      Author: abdullah-al-mamun
 *
 *  Here we use constant threshold value
 */

#include <iostream>
#include <cstdlib>
#include "RLA_CL.h"

using namespace std;

// driver function
int main(int argc, char** argv)
{
	if (argc < 3) {
		cerr << "rlacl config_file_name.xml number_of_records" << endl;
		return EXIT_FAILURE;
	}
		
	doRLA(argv[1], atoi(argv[2])); 
	return 0;
}



