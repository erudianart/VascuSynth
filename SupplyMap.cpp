/*=========================================================================
 
	Program: VascuSynth
	Module: $RCSfile: SupplyMap.cpp,v $
	Language: C++
	Date: $Date: 2011/02/08 10:43:00 $
	Version: $Revision: 1.0 $

Copyright (c) 2011 Medical Imaging Analysis Lab, Simon Fraser University, 
British Columbia, Canada.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

* The name of the Insight Consortium, nor the names of any consortium members,
nor of any contributors, may be used to endorse or promote products derived
from this software without specific prior written permission.

* Modified source versions must be plainly marked as such, and must not be
misrepresented as being the original software.
 
* Free for non-commercial use only.  For commercial use, explicit approval 
must be requested by contacting the Authors.

* If you use the code in your work, you must acknowledge it

* Modifications of the source code must also be released as open source 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <string>
#include <cmath>

#include "SupplyMap.h"

using namespace std;

// A parametric macro to allow the use of dynamic multi-dimensional arrays
#define arr(arr,w,x,y,z,dim) *(arr + (z + dim[3]*(y + dim[2]*(x + dim[1]*(w)))))


/**
 * calculate the reduction to target if supplied becomes a terminal node
 */
double SupplyMap::reduction(int supplied[], int target[]){
	if(supplied[0] == target[0] && supplied[1] == target[1] && supplied[2] == target[2])
		return 0;

	double dist = sqrt(pow((double)(supplied[0] - target[0]), 2) + pow((double)(supplied[1] - target[1]),2 )+ pow((double)(supplied[2] - target[2]), 2));
	double acc = 0;

	if(arr(map_d, target[0], target[1], target[2], dim[3]-1, dim) <= dist)
		return 1;
	
	for(int i = 0; i < dim[3]; i++){
		acc += pow(dist, i)*arr(map_d, target[0], target[1], target[2], i, dim);
	}
	
	acc = 1.0 - (1.0/acc);
	acc = (acc < 0 ? 0 : acc);
	
	
	return acc;
}

/**
 * load a file containing the description of an supply map
 */
void SupplyMap::loadMap(string filename){
	
	char * tempFilename;
	
	tempFilename = new char[filename.size()+1];
	strcpy(tempFilename, filename.c_str());
	
	int lastIndex = strlen(tempFilename)-1;
	
	while ((tempFilename[lastIndex] == '\r') || (tempFilename[lastIndex] == '\n') || (tempFilename[lastIndex] == ' ')) {
		tempFilename[lastIndex] = '\0';
		lastIndex--;
	}
	
	ifstream mapFile;
	mapFile.open(tempFilename, ios::in);
	string line;
		
	if (mapFile.is_open()) {
        
		// defining i,j,k,l outside of the for stament 
        // so this is compatable with both VC++ 6.0 and later versions
		int i;
        int j;
        int k;
        int l; 

		getline(mapFile, line);
		char * tok = new char[line.size()+1];
		strcpy(tok, line.c_str());
		
		for(i = 0; i < 4; i++){
			char * temp = strtok((i == 0? tok : NULL), " ");
			if(temp == NULL) {
				break;
			}
			dim[i] = atoi(temp);
		}
		
		//get the dimensions
		map_d = new double[dim[0]*dim[1]*dim[2]*dim[3]];

		for(i = 0; i < dim[0]; i++) {
			for(j = 0; j < dim[1]; j++) {
				for(k = 0; k < dim[2]; k++) {
					arr(map_d,i,j,k,dim[3]-1,dim) = 0;
				}
			}
		}

		//get the regions
		
		while(!mapFile.eof()){
			
			getline(mapFile, line);
			tok = new char[line.size()+1];
			strcpy(tok, line.c_str());
			
			int region[6];
			for(i = 0; i < 6; i++){
				char * temp = strtok((i == 0? tok : NULL), " ");
				if(temp == NULL) {
					break;
				}
				region[i] = atoi(temp);
			}
				
			if(mapFile.eof()) {
				break;
			}

			getline(mapFile, line);
			tok = new char[line.size()+1];
			strcpy(tok, line.c_str());

			double *values = new double[dim[3]];
			for(i = 0; i < dim[3]; i++) {
				values[i] = atof(strtok((i == 0? tok : NULL), " "));
			}

			for(i = region[0]; i < region[3]; i++) {
				for(j = region[1]; j< region[4]; j++) {
					for(k = region[2]; k < region[5]; k++) {
						for(l = 0; l < dim[3]; l++) {
							arr(map_d, i, j, k, l, dim) = values[l];
						}
					}
				}
			}
			 
		}
		 
	
		mapFile.close();
	
    } else {
        
        throw "Could not read the SupplyMap file";
        
    }
	 
}
