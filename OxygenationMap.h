/*=========================================================================
 
	Program: VascuSynth
	Module: $RCSfile: OxygenationMap.h,v $
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

#ifndef _oxygenationmap_h
#define _oxygenationmap_h

#include <string>
#include <vector>

#include "SupplyMap.h"
#include "MersenneTwister.h"

using namespace std;

/** \class OxygenationMap
 * \brief Class for the Oxygenation Demand Map
 *
 * Contains a map of information about the demand for oxygen in the volume
 * as well as methods to update the map based on the supply map
 */
class OxygenationMap {
public:
	SupplyMap *supply;

	double * map_d;
	double * effectiveMap_d;
	MTRand rand;

	int dim[3];

	OxygenationMap(SupplyMap *sMap, int randSeed);

	//load a file containing the description of an oxygenation map
	void loadMap(string filename);
	//calculate the sum of the current effective map
	double sum();
	//selecte a candidate terminal node
	void candidate(double sum, int *cand);
	//update the effective map based on using cand as a new terminal node
	void applyCandidate(int cand[]);
	//determine if source is visible from target with respect to the oxygenation map (uses original map - not effective map)
	bool visible(double source[], double target[]);
};

#endif
