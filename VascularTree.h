/*=========================================================================
 
	Program: VascuSynth
	Module: $RCSfile: VascularTree.h,v $
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

#ifndef _vasculartree_h
#define _vasculartree_h

#include "NodeTable.h"
#include "OxygenationMap.h"


/** \class VascularTree
 * \brief Class that iteratively builds the vascular tree 
 *
 * Iteratively builds the vascular structure by selecting a new candidate
 * node and connecting it to an existing edge, creating a bifurcation location.
 * Updates the node table with the new candidate node.
 */
class VascularTree{
public:
	NodeTable nt;
	OxygenationMap *oxMap;

	double Pperf;
	double Pterm;
	double Qperf;
	double Qterm;
	
	double rho;
	double gamma;
	double lambda;
	double mu;
	
	double minDistance;
	int numNodes;
	
	double* perf;
	
	int closestNeighbours;
	
	double mapVoxelWidth; //in cm
	
	VascularTree(OxygenationMap * oxMap, double* perf, double Pperf, double Pterm, double Qperf, double rho, double gamma, double lambda, double mu, double minDistance, int numNodes, double voxelWidth, int closestNeighbours);
		
	//calculate the distance between to nodes in the node table
	double distance(int from, int to);
	//calculate the reduced resistence of segment @ id
	void calculateReducedResistence(int id);	
	//calculates the ratio of radii of the segment @ id
	void calculateRatios(int id);	
	//update the tree @ the bifurication point @ id
	void updateAtBifurication(int id, int newChild);
	
	//calculate the radii throughout the tree
	void calculateRadius();
	void calculateRadius(int id);

	//calculate the fitness function
	double calculateFitness();

	//When used by local optimization, ignored is the segment to connect to
	//otherwise it should be -1;
	bool validateCandidate(double* x0, int ignored);
	
	//connect point to segment through bifPoint
	void connectPoint(double* point, int segment, double* bifPoint);	
	//update the flow throughout the tree
	void incrementFlow(int parent, double Qterm);	
	//the distance between a point and a segment
	double pointSegmentDistance(double* x0, int segment);
	//optimize the location of a bifurication point for terminal node point and segment 'segment'
	double* localOptimization(double * point, int segment, int steps);
	//determine if point is in the volume
	bool inVolume(double * point);
	//try to attach point to the tree
	bool connectCandidate(double * point, int steps);
	//build the tree
	void buildTree();
};

#endif