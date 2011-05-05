/*=========================================================================
 
	Program: VascuSynth
	Module: $RCSfile: TreeDrawer.h,v $
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

#ifndef _treedrawer_h
#define _treedrawer_h

#include "VascularTree.h"

#include <vector>

using namespace std;

/** \class TreeDrawer
 * \brief Constructs a volume of the vasculature that can be printed to an image
 *
 * Converts the vascular structure from the node table into a 3d matrix that can then
 * be saved as a 3d image.  Also adds noise to the matrix before it is saved as an image.
 */
class TreeDrawer {
public:
	double imageVoxelWidth;
	double mapVoxelWidth;
	VascularTree * vt;
	double *c1, *c2;
	unsigned char *image;
	int dim[3];
	
	TreeDrawer(VascularTree * vt, double width, double *corner1, double *corner2);

	/**
	 * maps a voxel to a point
	 * 8 subsections in a voxel
	 * if a voxel has corners (0,0,0) and (1,1,1)
	 * 	section 0 has corners (0.0, 0.0, 0.0) and (0.5, 0.5, 0.5) => center = voxel center + (-0.25, -0.25, -0.25)
	 * 	section 1 has corners (0.5, 0.0, 0.0) and (1.0, 0.5, 0.5) => center = voxel center + (+0.25, -0.25, -0.25)
	 * 	section 2 has corners (0.0, 0.5, 0.0) and (0.5, 1.0, 0.5) => center = voxel center + (-0.25, +0.25, -0.25)
	 * 	section 3 has corners (0.5, 0.5, 0.0) and (1.0, 1.0, 0.5) => center = voxel center + (+0.25, +0.25, -0.25)
	 * 	section 4 has corners (0.0, 0.0, 0.5) and (0.5, 0.5, 1.0) => center = voxel center + (-0.25, -0.25, +0.25)
	 * 	section 5 has corners (0.5, 0.0, 0.5) and (1.0, 0.5, 1.0) => center = voxel center + (+0.25, -0.25, +0.25)
	 * 	section 6 has corners (0.0, 0.5, 0.5) and (0.5, 1.0, 1.0) => center = voxel center + (-0.25, +0.25, +0.25)
	 * 	section 7 has corners (0.5, 0.5, 0.5) and (1.0, 1.0, 1.0) => center = voxel center + (+0.25, +0.25, +0.25)
	 */
	void voxelToPoint(int* voxel, int subSection, double*ret);	
	//checks if a point is inside a given tube
	bool inTube(double * point, double * p1, double * p2, double radius);	
	//checks if a point is in the end of a tube
	bool inEnd(double *point, double *p1, double radius);	
	//checks if a point is in the same region as a tube (used to improve performance)
	bool checkTube(double *point, double *to, double *from, double radius);	
	//determines the value at a specific voxel
	unsigned char valueAtVoxel(int* voxel);	

	//add some uniform noise to a voxel
	unsigned char addNoise_Uniform(unsigned char c, double lb, double ub);	
	//add some salt&pepper noise to a voxel
	unsigned char addNoise_saltPepper(unsigned char c, unsigned char valSalt, double probSalt, unsigned char valPepper, double probPepper);	
	//add some gaussian noise to a voxel
	unsigned char addNoise_gaussian(unsigned char c, double median, double sigma);

	//draw the tree to a matrix
	void drawImage();
	//add uniform noise to the matrix
	void addNoise_Uniform(double lb, double ub);
	//add salt&pepper noise to the matirx
	void addNoise_saltPepper(unsigned char valSalt, double probSalt, unsigned char valPepper, double probPepper);
	//add gaussian noise to the matrix
	void addNoise_gaussian(double median, double sigma);
	//add a single shadow @ (x,y,z)
	void addShadow(double x, double y, double z, double centerRatio, double r);	
	//add 'numShadows' shadows to the matix
	void addShadows(int numShadows);
	//get the value of the matrix at (x,y,z)
	unsigned char imageAt(int x, int y, int z);
	//copy the matrix
	TreeDrawer* copy();
};

#endif