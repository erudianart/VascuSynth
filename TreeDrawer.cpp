/*=========================================================================
 
	Program: VascuSynth
	Module: $RCSfile: TreeDrawer.cpp,v $
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

#include "TreeDrawer.h"
#include "VascularTree.h"
#include "MersenneTwister.h"

#include <vector>
#include <cmath>

using namespace std;

// A parametric macro to allow the use of dynamic multi-dimensional arrays
#define arr(arr,x,y,z,dim) *(arr + (z + dim[2]*(y + dim[1]*(x))))

/**
 * Constructor
 */
TreeDrawer::TreeDrawer(VascularTree * _vt, double width, double * corner1, double * corner2){
	vt = _vt;

	imageVoxelWidth = width;
	mapVoxelWidth = vt->mapVoxelWidth;
	
    c1 = new double[3];
	c1[0] = corner1[0]; c1[1] = corner1[1]; c1[2]= corner1[2];
	
    c2 = new double[3];
	c2[0] = corner2[0]; c2[1] = corner2[1]; c2[2]= corner2[2];
	
	dim[2] = (int) ((fabs(corner2[2] - corner1[2]) * mapVoxelWidth / width) + 1);
	dim[1] = (int) ((fabs(corner2[1] - corner1[1]) * mapVoxelWidth / width) + 1);
	dim[0] = (int) ((fabs(corner2[0] - corner1[0]) * mapVoxelWidth / width) + 1);
	
	image = new unsigned char[dim[2]*dim[1]*dim[0]];
	
}

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
void TreeDrawer::voxelToPoint(int* voxel, int subSection, double *ret){
	ret[0] = ((double)voxel[0])*imageVoxelWidth / mapVoxelWidth;
	ret[1] = ((double)voxel[1])*imageVoxelWidth / mapVoxelWidth;
	ret[2] = ((double)voxel[2])*imageVoxelWidth / mapVoxelWidth;
	
	if(subSection == 0 || subSection == 2 || subSection == 4 || subSection == 6)
		ret[0] -= 0.25*imageVoxelWidth/mapVoxelWidth;
	if(subSection == 1 || subSection == 3 || subSection == 5 || subSection == 7)
		ret[0] += 0.25*imageVoxelWidth/mapVoxelWidth;
	if(subSection == 0 || subSection == 1 || subSection == 4 || subSection == 5)
		ret[1] -= 0.25*imageVoxelWidth/mapVoxelWidth;
	if(subSection == 2 || subSection == 3 || subSection == 6 || subSection == 7)
		ret[1] += 0.25*imageVoxelWidth/mapVoxelWidth;
	if(subSection == 0 || subSection == 1 || subSection == 2 || subSection == 3)
		ret[2] -= 0.25*imageVoxelWidth/mapVoxelWidth;
	if(subSection == 4 || subSection == 5 || subSection == 6 || subSection == 7)
		ret[2] += 0.25*imageVoxelWidth/mapVoxelWidth;
}

/**
 * checks if a point is inside a given tube
 */
bool TreeDrawer::inTube(double * point, double * p1, double * p2, double radius){
	if(inEnd(point, p1, radius))
		return true;
	
	double t = -((p2[0] - p1[0])*(p1[0] - point[0]) + 
				 (p2[1] - p1[1])*(p1[1] - point[1]) + 
				 (p2[2] - p1[2])*(p1[2] - point[2])) / 
				((p2[0] - p1[0])*(p2[0] - p1[0]) + 
				 (p2[1] - p1[1])*(p2[1] - p1[1]) +
				 (p2[2] - p1[2])*(p2[2] - p1[2]));
	
	if (t < 0 || t > 1) {
		return false;
	} else {
		double distance = pow(pow((((p2[0] - p1[0])*t + p1[0]) - point[0]), 2) + 
						    	pow((((p2[1] - p1[1])*t + p1[1]) - point[1]), 2) + 
								pow((((p2[2] - p1[2])*t + p1[2]) - point[2]), 2), 0.5);
		distance *= mapVoxelWidth;
		
        // this should be return (distance <= radius);
        if (distance <= radius) {
			return true;
		} else { 
			return false;
        }
	}
}

/**
 * checks if a point is in the end of a tube
 */
bool TreeDrawer::inEnd(double * point, double * p1, double radius){
    
    //same as above, no need for if, just return
	if(pow(pow(p1[0] - point[0], 2) +
			pow(p1[1] - point[1], 2) +
			pow(p1[2] - point[2], 2), 0.5)*mapVoxelWidth < radius)
		return true;

	return false;
}

/**
 * checks if a point is in the same region as a tube (used to improve performance)
 */
bool TreeDrawer::checkTube(double *point, double *to, double *from, double radius){
	double c1[3];
	c1[0] = (to[0] < from[0]? to[0] : from[0]) - radius/mapVoxelWidth;
	c1[1] = (to[1] < from[1]? to[1] : from[1]) - radius/mapVoxelWidth;
	c1[2] = (to[2] < from[2]? to[2] : from[2]) - radius/mapVoxelWidth;
	
	double c2[3];
	c2[0] = (to[0] > from[0]? to[0] : from[0]) + radius/mapVoxelWidth;
	c2[1] = (to[1] > from[1]? to[1] : from[1]) + radius/mapVoxelWidth;
	c2[2] = (to[2] > from[2]? to[2] : from[2]) + radius/mapVoxelWidth;
	
    //can do the same thing, return blah
	if(point[0]  >= c1[0] && point[0] <= c2[0] &&
			point[1]  >= c1[1] && point[1] <= c2[1] &&
			point[2]  >= c1[2] && point[2] <= c2[2])
		return true;
	return false;
	
}

/**
 * determines the intensity value at a specific voxel
 * in the range of 255 to 0
 */
unsigned char TreeDrawer::valueAtVoxel(int* voxel){
	
    int value = 0;
	int segments[] = {0,1,2,3,4,5,6,7};
	double *point = new double[3];
    int size = (int) vt->nt.nodes.size();
    
	for(int j = 1; j < size; j++){
		double * p1 = vt->nt.getPos(j);
		int parent = vt->nt.getParent(j);
		double * p2 = vt->nt.getPos(parent);
		double radius = vt->nt.getRadius(j);
		
		voxelToPoint(voxel, -1, point);
		if(checkTube(point, p1, p2, radius)){
			//for each of the 8 offsets of a point - check if the point is in a tube
			for(int i = 0; i < 7; i++){
				if(segments[i] == -1)
					continue;
				
				voxelToPoint(voxel, segments[i], point);
				if(inTube(point, p1, p2, radius)){
					value += 32;
					segments[i] = -1;
				}
			}
		}
		if(value == 256)
			break;
	}
	
	value--;
	
	if (value >= 256) {
		value = 255;
    }
    
	if (value < 0) {
		value = 0;
    }
    
	delete point;
	return (unsigned char)value;
}

/**
 * add some uniform noise to a voxel
 */
unsigned char TreeDrawer::addNoise_Uniform(unsigned char c, double lb, double ub){
        
	MTRand rand;
	int x = c + (int) (rand.rand()*(ub-lb) + lb);
	if (x > 255) {
		return (unsigned char)255;
	} else if (x < 0) {
		return (unsigned char)0;
	} else {
		return (unsigned char) x;
    }
}

/**
 * add some salt&pepper noise to a voxel
 */
unsigned char TreeDrawer::addNoise_saltPepper(unsigned char c, unsigned char valSalt, double probSalt, unsigned char valPepper, double probPepper){
	
    MTRand rand;
	double r = rand.rand();
	if(r <= probSalt) {
		return valSalt;
	} else if(r <= probPepper+ probSalt) {
		return valPepper;
	} else {
		return c;
    }
}

/**
 * add some gaussian noise to a voxel
 */
unsigned char TreeDrawer::addNoise_gaussian(unsigned char c, double median, double sigma){
    
	MTRand rand;
	int x = c + (int)(rand.randNorm(median, sigma));
	if(x > 255) {
		return (unsigned char)255;
	} else if(x < 0) {
		return (unsigned char)0;
	} else {
		return (unsigned char) x;
    }
}

/** 
 * draw the tree to a matrix
 */
void TreeDrawer::drawImage(){
	int voxel[3];
	for(int i = 0; i < dim[0]; i++){
		for(int j = 0; j < dim[1]; j++) {
			for(int k = 0; k < dim[2]; k++){
				voxel[0] = i; voxel[1] = j; voxel[2] = k;
				arr(image, i, j, k, dim) = valueAtVoxel(voxel);
			}
        }
	}
}

/**
 * add uniform noise to the matrix
 */
void TreeDrawer::addNoise_Uniform(double lb, double ub){
    
	for(int i = 0; i < dim[0]; i++) {
		for(int j = 0; j < dim[1]; j++) {
			for(int k = 0; k < dim[2]; k++) {
				arr(image, i, j, k, dim) = this->addNoise_Uniform(arr(image, i, j, k, dim), lb, ub);
            }
        }
    }
}

/**
 * add salt&pepper noise to the matrix
 */
void TreeDrawer::addNoise_saltPepper(unsigned char valSalt, double probSalt, unsigned char valPepper, double probPepper){
    
	for(int i = 0; i < dim[0]; i++) {
		for(int j = 0; j < dim[1]; j++) {
			for(int k = 0; k < dim[2]; k++) {
				arr(image, i, j, k, dim) = this->addNoise_saltPepper(arr(image, i, j, k, dim), valSalt, probSalt, valPepper, probPepper);
            }
        }
    }
			
}

/**
 * add gaussian noise to the matrix
 */
void TreeDrawer::addNoise_gaussian(double median, double sigma){
    
	for(int i = 0; i < dim[0]; i++) {
		for(int j = 0; j < dim[1]; j++) {
			for(int k = 0; k < dim[2]; k++) {
				arr(image, i, j, k, dim) = this->addNoise_gaussian(arr(image, i, j, k, dim), median, sigma);
            }
        }
    }
    
}


/**
 * add a single shadow @ (x,y,z)
 */
void TreeDrawer::addShadow(double x, double y, double z, double centerRatio, double r){
    
	for(int i = 0; i < dim[0]; i++) {
		for(int j = 0; j < dim[1]; j++) {
			for(int k = 0; k < dim[2]; k++){
				
                int voxel[3]; voxel[0] = i; voxel[1] = j; voxel[2] = k;
				double point[3]; voxelToPoint(voxel, -1, point);
				double dist = sqrt(pow((point[0]-x), 2) + pow((point[1]-y), 2) + pow((point[2]-z), 2));
				
				if(dist > r) {
					dist = r;
                }

				arr(image, i, j, k, dim) *= (dist/r) + (1- dist/r)*centerRatio;
                
			}
        }
    }
    
}

/**
 * add 'numShadows' shadows to the matrix
 */
void TreeDrawer::addShadows(int numShadows){
    
	MTRand rand;
	for(int i = 0; i < numShadows; i++){
	
        int segment = (int)(rand.rand()*(vt->nt.nodes.size()-1))+1;
		double * endPoint = vt->nt.getPos(segment);
		double * startPoint = vt->nt.getPos(vt->nt.getParent(segment));
		double length = sqrt(pow((startPoint[0]-endPoint[0]), 2) + pow((startPoint[1]-endPoint[1]), 2) + pow((startPoint[2]-endPoint[2]), 2));

		double strength = rand.rand();
		this->addShadow((startPoint[0]+endPoint[0])/2, (startPoint[1]+endPoint[1])/2, (startPoint[2]+endPoint[2])/2, strength, length);
	}
    
}

/**
 * get the value of the matrix at (x,y,z)
 */
unsigned char TreeDrawer::imageAt(int x, int y, int z){
	return arr(image, x, y, z, dim);
}

/**
 * copy the matrix
 */
TreeDrawer* TreeDrawer::copy(){
	TreeDrawer * copy = new TreeDrawer(vt, imageVoxelWidth, c1, c2);

	for(int i = 0; i < dim[0]*dim[1]*dim[2]; i++) {
		copy->image[i] = image[i];
    }

	return copy;
}