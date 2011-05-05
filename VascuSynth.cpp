/*=========================================================================

	Program: VascuSynth
	Module: $RCSfile: VascuSynth.cpp,v $
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

//commands for mkdir
#ifdef _WIN32
#include "direct.h"
#else
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <iterator>
#include <cmath>

using namespace std;

/**
 * Utility function to read a text file and store the lines into a vector
 * that is returned.  This way we do not have to have a bunch of files open at
 * the same time.
 * @param const char* filename the filename to read
 * @return vector<string> a vector of file lines from the file
 * @throws string exception if the file cannot be read
 */
vector<string> * readFileLines(const char * filename){
	
	ifstream oFile;
	oFile.open(filename, ios::in);
	vector<string> * lines = new vector<string>;
	string line;
		
	if(oFile.is_open()){
	
		while(!oFile.eof()){
			getline(oFile, line);
			string copy = line;
			lines->push_back(copy);
		}
		
		oFile.close();
		
	} else {
		throw "Could not open file " + ( (string) filename);
	}
	
	return lines;
	
}


#include "OxygenationMap.h"
#include "SupplyMap.h"
#include "TreeDrawer.h"
#include "VascularTree.h"



/**
 * make dir that is cross platform
 */
int mmkdir(const char * dirname) {
#ifdef _WIN32
	return mkdir(dirname);
#else
    return mkdir(dirname, 0777);
#endif
}

/**
 * itoa is non standard so define it and use it,
 * converts an integer to a string
 */
string itoa(int value, int base) {
	
	string buf;
	
	// check that the base if valid
	if (base < 2 || base > 16) return buf;
	
	enum { kMaxDigits = 35 };
	buf.reserve( kMaxDigits ); // Pre-allocate enough space.
	
	int quotient = value;
	
	// Translating number to string with base:
	do {
		buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
		quotient /= base;
	} while ( quotient );
	
	// Append the negative sign
	if ( value < 0) buf += '-';
	
	reverse( buf.begin(), buf.end() );
	return buf;
}

/**
 *  Reads the parameters from the parameter file and then builds
 *  the vascular structure in the form of a tree
 * 
 *	Parameter File Entries:
 * 
 *	SUPPLY_MAP: supply map file
 *	OXYGENATION_MAP: oxygenation map file
 *	PERF_POINT: perf_x perf_y perf_z
 *	PERF_PRESSURE: perfussion pressure
 *	TERM_PRESSURE: termination pressure
 *	PERF_FLOW:	perfusion flow
 *	RHO: rho
 *	GAMMA: gamma
 *	LAMBDA: lambda
 *	MU: mu
 *	MIN_DISTANCE: minDistance
 *	NUM_NODES: numNodes 
 *	VOXEL_WIDTH: voxelWidth
 *	CLOSEST_NEIGHBOURS: closestNeighbours
 */
VascularTree * buildTree(const char * filename){
	
	SupplyMap * sm = NULL;
	OxygenationMap * om = NULL;
    
	double* perf = new double[3];
	double pperf;
	double pterm;
	double qperf;
	double rho;
	double gamma;
	double lambda;
	double mu;
	double minDistance;
	int numNodes;
	double voxelWidth;
	int closestNeighbours;
	int randomSeed = -1;
	string line;
    string supplyMapFileName;
    string oxygenMapFileName;
    
    //c++ doesn't allow us to check for undefined variables
    //so we need a boolean so that we know that the variables
    //are defined.
    bool perfSet = false;
    bool pperfSet = false;
    bool ptermSet = false;
    bool qperfSet = false;
    bool rhoSet = false;
    bool gammaSet = false;
    bool lambdaSet = false;
    bool muSet = false;
    bool minDistanceSet = false;
    bool numNodesSet = false;
    bool voxelWidthSet = false;
    bool closestNeighboursSet = false;
    bool supplyMapFileNameSet = false;
    bool oxygenMapFileNameSet = false;
	
	vector<string> *mapFilesLines = readFileLines(filename);
    int size = (int) mapFilesLines->size();
	
	for (int i=0; i < size; i++) {
		
		line = mapFilesLines->at(i);
		
		if (line.compare("") == 0) {
			break;
		}
		
		int colonPosition = line.find(":");
		string name = line.substr(0, colonPosition);
		string value = line.substr(colonPosition+2);
		
		if (name.compare("SUPPLY_MAP") == 0) {
			
            //store the supplyMapFileName for later
            //the ordering of things matters so just save the file name
            //then initialize the map later so that the user can type the fields
            //in whatever order
            supplyMapFileName = value;
            supplyMapFileNameSet = true;
			
		} else if (name.compare("OXYGENATION_MAP") == 0) {
            
            //same as supply map since we need random seed before
            //we init the oxygenation map
            oxygenMapFileName = value;
            oxygenMapFileNameSet = true;
             
		} else if (name.compare("PERF_POINT") == 0) {
			
			int spacePosition = value.find(" ");
			string pointvalue = value.substr(0, spacePosition);
			perf[0] = atof(pointvalue.c_str());
			
			int spacePosition2 = value.find(" ", spacePosition+1);
			pointvalue = value.substr(spacePosition+1, spacePosition2);
			perf[1] = atof(pointvalue.c_str());
			
			pointvalue = value.substr(spacePosition2+1);
			perf[2] = atof(pointvalue.c_str());
            
            perfSet = true;
						
		} else if (name.compare("PERF_PRESSURE") == 0){
			
			pperf = atof(value.c_str());
            pperfSet = true;
			
		} else if (name.compare("TERM_PRESSURE") == 0){
			
			pterm = atof(value.c_str());
            ptermSet = true;
			
		} else if (name.compare("PERF_FLOW") == 0){
			 
			qperf = atof(value.c_str());
            qperfSet = true;
			
		} else if (name.compare("RHO") == 0){
			
			rho = atof(value.c_str());
            rhoSet = true;
			
		} else if (name.compare("GAMMA") == 0){
			
			gamma = atof(value.c_str());
            gammaSet = true;
			
		} else if (name.compare("LAMBDA") == 0){
			 
			lambda = atof(value.c_str());
            lambdaSet = true;
			
		} else if (name.compare( "MU") == 0){
			 
			mu = atof(value.c_str());
            muSet = true;
			
		} else if (name.compare("MIN_DISTANCE") == 0){
			
			minDistance = atof(value.c_str());
            minDistanceSet = true;
			 
		} else if (name.compare("NUM_NODES") == 0){
			 
			numNodes = atoi(value.c_str());
            numNodesSet = true;
			
		} else if (name.compare("VOXEL_WIDTH") == 0){
			
			voxelWidth = atof(value.c_str());
            voxelWidthSet = true;
			
		} else if (name.compare("CLOSEST_NEIGHBOURS") == 0){
			 
			closestNeighbours = atoi(value.c_str());
            closestNeighboursSet = true;
			
		} else if (name.compare("RANDOM_SEED") == 0){
			 
			randomSeed = atoi(value.c_str());
			
		} else {
			
		}
		
		
	}
        
    //make sure that we have everything defined
    if (perfSet && pperfSet && ptermSet && qperfSet && rhoSet && gammaSet && lambdaSet && muSet && minDistanceSet && numNodesSet && voxelWidthSet && closestNeighboursSet && supplyMapFileNameSet && oxygenMapFileNameSet) {
    
        //load the supply map
        sm = new SupplyMap();
        
        try {
            sm->loadMap(supplyMapFileName);
        } catch (char * str) {
            throw (string) str;
        }
        
        //load the oxygenation map, rand seed will be -1 if there
        //is no randomSeed specified or it will be whatever the user specifies
        om = new OxygenationMap(sm, randomSeed);
        
        try {
            om->loadMap(oxygenMapFileName);
        } catch (char * str) {
            throw (string) str;
        }
        om->supply = sm;

        //TODO: should probably check and make sure everything is defined
        //and throw an error if it is not (and catch the error in the main function)
        VascularTree *vt = new VascularTree(om, perf, pperf, pterm, qperf, rho, gamma, lambda, mu, minDistance, numNodes, voxelWidth, closestNeighbours);
        vt->buildTree();
        
        return vt;
        
    } else {
     
        string errorStr = "Error while parsing the parameter file, not all parameters have been defined.";
        throw errorStr;
        
    }
}

/**
 * draws the tree to a matrix
 */
TreeDrawer *drawTree(VascularTree * vt, double * corner1, double * corner2, double imageVoxelWidth){
	TreeDrawer * td = new TreeDrawer(vt, imageVoxelWidth, corner1, corner2);
	td->drawImage();
	return td;
}

/**
 * draws the tree from the TreeDrawer into a volumetric 3D
 * image as a series of 2D png slices
 */
void drawImage(TreeDrawer * td, const char* rootName){
	typedef unsigned char PixelType;
	const unsigned int Dimension = 3;
	typedef itk::Image< PixelType, Dimension > ImageType;

	ImageType::Pointer image = ImageType::New();

	ImageType::SizeType size;
	size[0] = td->dim[0]; // size along X
	size[1] = td->dim[1]; // size along Y
	size[2] = td->dim[2]; // size along Z

	ImageType::IndexType start;
	start[0] = 0; // first index on X
	start[1] = 0; // first index on Y
	start[2] = 0; // first index on Z

	ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	image->SetRegions( region );
	image->Allocate();
	
	ImageType::IndexType pixelIndex;
	pixelIndex[0] = 0; // x position
	pixelIndex[1] = 0; // y position
	pixelIndex[2] = 0; // z position

	for(int i = 0; i < td->dim[0]; i++){
		for(int j = 0; j < td->dim[1]; j++){
			for(int k = 0 ; k < td->dim[2]; k++){
				pixelIndex[0] = i;
				pixelIndex[1] = j;
				pixelIndex[2] = k;
				
				image->SetPixel(pixelIndex, td->imageAt(i, j, k));
			}
		}
	}


	typedef itk::Image< unsigned char, 2 > Image2DType;
	typedef itk::ImageSeriesWriter< ImageType, Image2DType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput( image);

	typedef itk::NumericSeriesFileNames NameGeneratorType;
	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();

	std::string format = rootName;
	format += "%03d";
	format += ".jpg";
	nameGenerator->SetSeriesFormat( format.c_str() );

	const unsigned int firstSlice = start[2];
	const unsigned int lastSlice = start[2] + size[2] - 1;
	nameGenerator->SetStartIndex( firstSlice );
	nameGenerator->SetEndIndex( lastSlice );
	nameGenerator->SetIncrementIndex( 1 );

	writer->SetFileNames( nameGenerator->GetFileNames() );

	try{
		writer->Update();
	}catch( itk::ExceptionObject & excp ){
        
        throw "Exception thrown while reading the image";
        
	}

	return;
}

/** 
 *  applies noise to the volumetric image
 * 
 * noise.txt format:
 * 
 * SHADOW: numShadows
 * GAUSSIAN: median sigma
 * UNIFORM: lb ub
 * SALTPEPER: valSalt probsalt valpepper probpepper
 * 
 * noise will be added in the order specified by the file
 * 
 * one image for each noise file will be generated
 */
void applyNoise(TreeDrawer *td, const char* noiseFile){
    
	ifstream mapFile;
	mapFile.open(noiseFile);
	string line;

	double lb, ub, median, sigma, probSalt, probPepper;
	int numShadows;
	char valSalt, valPepper;

	if(mapFile.is_open()){		
		while(!mapFile.eof()){
			getline(mapFile, line);
			char* tok = new char[line.size()];
			strcpy(tok, line.c_str());
			
			if(line.length() == 0)
				continue;

			char * field = strtok(tok, ":");
					
			if(strcmp(field, "SHADOW") == 0){
				
				//apply shadow noise
				numShadows = atoi(strtok(NULL, " "));
				
                td->addShadows(numShadows);
			
			} else if(strcmp(field, "GAUSSIAN") == 0){
				
				//apply gaussian noise
				median = atof(strtok(NULL, " "));
				sigma = atof(strtok(NULL, " "));
				                
				td->addNoise_gaussian(median, sigma);
				
			} else if(strcmp(field, "UNIFORM") == 0){
				
				//applying uniform noise
				lb = atof(strtok(NULL, " "));
				ub = atof(strtok(NULL, " "));
				
                td->addNoise_Uniform(lb, ub);
				
			} else {
				
				if (strcmp(field, "SALTPEPPER") == 0) {
				
					//apply salt and pepper noise
					valSalt = (char)atoi(strtok(NULL, " "));
					probSalt = atof(strtok(NULL, " "));
					valPepper = (char)atoi(strtok(NULL, " "));
					probPepper = atof(strtok(NULL, " "));
					valPepper = (char)valPepper;
					
					td->addNoise_saltPepper(valSalt, probSalt, valPepper, probPepper);
				
				}
				
			} 
		}
	
    } else {
    
        throw "Could not read the noise file";
        
    }
	
	mapFile.close();
}


/**
 * prints a node into XML/GXL format from the NodeTable
 */
void subPrint_node(NodeTable *nt, int segment, ofstream &os){
	os<<"  <node id=\"n"<<segment<<"\">"<<endl;
	os<<"    <attr name=\" nodeType\">"<<endl;
	if(nt->getType(segment) == NodeTable::ROOT){
		os<<"      <string> root node </string>"<<endl;
	} else if(nt->getType(segment) == NodeTable::TERM){
		os<<"      <string> terminal node </string>"<<endl;
	} else if(nt->getType(segment) == NodeTable::BIF){
		os<<"      <string> bifurication </string>"<<endl;
	} else {
		os<<"      <string> unknown type </string>"<<endl;
	}
	os<<"    </attr>"<<endl;

	os<<"    <attr name=\" position\">"<<endl;
	os<<"      <tup>"<<endl;
	double *pos = nt->getPos(segment);
	os<<"        <float>"<<pos[0]<<"</float>"<<endl;
	os<<"        <float>"<<pos[1]<<"</float>"<<endl;
	os<<"        <float>"<<pos[2]<<"</float>"<<endl;
	os<<"      </tup>"<<endl;
	os<<"    </attr>"<<endl;	
	os<<"  </node>"<<endl;

	if(nt->getType(segment) != NodeTable::TERM){
		subPrint_node(nt, nt->getLeftChild(segment), os);
		if(nt->getType(segment) != NodeTable::ROOT)
			subPrint_node(nt, nt->getRightChild(segment), os);
	}
}

/**
 * prints an edge into XML/GXL format from a node table
 */
void subPrint_edge(NodeTable *nt, int segment, ofstream &os){

	if(nt->getType(segment) != NodeTable::ROOT){
		os<<"  <edge id=\"e"<<segment<<"\" to=\"n"<<segment<<"\" from=\"n"<<nt->getParent(segment)<<"\">"<<endl;
		os<<"    <attr name=\" flow\">"<<endl;
		os<<"      <float>"<<nt->getFlow(segment)<<"</float>"<<endl;
		os<<"    </attr>"<<endl;

		os<<"    <attr name=\" radius\">"<<endl;
		os<<"      <float>"<<nt->getRadius(segment)<<"</float>"<<endl;
		os<<"    </attr>"<<endl;

		os<<"  </edge>"<<endl;
	}

	if(nt->getType(segment) != NodeTable::TERM){
		subPrint_edge(nt, nt->getLeftChild(segment), os);
		if(nt->getType(segment) != NodeTable::ROOT)
			subPrint_edge(nt, nt->getRightChild(segment), os);
	}
}

/**
 * Writes the tree structure to a GXL file and stores information
 * about the nodes, edges, their heirarchy, radii and bifurcation locations
 */
void printTreeStructure(VascularTree * vt, const char * filePath){

	ofstream output;

	//writing the tree structure as GXL to the filePath specified
	output.open(filePath);
	output<<"<gxl><graph id=\""<<filePath<<"\" edgeids=\" true\" edgemode=\" directed\" hypergraph=\" false\">"<<endl;
	
	//this seems really really stupid to do, why would we
	//loop through the entire structure to find the root, and then
	//recursively the nodes and edges?
	//TODO: fix this - should have the root always stored and easily
	//accessible and then recursively output the node/edges
	for(int i = 0; i < vt->nt.length; i++){
		if(vt->nt.getType(i) == NodeTable::ROOT){
			subPrint_node(&vt->nt, i, output);
			subPrint_edge(&vt->nt, i, output);
			output<<"</graph></gxl>"<<endl;
			output.close();
			return;
		}
	}

    output.close();
    
	throw "Unable to find root node.  The GXL file has not been generated.";
	
}


/**
 * VascuSynth: takes a series of parameter files, image names and (optionally) noise files
 * and generates a vascular structure based on the parameters.  The 3d volume is saved
 * as a series of 2D png slices.  Information about the vascular structure is saved as
 * a GXL file that can be visualized using software such as GraphViz.
 *
 * Arguments: paramFile.txt imageNames.txt voxelWidth noiseFiles.txt
 * noiseFiles.txt is an optional parameter.
 *
 * For each imageName/parameterFile/noiseFile, a folder is generated with the name specified
 * which will contain the 2D slices and the GXL file.
 */
int main(int argc, char** argv){

	if (argc < 4) {
		//not enough parameters specified
		cout << "An error has occured: incorrect number of arguments" << endl;
		cout << "Usage: VascuSynth [paramFile] [imageNameFile] [voxelWidth]" << endl;
		return 0;
	}
	
    try {
    
        //read the param files and image name files
        vector<string> *paramFiles = readFileLines(argv[1]);
        vector<string> *imageNameFiles = readFileLines(argv[2]);

        //voxel widths
        string voxelWidth = argv[3];
        string *noiseFiles = new string[argc-4];
        
        int paramFilesSize = (int) paramFiles->size();
        
        //go through each param file and build tree an spit it out
        for(int m = 0; m < paramFilesSize; m++){
            
            string paramFile = paramFiles->at(m);
            string rootDirectory = imageNameFiles->at(m);

            for(int i = 4; i < argc; i++) {
                noiseFiles[i-4] = argv[i];
            }
            
            int numNoise = argc-4;
            
            cout << "Reading parameters and building the tree..." << endl;
            
            //build the tree
            VascularTree * vt = buildTree(paramFile.c_str());
            
            cout << "The vascular tree has been built sucessfully..." << endl;
            
            //filter out the /r that appears at times and messes up the directory name
            if (rootDirectory[ rootDirectory.length() - 1 ] == '\r') {
                rootDirectory = rootDirectory.substr(0, rootDirectory.length()-1);
            }
            
            //create the directory for the output
            rootDirectory = "./"+rootDirectory;
            mmkdir(rootDirectory.c_str());
            
            //create the GXL file
            string treeStructure = rootDirectory+"/tree_structure.xml";
            printTreeStructure(vt, treeStructure.c_str());
            
            cout << "The directory for the image has been created..." << endl;
            cout << "Information about the vascular structure has been saved in the gxl file tree_structure.xml..." << endl;

            double corner1[] = {0,0,0};
            double corner2[] = {vt->oxMap->dim[0], vt->oxMap->dim[1], vt->oxMap->dim[2]};
            TreeDrawer * td = drawTree(vt, corner1, corner2, atof(voxelWidth.c_str()));
                
            //create the subdirectory for the images
            string imageName = rootDirectory+"/original_image";
            mmkdir(imageName.c_str());

            //output the images
            imageName = imageName + "/image";
            drawImage(td, imageName.c_str());
            
            cout << "The volumetric image has been saved..." << endl;

            char *buff = new char[20];
            
            if (numNoise > 0) {
                cout << "The images are being degraded by noise..." << endl;
            }
            
            //apply noise to the images - creating niose_images
            for(int i = 0; i < numNoise; i++){
                
                TreeDrawer * td_c = td->copy();
                applyNoise(td_c, noiseFiles[i].c_str());
                
                string noiseImage = rootDirectory+"/noise_image_"+itoa(i, 10);
                mmkdir(noiseImage.c_str());
                
                noiseImage = noiseImage+"/image";
                drawImage(td_c, noiseImage.c_str());
            
            }
            
            if (numNoise > 0) {
                cout << "Images have been succesfully degraded by noise and saved..." << endl;
            }

            //clean up
            delete buff;
            delete td;
            delete vt;
        }
        
    } catch (string str) {
        cout << "ERROR: " << str << endl;
        cout << "Exiting VascuSynth" << endl;
    }

	return 0;
}

