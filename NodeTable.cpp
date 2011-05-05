/*=========================================================================
 
	Program: VascuSynth
	Module: $RCSfile: NodeTable.cpp,v $
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

#include "NodeTable.h"

#include <vector>
#include <map>
#include <iostream>

using namespace std;

int NodeTable::TERM = 0;
int NodeTable::ROOT = 1;
int NodeTable::BIF = 2;
int NodeTable::FIELDS = 12;

/**
 * constructor
 */
NodeTable::NodeTable(){
	prepareUndo = false;
}

/**
 * start recording entries for undos
 */
void NodeTable::startUndo(){
	prepareUndo = true;
	length = nodes.size();
}

/**
 * stop recording entries for undo
 */
void NodeTable::stopUndo(){
	prepareUndo = false;
}

/**
 * clear the undo table
 */
void NodeTable::clearUndo(){
    
    int size = (int) originalEntries.size();
	for(int i = 0; i < size; i++) {
		delete originalEntries[i];
    }

	originalEntries.clear();
}

/**
 * apply all accumulate undos
 */
void NodeTable::applyUndo(){
	stopUndo();
	map<int, double*>::iterator itr;
	for(itr = originalEntries.begin(); itr != originalEntries.end(); itr++){
		double *newNode = new double[FIELDS];
		for(int i = 0; i < FIELDS; i++)
			newNode[i] = itr->second[i];
		setNode(itr->first, newNode);
	}

	if( (int) nodes.size() > length) {
		nodes.resize(length);
    }
    
	clearUndo();
	startUndo();
}

/**
 * return the type of the node at index
 */
int NodeTable::getType(int index){
	if(index > (int) nodes.size()) {
		return -1;
    }
	return (int)nodes[index][0];
}

/**
 * set the type of a node at an index
 */
void NodeTable::setType(int index, int type){
    
	if(index < (int) nodes.size()) {

        //if this is the first time this entry has been modified - adds its original state to the undo table
        if(prepareUndo == true && originalEntries[index] == NULL){
            double * d = new double[FIELDS];
            double * ent = nodes[index];
            for(int i = 0 ; i < FIELDS; i++)
                d[i] = ent[i];
            originalEntries[index] = d;
        }

        nodes[index][0] = type;
        
    }
}

/**
 * get the position
 */
double* NodeTable::getPos(int index){
    
	if(index > (int) nodes.size()) {
		return NULL;
    }
	
	return nodes[index]+1;
}

/**
 * set the position
 */
void NodeTable::setPos(int index, double* source){
    
	if (index < (int) nodes.size()) {

        //if this is the first time this entry has been modified - adds its original state to the undo table
        if(prepareUndo == true && originalEntries[index] == NULL){
            double * d = new double[FIELDS];
            double * ent = nodes[index];
            for(int i = 0 ; i < FIELDS; i++)
                d[i] = ent[i];
            originalEntries[index] = d;
        }

        nodes[index][1] = source[0];
        nodes[index][2] = source[1];
        nodes[index][3] = source[2];
        
    }
}

/**
 * get the parent of the node at index
 */
int NodeTable::getParent(int index){
    
	if(index > (int) nodes.size()) {
		return NULL;
    }
	
	return (int) nodes[index][4];
}

/**
 * set the parent of the node at index
 */
void NodeTable::setParent(int index, int parent){
    
	if (index < (int) nodes.size()) {
		
        //if this is the first time this entry has been modified - adds its original state to the undo table
        if(prepareUndo == true && originalEntries[index] == NULL){
            double * d = new double[FIELDS];
            double * ent = nodes[index];
            for(int i = 0 ; i < FIELDS; i++)
                d[i] = ent[i];
            originalEntries[index] = d;
        }

        nodes[index][4] = parent;
        
    }
}

/**
 * get the ratio which is used to determine the radii of the branches
 * for the left segment
 */
double NodeTable::getLeftRatio(int index){
    
	if (index > (int) nodes.size()) {
		return NULL;
    }
    
	return nodes[index][5];
}

/**
 * set the left ratio for a node, ratio used for radii calculations
 */
void NodeTable::setLeftRatio(int index, double ratio){
    
	if (index < (int) nodes.size()) {

        //if this is the first time this entry has been modified - adds its original state to the undo table
        if(prepareUndo == true && originalEntries[index] == NULL){
            double * d = new double[FIELDS];
            double * ent = nodes[index];
            for(int i = 0 ; i < FIELDS; i++)
                d[i] = ent[i];
            originalEntries[index] = d;
        }
        nodes[index][5] = ratio;
        
    }
}

/**
 * get the ratio which is used to determine the radii of the branches
 * for the right segment
 */
double NodeTable::getRightRatio(int index){
    
	if (index > (int) nodes.size()) {
		return NULL;
    }
    
	return nodes[index][6];
}

/**
 * set the right ratio for a node, ratio used for radii calculations
 */
void NodeTable::setRightRatio(int index, double ratio){
    
	if (index < (int) nodes.size()) {

        //if this is the first time this entry has been modified - adds its original state to the undo table
        if(prepareUndo == true && originalEntries[index] == NULL){
            double * d = new double[FIELDS];
            double * ent = nodes[index];
            for(int i = 0 ; i < FIELDS; i++)
                d[i] = ent[i];
            originalEntries[index] = d;
        }
        nodes[index][6] = ratio;
        
    }
}

/**
 * get the flow at the node position
 */
double NodeTable::getFlow(int index){
    
	if (index > (int) nodes.size()) {
		return NULL;
    }
    
	return nodes[index][7];
}

/**
 * set the flow at the position
 */
void NodeTable::setFlow(int index, double flow){
    
	if (index < (int) nodes.size()) {

        //if this is the first time this entry has been modified - adds its original state to the undo table
        if(prepareUndo == true && originalEntries[index] == NULL){
            double * d = new double[FIELDS];
            double * ent = nodes[index];
            for(int i = 0 ; i < FIELDS; i++)
                d[i] = ent[i];
            originalEntries[index] = d;
        }
        nodes[index][7] = flow;
    
    }
}

/**
 * get the left child at a node
 */
int NodeTable::getLeftChild(int index){
    
	if(index > (int) nodes.size()) {
		return NULL;
    }
    
	return (int) nodes[index][8];
}

/**
 * set the left child at a node
 */
void NodeTable::setLeftChild(int index, int id){
    
	if (index < (int) nodes.size()) {

        //if this is the first time this entry has been modified - adds its original state to the undo table
        if(prepareUndo == true && originalEntries[index] == NULL){
            double * d = new double[FIELDS];
            double * ent = nodes[index];
            for(int i = 0 ; i < FIELDS; i++)
                d[i] = ent[i];
            originalEntries[index] = d;
        }
        nodes[index][8] = id;
    
    }
    
}

/**
 * get the right child at a node
 */
int NodeTable::getRightChild(int index){
    
	if (index > (int) nodes.size()) {
		return NULL;
    }
    
	return (int) nodes[index][9];
}

/**
 * set the right child at a node
 */
void NodeTable::setRightChild(int index, int id){
    
	if(index < (int) nodes.size()) {
        
        //if this is the first time this entry has been modified - adds its original state to the undo table
        if(prepareUndo == true && originalEntries[index] == NULL){
            double * d = new double[FIELDS];
            double * ent = nodes[index];
            for(int i = 0 ; i < FIELDS; i++)
                d[i] = ent[i];
            originalEntries[index] = d;
        }
        nodes[index][9] = id;
        
    }
}

/**
 * get the radius at a node position
 */
double NodeTable::getRadius(int index){
    
	if(index > (int) nodes.size()) {
		return NULL;
    }
    
	return nodes[index][11];
    
}

/**
 * set the radius at a node position
 */
void NodeTable::setRadius(int index, double radius){
    
	if (index < (int) nodes.size()) {
		
        //if this is the first time this entry has been modified - adds its original state to the undo table
        if(prepareUndo == true && originalEntries[index] == NULL){
            double * d = new double[FIELDS];
            double * ent = nodes[index];
            for(int i = 0 ; i < FIELDS; i++)
                d[i] = ent[i];
            originalEntries[index] = d;
        }
        nodes[index][11] = radius;
        
    }
}

/**
 * get the reduced resistence for a node at the index
 */
double NodeTable::getReducedResistance(int index){
    
	if (index > (int) nodes.size()) {
		return NULL;
    }
    
	return nodes[index][10];
    
}

/**
 * set the reduced resistence for a node at the index
 */
void NodeTable::setReducedResistence(int index, double resistance){
    
	if (index < (int) nodes.size()) {

        //if this is the first time this entry has been modified - adds its original state to the undo table
        if(prepareUndo == true && originalEntries[index] == NULL){
            double * d = new double[FIELDS];
            double * ent = nodes[index];
            for(int i = 0 ; i < FIELDS; i++) {
                d[i] = ent[i];
            }
            originalEntries[index] = d;
        }
        nodes[index][10] = resistance;
    
    }
    
}

/**
 * add a new node to the node table
 */
void NodeTable::addNode(double* node){
	double *row = new double[FIELDS];
	for(int i = 0; i < FIELDS; i++)
		row[i] = node[i];

	nodes.push_back(row);
}

/**
 * set a node in the node table
 */
void NodeTable::setNode(int index, double* node){
	//if this is the first time this entry has been modified - adds its original state to the undo table
	if(prepareUndo == true && originalEntries[index] == NULL){
		double * d = new double[FIELDS];
		double * ent = nodes[index];
		for(int i = 0 ; i < FIELDS; i++)
			d[i] = ent[i];
		originalEntries[index] = d;
	}
	
	if(index >= (int) nodes.size()) {
		nodes.resize(index+1);
    }
        
	delete nodes[index];
	nodes[index] = node;
}

/**
 * add a node to the node table and set the values of that node as well
 */
void NodeTable::addNode(int type, double* pos, int parent, double leftRatio, double rightRatio, double flow, int leftChild, int rightChild){
	double * row = new double[FIELDS];
	row[0] = type;
	row[1] = pos[0];
	row[2] = pos[1];
	row[3] = pos[2];
	row[4] = parent;
	row[5] = leftRatio;
	row[6] = rightRatio;
	row[7] = flow;
	row[8] = leftChild;
	row[9] = rightChild;
	
	nodes.push_back(row);
}

/**
 * copy the node table
 */
NodeTable NodeTable::copy(){
	NodeTable * nt = new NodeTable();
	vector<double*>::iterator itr;

	for(itr = nodes.begin(); itr != nodes.end(); itr++)
		nt->addNode(*itr);
	
	return *nt;
}
