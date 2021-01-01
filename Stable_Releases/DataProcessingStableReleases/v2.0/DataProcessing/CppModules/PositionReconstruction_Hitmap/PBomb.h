//
//  PBomb.h
//  PositionReconstruction
//
//  Created by Chang Lee on 3/20/13.
//  Copyright (c) 2013 Chang Lee. All rights reserved.
//

#ifndef __PositionReconstruction__PBomb__
#define __PositionReconstruction__PBomb__

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <ctime>
#include <cmath>
#include <string>
#include <cstddef>

#include "../Utilities/RQFile_IO.hh"

#define nTopPMT	61		// number of top PMTs in LUX
#define OFFSET 4
#define VERBOSITY 0

using namespace std;

//______________________________________________________________________________
class PBomb;
class HitLib;
typedef std::multimap<size_t,PBomb*> HitMapCon;
typedef HitMapCon::iterator HitMapItr;
//______________________________________________________________________________

class HitLib{
	
public:
	
	// constructor
	HitLib();
	
	// Build a hit map with bin file path and the disabled PMTs
	bool BuildFrom(const string&, const vector<size_t>&, const size_t);
	
	// accsessors
	vector<PBomb*> SearchKey(size_t);
	pair<HitMapItr,HitMapItr> ScarceMapRange();
	
	// mutators
	void Delete();
	
	// a vector containg top PMT serial numbers
	vector<size_t> liveTopPMT;
	map<size_t,size_t> pmtIdx;

private:
	
	// the actual hit pattern library
	HitMapCon HitMap;
	// scarce map for fast key search
	HitMapCon ScarceMap;
};

//______________________________________________________________________________

class PBomb{
    
public:
    
	// default constructor
	PBomb();
    
	// overload constructor with x, y coordinates
	PBomb(size_t, size_t, float, float);
    
	// accessors
	float ReadX() const;
	float ReadY() const;
	size_t ReadPMT(size_t);
	float ReadNrmPMT(size_t);
	size_t ReadTotalPhe();
	float ReadChiSq();
	float ReadNLL();
	size_t ReadScarceMapKey();
	size_t ReadPMTKey();
    
	// mutators
	void WritePMTHit(size_t, size_t);		// index, hit count
	void AddTotalPhe(size_t);
	void AddChiSq(float);
	void AddNLL(float);
	void ResetStat();

private:
    
	// values to be imported
	size_t PMT_key;
	size_t scarce_map_key;
	float xcoordinate;
	float ycoordinate;
	size_t PMT[nTopPMT]; // stores PMT hits. Index is PMTnumber-1, but the last element is PMT121
	
	// values to be computed
	float chiSqSum;
	float likelihood;
	size_t totalPhe;
	
};

#endif /* defined(__PositionReconstruction__PBomb__) */
