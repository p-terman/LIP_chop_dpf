//
//  PBomb.cpp
//  PositionReconstruction
//
//  Created by Chang Lee on 3/20/13.
//  Copyright (c) 2013 Chang Lee. All rights reserved.
//

#include "PBomb.h"
#include <algorithm>

//_________________________________________________________________________________
HitLib::HitLib()
{
	
}
//_________________________________________________________________________________
bool HitLib::BuildFrom(const string& filepath, const vector<size_t>& disabledPMTs, const size_t keyType)
{
	size_t numOfRow = OFFSET+nTopPMT;	// +2 for PMT key, x, y coordinates
	size_t size;
	int key;
	float * memblock;			// declare buffer
	
	clock_t	t;
	t = clock(); //tic

	ifstream rawBin;
	rawBin.open (filepath.c_str(),ios::in|ios::binary|ios::ate);
	if (rawBin.fail()){cout<<"Could not open the Hitmap file"<<endl; return 1;}
	
	rawBin.seekg (0, ios::end);
	size = (size_t)rawBin.tellg() / sizeof(float);
    
	memblock = new float [size];
    
	// rewind the pointer, read the file, and close the file
	rawBin.seekg (0, ios::beg);
	rawBin.read (reinterpret_cast<char*>(memblock), size*sizeof(float));
	rawBin.close();
	
#if VERBOSITY>0
	// print the number of top PMT + key + x + y coordinates
	cout << "The file size is " << numOfRow << " x " << size/numOfRow<< endl;
#endif

	// read in and store PMT numbers
	pmtIdx.clear();
	liveTopPMT.clear();
	size_t PMTNum;
	for(size_t it=OFFSET;it<numOfRow;++it){
		// make PMT serial index map
		PMTNum=memblock[it];
		pmtIdx.insert(pair<size_t,size_t>(PMTNum,it-OFFSET));
		// collect the live PMTs
		if (find(disabledPMTs.begin(),disabledPMTs.end(),PMTNum)==disabledPMTs.end())
			liveTopPMT.push_back(PMTNum);
	}

	for(size_t it=numOfRow;it<size; ++it){
		
		// at the beginning of each column
		if(it % numOfRow==0){
			
			// Make a PBomb with given info
			PBomb* tmpPBomb =
			new PBomb(fabs(memblock[it]),memblock[it+1],100*memblock[it+2],100*memblock[it+3]);
			
			// store the PMT hit info
			for (size_t pt=0; pt<nTopPMT; ++pt){	
				tmpPBomb->WritePMTHit(pt, memblock[it+pt+OFFSET]);
				if (find(disabledPMTs.begin(),disabledPMTs.end(),
						memblock[pt+OFFSET])==disabledPMTs.end())
					tmpPBomb->AddTotalPhe(memblock[it+pt+OFFSET]);
			}
			// Scarce map key
			if (keyType==0) {
				key = memblock[it];
				if (key<0)
					ScarceMap.insert(pair<size_t,PBomb*>(-1*key,tmpPBomb));
				HitMap.insert(pair<size_t,PBomb*>(fabs(key),tmpPBomb));
			}
			// hottest PMT key
			if (keyType==1)
				HitMap.insert(pair<size_t,PBomb*>(memblock[it+1],tmpPBomb));
		}
	}
		
	// garbage collection
	delete[] memblock;
	memblock = NULL;
    
	t = clock()-t; //toc
#if VERBOSITY>0
	cout<< "Hitmap is successfully imported in "
		<< ((float)t)/CLOCKS_PER_SEC<< " seconds!"<<endl;
#endif
	
	return 0;
}
//_________________________________________________________________________________
vector<PBomb*> HitLib::SearchKey(size_t key)
{
	pair<HitMapItr, HitMapItr> searchRange = HitMap.equal_range(key);

	vector<PBomb*> hitMapstoScan;
	for (HitMapItr ki=searchRange.first; ki!=searchRange.second; ++ki)
		hitMapstoScan.push_back(ki->second);
	
	return hitMapstoScan;
}
//_________________________________________________________________________________
pair<HitMapItr,HitMapItr> HitLib::ScarceMapRange()
{
	return pair<HitMapItr,HitMapItr> (ScarceMap.begin(),ScarceMap.end());
}
//_________________________________________________________________________________
void HitLib::Delete()
{
	for (HitMapItr ci=HitMap.begin(); ci!=HitMap.end();++ci)
	{delete ci->second; ci->second=NULL;}
	
	HitMap.clear();
}
//_________________________________________________________________________________
PBomb::PBomb()
{
	xcoordinate = 0;
	ycoordinate = 0;
	totalPhe = 0;
	chiSqSum = 0;
	likelihood = 0;
}
//_________________________________________________________________________________
PBomb::PBomb(size_t scarce, size_t pmt, float x, float y)
{
	PMT_key = pmt;
	scarce_map_key = scarce;
	xcoordinate = x;
	ycoordinate = y;
	totalPhe = 0;
	chiSqSum = 0;
	likelihood = 0;
}
//_________________________________________________________________________________
float PBomb::ReadX() const
{
	return xcoordinate;
}
//_________________________________________________________________________________
float PBomb::ReadY() const
{
	return ycoordinate;
}
//_________________________________________________________________________________
size_t PBomb::ReadPMT(size_t index)
{
	return PMT[index];
}
//_________________________________________________________________________________
float PBomb::ReadNrmPMT(size_t index)
{
	return (float)ReadPMT(index)/(float)totalPhe;
}
//_________________________________________________________________________________
size_t PBomb::ReadTotalPhe()
{
	return totalPhe;
}
//_________________________________________________________________________________
float PBomb::ReadChiSq()
{
	return chiSqSum;
}
//_________________________________________________________________________________
float PBomb::ReadNLL()
{
	return likelihood;
}
//_________________________________________________________________________________
size_t PBomb::ReadScarceMapKey()
{
	return scarce_map_key;
}
//_________________________________________________________________________________
size_t PBomb::ReadPMTKey()
{
	return PMT_key;
}
//_________________________________________________________________________________
void PBomb::WritePMTHit(size_t index, size_t hit)
{
	PMT[index] = hit;
}
//_________________________________________________________________________________
void PBomb::AddTotalPhe(size_t phe)
{
	totalPhe += phe;
}//_________________________________________________________________________________
void PBomb::AddChiSq(float chi)
{
	chiSqSum += chi;
}
//_________________________________________________________________________________
void PBomb::AddNLL(float NLL)
{
	likelihood += NLL;
}
//_________________________________________________________________________________
void PBomb::ResetStat()
{
	chiSqSum = 0;
	likelihood = 0;
}
