//
//  Pulse.h
//  PositionReconstruction
//
//  Created by Chang Lee on 6/3/13.
//  Copyright (c) 2013 Chang Lee. All rights reserved.
//

#ifndef __PositionReconstruction__Pulse__
#define __PositionReconstruction__Pulse__
#define BOOST_UBLAS_TYPE_CHECK 0 

#include <iostream>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "PBomb.h"
#include "../Utilities/RQFile_IO.hh"

namespace ublas=boost::numeric::ublas;
using namespace std;
using boost::math::chi_squared;

//__________________________________________________________________________
/* Matrix inversion routine.
    Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix(const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
 	typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;
 	// create a working copy of the input
 	ublas::matrix<T> A(input);
 	// create a permutation matrix for the LU-factorization
 	pmatrix pm(A.size1());
 	// perform LU-factorization
 	int res = lu_factorize(A,pm);
        if( res != 0 ) return false;
 	// create identity matrix of "inverse"
 	inverse.assign(ublas::identity_matrix<T>(A.size1()));
 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);
 	return true;
}
//__________________________________________________________________________
struct sort_by_ChiSq
{
    bool operator()(PBomb* a, PBomb* b)	{return a->ReadChiSq() < b->ReadChiSq();}
};
//__________________________________________________________________________
struct sort_by_NLL
{
    bool operator()(PBomb* a, PBomb* b)	{return a->ReadNLL() < b->ReadNLL();}
};
//__________________________________________________________________________
class Pulse{
    
public:
    
	// default constructor
	Pulse();
    
	// overload constructor with p pulse index
	Pulse( size_t, size_t, map<size_t,size_t>&, vector<size_t>&,size_t);
    
	// accessors &  mutators
	long int WidthCheck(Rq*, Rq*);
	bool ReadTopAreaAndCheck(Rq*,size_t);
	void FillNaN(Rq*, Rq*);
	void GeneratePMTKey(size_t);
	void GenerateScarceMapKey(size_t, HitLib*, vector<float>&);
	void SelectPBombsFromKey(HitLib*);
	void PBombLoop(vector<float>&);
	void InterpolateXY(size_t);
	double ReadFinalX();
	double ReadFinalY();
	double ReadSigmaRad();
	double ReadChiSqP();
	void Clear();

private:

	// internal methods
	double ChiSq(PBomb*,vector<float>&);
	double Likelihood(PBomb*);
	//float Resolution();
	double LikelihoodConicalFit(size_t&);
	void PrintLikelihood(size_t);
	
	size_t statMethod;	// 0 for chi-squared, 1 for likelihood
	size_t event_num,pulse_num;
	float logGamma;
	int ndf;
	float top_area_sum, minNLL;
	float xFinal,yFinal, sigmaRadius, chisqP;
	vector<size_t> liveTopPMT,bestKeys;
	vector<float> topHit;
	vector<PBomb*> sortedPBombs,hitMapstoScan,scarcePBombs;
	map<size_t,size_t> pmtIdx;

};

#endif /* defined(__PositionReconstruction__Pulse__) */
