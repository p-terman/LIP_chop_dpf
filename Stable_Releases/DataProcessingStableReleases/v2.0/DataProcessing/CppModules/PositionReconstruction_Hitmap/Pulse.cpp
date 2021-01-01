//
//  Pulse.cpp
//  PositionReconstruction
//
//  Created by Chang Lee on 6/3/13.
//  Copyright (c) 2013 Chang Lee. All rights reserved.
//

#include "Pulse.h"
//__________________________________________________________________________
double chisqr(int Dof, double Cv)
{
	boost::math::chi_squared dist(Dof);
	return boost::math::cdf(dist,Cv);
}
//__________________________________________________________________________
Pulse::Pulse()
{
	
}
//__________________________________________________________________________
Pulse::Pulse( size_t e, size_t p, map<size_t,size_t>&idx, vector<size_t>& v, size_t method)
{
	event_num = e;
	pulse_num = p;
	liveTopPMT = v;
	statMethod = method;
	pmtIdx=idx;
#if VERBOSITY>0
	cout <<endl<<"Pulse "<<p<<": ";
#endif
}
//__________________________________________________________________________
long int Pulse::WidthCheck(Rq* pulse_end_samples, Rq* pulse_start_samples)
{
	long int pulse_width = pulse_end_samples->GetInt(pulse_num)
							- pulse_start_samples->GetInt(pulse_num);
#if VERBOSITY>0
	cout <<pulse_width<<" t-chnls, ";
#endif
	
	return pulse_width;
}
//__________________________________________________________________________
bool Pulse::ReadTopAreaAndCheck(Rq* areas, size_t threshold)
{
	float hit;
	topHit.clear();
	topHit.resize(liveTopPMT.size());
	top_area_sum=0;
	ndf=0;
	logGamma=0;
	
	for (size_t tp=0; tp<liveTopPMT.size(); ++tp){
		hit=areas->GetDouble(pulse_num, liveTopPMT[tp]-1);
        if (hit<0) {topHit[tp]=0; continue;}
		if (hit>0) ++ndf;
		logGamma+=lgamma(hit+1);
		topHit[tp]=hit;
		top_area_sum+=hit;
	}
	
#if VERBOSITY>0
	cout <<top_area_sum<<" phe hit top, ";
#endif
	
	// If there is too few phe, use likelihood method
	//if (top_area_sum < threshold) statMethod=1;
	
	return ndf==0;
}
//__________________________________________________________________________
void Pulse::FillNaN(Rq* x_cm_tmplt, Rq* y_cm_tmplt)
{	
	x_cm_tmplt->SetDouble(-100,pulse_num);
	y_cm_tmplt->SetDouble(-100,pulse_num);
	
	Clear();
}
//__________________________________________________________________________
void Pulse::GeneratePMTKey(size_t nKey)
{
	bestKeys.clear();
	bestKeys.resize(nKey);
	
	size_t maxIdx;
	vector<float> topHitCopy(topHit);
	float minHit=*min_element(topHitCopy.begin(),topHitCopy.end());
	for (size_t m=0;m<nKey;++m){
		maxIdx=max_element(topHitCopy.begin(),topHitCopy.end())-topHitCopy.begin();
		bestKeys[m]=liveTopPMT[maxIdx];
		topHitCopy[maxIdx]=minHit;
	}

#if VERBOSITY>0
	cout<<"max hit at PMT"<<bestKeys.front();
#endif
}
//__________________________________________________________________________
void Pulse::GenerateScarceMapKey(size_t nKey, HitLib* Hitmap, vector<float>& ChiUnc)
{
	hitMapstoScan.clear();
	pair<HitMapItr,HitMapItr> range=Hitmap->ScarceMapRange();
	hitMapstoScan.reserve(distance(range.first, range.second));
	for (HitMapItr sit=range.first;sit!=range.second;++sit)
		hitMapstoScan.push_back(sit->second);
	
	PBombLoop(ChiUnc);

	bestKeys.clear();
	bestKeys.resize(nKey);
	for (size_t bp=0;bp<nKey;++bp)
		bestKeys[bp]=sortedPBombs[bp]->ReadScarceMapKey();		

#if VERBOSITY>0
	cout<<"max hit at Key"<<bestKeys.front();
#endif
	
}
//__________________________________________________________________________
void Pulse::SelectPBombsFromKey(HitLib* Hitmap)
{
	//hitMapstoScan.clear();
	vector<PBomb*> tempSearch;
	// based on the key, collect the nearby PBombs
	for(size_t key=0;key<bestKeys.size();++key){
		tempSearch.clear();
		tempSearch = Hitmap->SearchKey(bestKeys[key]);
		hitMapstoScan.insert(hitMapstoScan.end(),tempSearch.begin(),tempSearch.end());
	}
}
//__________________________________________________________________________
void Pulse::PBombLoop(vector<float>&ChiUnc)
{
	sortedPBombs.clear();
	sortedPBombs=hitMapstoScan;

	if(statMethod==0){
		for (size_t pit=0;pit<hitMapstoScan.size();++pit)
			ChiSq(hitMapstoScan[pit],ChiUnc);
		sort(sortedPBombs.begin(),sortedPBombs.end(),sort_by_ChiSq());		
	}
	if(statMethod==1){
		for (size_t pit=0;pit<hitMapstoScan.size();++pit)
			Likelihood(hitMapstoScan[pit]);
		sort(sortedPBombs.begin(),sortedPBombs.end(),sort_by_NLL());
		minNLL = sortedPBombs[0]->ReadNLL();		
	}
}
//__________________________________________________________________________
double Pulse::ChiSq(PBomb* pit, vector<float>&ChiUnc)
{
	double diff, obs;
	pit->ResetStat();

	for (size_t pp=0; pp<liveTopPMT.size();++pp) {
			
		// filter the observed value
		obs = topHit[pp];
		if (topHit[pp]<=0) obs=0.1;
		
		diff = obs-top_area_sum*pit->ReadNrmPMT(pmtIdx[liveTopPMT[pp]]);
		pit->AddChiSq(diff*diff/obs/ChiUnc[liveTopPMT[pp]-1]);
	}
	return pit->ReadChiSq();
}
//__________________________________________________________________________
double Pulse::Likelihood(PBomb* pit)
{
	// Based on Poisson distribution, compute -log L,
	// and stores the value in each PBombs.

	float expected;
	pit->ResetStat();
	pit->AddNLL(top_area_sum+logGamma);
	for (size_t pp=0; pp<liveTopPMT.size();++pp){
		if (topHit[pp]<=0) continue;
		expected= top_area_sum*pit->ReadNrmPMT(pmtIdx[liveTopPMT[pp]]);
		if (expected==0) expected=0.01;
		pit->AddNLL(-topHit[pp]*log(expected));
	}
	return pit->ReadNLL();
}
//__________________________________________________________________________
void Pulse::InterpolateXY(size_t nItp)
{
	float p;
	float xpsum=0;
	float ypsum=0;
	float pValSum=0;
	
#if VERBOSITY>0
	cout<<", "<<hitMapstoScan.size()<<" points scanned";
#endif
	
	if (statMethod==0){
		for (size_t bnp=0;bnp<sortedPBombs.size() & bnp<nItp ;++bnp){
			p=chisqr(ndf,sortedPBombs[bnp]->ReadChiSq());
			xpsum+=sortedPBombs[bnp]->ReadX()*p;
			ypsum+=sortedPBombs[bnp]->ReadY()*p;
			pValSum+=p;
		}
		xFinal = xpsum/pValSum;
		yFinal = ypsum/pValSum;
		chisqP=chisqr(ndf,sortedPBombs[0]->ReadChiSq());
#if VERBOSITY>0
		cout<<", min Chi2/ndf is "<<sortedPBombs[0]->ReadChiSq()<<'/'<<ndf
			<<", P-value is "<<chisqP;
#endif
	}
	
	if (statMethod==1){
		sigmaRadius=LikelihoodConicalFit(nItp);
#if VERBOSITY>0
		cout<<", min -Log Likelihood is "<<sortedPBombs[0]->ReadNLL()
			<<", 1 sigma radius is "<<sigmaRadius<<"cm";
#endif
	}

#if VERBOSITY>0
	cout<<", X : "<<xFinal<<", Y : "<<yFinal;
	size_t scanEvent[] = {877,908,919,922,979,985,996,1000};
	size_t s2Idx;
	size_t scanIdx[]={1,1,1,1,1,1,1,1};

	/*string line;
	ifstream scanlist;
	scanlist.open("scanlist.txt");
	if (scanlist.is_open())
	{
		while ( getline(scanlist,line) )
		{
			cout << line << endl;
		}	
		scanlist.close();
	}*/
	
	if (find(scanEvent,scanEvent+8,event_num)!=scanEvent+8){
		s2Idx=find(scanEvent,scanEvent+8,event_num)-scanEvent;
		cout<<"***Event matched!***"<<s2Idx<<":"<<scanIdx[s2Idx]<<endl;
		if (pulse_num==scanIdx[s2Idx]) 
			PrintLikelihood(event_num);
	}

#endif
}
//__________________________________________________________________________
/*float Pulse::Resolution(){
	size_t i=0;
	float variance = 0;
	float threshold = 2.31; // 1 sigma for 2D
	
	// move farther out until the variance is 1 sigma away in 2d (2.31)
	while ((variance<threshold) & (i<=sortedPBombs.size())) {
		if (i==sortedPBombs.size()) return -999;
		variance=2*(sortedPBombs[i]->ReadNLL()-minNLL);
		++i;
	};
	
	// calculate the 1 sigma radius, and parabolically interpolate at the threshold
	float xdiff=sortedPBombs[i-1]->ReadX()-xFinal;
	float ydiff=sortedPBombs[i-1]->ReadY()-yFinal;
	float radSq=xdiff*xdiff+ydiff*ydiff;
	if (radSq>20*20) return -999;
	return sqrt(threshold/variance*(radSq));
}*/
//__________________________________________________________________________
double Pulse::ReadFinalX(){
	return xFinal;
}
//__________________________________________________________________________
double Pulse::ReadFinalY(){
	return yFinal;
}
//__________________________________________________________________________
double Pulse::ReadSigmaRad(){
	return sigmaRadius;
}
//__________________________________________________________________________
double Pulse::ReadChiSqP(){
	return chisqP;
}
//__________________________________________________________________________
void Pulse::Clear(){

	pulse_num=NAN;
	liveTopPMT.clear();
	
	top_area_sum=0;
	bestKeys.clear();
	hitMapstoScan.clear();
	sortedPBombs.clear();
}
//__________________________________________________________________________
double Pulse::LikelihoodConicalFit(size_t& nItp){
	if (nItp>=sortedPBombs.size()) return -999;

	const double radius = 23.485/cos(M_PI/12);
	const double threshold = 2.31; // 1 sigma for 2D
	
	double x,y,L;
	double xsum=0;	double ysum=0;
	double x2sum=0;	double xysum=0;		double y2sum=0;
	double x3sum=0;	double x2ysum=0;	double xy2sum=0;	double y3sum=0;
	double x4sum=0;	double x3ysum=0;	double x2y2sum=0;	double xy3sum=0;
	double y4sum=0;
	double Lsum=0;
	double xLsum=0;	double yLsum=0;
	double x2Lsum=0;double y2Lsum=0;	double xyLsum=0;

	for (size_t i=0; i<nItp;++i){
		x=sortedPBombs[i]->ReadX();
		y=sortedPBombs[i]->ReadY();
		L=sortedPBombs[i]->ReadNLL();

		xsum+=x;		ysum+=y;
		x2sum+=x*x;		y2sum+=y*y;	xysum+=x*y;
		x3sum+=x*x*x;	x2ysum+=x*x*y;	xy2sum+=x*y*y;	y3sum+=y*y*y;
		x4sum+=x*x*x*x;	x3ysum+=x*x*x*y;	x2y2sum+=x*x*y*y;
		xy3sum+=x*y*y*y;y4sum+=y*y*y*y;

		Lsum+=L;	xLsum+=x*L;	yLsum+=y*L;
		x2Lsum+=x*x*L;	y2Lsum+=y*y*L;	xyLsum+=x*y*L;		
	}

	size_t nParm=6;
	ublas::matrix<double> AL (nParm,nParm);
	ublas::vector<double> Lx(nParm);
	ublas::vector<double> Parms(nParm);	// vector to hold fit parmeters
	ublas::matrix<double> invAL (nParm,nParm);

	AL(0,0)=nItp;	AL(0,1)=xsum;	AL(0,2)=ysum;	
	AL(0,3)=x2sum;	AL(0,4)=y2sum;	AL(0,5)=xysum;
	AL(1,0)=xsum;	AL(1,1)=x2sum;	AL(1,2)=xysum;
	AL(1,3)=x3sum;	AL(1,4)=xy2sum;	AL(1,5)=x2ysum;
	AL(2,0)=ysum;	AL(2,1)=xysum;	AL(2,2)=y2sum;
	AL(2,3)=x2ysum; AL(2,4)=y3sum;	AL(2,5)=xy2sum;
	AL(3,0)=x2sum;	AL(3,1)=x3sum;	AL(3,2)=x2ysum;
	AL(3,3)=x4sum;	AL(3,4)=x2y2sum;AL(3,5)=x3ysum;
	AL(4,0)=y2sum;	AL(4,1)=xy2sum;	AL(4,2)=y3sum;
	AL(4,3)=x2y2sum;AL(4,4)=y4sum;	AL(4,5)=xy3sum;
	AL(5,0)=xysum;	AL(5,1)=x2ysum;	AL(5,2)=xy2sum;
	AL(5,3)=x3ysum;	AL(5,4)=xy3sum;	AL(5,5)=x2y2sum;

	Lx(0)=Lsum;		Lx(1)=xLsum;	Lx(2)=yLsum;
	Lx(3)=x2Lsum;	Lx(4)=y2Lsum;	Lx(5)=xyLsum;

	InvertMatrix(AL,invAL);
	Parms = prod(invAL,Lx);	

	double A=Parms(3);
	double B=Parms(5);
	double C=Parms(4);
	double D=Parms(1);
	double E=Parms(2);
	double F=Parms(0);
	
	double det=4*A*C-B*B;	// determinent of the A33 matrix
	xFinal = sortedPBombs[0]->ReadX();	// added for AC1.1
	yFinal = sortedPBombs[0]->ReadY();	// default xy if fit is not good
	if (det<=0) return -999;	// if not paraboloid, return bad fit
	
	// if paraboloid caculate the minimum
	xFinal=(B*E-2*C*D)/det;
	yFinal=(B*D-2*A*E)/det;
	
	//tangent of eigenvectors of A33 matrix
	double tan1=(C-A+sqrt(B*B+(A-C)*(A-C)))/B; 
	double tan2=(C-A-sqrt(B*B+(A-C)*(A-C)))/B;
	
	// calculate the radii of the confidence interval
	double quadA1=A+B*tan1+C*tan1*tan1;
	double quadB1=2*A*xFinal+2*C*tan1+D+E*tan1+B*yFinal+B*tan1*xFinal;
	double x1dist=(-fabs(quadB1)+sqrt(quadB1*quadB1+4*quadA1*threshold))
					/2/quadA1;
	double r1=fabs(x1dist)*sqrt(1+tan1*tan1);

	double quadA2=A+B*tan2+C*tan2*tan2;
	double quadB2=2*A*xFinal+2*C*tan2+D+E*tan2+B*yFinal+B*tan2*xFinal;
	double x2dist=(-fabs(quadB2)+sqrt(quadB2*quadB2+4*quadA2*threshold))
					/2/quadA2;
	double r2=fabs(x2dist)*sqrt(1+tan2*tan2);
	
	// project the radii vectors to the r unit vector
	ublas::vector<double> r1vec(2);
	r1vec(0)=r1/sqrt(1+tan1*tan1);
	r1vec(1)=r1vec(0)*tan1;
	ublas::vector<double> r2vec(2);
	r2vec(0)=r2/sqrt(1+tan2*tan2);
	r2vec(1)=r1vec(0)*tan2;
	ublas::vector<double> rUnitvec(2);
	rUnitvec(0)=xFinal/sqrt(xFinal*xFinal+yFinal*yFinal);
	rUnitvec(1)=yFinal/sqrt(xFinal*xFinal+yFinal*yFinal);

	// calculate the square sum of the projections and return
	double r1InnerProd=inner_prod(r1vec,rUnitvec);
	double r2InnerProd=inner_prod(r2vec,rUnitvec);

	// code below abolished since AC1.1
	/*// if the reconstructed position is outside the detector wall,
	double angle, minEdgeNLL,th,edgeX,edgeY,edgeNLL;
	double slice = M_PI/12.;
	double modangle = fmod((double)atan2(yFinal,xFinal),2*slice);
	if (modangle<0) modangle+=2*slice;

	double wallR = radius*cos(slice)/cos(modangle-slice);
	// if the reconstructed position is outside the detector wall, 
	if (xFinal*xFinal+yFinal*yFinal > wallR*wallR){
		angle = -M_PI/2+atan2(yFinal,xFinal);
		modangle = fmod(angle,2*slice);
		if (modangle<0) modangle+=2*slice;
		wallR = radius*cos(slice)/cos(modangle-slice);
		xFinal=wallR*cos(angle);
		yFinal=wallR*sin(angle);
		minEdgeNLL=A*xFinal*xFinal+B*xFinal*yFinal+C*yFinal*yFinal+D*xFinal+E*yFinal+F;
		// find the wall position with lowest NLL
		for (size_t t=0;t<1000;++t){
			th = t*M_PI/1000-M_PI/2+atan2(yFinal,xFinal);
			modangle = fmod(th,2*slice);
			if (modangle<0) modangle+=2*slice;
			wallR = radius*cos(slice)/cos(modangle-slice);
			edgeX = wallR*cos(th);
			edgeY = wallR*sin(th);
			edgeNLL=A*edgeX*edgeX+B*edgeX*edgeY+C*edgeY*edgeY+D*edgeX+E*edgeY+F;
			if (edgeNLL<minEdgeNLL){
				xFinal=edgeX;
				yFinal=edgeY;
				minEdgeNLL=edgeNLL;
			}
		}
	}*/	

	return sqrt(r1InnerProd*r1InnerProd+r2InnerProd*r2InnerProd);
}
//__________________________________________________________________________
void Pulse::PrintLikelihood(size_t event_num){
	ofstream LikelihoodLog;
	LikelihoodLog.open("LikelihoodMap.txt",ios::out | ios::app);
	LikelihoodLog<<"event/pulse"<<endl;
	LikelihoodLog<<event_num<<' '<<pulse_num<<endl;
	LikelihoodLog<<"xy"<<endl;
	LikelihoodLog<<xFinal<<' '<<yFinal<<endl;
	LikelihoodLog<<"PMT"<<endl;
	for (size_t pix=0;pix<liveTopPMT.size();++pix)
		LikelihoodLog<<liveTopPMT[pix]<<' '<<topHit[pix]<<endl;
	/*LikelihoodLog<<"Keys"<<endl;
	for (size_t pik=0;pik<bestKeys.size();++pik)
		LikelihoodLog<<bestKeys[pik]<<' '<<topHit[pik]<<endl;*/
	LikelihoodLog<<"Likelihood"<<endl;
	for (size_t pit=0;pit<sortedPBombs.size();++pit)
		LikelihoodLog<<sortedPBombs[pit]->ReadX()<<' '
					<<sortedPBombs[pit]->ReadY()<<' '
					<<sortedPBombs[pit]->ReadNLL()<<endl;
	LikelihoodLog<<endl;
	LikelihoodLog.close(); 
}

