/*****************************************************************************
 * Project: ProfileLikelihood                                                *
 *                                                                           *
 * Markus Horn - 11 July 2013                                                *
 *****************************************************************************/

#ifndef ROOS1RESPONSE
#define ROOS1RESPONSE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooS1Response : public RooAbsPdf {
    
public:
    RooS1Response() {};
    RooS1Response(const char *name, const char *title,
                  RooAbsReal& _S1,
                  RooAbsReal& _nPhot,
                  RooAbsReal& _r,
                  RooAbsReal& _z);
    RooS1Response(const RooS1Response& other, const char *name=0); // copy constructor?
    virtual TObject* clone(const char *newname) const { return new RooS1Response(*this, newname); }
    inline virtual ~RooS1Response() { }; // destructor
    
protected:

    RooRealProxy S1;
    RooRealProxy nPhot ;
    RooRealProxy r;
    RooRealProxy z;
 
    Double_t evaluate() const ;
    
private:
    
    ClassDef(RooS1Response,1); // to integrate class to ROOT
    double LightCollectionEfficiency(double _r, double _z) const;
    
};

#endif