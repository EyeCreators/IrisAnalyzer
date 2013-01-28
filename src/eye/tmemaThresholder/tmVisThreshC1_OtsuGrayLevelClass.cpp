#include "tmVisThreshC1_OtsuGrayLevelClass.h"



//==============================  
tmVisThreshC1_OtsuGrayLevelClass::tmVisThreshC1_OtsuGrayLevelClass (int nPixels, bool f) {
	
	mu 		= 0;
	omega	= 0;
	first	= f;
	N 		= nPixels;
	
	probabilityHistogram = new float[256];
	for (int i=0; i<256; i++){
		probabilityHistogram[i] = 0;
	}
	
}

//============================== 
void tmVisThreshC1_OtsuGrayLevelClass::initialize (int* srcHistogram) {
	
    for (int i=0; i<256; i++) {
		probabilityHistogram[i] = ((float) srcHistogram[i])/((float) N);
    }

	if (first) {
	    index = 1;
	    omega = probabilityHistogram[index-1];
	    if (omega == 0){
			mu = 0;
		} else {
			mu =  1*probabilityHistogram[index-1]/omega;
		}
	}
	else {
	    index = 2;
	    omega = 0;
	    mu = 0;
	    for (int i=index; i<256; i++) {
			omega +=  probabilityHistogram[i-1];
			mu    +=  probabilityHistogram[i-1]*i;
	    }
	    if(omega == 0){
			mu = 0;
		} else {
			mu /= omega;
		}
	}
}
 
 
//==============================   
void tmVisThreshC1_OtsuGrayLevelClass::removeFromBeginning() {
	index++;
	mu    = 0;
	omega = 0;

	for (int i=index; i<256 ; i++) {
	    omega +=   probabilityHistogram[i-1];
	    mu    += i*probabilityHistogram[i-1];
	}
	if(omega == 0){
	    mu = 0;
	} else {
	    mu /= omega;
    }
}

//==============================  
void tmVisThreshC1_OtsuGrayLevelClass::addToEnd() {
	index++;
	mu 	  = 0;
	omega = 0;
	for (int i=1; i<index ; i++) {
	    omega 	+=   probabilityHistogram[i-1];
	    mu 		+= i*probabilityHistogram[i-1];
	}
	if (omega == 0){
	    mu = 0;
	} else {
	    mu /= omega;
	}
}

//==============================  
float tmVisThreshC1_OtsuGrayLevelClass::getMu() 		{ return mu; }
float tmVisThreshC1_OtsuGrayLevelClass::getOmega() 		{ return omega; }
int   tmVisThreshC1_OtsuGrayLevelClass::getThreshold() 	{ return index; }
