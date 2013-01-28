#include "tmVisThreshC1_MixtureModeling2.h"



//----------------------   
tmVisThreshC1_MixtureModeling2::tmVisThreshC1_MixtureModeling2 (int w, int h) {
    INDEX_MIN 	= 1;
    INDEX_MAX 	= 253;
    MIN_VAL_8u 	= 0;
    MAX_VAL_8u 	= 255;
    
	cardinal 	= w*h;
    histogram 	= NULL;
	index 		= INDEX_MIN-1;

}

//----------------------
void tmVisThreshC1_MixtureModeling2::setHistogram (int *srcHistogram){
	histogram = srcHistogram;
}


//----------------------
bool tmVisThreshC1_MixtureModeling2::addToIndex() {
	index++;
	if(!(index<=INDEX_MAX)){
	    return false;	
	}

	setValues();
	return true;
}

//----------------------
float tmVisThreshC1_MixtureModeling2::calculateMax(int index) {
	float sum = histogram[index];
	float num = 1;
	if(index-1>=0) {
	    sum += histogram[index-1];
	    num++;
	}
	if(index+1<MAX_VAL_8u) {
	    sum += histogram[index+1];
	    num++;
	}
	return sum/num;
}


	

//----------------------
float tmVisThreshC1_MixtureModeling2::getCardinal() {
	return cardinal;
}

//----------------------
float tmVisThreshC1_MixtureModeling2::getMu1() {
	return mu1;
}

//----------------------
float tmVisThreshC1_MixtureModeling2::getMu2() {
	return mu2;
}

//----------------------
float tmVisThreshC1_MixtureModeling2::getMax1() {
	return max1;
}

//----------------------
float tmVisThreshC1_MixtureModeling2::getMax2() {
	return max2;
}

//----------------------
float tmVisThreshC1_MixtureModeling2::getVariance1() {
	return sigma2_1;
}

//----------------------
float tmVisThreshC1_MixtureModeling2::getVariance2() {
	return sigma2_2;
}

//----------------------
float tmVisThreshC1_MixtureModeling2::getCardinal1() {
	return cardinal1;
}

//----------------------
float tmVisThreshC1_MixtureModeling2::getCardinal2() {
	return cardinal2;
}

//----------------------
int tmVisThreshC1_MixtureModeling2::getThreshold() {
	return index;
}

//----------------------
void tmVisThreshC1_MixtureModeling2::setIndex(int i) {
	index = i;
	setValues();
}

//----------------------
void tmVisThreshC1_MixtureModeling2::setValues() {	
	mu1 = 0; 
	mu2 = 0;
	sigma2_1 = 0; 
	sigma2_2 = 0;
	max1 = 0; 
	max2 = 0;
	cardinal1 = 0; 
	cardinal2 = 0;

	for (int i=MIN_VAL_8u; i<=index ; i++) {
	    cardinal1 +=  histogram[i];
	    mu1 += i*histogram[i];
	}

	for (int i=index+1; i<=MAX_VAL_8u ; i++) {
	    cardinal2 +=  histogram[i];
	    mu2 += i*histogram[i];
	}

	if(cardinal1 == 0) {
	    mu1 = 0;
	    sigma2_1 = 0;
	} else {
   		mu1 /= (float)cardinal1; 
   	}

	if(cardinal2 == 0) {
	    mu2 = 0;
	    sigma2_2 = 0;
	} else {
    	mu2 /= (float)cardinal2; 
    }

	if( mu1 != 0 ) {
	    for(int i=MIN_VAL_8u; i<=index ; i++) {
			sigma2_1 += histogram[i]*pow(i-mu1,2);
		}
	    sigma2_1 /= (float)cardinal1;

	    max1 = calculateMax((int) mu1);
	    mult1 = (float) max1;
	    twoVariance1 = 2*sigma2_1; 
	}
	
	if( mu2 != 0 ) {
	    for(int i=index+1; i<=MAX_VAL_8u ; i++) {
			sigma2_2 += histogram[i]*pow(i-mu2,2);
		}
	    sigma2_2 /= (float)cardinal2;

	    max2 = calculateMax((int) mu2);
	    mult2 = (float) max2;
	    twoVariance2 = 2*sigma2_2; 
	}
}

//----------------------
float tmVisThreshC1_MixtureModeling2::gamma1(int i) {
	if (sigma2_1 == 0) 
    	return 0;
	return (float)(mult1*expf(-(powf((float)i-mu1,2.0f))/twoVariance1));
}

//----------------------
float tmVisThreshC1_MixtureModeling2::gamma2(int i) {
	if (sigma2_2 == 0) 
    	return 0;
	return (float)(mult2*expf(-(powf((float)i-mu2,2.0f))/twoVariance2));
}

//----------------------
float tmVisThreshC1_MixtureModeling2::gamma (int i) {
	return gamma1(i)+gamma2(i);
}

//----------------------
float tmVisThreshC1_MixtureModeling2::differenceGamma (int i) {
	return gamma1(i)-gamma2(i);
}

//----------------------
int	tmVisThreshC1_MixtureModeling2::getHistogram(int i) {
	return histogram[i];
}


