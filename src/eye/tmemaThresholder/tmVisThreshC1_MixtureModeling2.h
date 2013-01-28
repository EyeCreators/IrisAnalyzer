#ifndef _TM_VIS_THRESHOLDERC1_MMOD2
#define _TM_VIS_THRESHOLDERC1_MMOD2

#include "ofMain.h"

class tmVisThreshC1_MixtureModeling2 {
	

  public:
	

	int INDEX_MIN;
	int INDEX_MAX;
	int MIN_VAL_8u;
	int MAX_VAL_8u;

	int 	index;
	float 	mu1, mu2;
	float 	sigma2_1, sigma2_2;
	float 	mult1, mult2;
	float 	twoVariance1, twoVariance2;
	float 	max1, max2;
	int 	cardinal1, cardinal2;
	int 	cardinal;
	
	int		*histogram;


	tmVisThreshC1_MixtureModeling2 (int w, int h);
	void 	setHistogram (int *srcHistogram);
	bool 	addToIndex();
	float 	calculateMax (int index);

	float 	getCardinal();
	float 	getMu1();
	float 	getMu2();
	float 	getMax1();
	float 	getMax2();
	float 	getVariance1();
	float 	getVariance2();
	float 	getCardinal1();
	float 	getCardinal2();
	int 	getThreshold();
	void 	setIndex(int i);
	void 	setValues();
	float 	gamma1(int i);
	float 	gamma2(int i);
	float 	gamma (int i);
	float 	differenceGamma (int i);
	int 	getHistogram(int i);

};

#endif /* _TM_VIS_THRESHOLDERC1_MMOD2 */
