#include "PupilContourAnalyzer.h"

//==========================================================
PupilContourAnalyzer::PupilContourAnalyzer( ){
	
	slice_w = 64;
	slice_h = 1;
	slice_s = 64;
	ipp8u_C1_slice     = ippiMalloc_8u_C1 (slice_w, slice_h, &slice_s);
    ipp8u_C1_sliceCopy = ippiMalloc_8u_C1 (slice_w, slice_h, &slice_s);
	
	sliceTex.allocate (ofNextPow2(slice_w), ofNextPow2(slice_h), GL_LUMINANCE);
	
	dst_quad[0][0] = 0;
	dst_quad[0][1] = 0;
	dst_quad[1][0] = 0;
	dst_quad[1][1] = slice_h;
	dst_quad[2][0] = slice_w;
	dst_quad[2][1] = slice_h;
	dst_quad[3][0] = slice_w;
	dst_quad[3][1] = 0;

	handyVec3f0 = new ofVec3f();
	initialPupilEstimate = new ofVec3f();
	
	pupilGeometries.bValid = false;
    medianIrisRadius = 0;
}


//==========================================================
void PupilContourAnalyzer::fitEllipseToProvidedBlob (ofxCvBlob pupilBlob){

	pupilReport.set(0,0,0); 
	ellipseNorCent.set(0,0,0); 
	for (int i=0; i<N_ELLIPSE_SAMPLES; i++){
		ellipsePoints[i].set(0,0,0); 
		ellipsePtsNorm[i].set(0,0,0); 
	}
	
	//----------------------------------
	// fetch raw ellipse points from detected pupil ofxCvBlob
	
	int nPupilBlobPts = pupilBlob.nPts;
	if (nPupilBlobPts >= N_ELLIPSE_SAMPLES){
		// thin the incoming number of blob points to just N_ELLIPSE_SAMPLES
		float thinningFactorf = ((float)nPupilBlobPts / (float)N_ELLIPSE_SAMPLES);
	
		int count = 0; 
		for (float fi=0; fi<nPupilBlobPts; fi+=thinningFactorf){ 
			if (count < N_ELLIPSE_SAMPLES){
				int i = min((int)fi, nPupilBlobPts-1);
				handyVec3f0->x = pupilBlob.pts.at(i).x;
				handyVec3f0->y = pupilBlob.pts.at(i).y;
				ellipsePoints[count].set(*handyVec3f0);
				count++;
			}
		}
		
		//----------------------------------
		// normalize ellipse points
		float dis_scale = normalize_edge_points();
		pupil_fitting_inliers (dis_scale, 0.825);
	}
	
	
	
}

//==========================================================
ofVec3f PupilContourAnalyzer::getCircleReport(){
	
	// pupil_ellipse_param [5] is the parameters of an ellipse {ellipse_a, ellipse_b, cx, cy, theta}; 
	// a & b is the major or minor axis; 
	// cx & cy is the ellipse center; 
	// theta is the ellipse orientation
	
	pupilReport.set(0,0,0); 
	if (pupilGeometries.bValid){
	
		float *params_ellipse = pupilGeometries.params_ellipse;
		float axis_a = params_ellipse[0];
		float axis_b = params_ellipse[1];
		float cx 	 = params_ellipse[2];
		float cy 	 = params_ellipse[3];
		float theta	 = params_ellipse[4] * RAD_TO_DEG;
		float aspect = axis_b/axis_a;
		float radius = (axis_a + axis_b)/2.0; 
		
		pupilReport.x = cx;
		pupilReport.y = cy; 
		pupilReport.z = radius;  // encode radius in z, as a convenience.
	}
	
	return pupilReport;
}



//==========================================================
void PupilContourAnalyzer::findIrisEllipse (
                                    Ipp8u *srcImage,
                                    int srcw, int srch, int srcs,
                                    float initialCx,
                                    float initialCy,
                                    float minRadius,
                                    float maxRadius,
                                    float radialNonlinearity){

	providedIrisCx = initialCx;
    providedIrisCy = initialCy;
    
    pupilReport.set(0,0,0);
	ellipseNorCent.set(0,0,0);
	for (int i=0; i<N_ELLIPSE_SAMPLES; i++){
		ellipsePoints[i].set(0,0,0);
		ellipsePtsNorm[i].set(0,0,0);
	}
	
	//----------------------------------
	initialPupilEstimate->set(initialCx,initialCy,0);
    float search_Radius = maxRadius;
    
    for (int j=0; j< N_ELLIPSE_SAMPLES; j++){
        
        float jfrac = ((float)j/(float)N_ELLIPSE_SAMPLES);
		float jfracShaped = function_DoubleExponentialSeat (jfrac, radialNonlinearity);
		float jangle = (jfracShaped * PI - HALF_PI) * 0.75;
		if (j%2 == 1){ jangle += PI;}
		float search_Angle = jangle;
        
        ofVec3f ellipsePoint = getEllipsePoint (initialPupilEstimate, minRadius, maxRadius, search_Angle, srcImage, srcw, srch, srcs);
        handyVec3f0->x = ellipsePoint.x;
        handyVec3f0->y = ellipsePoint.y;
        ellipsePoints[j].set(*handyVec3f0);
        
        
        float dx = ellipsePoint.x - initialCx;
        float dy = ellipsePoint.y - initialCy;
        float dh = sqrtf(dx*dx + dy*dy);
        radiusArray[j] = dh;
    }
    
    
	
    
    //----------------------------------
    // normalize ellipse points
    bool bDoIrisEllipseFitting = true; //not actually used.
    if (bDoIrisEllipseFitting){
        float dis_scale = normalize_edge_points();
        pupil_fitting_inliers (dis_scale, 0.90);
    }
    
    
    // bogus it up a bit with heuristics
    //float *params_ellipse = pupilGeometries.params_ellipse;
    //float axis_a = params_ellipse[0];
    //float axis_b = params_ellipse[1];
    float medRad = computeMedianRadius();
    medianIrisRadius = medRad;
    

}

float PupilContourAnalyzer::getMedianRadius(){
    return medianIrisRadius;
}

//==========================================================
float PupilContourAnalyzer::computeMedianRadius() {
    
 
    for (int i = 0; i < N_ELLIPSE_SAMPLES; ++i) {
        radiusArrayCopy[i] = radiusArray[i];
    }
    for (int i = N_ELLIPSE_SAMPLES - 1; i > 0; --i) {
        for (int j = 0; j < i; ++j) {
            if (radiusArrayCopy[j] > radiusArrayCopy[j+1]) {
                float dTemp = radiusArrayCopy[j];
                radiusArrayCopy[j] = radiusArrayCopy[j+1];
                radiusArrayCopy[j+1] = dTemp;
            }
        }
    }
    
    // Middle or average of middle values in the sorted array.
    float dMedian = 0.0;
    if ((N_ELLIPSE_SAMPLES % 2) == 0) {
        dMedian = (radiusArrayCopy[N_ELLIPSE_SAMPLES/2] + radiusArrayCopy[(N_ELLIPSE_SAMPLES/2) - 1])/2.0;
    } else {
        dMedian = radiusArrayCopy[N_ELLIPSE_SAMPLES/2];
    }
    
    /*
    // or how about a winsorized mean instead. nope, it makes it worse. 
    // eliminate bottom 40 and top 40 percent.
    int winsor = 0.40;
    int quarterW0 = (int)(N_ELLIPSE_SAMPLES * winsor);
    int quarterW1 = (int)(N_ELLIPSE_SAMPLES * (1.0 - winsor));
    int count = 0;
    float sum = 0; 
    for (int i=quarterW0; i<quarterW1; i++){
        sum += radiusArrayCopy[i];
        count++;
    }
    dMedian = sum / (float)count; 
     */
    
    
    // or how about a heuristic in which we take an average
    // of those radii which are within some range of tolerance of the median? ... not bad.
    float lowR  = dMedian * 0.90;
    float highR = dMedian * 1.10;
    float sum = 0;
    int count = 0;
    for (int i=0; i<N_ELLIPSE_SAMPLES; i++){
        if ((radiusArrayCopy[i] > lowR) && (radiusArrayCopy[i] < highR)){
            sum += radiusArrayCopy[i];
            count ++; 
        }
    }
    dMedian = sum / (float)count;
    
    
    return dMedian;
}



//==========================================================
void  PupilContourAnalyzer::pupil_fitting_inliers (float dis_scale, float minAspectRatio){

	long then = ofGetElapsedTimeMillis();
	int nEllipsePts = N_ELLIPSE_SAMPLES;
		
	
	double dis_threshold = 0.1 *dis_scale;
	memset(inliers_index, int(0),     sizeof(int)*nEllipsePts);
	memset(max_inliers_index, int(0), sizeof(int)*nEllipsePts);
	
	int rand_index[5];
	double A[6][6];
	int M = 6; 
	int N = 6; //M is row; N is column
	for (int j=0; j<N; j++) {
		A[j][5] = 1;
		A[5][j] = 0;
	}

	//-----------------------
	// MALLOCs (matched by FREEs below)
	double **ppa = (double**)malloc(sizeof(double*)*M);
	double **ppu = (double**)malloc(sizeof(double*)*M);
	double **ppv = (double**)malloc(sizeof(double*)*N);
	for (int i=0; i<M; i++) {
		ppa[i] = A[i];
		ppu[i] = (double*)malloc(sizeof(double)*N);
	}
	for (int i=0; i<N; i++) {
		ppv[i] = (double*)malloc(sizeof(double)*N);
	}
	double 	pd[6]; 
	int		min_d_index;
	double 	conic_par[6] = {0};
	double 	ellipse_par[5] = {0};
	double 	best_conic_par[6] = {0};
	double 	best_ellipse_par[5] = {0};
	double 	ratio;
	
	
	//-----------------------
	int max_inliers = 0;
	int ransac_count = 0;
	int best_guess_count = 0;
	int sample_num = 1600;
	int ninliers = 0;
	double dis_error;
	double cx,cy;
	double theta;
	double axis_a, axis_b;
	
	double minAspect = minAspectRatio; //0.825 for pupils, 0.9... for irises..
	double leastPreviousRatio = 0;
	double leastPreviousError = 99999;
	
	while (sample_num > ransac_count) {
		get_5_random_num((nEllipsePts-1), rand_index);

		//svd decomposition to solve the ellipse parameter
		for (int i=0; i<5; i++) {
			int 		randomIndex = rand_index[i];
			ofVec3f	randomEllipsePt = ellipsePtsNorm[randomIndex];
			A[i][0] = randomEllipsePt.x * randomEllipsePt.x;
			A[i][1] = randomEllipsePt.x * randomEllipsePt.y;
			A[i][2] = randomEllipsePt.y * randomEllipsePt.y;
			A[i][3] = randomEllipsePt.x;
			A[i][4] = randomEllipsePt.y;
		}

		svd (M, N, ppa, ppu, pd, ppv);
		min_d_index = 0;
		for (int i=1; i<N; i++) {
			if (pd[i] < pd[min_d_index]){
				min_d_index = i;
			}
		}

		// the column of v that corresponds to the smallest singular value, 
        // which is the solution of the equations
		for (int i=0; i<N; i++){
			conic_par[i] = ppv[i][min_d_index];
		}	
      
    	ninliers = 0;
    	memset(inliers_index, 0, sizeof(int)*nEllipsePts);
    
		for (int i=0; i<nEllipsePts; i++) {
			ofVec3f ithEllipsePt = ellipsePtsNorm[i];
			dis_error =     conic_par[0] * ithEllipsePt.x * ithEllipsePt.x + 
							conic_par[1] * ithEllipsePt.x * ithEllipsePt.y +
							conic_par[2] * ithEllipsePt.y * ithEllipsePt.y + 
							conic_par[3] * ithEllipsePt.x + 
							conic_par[4] * ithEllipsePt.y + 
							conic_par[5];
			if (fabs(dis_error) < dis_threshold) {
				inliers_index[ninliers] = i;
				ninliers++;
			}
		}
    
    	//--------------------
    	// filtering of the ellipses
		if (ninliers > max_inliers) {
			if (solve_ellipse(conic_par, ellipse_par)) {
				double norm_cx = (double) ellipseNorCent.x;
				double norm_cy = (double) ellipseNorCent.y;
				denormalize_ellipse_param (ellipse_par, ellipse_par, dis_scale, norm_cx, norm_cy);

				axis_a = ellipse_par[0];
				axis_b = ellipse_par[1];
				cx 	   = ellipse_par[2];
				cy 	   = ellipse_par[3];
				ratio  = axis_a / axis_b;
	
				if ( (cx > 0) && (cx < 10000) && 
					 (cy > 0) && (cy < 10000) &&
					 (ratio > minAspect) && (ratio < (1.0/minAspect))) {
				
					double currentRatio = (ratio > 1) ? (1.0/ratio) : ratio;
					double currentError = fabs(dis_error);
					if ((currentError < leastPreviousError) && (currentRatio > leastPreviousRatio)) {
					//if (true){ 
						leastPreviousRatio = currentRatio;
						leastPreviousError = currentError;
					
						memcpy(max_inliers_index, inliers_index, sizeof(int)*nEllipsePts);
						for (int i=0; i<5; i++) {
							best_ellipse_par[i] = ellipse_par[i];
						}
						for (int i=0; i<6; i++){
							best_conic_par[i] = conic_par[i];
						}
						max_inliers = ninliers;
						sample_num = (int)(log((double)(1-0.99))/log(1.0-pow(ninliers/(float)nEllipsePts, 5)));
						best_guess_count++;
					}
					
				}
			}
		}
		ransac_count++;
		if (ransac_count > 2400) { // about 10 milliseconds!
			//printf("Error! ransac_count exceed! ransac break! sample_num=%d, ransac_count=%d\n", sample_num, ransac_count);
			break;
		}
		
	}
	
	// end of Ransac. Fill the pupilGeometries.
	if (best_ellipse_par[0] > 0 && best_ellipse_par[1] > 0) {
		pupilGeometries.bValid = true;
		for (int i=0; i<5; i++) {
			pupilGeometries.params_ellipse[i] = best_ellipse_par[i];
		}
		for (int i=0; i<6; i++){
			pupilGeometries.params_conic[i] = best_conic_par[i];
		}
		
	} else {
		pupilGeometries.bValid = false;
		max_inliers = 0;
	}

	
	//-----------------------
	// FREEs
	for (int i=0; i<M; i++) {
		free(ppu[i]);
		free(ppv[i]);
	}
	free(ppu);
	free(ppv);
	free(ppa);
	

	//-----------------------
	long now = ofGetElapsedTimeMillis();
	// printf("---	count %d,	gus %d,	took: %d\n", ransac_count, best_guess_count, (int)(now-then)); 

}



//==========================================================
// solve_ellipse
// conic_param[6] is the parameters of a conic {a, b, c, d, e, f}; conic equation: ax^2 + bxy + cy^2 + dx + ey + f = 0;
// ellipse_param[5] is the parameters of an ellipse {ellipse_a, ellipse_b, cx, cy, theta}; a & b is the major or minor axis; 
// cx & cy is the ellipse center; theta is the ellipse orientation
bool PupilContourAnalyzer::solve_ellipse (double* conic_param, double* ellipse_param){

  double a = conic_param[0];
  double b = conic_param[1];
  double c = conic_param[2];
  double d = conic_param[3];
  double e = conic_param[4];
  double f = conic_param[5];
  
  //get ellipse orientation
  double theta = atan2(b, a-c)/2;

  //get scaled major/minor axes
  double ct = cos(theta);
  double st = sin(theta);
  double ap = a*ct*ct + b*ct*st + c*st*st;
  double cp = a*st*st - b*ct*st + c*ct*ct;

  //get translations
  double cx = (2*c*d - b*e) / (b*b - 4*a*c);
  double cy = (2*a*e - b*d) / (b*b - 4*a*c);

  //get scale factor
  double val = a*cx*cx + b*cx*cy + c*cy*cy;
  double scale_inv = val - f;

  if (scale_inv/ap <= 0 || scale_inv/cp <= 0) {
    //printf("Error! ellipse parameters are imaginary a=sqrt(%lf), b=sqrt(%lf)\n", scale_inv/ap, scale_inv/cp);
    memset(ellipse_param, 0, sizeof(double)*5);
    return false;
  }

  ellipse_param[0] = sqrt(scale_inv / ap);
  ellipse_param[1] = sqrt(scale_inv / cp);
  ellipse_param[2] = cx;
  ellipse_param[3] = cy;
  ellipse_param[4] = theta;
  return true;
}

//==========================================================
// Randomly select 5 indices
void PupilContourAnalyzer::get_5_random_num (int max_num, int* rand_num) {
  int rand_index = 0;
  int r;
  int i;
  bool is_new = true;

  if (max_num == 4) {
    for (i = 0; i < 5; i++) {
      rand_num[i] = i;
    }
    return;
  }

  while (rand_index < 5) {
    is_new = true;
    r = (int)(ofRandomuf() * (float)max_num);
    for (i=0; i<rand_index; i++) {
      if (r == rand_num[i]) {
        is_new = false;
        break;
      }
    }
    if (is_new) {
      rand_num[rand_index] = r;
      rand_index++;
    }
  }
}


//==========================================================
float PupilContourAnalyzer::normalize_edge_points (){

	float sumx = 0; 
	float sumy = 0;
	float sumdis = 0;
	float dis_scale = 0;

	int nEllipsePts = N_ELLIPSE_SAMPLES; 
	for (int j=0; j<nEllipsePts; j++) {
		float ex = ellipsePoints[j].x;
		float ey = ellipsePoints[j].y;
		sumx += ex;
		sumy += ey;
		sumdis += sqrtf((ex*ex + ey*ey));
	}

	if (sumdis > 0){
		dis_scale =   sqrt(2.0)*(float)nEllipsePts/sumdis;
		ellipseNorCent.x = sumx/(float)nEllipsePts;
		ellipseNorCent.y = sumy/(float)nEllipsePts;

		for (int j=0; j<nEllipsePts; j++) {

			float ex = ellipsePoints[j].x;
			float ey = ellipsePoints[j].y;
			ofVec3f normedPt; 
			normedPt.x = (ex - ellipseNorCent.x)*dis_scale;
			normedPt.y = (ey - ellipseNorCent.y)*dis_scale;

			ellipsePtsNorm[j].x = (ex - ellipseNorCent.x)*dis_scale;
			ellipsePtsNorm[j].y = (ey - ellipseNorCent.y)*dis_scale;
		}
	} else {
		printf("Normalization failure!\n"); 
	}
	
	return dis_scale;
}



//==========================================================
ofVec3f	PupilContourAnalyzer::getEllipsePoint (
		ofVec3f *oldCentroid,  float min_R, float max_R, float search_A,
		Ipp8u *srcImage, int src_w, int src_h, int src_s){
	
	// from a raster image, search through for an edge 
	// that represents a point on the edge of an ellipse.
	ofVec3f outputVec3f;
	
	// compute remote location
	float cx = oldCentroid->x;
	float cy = oldCentroid->y;
	
	// safety constraint
	float limR = ceil(max_R);
	cx = min( src_w-1-limR, max(limR, cx)); 
	cy = min( src_h-1-limR, max(limR, cy)); 
	
	// compute endpoints
	float cosA = cos(search_A);
	float sinA = sin(search_A);
	float rx = cx + (max_R*cosA);
	float ry = cy + (max_R*sinA);
	float th = 0.125;
	
	// construct source quad
	src_quad[0][0] = cx - th*sinA;	src_quad[0][1] = cy + th*cosA;
	src_quad[1][0] = cx + th*sinA;	src_quad[1][1] = cy - th*cosA;
	src_quad[2][0] = rx + th*sinA;	src_quad[2][1] = ry - th*cosA;
	src_quad[3][0] = rx - th*sinA;	src_quad[3][1] = ry + th*cosA;

	// fetch interpolated slice of pixels
	IppiSize srcSize = {src_w, src_h};
	IppiRect srcRect = {0,0, src_w, src_h};
	IppiRect sliceRect = {0,0, slice_w, slice_h};
	int interpolation = IPPI_INTER_CUBIC;
	ippiWarpBilinearQuad_8u_C1R(
                                srcImage, srcSize, src_s, srcRect, src_quad, 
                                ipp8u_C1_slice, slice_s, sliceRect, dst_quad, 
                                interpolation);

    // blur it!
    int nBlurRepeats = 2;
    for (int r=0; r<nBlurRepeats; r++){
        for (int x=1; x<(slice_w-1); x++){
            int valA = ipp8u_C1_slice[x-1];
            int valB = ipp8u_C1_slice[x  ];
            int valC = ipp8u_C1_slice[x+1];
            ipp8u_C1_sliceCopy[x] = (valA+ (2*valB) +valC)/4;
        }
        for (int x=0; x<slice_w; x++){
            ipp8u_C1_slice[x] = ipp8u_C1_sliceCopy[x];
        }
    }
    
    
	// find brightest pixel in slice; 
	// convert its location to a coordinate in the eye image.
	IppiSize sliceRoi = {slice_w, slice_h};
	float sliceFrac;
	float point_R = 0;
	float point_X = cx;
	float point_Y = cy;
	
	const int searchMethod_forMaxBrightnessPoint 		= 0;
	const int searchMethod_forFirstMajorDrop			= 1;
	const int searchMethod_forFirstMajorStepUp			= 2;

	int searchMethod = searchMethod_forFirstMajorStepUp;

	Ipp8u pMax = 0;
	int pMaxIndexX = 0;
	int pMaxIndexY = 0;

    
    int compareStep = 3;
    int searchStartIndex =  (int) ofMap(min_R, 0,max_R, 1,slice_w-1);
    searchStartIndex = max(compareStep, searchStartIndex);
    
    
	int indexOfLargestDrop = 0;
	int largestDrop = 0;
    
    
	switch (searchMethod){
	
		case searchMethod_forMaxBrightnessPoint:
			ippiMaxIndx_8u_C1R (ipp8u_C1_slice, slice_s, sliceRoi, &pMax, &pMaxIndexX, &pMaxIndexY);
			sliceFrac = (float)pMaxIndexX/(float)slice_w;
			break;
			
		case searchMethod_forFirstMajorDrop: 
			for (int i=1; i<slice_w; i++){
				Ipp8u v0 = ipp8u_C1_slice[i-1];
				Ipp8u v1 = ipp8u_C1_slice[i];
				if ((v0 - v1) > largestDrop){
					indexOfLargestDrop = i;
					largestDrop = (v0 - v1);
				}
			}
			sliceFrac = (float)indexOfLargestDrop/(float)slice_w;
			break;	
		
		case searchMethod_forFirstMajorStepUp: 
			int indexOfLargestStepUp = 0;
			int largestStepUp = 0;
			for (int i=searchStartIndex; i<slice_w; i++){
				Ipp8u v0 = ipp8u_C1_slice[i-compareStep];
				Ipp8u v1 = ipp8u_C1_slice[i];
				if ((v1 - v0) > largestStepUp){
					indexOfLargestStepUp = i;
					largestStepUp = (v1 - v0);
				}
			}
			sliceFrac = (float)(indexOfLargestStepUp - compareStep)/(float)slice_w;
			break;	
			
	}					
	
	point_R = sliceFrac * max_R;
	point_X = cx + (point_R * cosA) + 0.5;
	point_Y = cy + (point_R * sinA) + 0.5;		
	outputVec3f.set(point_X, point_Y, 0); 	
	return outputVec3f;		
}

//==========================================================
void PupilContourAnalyzer::drawEllipse(){
	
	// pupil_ellipse_param [5] is the parameters of an ellipse {ellipse_a, ellipse_b, cx, cy, theta}; a & b is the major or minor axis; 
	// cx & cy is the ellipse center; theta is the ellipse orientation
	
	if (pupilGeometries.bValid){
	
		float *params_ellipse = pupilGeometries.params_ellipse;
		float axis_a = params_ellipse[0];
		float axis_b = params_ellipse[1];
		float cx 	 = params_ellipse[2];
		float cy 	 = params_ellipse[3];
		float theta	 = params_ellipse[4] * RAD_TO_DEG;
		float aspect = axis_b/axis_a;

		int resolution = 36;
		glPushMatrix();
		glTranslatef(cx,cy,0); 
		glRotatef(theta, 0,0,1);
		
		glColor3f(0,0.2,1.0);
		glBegin(GL_LINE_LOOP);
		for (int i=0; i<resolution; i++){
			float t = TWO_PI * (float)i/(float)resolution;
			float ex = (axis_a * cos(t));
			float ey = (axis_b * sin(t)); 
			glVertex2f (ex,ey); 
		}
		
		glEnd();
		glBegin(GL_LINES);
		glVertex2f(0,0); 
		glVertex2f(axis_a,0);
		glEnd();
		
		glPointSize(3.0);
		glBegin(GL_POINTS);
		glVertex2f(0,0);
		glEnd();
		
		glPopMatrix();
	}
}


//==========================================================
void PupilContourAnalyzer::draw(float px, float py){
	
	glPushMatrix();
	glTranslatef(px,py,0);
	
	bool bDrawEllipsePoints   = true;
	bool bDrawEllipseSolution = true;
	
	// draw detected ellipse points
	if (bDrawEllipsePoints){
		glPointSize(1.0);
        ofEnableAlphaBlending();
		glColor4f(0,0.8,0,  0.5);
        
		
        glPointSize(4.0);
		glBegin(GL_POINTS);
		int nPts = N_ELLIPSE_SAMPLES;
		for (int j=0; j<nPts; j++){
			float ex = ellipsePoints[j].x;
			float ey = ellipsePoints[j].y;
			glVertex2f(ex,ey); 
		}
		glEnd();
        glPointSize(1.0);
	}
	
	// draw computed ellipse geometry
	if (bDrawEllipseSolution){
		drawEllipse();
	}
    
    ofNoFill();
    ofSetColor(0,255,100);
    ofEllipse(providedIrisCx, providedIrisCy, medianIrisRadius*2.0, medianIrisRadius*2.0);

	glPopMatrix();
}

//==========================================================
void PupilContourAnalyzer::drawFromPixelDetectionMethod (float px, float py){
	
	
	glPushMatrix();
	glTranslatef(px,py,0); ; 
	
	bool bDrawEllipseSolution = true;
	bool bDrawEllipsePoints = true;
	bool bDrawSlicePixels = false;
	bool bDrawOriginalCentroid = true;
	bool bDrawSliceSelectionBoundary = false;
		
	// draw centroid of original contour
	if (bDrawOriginalCentroid){
		glColor3f(0,0,1);
		float cx = initialPupilEstimate->x;
		float cy = initialPupilEstimate->y;
		glBegin(GL_POINTS);
		glVertex2f(cx,cy); 
		glEnd();
	}
	
	// draw detected ellipse points
	if (bDrawEllipsePoints){
		glPointSize(2.0); 
		glColor3f(1,0,0); 
		
		glBegin(GL_POINTS);
		for (int j=0; j<N_ELLIPSE_SAMPLES; j++){
			float ex = ellipsePoints[j].x;
			float ey = ellipsePoints[j].y;
			glVertex2f(ex,ey); 
		}
		glEnd();
	}
	
	// draw computed ellipse geometry
	if (bDrawEllipseSolution){
		drawEllipse();
	}
	
	// draw slice selection boundary
	if (bDrawSliceSelectionBoundary){
		glColor3f(1,0,0);
		glBegin(GL_LINE_LOOP);
		for (int j=0; j<4; j++){
			glVertex2f(src_quad[j][0], src_quad[j][1]);
		}
		glEnd();
	}
	
	// draw actual slice
	if (bDrawSlicePixels){
		glColor3f(1,1,1);
		sliceTex.loadData(ipp8u_C1_slice, slice_w, slice_h, GL_LUMINANCE);
		sliceTex.draw(0, 300, slice_w, slice_h);
	}
		

	
	glPopMatrix();
}



//==========================================================
void PupilContourAnalyzer::denormalize_ellipse_param(
	double* par, double* normalized_par, double dis_scale, double norm_cx, double norm_cy){
	par[0] = normalized_par[0] / dis_scale;				//major or minor axis
	par[1] = normalized_par[1] / dis_scale;
	par[2] = normalized_par[2] / dis_scale + norm_cx;	//ellipse center
	par[3] = normalized_par[3] / dis_scale + norm_cy;
}

//------------------------------------------------------------------
float PupilContourAnalyzer::function_DoubleExponentialSeat (float x, float a){

  float epsilon = 0.00001;
  float min_param_a = 0.0 + epsilon;
  float max_param_a = 1.0 - epsilon;
  a = min(max_param_a, max(min_param_a, a));

  float y = 0;
  if (x<=0.5){
    y = (powf(2.0*x, 1-a))/2.0;
  } 
  else {
    y = 1.0 - (powf(2.0*(1.0-x), 1-a))/2.0;
  }
  return y;
}

//------------------------------------------------------------------
float PupilContourAnalyzer::function_DoubleExponentialSigmoid (float x, float a){

  float epsilon = 0.00001;
  float min_param_a = 0.0 + epsilon;
  float max_param_a = 1.0 - epsilon;
  a = min(max_param_a, max(min_param_a, a)); 
  a = 1-a;
  
  float y = 0;
  if (x<=0.5){
    y = (pow(2.0*x, 1.0/a))/2.0;
  } 
  else {
    y = 1.0 - (pow(2.0*(1.0-x), 1.0/a))/2.0;
  }
  return y;
}


