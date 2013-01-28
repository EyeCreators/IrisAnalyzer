// DarkPupilContourAnalyzer.h
#ifndef _PUPIL_CONTOUR_ANALYZER
#define _PUPIL_CONTOUR_ANALYZER

#include "ofMain.h"
#include "ipp.h"
#include "cv.h"
#include "PupilStructs.h"
#include "ofxOpenCv.h" // for ofxCvBlob

extern "C" {
	#include "svd.h"
}


#define N_ELLIPSE_SAMPLES	128

//================================================
class PupilContourAnalyzer {

public:
	PupilContourAnalyzer ();

	void	fitEllipseToProvidedBlob (ofxCvBlob pupilBlob);
    
    void    findIrisEllipse (Ipp8u *eyeImageC1,
                            int srcw, int srch, int srcs,
                            float initialCx,
                            float initialCy,
                            float minRadius,
                            float maxRadius,
                            float radialNonlinearity);
    
    float   getMedianRadius();
	
	void	drawFromPixelDetectionMethod (float px, float py); 
	void	draw (float px, float py);
	
	PupilGeometry	pupilGeometries;
	ofVec3f getCircleReport();
    
    
	
private:
    
    float providedIrisCx;
    float providedIrisCy;
    float medianIrisRadius;
    float computeMedianRadius();
    float radiusArray[N_ELLIPSE_SAMPLES];
    float radiusArrayCopy[N_ELLIPSE_SAMPLES];
    

	ofTexture	sliceTex;
	
	Ipp8u	*ipp8u_C1_slice;
    Ipp8u	*ipp8u_C1_sliceCopy;
	int slice_w;
	int slice_h;
	int slice_s;
	
	int SRC_W;
	int SRC_H;
	int SRC_S;
	
	ofVec3f					ellipsePoints	[N_ELLIPSE_SAMPLES]; // raw ellipse points
	ofVec3f					ellipsePtsNorm	[N_ELLIPSE_SAMPLES]; // "normalized" ellipse points
	ofVec3f					ellipseNorCent;
	ofVec3f					pupilReport;
    ofVec3f					irisReport;
	
	// for Ransac
	int						inliers_index		[N_ELLIPSE_SAMPLES];
	int						max_inliers_index	[N_ELLIPSE_SAMPLES];
	
	double 					pupil_ellipse_param [1][5];
	double 					pupil_conic_param   [1][6];
	
	ofVec3f	getEllipsePoint (
							ofVec3f *oldCentroid, float min_R, float max_R, float search_A, 
							Ipp8u *srcImage, int src_w, int src_h, int src_s);
	
	void	pupil_fitting_inliers (float dis_scale, float minAspectRatio);
	float	normalize_edge_points ();
	void 	get_5_random_num (int max_num, int* rand_num);
	bool 	solve_ellipse (double* conic_param, double* ellipse_param);
	void 	drawEllipse();
	
	float	getBaseAngleBetweenPupils();
	void 	denormalize_ellipse_param(double* par, double* normalized_par, double dis_scale,  double norm_cx, double norm_cy);
	float	function_DoubleExponentialSeat (float x, float a);
	float	function_DoubleExponentialSigmoid (float x, float a);
	
	
	
	
		
	// parameters for acquiring line samples
	double 					src_quad[4][2];
	double 					dst_quad[4][2];
	
	ofVec3f					*handyVec3f0;
	ofVec3f					*initialPupilEstimate;
    
	
};


#endif // _PUPIL_CONTOUR_ANALYZER
