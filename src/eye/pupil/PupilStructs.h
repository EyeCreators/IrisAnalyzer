// PupilStructs.h

#ifndef _PUPIL_STRUCTS
#define _PUPIL_STRUCTS

#include "ofMain.h"
#include "ipp.h"

#define INVALID_CID	-1


//----------------------------
typedef struct {
	Ipp8u	*pupilContourImage;
	int	img_w;
	int	img_h;
	int	img_s;
	int whichEye;
	int	contour_index;
	
	float cx, cy;	// centroid, from report, pixel units in eye image
	float bx, by;	// bounding rect, from report, pixel units in eye image
	float bw, bh;	// bounding rect, from report, pixel units in eye image
	
	ofVec3f	ncCentroid; // 0..1 in original source image, from REPORT
	
} PupilContourImagePacket;

//----------------------------
typedef struct {
	bool		bValid;
	
	float		params_ellipse [5];
	float		params_conic [6];
	ofVec3f		normals[2];
	ofVec3f		centroidInEyeImg;

} PupilGeometry;


//----------------------------
#endif // _PUPIL_STRUCTS