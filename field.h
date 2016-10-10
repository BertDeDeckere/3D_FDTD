#include <math.h>
#include <iostream>
#include <vector>

class WaveguideSourceXY;
class Field;

class Field {
	int cellsX, cellsY, cellsZ;
	int sizeOB;

	float *OBEx, *OBEy, *OBEz, *OBHx, *OBHy, *OBHz;
	std::vector<WaveguideSourceXY> waveguideSourcesXY;

public:
	void update(int , int , int , int, double, double , \
		float* , float* , float* , float* , float* , float* , \
		float* , float* , float* , float* , float* , float* , \
		float* , float* , float* , float* , float* , float* ,
		float* , float* , float*);

	void setOutputBuffer(int, float*, float*, float*, float*, float*, float*);
	void addWaveguideSourceXY(int, int, int, int, int, float*, float*, float*, float*, float*, float*);
	int a(int, int, int);	// Compute argument for array indexing
};

class WaveguideSourceXY {
public:
	int imin, imax, jmin, jmax, k, sizeY0, sizeY1;
	float *modE, *modH;
	float *Ex, *Ey, *Hx, *Hy;

	WaveguideSourceXY(int, int, int, int, int, float*, float*, float*, float*, float*, float*);
	int a(int, int);
	int b(int, int);
};