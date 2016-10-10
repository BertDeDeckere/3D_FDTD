#include "field.h"

void Field::update(int cellsX, int cellsY, int cellsZ, int sampleDistance, double nmax, double dt, \
	float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, \
	float* CExdy, float* CExdz, float* CEydx, float* CEydz, float* CEzdx, float* CEzdy, \
	float* CHxdy, float* CHxdz, float* CHydx, float* CHydz, float* CHzdx, float* CHzdy,
	float* CEx, float* CEy, float* CEz) {

	// std::cout << cellsX << " " << cellsY << " "  << cellsZ << " "  << sampleDistance << " "  << nmax << " " << dt << "\n" ;

	this->cellsX = cellsX;
	this->cellsY = cellsY;
	this->cellsZ = cellsZ;

	float *oldEx, *newEx, *oldEy, *newEy, *oldEz, *newEz;
	float *oldHx, *newHx, *oldHy, *newHy, *oldHz, *newHz;
	int pos, pos1, pos2;

	for(int n=0; n<nmax; n++) {	
		if(sampleDistance != 1) {
			if(n%sampleDistance == 0) {
				pos = n/sampleDistance*cellsX*cellsY*cellsZ;
				newEx = &OBEx[pos];		oldEx = Ex;
				newEy = &OBEy[pos];		oldEy = Ey;
				newEz = &OBEz[pos];		oldEz = Ez;
				newHx = &OBHx[pos];		oldHx = Hx;
				newHy = &OBHy[pos];		oldHy = Hy;
				newHz = &OBHz[pos];		oldHz = Hz;			
			}
			else if(n%sampleDistance == 1) {
				pos = (n-1)/sampleDistance*cellsX*cellsY*cellsZ;
				newEx = Ex;		oldEx = &OBEx[pos];
				newEy = Ey;		oldEy = &OBEy[pos];
				newEz = Ez;		oldEz = &OBEz[pos];
				newHx = Hx;		oldHx = &OBHx[pos];
				newHy = Hy;		oldHy = &OBHy[pos];
				newHz = Hz;		oldHz = &OBHz[pos];
			}
			else {
				newEx = Ex;		oldEx = Ex;
				newEy = Ey;		oldEy = Ey;
				newEz = Ez;		oldEz = Ez;
				newHx = Hx;		oldHx = Hx;
				newHy = Hy;		oldHy = Hy;
				newHz = Hz;		oldHz = Hz;
			}
		}
		else {
			if(n>0) {
				pos1 = n*cellsX*cellsY*cellsZ;
				pos2 = (n-1)*cellsX*cellsY*cellsZ;
			}
			else {
				pos1 = n*cellsX*cellsY*cellsZ;
				pos2 = n*cellsX*cellsY*cellsZ;
			}
			newEx = &OBEx[pos1];		oldEx = &OBEx[pos2];
			newEy = &OBEy[pos1];		oldEy = &OBEy[pos2];
			newEz = &OBEz[pos1];		oldEz = &OBEz[pos2];
			newHx = &OBHx[pos1];		oldHx = &OBHx[pos2];
			newHy = &OBHy[pos1];		oldHy = &OBHy[pos2];
			newHz = &OBHz[pos1];		oldHz = &OBHz[pos2];			
		}

		// Update E
		for(int i=0; i<cellsX-1; i++) {
			for(int j=0; j<cellsY-1; j++) {
				for(int k=0; k<cellsZ-1; k++) {
					newEx[a(i,j,k)] = CEx[a(i,j,k)]*oldEx[a(i,j,k)] + CExdy[a(i,j,k)]*(oldHz[a(i,j+1,k)] - oldHz[a(i,j,k)]) + CExdz[a(i,j,k)]*(oldHy[a(i,j,k+1)] - oldHy[a(i,j,k)]);
					newEy[a(i,j,k)] = CEy[a(i,j,k)]*oldEy[a(i,j,k)] + CEydx[a(i,j,k)]*(oldHz[a(i+1,j,k)] - oldHz[a(i,j,k)]) + CEydz[a(i,j,k)]*(oldHx[a(i,j,k+1)] - oldHx[a(i,j,k)]);
					newEz[a(i,j,k)] = CEz[a(i,j,k)]*oldEz[a(i,j,k)] + CEzdx[a(i,j,k)]*(oldHy[a(i+1,j,k)] - oldHy[a(i,j,k)]) + CEzdy[a(i,j,k)]*(oldHx[a(i,j+1,k)] - oldHx[a(i,j,k)]);
				}
			}
		}

		// Deal with the sources
		for(int m=0; m<waveguideSourcesXY.size(); m++) {
			WaveguideSourceXY WG = waveguideSourcesXY.at(m);

			for(int i=WG.imin; i<WG.imax; i++) {
				for(int j=WG.jmin; j<WG.jmax-1; j++) {
					newEx[a(i,j,WG.k)] += CExdz[a(i,j,WG.k)]*WG.Hy[WG.b(i-WG.imin,j-WG.jmin)]*WG.modE[n];
				}
			}

			for(int i=WG.imin; i<WG.imax-1; i++) {
				for(int j=WG.jmin; j<WG.jmax; j++) {
					newEy[a(i,j,WG.k)] += CEydz[a(i,j,WG.k)]*WG.Hx[WG.a(i-WG.imin,j-WG.jmin)]*WG.modE[n];
				}
			}
		}

		// newEx[a(cellsX/2, cellsY/2, cellsZ/2)] += CExdy[a(cellsX/2, cellsY/2, cellsZ/2)]*sin(2*3.1416*n*dt*1E9)*1E-3;
		// newEx[a(1,1,1)] += CExdy[a(1,1,1)]*sin(2*3.1416*n*dt*1E9)*1E-3;

		// Update H
		for(int i=1; i<cellsX; i++) {
			for(int j=1; j<cellsY; j++) {
				for(int k=1; k<cellsZ; k++) {
					newHx[a(i,j,k)] = oldHx[a(i,j,k)] + CHxdy[a(i,j,k)]*(newEz[a(i,j,k)] - newEz[a(i,j-1,k)]) + CHxdz[a(i,j,k)]*(newEy[a(i,j,k)] - newEy[a(i,j,k-1)]);
					newHy[a(i,j,k)] = oldHy[a(i,j,k)] + CHydx[a(i,j,k)]*(newEz[a(i,j,k)] - newEz[a(i-1,j,k)]) + CHydz[a(i,j,k)]*(newEx[a(i,j,k)] - newEx[a(i,j,k-1)]);
					newHz[a(i,j,k)] = oldHz[a(i,j,k)] + CHzdx[a(i,j,k)]*(newEy[a(i,j,k)] - newEy[a(i-1,j,k)]) + CHzdy[a(i,j,k)]*(newEx[a(i,j,k)] - newEx[a(i,j-1,k)]);
				}
			}
		}

		// Deal with the sources H-field
		for(int m=0; m<waveguideSourcesXY.size(); m++) {
			WaveguideSourceXY WG = waveguideSourcesXY.at(m);

			for(int i=WG.imin; i<WG.imax-1; i++) {
				for(int j=WG.jmin; j<WG.jmax; j++) {
					newHx[a(i,j,WG.k)] -= CHxdz[a(i,j,WG.k)]*WG.Ey[WG.a(i-WG.imin,j-WG.jmin)]*WG.modH[n];
				}
			}

			for(int i=WG.imin; i<WG.imax; i++) {
				for(int j=WG.jmin; j<WG.jmax-1; j++) {
					newHy[a(i,j,WG.k)] -= CHydz[a(i,j,WG.k)]*WG.Ex[WG.b(i-WG.imin,j-WG.jmin)]*WG.modH[n];
				}
			}
		}
	}
}

void Field::setOutputBuffer(int sizeOB, float* OBEx, float* OBEy, float* OBEz, float* OBHx, float* OBHy, float* OBHz) {
	this->sizeOB = sizeOB;
	this->OBEx = OBEx;
	this->OBEy = OBEy;
	this->OBEz = OBEz;
	this->OBHx = OBHx;
	this->OBHy = OBHy;
	this->OBHz = OBHz;
}

void Field::addWaveguideSourceXY(int imin, int imax, int jmin, int jmax, int k, float *modE, float* modH, float *Ex, float *Ey, float *Hx, float *Hy) {
	waveguideSourcesXY.push_back(WaveguideSourceXY(imin, imax, jmin, jmax, k, modE, modH, Ex, Ey, Hx, Hy));
}

int Field::a(int i, int j, int k) {
	return i*cellsY*cellsZ + j*cellsZ + k;
}

WaveguideSourceXY::WaveguideSourceXY(int imin, int imax, int jmin, int jmax, int k, float *modE, float *modH, float *Ex, float *Ey, float *Hx, float *Hy) {
	this->imin = imin;
	this->imax = imax;
	this->jmin = jmin;
	this->jmax = jmax;
	this->k = k;
	this->modE = modE;
	this->modH = modH;
	this->Ex = Ex;
	this->Ey = Ey;
	this->Hx = Hx;
	this->Hy = Hy;

	this->sizeY0 = imax-imin;
	this->sizeY1 = imax-imin-1;
}

int WaveguideSourceXY::a(int i, int j) {
	return i*sizeY0 + j;
}

int WaveguideSourceXY::b(int i, int j) {
	return i*sizeY1 + j;
}

