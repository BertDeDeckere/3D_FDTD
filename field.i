/*  Example of wrapping a C function that takes a C double array as input using
 *  numpy typemaps for SWIG. */

%module Field
%{
    #define SWIG_FILE_WITH_INIT     /* the resulting C file should be built as a python extension */
    #include "field.h"              /*  Includes the header in the wrapper code */
%}

%include "numpy.i"                  /*  include the numpy typemaps */
%init %{                            /*  need this for correct module initialization */
    import_array();
%}

/*  typemaps for the two arrays, the second will be modified in-place */
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* Ex, int cellsX01, int cellsY01, int cellsZ01)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* Ey, int cellsX02, int cellsY02, int cellsZ02)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* Ez, int cellsX03, int cellsY03, int cellsZ03)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* Hx, int cellsX04, int cellsY04, int cellsZ04)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* Hy, int cellsX05, int cellsY05, int cellsZ05)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* Hz, int cellsX06, int cellsY06, int cellsZ06)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CExdy, int cellsX07, int cellsY07, int cellsZ07)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CExdz, int cellsX08, int cellsY08, int cellsZ08)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CEydx, int cellsX09, int cellsY09, int cellsZ09)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CEydz, int cellsX10, int cellsY10, int cellsZ10)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CEzdx, int cellsX11, int cellsY11, int cellsZ11)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CEzdy, int cellsX12, int cellsY12, int cellsZ12)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CHxdy, int cellsX13, int cellsY13, int cellsZ13)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CHxdz, int cellsX14, int cellsY14, int cellsZ14)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CHydx, int cellsX15, int cellsY15, int cellsZ15)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CHydz, int cellsX16, int cellsY16, int cellsZ16)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CHzdx, int cellsX17, int cellsY17, int cellsZ17)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CHzdy, int cellsX18, int cellsY18, int cellsZ18)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CEx, int cellsX19, int cellsY19, int cellsZ19)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CEy, int cellsX20, int cellsY20, int cellsZ20)}
%apply (float* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(float* CEz, int cellsX21, int cellsY21, int cellsZ21)}

%apply (float* IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) {(float* OBEx, int n25, int cellsX25, int cellsY25, int cellsZ25)}
%apply (float* IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) {(float* OBEy, int n26, int cellsX26, int cellsY26, int cellsZ26)}
%apply (float* IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) {(float* OBEz, int n27, int cellsX27, int cellsY27, int cellsZ27)}
%apply (float* IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) {(float* OBHx, int n28, int cellsX28, int cellsY28, int cellsZ28)}
%apply (float* IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) {(float* OBHy, int n29, int cellsX29, int cellsY29, int cellsZ29)}
%apply (float* IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) {(float* OBHz, int n30, int cellsX30, int cellsY30, int cellsZ30)}

%apply (float* IN_ARRAY1, int DIM1) {(float* WGModE, int size1)}
%apply (float* IN_ARRAY1, int DIM1) {(float* WGModH, int size2)}
%apply (float* IN_ARRAY2, int DIM1, int DIM2) {(float* WGEx, int cellsX1, int cellsY1)}
%apply (float* IN_ARRAY2, int DIM1, int DIM2) {(float* WGEy, int cellsX2, int cellsY2)}
%apply (float* IN_ARRAY2, int DIM1, int DIM2) {(float* WGHx, int cellsX3, int cellsY3)}
%apply (float* IN_ARRAY2, int DIM1, int DIM2) {(float* WGHy, int cellsX4, int cellsY4)}

%include "field.h"

%extend Field {
    void update_cpp(int sampleDistance, double nmax, double dt, \
        float* Ex, int cellsX01, int cellsY01, int cellsZ01, float* Ey, int cellsX02, int cellsY02, int cellsZ02, \
        float* Ez, int cellsX03, int cellsY03, int cellsZ03, float* Hx, int cellsX04, int cellsY04, int cellsZ04, \
        float* Hy, int cellsX05, int cellsY05, int cellsZ05, float* Hz, int cellsX06, int cellsY06, int cellsZ06, \
        float* CExdy, int cellsX07, int cellsY07, int cellsZ07, float* CExdz, int cellsX08, int cellsY08, int cellsZ08, \
        float* CEydx, int cellsX09, int cellsY09, int cellsZ09, float* CEydz, int cellsX10, int cellsY10, int cellsZ10, \
        float* CEzdx, int cellsX11, int cellsY11, int cellsZ11, float* CEzdy, int cellsX12, int cellsY12, int cellsZ12, \
        float* CHxdy, int cellsX13, int cellsY13, int cellsZ13, float* CHxdz, int cellsX14, int cellsY14, int cellsZ14, \
        float* CHydx, int cellsX15, int cellsY15, int cellsZ15, float* CHydz, int cellsX16, int cellsY16, int cellsZ16, \
        float* CHzdx, int cellsX17, int cellsY17, int cellsZ17, float* CHzdy, int cellsX18, int cellsY18, int cellsZ18, \
        float* CEx, int cellsX19, int cellsY19, int cellsZ19, float* CEy, int cellsX20, int cellsY20, int cellsZ20, \
        float* CEz, int cellsX21, int cellsY21, int cellsZ21) {

        $self->update(cellsX01, cellsY01, cellsZ01, sampleDistance, nmax, dt, Ex, Ey, Ez, Hx, Hy, Hz, \
        CExdy, CExdz, CEydx, CEydz, CEzdx, CEzdy, CHxdy, CHxdz, CHydx, CHydz, CHzdx, CHzdy, CEx, CEy, CEz);
    }

    void setOutputBuffer_cpp(float* OBEx, int n25, int cellsX25, int cellsY25, int cellsZ25, \
        float* OBEy, int n26, int cellsX26, int cellsY26, int cellsZ26, float* OBEz, int n27, int cellsX27, int cellsY27, int cellsZ27, \
        float* OBHx, int n28, int cellsX28, int cellsY28, int cellsZ28, float* OBHy, int n29, int cellsX29, int cellsY29, int cellsZ29, \
        float* OBHz, int n30, int cellsX30, int cellsY30, int cellsZ30) {

        $self->setOutputBuffer(n25, OBEx, OBEy, OBEz, OBHx, OBHy, OBHz);
        }

    void addWaveguideSourceXY_cpp(int imin, int imax, int jmin, int jmax, int k, float *WGModE, int size1, float* WGModH, int size2, \
        float *WGEx, int cellsX1, int cellsY1, float *WGEy, int cellsX2, int cellsY2, float *WGHx, int cellsX3, int cellsY3, \
        float *WGHy, int cellsX4, int cellsY4) {

        $self->addWaveguideSourceXY(imin, imax, jmin, jmax, k, WGModE, WGModH, WGEx, WGEy, WGHx, WGHy);
        }
};

