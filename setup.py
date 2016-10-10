from distutils.core import setup, Extension
import numpy

setup(ext_modules=[Extension("_Field",
      swig_opts=['-c++'],
      sources=["field.cpp", "field.i"],
      extra_compile_args=["-O3"],
      include_dirs=[numpy.get_include()])]
      )

 	# run with build_ext --inplace