from setuptools import setup, Extension

setup(name='fourier',
      version='1.0',
      ext_modules=[Extension('fourier', ['exportFourier.cpp'], libraries=['fftw3'])], extra_compile_args=['-O3 -flto -march=native'])
