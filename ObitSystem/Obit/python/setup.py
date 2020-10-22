from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc
setup( name="Obit", version="1.0",
       ext_modules=[Extension("_Obit",
                              ['Obit_wrap.c'],
                              extra_compile_args=['-DPACKAGE_NAME="Obit"', '-DPACKAGE_TARNAME="obit"', '-DPACKAGE_VERSION="1.0"', '-DPACKAGE_STRING="Obit_1.0"', '-DPACKAGE_BUGREPORT="bcotton@nrao.edu"', '-DPACKAGE="obit"', '-DVERSION="1.0"', '-DSTDC_HEADERS=1', '-DHAVE_SYS_TYPES_H=1', '-DHAVE_SYS_STAT_H=1', '-DHAVE_STDLIB_H=1', '-DHAVE_STRING_H=1', '-DHAVE_MEMORY_H=1', '-DHAVE_STRINGS_H=1', '-DHAVE_INTTYPES_H=1', '-DHAVE_STDINT_H=1', '-DHAVE_UNISTD_H=1', '-DHAVE_STDLIB_H=1', '-DHAVE_SSE=1', '-DHAVE_AVX=1', '-DHAVE_LIBCFITSIO=1', '-DHAVE_CFITSIO=1', '-DHAVE_FFTW3=1', '-DHAVE_GSL=1', '-DHAVE_XMLRPC=1', '-DHAVE_ZLIB=1', '-DHAVE_PLPLOT=1', '-DHAVE_PYTHON=1', '-DHAVE_WVR=1', '-DHAVE_FSEEKO=1', '-DOBIT_THREADS_ENABLED=1'],
                              library_dirs=['../lib', '/users/bcotton/lib', '/export/ssd/bcotton/Git/Obit/trunk/other//lib', '/usr/lib/gcc/x86_64-redhat-linux/3.4.6', '/usr/lib/gcc/x86_64-redhat-linux/3.4.6/../../../../lib64', '/usr/lib/gcc/x86_64-redhat-linux/3.4.6/../../..', '/lib64', '/usr/local/cuda-10.1/lib64', '/export/ssd/bcotton/Git/Obit/trunk/other/lib'],
                              libraries=['Obit', 'Obit', 'plplot', 'cfitsio', 'fftw3f', 'glib-2.0', 'gsl', 'gslcblas', 'm', 'csirocsa', 'dl', 'freetype', 'pthread', 'gthread-2.0', 'rt', 'stdc++', 'cudart', 'z', 'xmlrpc_client', 'curl', 'xmlrpc_client', 'xmlrpc_server_abyss'],
                              runtime_library_dirs=['/export/ssd/bcotton/Git/Obit/trunk/other/lib', '/export/ssd/bcotton/Git/Obit/trunk/other/lib', '/export/ssd/bcotton/Git/Obit/trunk/other/lib'])],
       include_dirs=['get_python_inc', '/usr/include/glib-2.0', '/usr/lib64/glib-2.0/include', '/export/ssd/bcotton/Git/Obit/trunk/other/include', '/users/bcotton/include/plplot', '/export/ssd/bcotton/Git/Obit/trunk/ObitSystem/Obit/include', '/usr/include/cfitsio', '/export/ssd/bcotton/Git/Obit/trunk/other/include/', '/users/bcotton/include/plplot', '/export/ssd/bcotton/Git/Obit/trunk/other/include', '/export/ssd/bcotton/Git/Obit/trunk/other/include']
)