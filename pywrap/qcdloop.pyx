# distutils: language = c++
# distutils: extra_compile_args = [-std=c++11, -fext-numeric-literals, -fPIC]
# distutils: libraries = ['qcdloop','quadmath']

from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.vector cimport vector
cimport qcdloop_cxx as cpp

cdef class QCDLoop:
    cdef cpp.QCDLoop[double complex,double,double] *thisptr

    def __cinit__(self):
        self.thisptr = new cpp.QCDLoop[double complex,double,double]()

    def __dealloc__(self):
        del self.thisptr

    def integral(self, mu2, m, p = []):
        cdef vector[double complex] res
        self.thisptr.integral(res,mu2,m,p)
        out = [res[0], res[1], res[2]]
        return out

    def integral(self, mu2, m):
        cdef vector[double complex] res
        self.thisptr.integral(res,mu2,m,[])
        out = [res[0], res[1], res[2]]
        return out

