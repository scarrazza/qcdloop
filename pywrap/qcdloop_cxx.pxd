from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.vector cimport vector

cdef extern from "../src/qcdloop/qcdloop.h" namespace "ql":
    cdef cppclass QCDLoop[TOutput,TMass,TScale]:
        QCDLoop() except +
        void integral(vector[TOutput] & out,const TScale & mu2,const vector[TMass] & m,const vector[TMass] & p)
