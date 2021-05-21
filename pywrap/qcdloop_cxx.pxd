from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.vector cimport vector

cdef extern from "../src/qcdloop/qcdloop.h" namespace "ql":
    cdef cppclass QCDLoop[TOutput,TMass,TScale]:
        QCDLoop() except +
        void integral(vector[TOutput] & out,const TScale & mu2,const vector[TMass] & m,const vector[TScale] & p)

cdef extern from "../src/qcdloop/tadpole.h" namespace "ql":
    cdef cppclass TadPole[TOutput,TMass,TScale]:
        TadPole() except +
        void integral(vector[TOutput] & out,const TScale & mu2,const vector[TMass] & m,const vector[TScale] & p)

cdef extern from "../src/qcdloop/bubble.h" namespace "ql":
    cdef cppclass Bubble[TOutput,TMass,TScale]:
        Bubble() except +
        void integral(vector[TOutput] & out,const TScale & mu2,const vector[TMass] & m,const vector[TScale] & p)

cdef extern from "../src/qcdloop/triangle.h" namespace "ql":
    cdef cppclass Triangle[TOutput,TMass,TScale]:
        Triangle() except +
        void integral(vector[TOutput] & out,const TScale & mu2,const vector[TMass] & m,const vector[TScale] & p)

cdef extern from "../src/qcdloop/box.h" namespace "ql":
    cdef cppclass Box[TOutput,TMass,TScale]:
        Box() except +
        void integral(vector[TOutput] & out,const TScale & mu2,const vector[TMass] & m,const vector[TScale] & p)
