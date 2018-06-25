lib_LTLIBRARIES = libqcdloop.la
libqcdloop_la_SOURCES = qcdloop.cc \
                        tools.cc \
			cache.cc \
                        topology.cc \
                        tadpole.cc \
                        types.cc \
                        bubble.cc \
                        wrapper.cc \
                        triangle.cc \
			box.cc

otherincludedir = $(includedir)/qcdloop
otherinclude_HEADERS = qcdloop/qcdloop.h \
                       qcdloop/tools.h \
                       qcdloop/topology.h \
                       qcdloop/tadpole.h \
		       qcdloop/types.h \
                       qcdloop/maths.h \
                       qcdloop/bubble.h \
                       qcdloop/wrapper.h \
                       qcdloop/cache.h \
                       qcdloop/timer.h \
                       qcdloop/exceptions.h \
                       qcdloop/triangle.h \
		       qcdloop/box.h

libqcdloop_la_LDFLAGS = -version-info 0:0:0
libqcdloop_la_CPPFLAGS= -Wall -Wextra -fPIC $(MACROS)
