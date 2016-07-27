noinst_PROGRAMS = trigger_test cache_test cmass_test
trigger_test_SOURCES = trigger_test.cc
cache_test_SOURCES = cache_test.cc
cmass_test_SOURCES = cmass_test.cc

AM_CPPFLAGS = -I$(top_srcdir)/src -Wall -Wextra -fPIC
AM_LDFLAGS = -L$(top_builddir)/src/.libs
LIBS = -lqcdloop -lquadmath
