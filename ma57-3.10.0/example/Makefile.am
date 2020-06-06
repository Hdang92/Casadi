AUTOMAKE_OPTIONS = serial-tests
AM_FCFLAGS = -I $(top_builddir)/src

TESTS_ENVIRONMENT = $(srcdir)/runtest.sh
EXTRA_DIST = runtest.sh
TESTS = \
	ma57ss \
	ma57ds

check_PROGRAMS = \
	ma57ss \
	ma57ds

ma57ss_SOURCES = ma57ss.f
ma57ss_LDADD = $(top_builddir)/src/libma57.la $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += ma57ss.output ma57ss.data
ma57ds_SOURCES = ma57ds.f
ma57ds_LDADD = $(top_builddir)/src/libma57.la $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += ma57ds.output ma57ds.data
