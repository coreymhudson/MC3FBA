bin_PROGRAMS = fba
fba_SOURCES = main.cpp Biomass.cpp Biomass.h random.cpp random.h stoich_mat.cpp \
				stoich_mat.h Tokenizer.h
if INTERNAL_GLPK
  fba_CPPFLAGS = -I$(top_srcdir)/glpk
  fba_LDADD = ../glpk/libglpk.a
endif
if SYSTEM_GLPK
  AM_CPPFLAGS = $(DEPS_CFLAGS)
  fba_LDADD = $(DEPS_LIBS)
endif

