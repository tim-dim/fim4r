#-----------------------------------------------------------------------
# File    : Makefile
# Contents: build shared object for frequent item set mining in R
# Author  : Christian Borgelt
# History : 2016.03.30 file created
#-----------------------------------------------------------------------
SHELL    = /bin/bash
THISDIR  = ../src
UTILDIR  = ../util/src
MATHDIR  = ../math/src
TRACTDIR = ../tract/src
APRIDIR  = ../apriori/src
ECLATDIR = ../eclat/src
FPGDIR   = ../fpgrowth/src
SAMDIR   = ../sam/src
RELIMDIR = ../relim/src
CARPDIR  = ../carpenter/src
ISTADIR  = ../ista/src
ACCDIR   = ../accretion/src

CC       = gcc -std=c99
# CC       = g++
RINC     = -I/usr/share/R/include
INCS     = -I$(UTILDIR) -I$(MATHDIR)  -I$(TRACTDIR) \
           -I$(APRIDIR) -I$(ECLATDIR) -I$(FPGDIR)   \
           -I$(SAMDIR)  -I$(RELIMDIR) -I$(CARPDIR)  \
           -I$(ISTADIR) -I$(ACCDIR)
CFBASE   = -Wall -Wextra -Wconversion -Wno-unused-parameter \
           -pedantic -c -fPIC $(ADDFLAGS)
CFLAGS   = $(CFBASE) -DNDEBUG -DQUIET -O3 -funroll-loops
# CFLAGS   = $(CFBASE) -DNDEBUG -O3 -DALIGN8
# CFLAGS   = $(CFBASE) -g
# CFLAGS   = $(CFBASE) -g -DSTORAGE
# CFLAGS   = $(CFBASE) -g -DALIGN8

LD       = gcc
LDFLAGS  = $(ADDFLAGS)
LIBS     = -lm

HDRS     = $(UTILDIR)/arrays.h    $(UTILDIR)/memsys.h  \
           $(UTILDIR)/symtab.h    $(UTILDIR)/random.h  \
           $(UTILDIR)/sigint.h    $(MATHDIR)/gamma.h   \
           $(MATHDIR)/chi2.h      $(MATHDIR)/ruleval.h \
           $(TRACTDIR)/tract.h    $(TRACTDIR)/fim16.h  \
           $(TRACTDIR)/patspec.h  $(TRACTDIR)/clomax.h \
           $(TRACTDIR)/report.h   $(TRACTDIR)/patred.h \
           $(APRIDIR)/apriori.h   $(ECLATDIR)/eclat.h  \
           $(FPGDIR)/fpgrowth.h   $(FPGDIR)/fpgpsp.h   \
           $(SAMDIR)/sam.h        $(RELIMDIR)/relim.h  \
           $(CARPDIR)/carpenter.h $(ISTADIR)/ista.h    \
           $(ACCDIR)/accretion.h
OBJS     = arrays.o memsys.o idmap.o random.o chi2.o sigint.o \
           gamma.o ruleval.o \
           tatree.o fim16.o patspec.o clomax.o report.o \
           istree.o apriori.o eclat.o fpgrowth.o sam.o relim.o \
           repotree.o carpenter.o pfxtree.o pattree.o ista.o \
           accretion.o fpgpsp.o fim4r.o $(ADDOBJS)

#-----------------------------------------------------------------------
# Build Shared Library
#-----------------------------------------------------------------------
all:            fim4r.so

fim4r.so:       $(OBJS) Makefile
	R CMD SHLIB -o fim4r.so $(OBJS)

#-----------------------------------------------------------------------
# Array Operations
#-----------------------------------------------------------------------
arrays.o:     $(UTILDIR)/arrays.h $(UTILDIR)/fntypes.h
arrays.o:     $(UTILDIR)/arrays.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(UTILDIR)/arrays.c -o $@

#-----------------------------------------------------------------------
# Memory Management System for Objects of Equal Size
#-----------------------------------------------------------------------
memsys.o:     $(UTILDIR)/memsys.h
memsys.o:     $(UTILDIR)/memsys.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(UTILDIR)/memsys.c -o $@

#-----------------------------------------------------------------------
# Symbol Table Management
#-----------------------------------------------------------------------
idmap.o:      $(UTILDIR)/symtab.h $(UTILDIR)/fntypes.h \
              $(UTILDIR)/arrays.h
idmap.o:      $(UTILDIR)/symtab.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DIDMAPFN $(UTILDIR)/symtab.c -o $@

#-----------------------------------------------------------------------
# Random Number Generator Management
#-----------------------------------------------------------------------
random.o:     $(UTILDIR)/random.h
random.o:     $(UTILDIR)/random.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(UTILDIR)/random.c -o $@

#-----------------------------------------------------------------------
# Interrupt Signal Handling
#-----------------------------------------------------------------------
sigint.o:     $(UTILDIR)/sigint.h
sigint.o:     $(UTILDIR)/sigint.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(UTILDIR)/sigint.c -o $@

#-----------------------------------------------------------------------
# Gamma Function
#-----------------------------------------------------------------------
gamma.o:      $(MATHDIR)/gamma.h
gamma.o:      $(MATHDIR)/gamma.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(MATHDIR)/gamma.c -o $@

chi2.o:       $(MATHDIR)/chi2.h
chi2.o:       $(MATHDIR)/chi2.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(MATHDIR)/chi2.c -o $@

ruleval.o:    $(MATHDIR)/ruleval.h
ruleval.o:    $(MATHDIR)/ruleval.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(MATHDIR)/ruleval.c -o $@

#-----------------------------------------------------------------------
# 16 Items Machine
#-----------------------------------------------------------------------
fim16.o:      $(TRACTDIR)/tract.h $(TRACTDIR)/report.h \
              $(UTILDIR)/arrays.h  $(UTILDIR)/symtab.h
fim16.o:      $(TRACTDIR)/fim16.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(TRACTDIR)/fim16.c -o $@

#-----------------------------------------------------------------------
# Item and Transaction Management
#-----------------------------------------------------------------------
tatree.o:     $(TRACTDIR)/tract.h $(UTILDIR)/arrays.h \
              $(UTILDIR)/symtab.h
tatree.o:     $(TRACTDIR)/tract.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DTATREEFN -DTA_SURR \
              $(TRACTDIR)/tract.c -o $@

#-----------------------------------------------------------------------
# Item Set Reporter Management
#-----------------------------------------------------------------------
patspec.o:    $(TRACTDIR)/patspec.h $(TRACTDIR)/tract.h
patspec.o:    $(TRACTDIR)/patspec.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DPSP_ESTIM \
              $(TRACTDIR)/patspec.c -o $@

clomax.o:     $(TRACTDIR)/clomax.h $(TRACTDIR)/tract.h \
              $(UTILDIR)/arrays.h
clomax.o:     $(TRACTDIR)/clomax.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(TRACTDIR)/clomax.c -o $@

report.o:     $(TRACTDIR)/report.h $(TRACTDIR)/clomax.h \
              $(TRACTDIR)/tract.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/symtab.h
report.o:     $(TRACTDIR)/report.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DISR_PATSPEC -DISR_CLOMAX \
              -DISR_NONAMES $(TRACTDIR)/report.c -o $@

#-----------------------------------------------------------------------
# Pattern Set Reduction Functions
#-----------------------------------------------------------------------
patred.o:     $(UTILDIR)/arrays.h  $(TRACTDIR)/tract.h \
              $(TRACTDIR)/report.h $(TRACTDIR)/patred.h
patred.o:     $(TRACTDIR)/patred.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(TRACTDIR)/patred.c -o $@

#-----------------------------------------------------------------------
# Apriori
#-----------------------------------------------------------------------
istree.o:     $(HDRS) $(APRIDIR)/istree.h
istree.o:     $(APRIDIR)/istree.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DTATREEFN $(APRIDIR)/istree.c -o $@

apriori.o:    $(HDRS) $(APRIDIR)/istree.h $(APRIDIR)/apriori.h \
              $(UTILDIR)/fntypes.h
apriori.o:    $(APRIDIR)/apriori.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DISR_CLOMAX -DAPR_ABORT \
              $(APRIDIR)/apriori.c -o $@

#-----------------------------------------------------------------------
# Eclat
#-----------------------------------------------------------------------
eclat.o:      $(HDRS) $(ECLATDIR)/eclat.h $(UTILDIR)/fntypes.h
eclat.o:      $(ECLATDIR)/eclat.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DECL_ABORT $(ECLATDIR)/eclat.c -o $@

#-----------------------------------------------------------------------
# FP-growth
#-----------------------------------------------------------------------
fpgrowth.o:   $(HDRS) $(FPGDIR)/fpgrowth.h $(UTILDIR)/fntypes.h
fpgrowth.o:   $(FPGDIR)/fpgrowth.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DFPG_ABORT $(FPGDIR)/fpgrowth.c -o $@

fpgpsp.o:     $(HDRS) $(FPGDIR)/fpgrowth.h $(FPGDIR)/fpgpsp.h \
              $(UTILDIR)/fntypes.h
fpgpsp.o:     $(FPGDIR)/fpgrowth.c $(FPGDIR)/fpgpsp.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DFPG_ABORT $(FPGDIR)/fpgpsp.c -o $@

#-----------------------------------------------------------------------
# SaM
#-----------------------------------------------------------------------
sam.o:        $(HDRS) $(SAMDIR)/sam.h $(UTILDIR)/fntypes.h
sam.o:        $(SAMDIR)/sam.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DSAM_ABORT $(SAMDIR)/sam.c -o $@

#-----------------------------------------------------------------------
# RElim
#-----------------------------------------------------------------------
relim.o:      $(HDRS) $(RELIMDIR)/relim.h $(UTILDIR)/fntypes.h
relim.o:      $(RELIMDIR)/relim.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DRELIM_ABORT $(RELIMDIR)/relim.c -o $@

#-----------------------------------------------------------------------
# Carpenter
#-----------------------------------------------------------------------
repotree.o:   $(HDRS) $(CARPDIR)/repotree.h
repotree.o:   $(CARPDIR)/repotree.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(CARPDIR)/repotree.c -o $@

carpenter.o:  $(HDRS) $(CARPDIR)/carpenter.h $(CARPDIR)/repotree.h \
              $(UTILDIR)/fntypes.h
carpenter.o:  $(CARPDIR)/carpenter.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DCARP_ABORT \
              $(CARPDIR)/carpenter.c -o $@

#-----------------------------------------------------------------------
# IsTa
#-----------------------------------------------------------------------
pfxtree.o:    $(HDRS) $(ISTADIR)/pfxtree.h $(UTILDIR)/memsys.h
pfxtree.o:    $(ISTADIR)/pfxtree.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(ISTADIR)/pfxtree.c -o $@

pattree.o:    $(HDRS) $(ISTADIR)/pattree.h
pattree.o:    $(ISTADIR)/pattree.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(ISTADIR)/pattree.c -o $@

ista.o:       $(HDRS) $(ISTADIR)/pattree.h $(ISTADIR)/pfxtree.h \
              $(ISTADIR)/ista.h
ista.o:       $(ISTADIR)/ista.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DISTA_ABORT $(ISTADIR)/ista.c -o $@

#-----------------------------------------------------------------------
# Accretion
#-----------------------------------------------------------------------
accretion.o:  $(HDRS) $(ACCDIR)/accretion.h $(UTILDIR)/fntypes.h
accretion.o:  $(ACCDIR)/accretion.c Makefile
	$(CC) $(CFLAGS) $(INCS) -DISR_CLOMAX -DACC_ABORT \
              $(ACCDIR)/accretion.c -o $@

#-----------------------------------------------------------------------
# R Stuff
#-----------------------------------------------------------------------
fim4r.o:      fim4r.c Makefile
	$(CC) $(CFLAGS) $(INCS) $(RINC) fim4r.c -o $@

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ..; rm -f fim4r.zip fim4r.tar.gz; \
        zip -rq fim4r.zip \
          fim4r/{R,src,data,man,tests,DESCRIPTION,NAMESPACE} \
          fim4r/{apriori,eclat,fpgrowth,sam,relim} \
	  fim4r/{carpenter,ista,accretion} \
	  fim4r/{tract,math,util}; \
        tar cfz fim4r.tar.gz \
          fim4r/{R,src,data,man,tests,DESCRIPTION,NAMESPACE} \
          fim4r/{fimnad,tract,math,util}

#-----------------------------------------------------------------------
# Install
#-----------------------------------------------------------------------
install:
	cd ../..; R CMD INSTALL fim4r

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
clean:
	-rm -f *.o *.so *~ *.flc *.log core
