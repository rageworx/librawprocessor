# Makefile for librawprocessor
# (C)2016 Raphael Kim / rageworx
#

# To enable build for embedded linux, you may encomment next 2 lines.
# CCPREPATH = ${ARM_LINUX_GCC_PATH}
# CCPREFIX  = arm-linux-

# To enable build for embedded linux, change following line.
# CCPATH    = ${CCPREPATH}/${CCPREFIX}
CCPATH =

# Compiler configure.
GCC = ${CCPATH}gcc
GPP = ${CCPATH}g++
AR  = ${CCPATH}ar

SOURCEDIR = ./src
OBJDIR    = ./obj/Release
OUTBIN    = librawprocessor.a
OUTDIR    = ./lib
DEFINEOPT = -DRAWPROCESSOR_USE_LOCALTCHAR
OPTIMIZEOPT = -O3 -s
CFLAGS    = -I$(SOURCEDIR) $(DEFINEOPT) $(OPTIMIZEOPT)

all: prepare clean ${OUTDIR}/${OUTBIN}

prepare:
	@mkdir -p ${OBJDIR}
	@mkdir -p ${OUTDIR}

${OBJDIR}/stdunicode.o:
	$(GPP) -c ${SOURCEDIR}/stdunicode.cpp ${CFLAGS} -o $@

${OBJDIR}/rawscale.o:
	$(GPP) -c ${SOURCEDIR}/rawscale.cpp ${CFLAGS} -o $@

${OBJDIR}/rawprocessor.o:
	$(GPP) -c ${SOURCEDIR}/rawprocessor.cpp ${CFLAGS} -o $@

${OUTDIR}/${OUTBIN}: ${OBJDIR}/rawscale.o ${OBJDIR}/stdunicode.o ${OBJDIR}/rawprocessor.o
	$(AR) -q $@ ${OBJDIR}/*.o

clean:
	@rm -rf ${OBJDIR}/*
	@rm -rf ${OUTDIR}/${OUTBIN}
