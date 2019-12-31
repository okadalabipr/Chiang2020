#----------------------------------------------------------------------------
#                                CVODE
#----------------------------------------------------------------------------
# This file will compile all the CVODE modules in order to make the Unix
# library cvodelib.a. All object (.o) files are removed after the library
# has been created. The library cvodelib.a is placed in CVODE/lib.
#
# The following variables are used:
#
# COMPILER - set to use the gcc compiler
# OPTS     - list of compiler options
# OBJS     - list of object files in cvodelib.a
#
# Modify the COMPILER and OPTS variables as needed.
#
# If the ranlib utility is not available on your system, then remove the call
# to ranlib from the commands to create lib/libcvode.a.


COMPILER = cc

INC=../include
LIB=../../lib

OPTS = -I$(INC) -O3

OBJS = cvode.o cvode_band.o cvode_bandpre.o cvode_bbdpre.o cvode_dense.o cvode_diag.o cvode_direct.o cvode_io.o \
       cvode_lapack.o cvode_sparse.o cvode_spbcgs.o cvode_spgmr.o cvode_spils.o cvode_sptfqmr.o \
       sundials_band.o sundials_dense.o sundials_direct.o sundials_iterative.o sundials_math.o sundials_nvector.o \
       sundials_pcg.o sundials_sparse.o sundials_spbcgs.o sundials_spfgmr.o sundials_spgmr.o sundials_sptfqmr.o nvector_serial.o

lib $(LIB)/libsundials_cvode.a : $(OBJS)
	ar rcv $(LIB)/libsundials_cvode.a  $(OBJS)
	ranlib $(LIB)/libsundials_cvode.a

cvode.o: cvode.c $(INC)/cvode_impl.h $(INC)/sundials_types.h \
                $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c cvode.c

cvode_band.o: cvode_band.c $(INC)/cvode_band.h \
                 $(INC)/cvode_direct_impl.h $(INC)/cvode_impl.h \
                 $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c cvode_band.c

cvode_bandpre.o: cvode_bandpre.c $(INC)/cvode_bandpre_impl.h $(INC)/cvode_impl.h \
                    $(INC)/cvode_spils_impl.h $(INC)/cvode_sptfqmr.h $(INC)/cvode_spbcgs.h \
                    $(INC)/cvode_spgmr.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c cvode_bandpre.c

cvode_bbdpre.o: cvode_bbdpre.c $(INC)/cvode_impl.h $(INC)/cvode_bbdpre_impl.h $(INC)/cvode_spils_impl.h \
                $(INC)/cvode_sptfqmr.h $(INC)/cvode_spbcgs.h $(INC)/cvode_spgmr.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c cvode_bbdpre.c

cvode_dense.o: cvode_dense.c $(INC)/cvode_dense.h $(INC)/cvode_direct_impl.h \
                  $(INC)/cvode_impl.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c cvode_dense.c

cvode_diag.o: cvode_diag.c $(INC)/cvode_diag_impl.h  $(INC)/cvode_impl.h
	$(COMPILER) $(OPTS) -c cvode_diag.c

cvode_direct.o: cvode_direct.c $(INC)/cvode_impl.h  $(INC)/cvode_direct_impl.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c cvode_direct.c

cvode_io.o: cvode_io.c $(INC)/cvode_impl.h $(INC)/sundials_types.h
	$(COMPILER) $(OPTS) -c cvode_io.c

cvode_lapack.o: cvode_lapack.c $(INC)/cvode_lapack.h $(INC)/cvode_direct_impl.h $(INC)/cvode_impl.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c cvode_lapack.c

cvode_sparse.o: cvode_sparse.c $(INC)/cvode_impl.h $(INC)/cvode_sparse_impl.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c cvode_sparse.c

cvode_spbcgs.o: cvode_spbcgs.c $(INC)/cvode_spbcgs.h $(INC)/cvode_spils_impl.h $(INC)/cvode_impl.h \
                $(INC)/sundials_spbcgs.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c cvode_spbcgs.c

cvode_spgmr.o: cvode_spgmr.c $(INC)/cvode_spgmr.h $(INC)/cvode_impl.h \
                  $(INC)/cvode_spils_impl.h $(INC)/sundials_spgmr.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c cvode_spgmr.c

cvode_spils.o: cvode_spils.c $(INC)/cvode_impl.h $(INC)/cvode_spils_impl.h
	$(COMPILER) $(OPTS) -c cvode_spils.c

cvode_sptfqmr.o: cvode_sptfqmr.c $(INC)/cvode_sptfqmr.h $(INC)/cvode_spils_impl.h $(INC)/cvode_impl.h \
                 $(INC)/sundials_sptfqmr.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c cvode_sptfqmr.c


sundials_band.o: sundials_band.c $(INC)/sundials_band.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c sundials_band.c

sundials_dense.o: sundials_dense.c $(INC)/sundials_dense.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c sundials_dense.c

sundials_direct.o: sundials_direct.c $(INC)/sundials_direct.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c sundials_direct.c

sundials_iterative.o: sundials_iterative.c $(INC)/sundials_iterative.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c sundials_iterative.c

sundials_math.o: sundials_math.c $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c sundials_math.c

sundials_nvector.o: sundials_nvector.c $(INC)/sundials_nvector.h
	$(COMPILER) $(OPTS) -c sundials_nvector.c

sundials_pcg.o: sundials_pcg.c $(INC)/sundials_pcg.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c sundials_pcg.c

sundials_sparse.o: sundials_sparse.c $(INC)/sundials_sparse.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c sundials_sparse.c

sundials_spbcgs.o: sundials_spbcgs.c $(INC)/sundials_spbcgs.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c sundials_spbcgs.c

sundials_spfgmr.o: sundials_spfgmr.c $(INC)/sundials_spfgmr.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c sundials_spfgmr.c

sundials_spgmr.o: sundials_spgmr.c $(INC)/sundials_spgmr.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c sundials_spgmr.c

sundials_sptfqmr.o: sundials_sptfqmr.c $(INC)/sundials_sptfqmr.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c sundials_sptfqmr.c

nvector_serial.o: nvector_serial.c $(INC)/nvector_serial.h $(INC)/sundials_math.h
	$(COMPILER) $(OPTS) -c nvector_serial.c

clean:
	rm $(OBJS)
