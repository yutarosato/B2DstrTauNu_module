#####################################
# Makefile (for Linux)              #
#  by S.Nishida	 (update: 05/05/22) #
#  by T.Higuchi	 (update: 07/01/27) #
#####################################

# Generated automatically from Makefile.in by configure.
#+
# File: Makefile.in
# Description : Makefile Template for fsim_analysis
#-

#---------------------------------------------------------------
# 1. User Specifications
#      MODULE1   : module name (shared object file)
#      OBJS1     : additonal obj files
#      DEPFILES1 : additonal files on which $(MODULE1) depends
#      LIBS1     : libraries in addition to $(LIBS)
#---------------------------------------------------------------

MODULE1 = DSTRTAUNU.so
OBJS1 = DSTRTAUNU.o DSTRTAUNU_event.o UserInfo.o Util.o Util_Gen.o Gen_info.o NP_Weight.o Ks.o Kaon.o Gamma.o Mu.o EE.o Tracks.o Pion.o Pi0.o ksfwmoments.o Rec_D_2body.o Rec_D_3body.o Rec_D_4body.o Rec_Dstr.o Rec_Dtau.o Rec_CC.o FullRecon.o Ana.o Ana_single.o Ana_hadtag.o Util_dstrtaunu.o
DEPFILES1 =
LIBS1 = $(LIBS)

MODULE2 = 
OBJS2 = 
DEPFILES2 =
LIBS2 = $(LIBS)

MODULE3 = 
OBJS3 = 
DEPFILES3 =
LIBS3 = $(LIBS)

MODULE4 = 
OBJS4 = 
DEPFILES4 =
LIBS4 = $(LIBS)

MODULE5 = 
OBJS5 = 
DEPFILES5 =
LIBS5 = $(LIBS)

MODULE6 = 
OBJS6 = 
DEPFILES6 =
LIBS6 = $(LIBS)

MODULE7 =
OBJS7 =
DEPFILES7 =
LIBS7 = $(LIBS)

LIBS = -L$(MY_RUN_DIR)/lib/so -L$(MY_RUN_DIR)/lib -L$(BELLE_RUN_DIR)/lib/so \
	-ltuple -lparticle -lpanther -lmdst -lip -leid -lhamlet -lkfitter -ltagv -lbenergy

CERNLIB = -L$(CERN)/$(CERN_LEVEL)/lib64 \
	-lmathlib

#---------------------------------------------------------------
# 2. System Specifications
#    --- Do not change without knowledge
#---------------------------------------------------------------

# Compiler Setup with machine dependence

FC =
CC = gcc
CXX = c++

DEFS = -DHAVE_LIBCURSES=1 -DHAVE_LIBREADLINE=1 -DHAVE_POSTGRES=1 -DHAVE_LIBCURSES=1 -DHAVE_LIBTERMCAP=1 -DHAVE_LIBHISTORY=1 -DHAVE_LIBREADLINE=1 -DHAVE_HISTORY=1 -DHAVE_LIBBSD=1 -DHAVE_LIBM=1 -DHAVE_LIBDL=1 -DHAVE_LIBNSL=1 -DHAVE_LIBCRYPT=1 -DHAVE_LIBNSL=1 -DHAVE_LIBDL=1 -DHAVE_LIBCRYPT=1  -DCERNLIB_TYPE
CPPFLAGS =
DEPCPPFLAGS = -MM

FFLAGS = -O -g -DBELLE_TARGET_H=\"belle-x86_64-unknown-linux-gnu-g++.h\"
CFLAGS = -O -g -DBELLE_TARGET_H=\"belle-x86_64-unknown-linux-gnu-g++.h\" -fPIC
CXXFLAGS = -g -DHEP_SHORT_NAMES -DBELLE_SHORT_NAMES -DDSTXX_NOINLINE -DBELLE_TARGET_H=\"belle-x86_64-unknown-linux-gnu-g++.h\" -fPIC
SOFLAGS = -shared -Wl,-export-dynamic
LDFLAGS =

SYSLIB = -lcrypt

CLHEPLIB = -lbelleCLHEP



# Include directories

INCLUDES_C = $(MY_TOP_DIR)/include $(BELLE_TOP_DIR)/include
INCLUDES_FORTRAN = $(BELLE_TOP_DIR)/inc

# Dependence description

include $(BELLE_RUN_DIR)/src/config/Makefile.panther

COMPILE_FCPP := $(FC) -c $(PANTHER_FMACROS) $(INCLUDES_FORTRAN:%=-I%) $(CPPFLAGS) $(FFLAGS)
COMPILE_FC := $(FC) -c $(INCLUDES_FORTRAN:%=-I%) $(FFLAGS)
COMPILE_CC := $(CC) -c $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CFLAGS)
COMPILE_CXX := $(CXX) -c $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CXXFLAGS)

LINK_FCPP := $(FC)
LINK_FC := $(FC)
LINK_CC := $(CC)
LINK_CXX := $(CXX)

DEPEND_FCPP :=  -M $(DEFS) $(PANTHER_FMACROS) $(INCLUDES_FORTRAN:%=-I%) $(CPPFLAGS) $(FFLAGS)
DEPEND_CC := $(CC) $(DEPCPPFLAGS) $(DEFS) $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CFLAGS)
DEPEND_CXX := $(CXX) $(DEPCPPFLAGS) $(DEFS) $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CXXFLAGS)

LINUX_G77_BUG = @LINUX_G77_BUG@

%.o:%.c
	$(COMPILE_CC) $<

%.d:%.c
	$(SHELL) -ec '$(DEPEND_CC) $< | sed -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.h/.tdf/g' > $@'

%.o:%.cc
	$(COMPILE_CXX) $<

%.d:%.cc
	$(SHELL) -ec '$(DEPEND_CXX) $< | sed -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.h/.tdf/g'> $@'


%.o:%.F
	$(COMPILE_FCPP) $<

%.d:%.F
	$(SHELL) -ec '$(DEPEND_FCPP) $< | sed $(LINUX_G77_BUG) -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.inc/.tdf/g'> $@'


# Dependencies
all::	$(MODULE1) $(MODULE2) $(MODULE3) $(MODULE4) $(MODULE5) $(MODULE6)

1::	$(MODULE1)
2::	$(MODULE2)
3::	$(MODULE3)
4::	$(MODULE4)
5::	$(MODULE5)
6::	$(MODULE6)

$(OBJS1): $(DEPFILES1)
$(OBJS2): $(DEPFILES2)
$(OBJS3): $(DEPFILES3)
$(OBJS4): $(DEPFILES4)
$(OBJS5): $(DEPFILES5)
$(OBJS6): $(DEPFILES6)

$(MODULE1): $(OBJS1) $(DEPFILES1)
	$(LINK_CXX) -o $(MODULE1) $(SOFLAGS) \
	$(OBJS1) $(LIBS1) $(CLHEPLIB) $(CERNLIB) $(SYSLIB)

$(MODULE2): $(OBJS2) $(DEPFILES2)
	$(LINK_CXX) -o $(MODULE2) $(SOFLAGS) \
	$(OBJS2) $(LIBS2) $(CLHEPLIB) $(CERNLIB) $(SYSLIB)

$(MODULE3): $(OBJS3) $(DEPFILES3)
	$(LINK_CXX) -o $(MODULE3) $(SOFLAGS) \
	$(OBJS3) $(LIBS3) $(CLHEPLIB) $(CERNLIB) $(SYSLIB)

$(MODULE4): $(OBJS4) $(DEPFILES4)
	$(LINK_CXX) -o $(MODULE4) $(SOFLAGS) \
	$(OBJS4) $(LIBS4) $(CLHEPLIB) $(CERNLIB) $(SYSLIB)

$(MODULE5): $(OBJS5) $(DEPFILES5)
	$(LINK_CXX) -o $(MODULE5) $(SOFLAGS) \
	$(OBJS5) $(LIBS5) $(CLHEPLIB) $(CERNLIB) $(SYSLIB)

$(MODULE6): $(OBJS6) $(DEPFILES6)
	$(LINK_CXX) -o $(MODULE6) $(SOFLAGS) \
	$(OBJS6) $(LIBS6) $(CLHEPLIB) $(CERNLIB) $(SYSLIB)

$(MODULE7): $(OBJS7) $(DEPFILES7)
	$(LINK_CXX) -o $(MODULE7) $(SOFLAGS) \
	$(OBJS7) $(LIBS7) $(CLHEPLIB) $(CERNLIB) $(SYSLIB)

clean::
	rm -f $(MODULE1) $(OBJS1) \
		$(MODULE2) $(OBJS2) \
		$(MODULE3) $(OBJS3) \
		$(MODULE4) $(OBJS4) \
		$(MODULE5) $(OBJS5) \
		$(MODULE6) $(OBJS6) \
		$(MODULE7) $(OBJS7) \
		fort.* fpda_pid.* core core.*

