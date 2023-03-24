XSQUEEZEITPATH = ../xSqueezeIt

HTSLIB_PATH = ${XSQUEEZEITPATH}/htslib
ZSTD_PATH := ${XSQUEEZEITPATH}/zstd/lib

# C++ Compiler
CXX=g++
INCLUDE_DIRS=-I include -I ../include -I $(HTSLIB_PATH)/htslib -I $(ZSTD_PATH) -I ${XSQUEEZEITPATH}/include
#EXTRA_FLAGS=-fsanitize=address -fsanitize=undefined -fsanitize=pointer-subtract -fsanitize=pointer-compare -fno-omit-frame-pointer -fstack-protector-all -fcf-protection
ifeq ($(ADD_EXTRA),y)
EXTRA_FLAGS=-fsanitize=address -fsanitize=undefined -fsanitize=pointer-subtract -fsanitize=pointer-compare -fno-omit-frame-pointer -fstack-protector-all -fcf-protection
endif
ifeq ($(OLEVEL),)
OLEVEL=3
endif
CXXFLAGS=-O$(OLEVEL) -g -Wall -std=c++17 $(INCLUDE_DIRS) $(CXXEXTRAFLAGS) $(EXTRA_FLAGS)
# Linker
CC=g++
LDLIBS+=-lhts -lzstd
LDFLAGS+=-O$(OLEVEL) $(EXTRA_FLAGS) -L $(HTSLIB_PATH) -L $(ZSTD_PATH)

# Project specific :
TARGETS := analyze_bin bin_diff bin_splitter bin_merger
XOBJS := ${XSQUEEZEITPATH}/xcf.o ${XSQUEEZEITPATH}/bcf_traversal.o
CPP_SOURCES := $(wildcard *.cpp)
CPP_OBJS := $(CPP_SOURCES:.cpp=.o)
DEPENDENCIES := $(CPP_SOURCES:.cpp=.d)

NOT_IN_GIT_REPO := $(shell git rev-parse --short HEAD >/dev/null 2>/dev/null; echo $$?)
ifneq ($(NOT_IN_GIT_REPO),0)
	GEN_GIT_REV := gen_git_rev
	GIT_REVISION := 0
else
	GIT_REVISION := $(shell git rev-parse --short HEAD)
	REV_NOT_UP_TO_DATE := $(shell grep $(GIT_REVISION) include/git_rev.h  >/dev/null 2>/dev/null; echo $$?)

	ifneq ($(REV_NOT_UP_TO_DATE),0)
		GEN_GIT_REV := gen_git_rev
	endif
endif

# Rules
all : $(GEN_GIT_REV) $(DEPENDENCIES) $(TARGETS)

datetime :
	date | figlet 2> /dev/null

gen_git_rev :
	echo "#define GIT_REVISION 0x$(GIT_REVISION)" > include/git_rev.h

# Link the targets (use implicit built-in rule)
$(TARGETS) : % : $(XOBJS) %.o

# Do not include the depency rules for "clean"
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPENDENCIES)
endif

# Rule to generate the dependency files
%.d : %.cpp
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@
%.d : %.c
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@

# Remove artifacts
clean :
	rm -f $(TARGETS) $(DEPENDENCIES) $(CPP_OBJS)

# Rules that don't generate artifacts
.PHONY : all clean datetime gen_git_rev