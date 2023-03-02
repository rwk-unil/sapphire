XSQUEEZEITPATH = ../xSqueezeIt

HTSLIB_PATH = ${XSQUEEZEITPATH}/htslib
ZSTD_PATH := ${XSQUEEZEITPATH}/zstd/lib

# C++ Compiler
CXX=g++
INCLUDE_DIRS=-I include -I $(HTSLIB_PATH)/htslib -I $(ZSTD_PATH) -I ${XSQUEEZEITPATH}/include
#EXTRA_FLAGS=-fsanitize=address -fsanitize=undefined -fsanitize=pointer-subtract -fsanitize=pointer-compare -fno-omit-frame-pointer -fstack-protector-all -fcf-protection
ifeq ($(ADD_EXTRA),y)
EXTRA_FLAGS=-fsanitize=address -fsanitize=undefined -fsanitize=pointer-subtract -fsanitize=pointer-compare -fno-omit-frame-pointer -fstack-protector-all -fcf-protection
endif
ifeq ($(OLEVEL),)
OLEVEL=3
endif
CXXFLAGS=-O$(OLEVEL) -g -Wall -std=c++17 $(INCLUDE_DIRS) $(CXXEXTRAFLAGS) $(EXTRA_FLAGS)
# Linker
LD=g++
LIBS=-lpthread -lhts -lzstd
LDFLAGS=-O$(OLEVEL) $(EXTRA_FLAGS) -L $(HTSLIB_PATH) -L $(ZSTD_PATH)

# Project specific :
TARGET := pp_extract
#OBJS := xcf.o bcf_traversal.o accessor.o c_api.o xsi_mixed_vcf.o $(OBJ)
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
all : $(GEN_GIT_REV) $(DEPENDENCIES) $(TARGET)

datetime :
	date | figlet 2> /dev/null

gen_git_rev :
	echo "#define GIT_REVISION 0x$(GIT_REVISION)" > include/git_rev.h

# Link the target
$(TARGET) : $(XOBJS) pp_extractor.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

pp_show : $(XOBJS) pp_show.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

# Do not include the depency rules for "clean"
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPENDENCIES)
endif

# Compile
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
%.o : %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to generate the dependency files
%.d : %.cpp
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@
%.d : %.c
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@

# Remove artifacts
clean :
	rm -f $(TARGET) $(DEPENDENCIES) $(CPP_OBJS) pp_show

# Rules that don't generate artifacts
.PHONY : all clean