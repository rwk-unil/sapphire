HTSLIB_PATH = ../xSqueezeIt/htslib

# C++ Compiler
CXX=g++
INCLUDE_DIRS=-I include -I ../include -I $(HTSLIB_PATH)/htslib
#EXTRA_FLAGS=-fsanitize=address -fsanitize=undefined -fsanitize=pointer-subtract -fsanitize=pointer-compare -fno-omit-frame-pointer -fstack-protector-all -fcf-protection
ifeq ($(ADD_EXTRA),y)
EXTRA_FLAGS=-fsanitize=address -fsanitize=undefined -fsanitize=pointer-subtract -fsanitize=pointer-compare -fno-omit-frame-pointer -fstack-protector-all -fcf-protection
endif
ifeq ($(OLEVEL),)
OLEVEL=3
endif
CXXFLAGS=-O$(OLEVEL) -g -Wall -pthread -std=c++17 $(INCLUDE_DIRS) $(CXXEXTRAFLAGS) $(EXTRA_FLAGS)
# Linker
LD=g++
LIBS=-pthread -lhts
LDFLAGS=-O$(OLEVEL) $(EXTRA_FLAGS) -L $(HTSLIB_PATH)

# Project specific :
TARGET := cram_accessor
SOURCE := main.cpp
OBJ := $(SOURCE:.cpp=.o)
CPP_SOURCES := $(wildcard *.cpp)
CPP_OBJS := $(CPP_SOURCES:.cpp=.o)
CPP_OBJS := $(CPP_OBJS:.c=.o)
OBJS = $(OBJ) $(CPP_OBJS)
DEPENDENCIES := $(CPP_SOURCES:.cpp=.d)
DEPENDENCIES := $(DEPENDENCIES:.c=.d)

# Rules
all : $(DEPENDENCIES) $(TARGET)

# Link the target
$(TARGET) : $(OBJS) $(XOBJS)
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
	rm -f $(OBJS) $(TARGET) $(DEPENDENCIES)

# Rules that don't generate artifacts
.PHONY : all clean