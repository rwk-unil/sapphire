PP_MAIN_ABSPATH := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PP_MAIN_RELPATH := $(dir $(lastword $(MAKEFILE_LIST)))

CPP_SOURCES := $(wildcard *.cpp)
CPP_OBJS := $(CPP_SOURCES:.cpp=.o)
DEPENDENCIES := $(CPP_SOURCES:.cpp=.d)

XSQUEEZEITPATH := $(PP_MAIN_RELPATH)/xSqueezeIt/

HTSLIB_PATH := $(XSQUEEZEITPATH)/htslib
ZSTD_PATH := $(XSQUEEZEITPATH)/zstd/lib

ifeq ($(ADD_EXTRA),y)
EXTRA_FLAGS := -fsanitize=address -fsanitize=undefined -fsanitize=pointer-subtract -fsanitize=pointer-compare -fno-omit-frame-pointer -fstack-protector-all -fcf-protection
endif

ifeq ($(OLEVEL),)
OLEVEL := 3
endif

# Use g++
CXX=g++
CC=g++

INCLUDE_DIRS=-I include -I $(PP_MAIN_RELPATH)/include -I $(HTSLIB_PATH) -I $(HTSLIB_PATH)/htslib -I $(ZSTD_PATH) -I ${XSQUEEZEITPATH}/include
CXXFLAGS=-O$(OLEVEL) -g -Wall -std=c++17 $(INCLUDE_DIRS) $(CXXEXTRAFLAGS) $(EXTRA_FLAGS)

ifeq ($(STATIC_BINS),y)
A_LIBS := $(HTSLIB_PATH)/libhts.a $(ZSTD_PATH)/libzstd.a
LDLIBS+=-llzma -lbz2 -lz -lm -lcurl -lcrypto -ldeflate -pthread
else
LDLIBS+=-lhts -lzstd -lcrypto -ldeflate -pthread
endif

LDFLAGS+=-O$(OLEVEL) $(EXTRA_FLAGS)
