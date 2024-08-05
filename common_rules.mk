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

BINARY_PATH=../bin
TARGET_BINARIES=$(foreach binary,$(TARGETS),$(BINARY_PATH)/$(binary))

# Rules
all : $(GEN_GIT_REV) $(DEPENDENCIES) $(TARGETS) $(TARGET_BINARIES)
	$(info Finished building $(TARGETS))

gen_git_rev :
	echo "#define GIT_REVISION 0x$(GIT_REVISION)" > include/git_rev.h

# Link the targets (use implicit built-in rule)
$(TARGETS) : % : $(XOBJS) %.o $(A_LIBS)

# Copy the binaries to the binary directory
$(BINARY_PATH)/% : %
	cp $< $@

# Do not include the depency rules for "clean"
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPENDENCIES)
endif

# Rules to generate the dependency files
%.d : %.cpp
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@
%.d : %.c
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@

# Remove artifacts
clean :
	rm -f $(TARGETS) $(DEPENDENCIES) $(CPP_OBJS)

# Rules that don't generate artifacts
.PHONY : all clean gen_git_rev
