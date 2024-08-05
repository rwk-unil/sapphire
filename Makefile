all : .dependencies_ready binaries

.dependencies_ready :
	./setup.sh

binaries :
	make -C bin_tools
	make -C pp_extractor
	make -C phase_caller
	make -C pp_update

clean :
	make -C bin_tools clean
	make -C pp_extractor clean
	make -C phase_caller clean
	make -C pp_update clean

.PHONY : all binaries clean