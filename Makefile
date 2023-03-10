all :
	make -C pp_extractor
	make -C phase_caller
	make -C pp_update

clean :
	make -C pp_extractor clean
	make -C phase_caller clean
	make -C pp_update clean

.PHONY : all clean