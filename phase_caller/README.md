# Phase Caller

This is a phase caller tool that rephases data given sequencing reads (BAM/CRAM).


## Machinery

- The `PhaseCaller` class will orchestrate the rephasing of samples
- It has a `rephase_orchestrator_multi_thread()` function that will launch a thread per sample
- The threaded function will call `rephase_sample()`
- `rephase_sample()` will get het genotypes from memory mapped file, create a linked list (trios)
- Then a `Rephaser` class is instanciated to rephase the "trios"
- The `Rephaser` does the following
        - Pile-up reads for each and every SNV
        - Go through the trios and rephase a low phased het genotype according to its neighbors