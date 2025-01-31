#include "CLI11.hpp"
#include "het_info_loader.hpp"

int main(int argc, char**argv) {
    CLI::App app{"Binary file switch utility"};
    std::string bin_fname = "-";
    app.add_option("-b,--input", bin_fname, "Binary file to switch");

    CLI11_PARSE(app, argc, argv);

    if (bin_fname.compare("-") == 0) {
        std::cerr << "Requires binary filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    {
        HetInfoMemoryMap himm(bin_fname, PROT_READ | PROT_WRITE);
        if (!himm.integrity_check_pass()) {
            std::cerr << "File " << bin_fname << " doesn't pass integrity checks" << std::endl;
        }

        /* If the program crashes here the file will be corrupted as switching
           is done "in-place", just a warning... */
        /** @todo maybe switch phase in new file (requires twice the space...) */
        for (uint32_t i = 0; i < himm.num_samples; ++i) {
            HetInfoMemoryMap::HetInfoPtrContainer hipc(himm, i);
            hipc.switch_phase();
        }
    } /* File is synced/closed here in the destructor of HetInfoMemoryMap */

    std::cout << "Done switching phase of file " << bin_fname << std::endl;
    return 0;
}
