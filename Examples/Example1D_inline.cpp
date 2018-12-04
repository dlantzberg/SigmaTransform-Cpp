// for std::cerr,std::cout
#include <iostream>
// for std::vector
#include <vector>
// for std::complex
#include <complex>
// the class-templace
#include "SigmaTransformN.h"

namespace st = SigmaTransform;
using cxVec = std::vector<std::complex<double>>;
// Chronometer, for benchmarking purposes
st::Chronometer    Chrono;

int main( int argc, char** argv ) {
    try {
        // start chronometer
        Chrono.tic();
        // perform an inline transform,mask,synthesize and save reconstruction
        st::save2file_bin("bat_out_inline.bin",                                 // save as file
            st::SigmaTransform<1>([](const st::point<1> &p)->st::point<1> {
                return p;
            }, (st::point<1>)0,143000, 400,{400*100})                           // construct (STF) transform
            .analyze(st::loadAscii1D( "Signals/bat.asc" ))                      // analyze
            //.applyMask([](){...} maskFunc)                                    // mask coefficients
            .synthesize()                                                       // synthesize
            .getReconstruction()                                                // return reconstruction
        );                                                                      // destroy
        // stop chronometer
        Chrono.toc("Inline WaveletTrnsf. with 16000000 Coeffs done.");
    } catch( std::exception &e ) {
        // error?
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
