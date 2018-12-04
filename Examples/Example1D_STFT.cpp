// for std::cerr,std::cout
#include <iostream>
// for std::vector
#include <vector>
// for std::complex
#include <complex>
// for to_string-conversion
#include <sstream>
// the class-templace
#include "SigmaTransformN.h"
// specific implementations, like STFT, WaveletTransform, etc.
#include "SigmaTransform1D.h"

namespace sigma  = SigmaTransform;
using cxVec = std::vector<std::complex<double>>;

int main( int argc, char** argv ) {
    try {
        // Chronometer, for benchmarking purposes
        sigma::Chronometer    Chrono;

        // start chronometer
        Chrono.tic();

        // load bat signal
        cxVec bat_signal = sigma::loadAscii1D( "Signals/bat.asc" );

        // setup
        double Fs = 143000, len = bat_signal.size(), numsteps = len * 200;
        std::vector<sigma::point<1>> chans = sigma::meshgridN<1>( sigma::linspace( log2(Fs*0.005) , log2(Fs/2*1.1) , numsteps ) );

        //construct 1D STFT transform
        sigma::STFT1D    Stft1D(
            (sigma::point<1>)4.0,          // window or: width (in steps) of a warped Gaussian window
            Fs ,                           // spatial/temporal sampling rate  ( point<N> )
            len ,                          // signal length ( point<N> )
            chans                         // the channels in the warped Fourier domain
        );
        std::cout <<"Wavelet transform with " << std::fixed << (int)(numsteps * len / 1000000) << " million coefficients.\n";

        // stop and restart chronometer
        Chrono.toc("construct").tic();

        // analyze
        Stft1D.analyze( bat_signal );
        Chrono.toc("analyze").tic();

        // save windows
        // sigma::save2file_bin("bat_windows.bin",Stft1D.getWindows() );
        // Chrono.toc("saveWindows").tic();

        // save coefficients
        // sigma::save2file_bin("bat_coeff.bin",Stft1D.getCoeffs() );
        // Chrono.toc("saveCoeffs").tic();

        // mask coefficients
        Stft1D.applyMask( [&]( sigma::point<1>const& x, sigma::point<1>const& step )->sigma::cmpx {
            return ( x > .0005 && x < .002 ) && (step < 16 && step > 12);
        } );
        Chrono.toc("mask").tic();

        // save masked coefficients
        // sigma::save2file_bin("bat_coeff_masked.bin",Stft1D.getCoeffs() );
        // Chrono.toc("saveMaskedCoeffs").tic();

        // synthesize
        Stft1D.synthesize();
        Chrono.toc("synthesize").tic();

        // save reconstruction
        // sigma::save2file_bin("bat_rec.bin",Stft1D.getReconstruction() );
        // Chrono.toc("saveRecon");
    } catch( std::exception &e ) {
        // error?
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
