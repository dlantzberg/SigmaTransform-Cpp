// for std::cerr,std::cout
#include <iostream>
// for std::vector
#include <vector>
// for std::complex
#include <complex>
// the class-template
#include "SigmaTransformN.h"
// specific implementations, like STFT2D, SIM2Transform, ShearletTransform, etc.
#include "SigmaTransform2D.h"

namespace sigma  = SigmaTransform;
using cxVec = std::vector<std::complex<double>>;

int main( int argc, char** argv ) {
    try {
        // Chronometer, for benchmarking purposes
        sigma::Chronometer    Chrono;
        // load lena
        int x, y;
        cxVec lena = sigma::loadAscii2D( "Signals/lena.asc", x, y );
        sigma::point<2>    sz( std::array<double,2>{ (double)x , (double)y } ), Fs( sz ), num_chans(6);

        // start chronometer
        Chrono.tic();

        // make rectangular window in warped Fourier domain in (0,1] x (0,1]
        //auto rect2D = []( sigma::point<2> const& p ) { return p>0 && p <= 1; };

        //construct 2D WaveletTransform transform
        sigma::WaveletTransform2D    WT2D(
            []( sigma::point<2> const& p ) { return p>0 && p<=1; } ,   // window (point<N>->cmpx (NULL for warped Gaussian))
            Fs ,                                                    // spatial/temporal sampling rate  ( point<N> )
            sz ,                                                    // signal length ( point<N> )
            sigma::meshgridN<2>( sigma::linspace( log2(1), log2(32), 6 ) )// shift-steps in warped Fourier domain ( point<N>-vec )
        );

        // stop and restart chronometer
        Chrono.toc("construct").tic();

        // analyze
        WT2D( lena ); // == WT2D.analyze( lena );
        Chrono.toc("analyze").tic();

        // synthesize
        WT2D.synthesize();
        Chrono.toc("synthesize").tic();

        // save coefficients
        // sigma::save2file( "lena_coeff.bin", WT2D.getCoeffs() , sigma::point<3>( { sz[0], sz[1], num_chans.prod() } ) );
        // Chrono.toc("saveCoeffs").tic();

        // save windows
        //sigma::save2file_bin( "lena_windows.bin", WT2D.getWindows() );
        //Chrono.toc("saveWindows").tic();

        // save reconstruction
        //sigma::save2file( "lena_rec.bin" , WT2D.getReconstruction() , sz );
        //Chrono.toc("saveRecon");
    } catch( std::exception &e ) {
        // error?
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
