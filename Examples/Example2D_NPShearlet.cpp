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

        // input vector
        cxVec lena      = sigma::loadAscii2D( "Signals/lena.asc", x, y );
        sigma::point<2> sz( std::array<double,2>{ (double)x , (double)y } ), Fs( sz );

        // start chronometer
        Chrono.tic();

        // make rectangular window in warped (polar) Fourier domain in (0,2] x (-1/2,1/2]
        auto rect2D = []( sigma::point<2> const& p ) {
            return (0<p[0] && p[0]<=2) && (-0.5<p[1] && p[1]<=0.5);
        };

        // make channels
        std::vector<sigma::point<2>> grid = sigma::meshgridN<2>(std::array<std::vector<double>,2>{ sigma::linspace(  0, 8, 5 ) ,
                                                                                                   sigma::linspace( -4, 4, 9 ) });
        //construct 2D WaveletTransform transform
        sigma::NPShearlet2D    NPShearlet2D(
            rect2D ,   // window (point<N>->cmpx (NULL for warped Gaussian))
            Fs ,       // spatial/temporal sampling rate  ( point<N> )
            sz ,       // signal length ( point<N> )
            grid       // shift-steps in warped Fourier domain ( point<N>-vec )
        );

        // stop and restart chronometer
        Chrono.toc("construct").tic();

        // analyze
        NPShearlet2D( lena ); // == NPShearlet2D.analyze( lena );
        Chrono.toc("analyze").tic();

        // synthesize
        NPShearlet2D.synthesize();
        Chrono.toc("synthesize").tic();

        // save coefficients
        // sigma::save2file_bin( "lena_coeff.bin", NPShearlet2D.getCoeffs() );
        // Chrono.toc("saveCoeffs").tic();

        // save windows
        //sigma::save2file_bin( "lena_windows.bin", NPShearlet2D.getWindows() );
        //Chrono.toc("saveWindows").tic();

        // save reconstruction
        //sigma::save2file( "lena_rec.bin" , NPShearlet2D.getReconstruction() , sz );
        //Chrono.toc("saveRecon");

        /************************************************************************/

    } catch( std::exception &e ) {
        // error?
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
