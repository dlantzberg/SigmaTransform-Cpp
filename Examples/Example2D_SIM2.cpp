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

int main( int argc , char** argv ) {
    try {
        // Chronometer, for benchmarking purposes
        sigma::Chronometer    Chrono;

        // input vector
        int x, y;
        // load lena
        cxVec lena = sigma::loadAscii2D( "Signals/lena.asc", x, y );
        sigma::point<2>  sz( std::array<double,2>{ (double)x , (double)y } ), Fs( sz ),
			             num_chans(std::array<double,2>{7,17});

        // start chronometer
        Chrono.tic();

        // make rectangular window in warped (polar) Fourier domain in (0,1] x (0,pi/16]
        auto rect2D = []( sigma::point<2> const& p ) {
            return (0<p[0] && p[0]<=1) && (0<=p[1] && p[1]<(M_PI/16));
        };

        // make channels
        std::vector<sigma::point<2>> grid = sigma::meshgridN<2>( std::array<std::vector<double>,2>{
                                                                    sigma::linspace( log2(1), log2(64)*1.1, 7 ) ,
                                                                    sigma::linspace( -M_PI/2*1.1, M_PI/2*1.1,  17 ) } );

        //construct 2D WaveletTransform transform
        sigma::SIM2D    SIM2T(
            sigma::point<2>(2) ,   // window (point<N>->cmpx (NULL for warped Gaussian))
            Fs ,       // spatial/temporal sampling rate  ( point<N> )
            sz ,       // signal length ( point<N> )
            grid       // shift-steps in warped Fourier domain ( point<N>-vec )
        );
        
        // stop and restart chronometer
        Chrono.toc("construct").tic();

        // analyze
        SIM2T( lena ); // == SIM2T.analyze( lena );
        Chrono.toc("analyze").tic();

        // synthesize
        SIM2T.synthesize();
        Chrono.toc("synthesize").tic();

        // save coefficients
        // sigma::save2file_bin( "lena_coeff.bin", SIM2T.getCoeffs() );
        // Chrono.toc("saveCoeffs").tic();

        // save windows
        //sigma::save2file_bin( "lena_windows.bin", SIM2T.getWindows() );
        //Chrono.toc("saveWindows").tic();

        // save reconstruction
        //sigma::save2file( "lena_rec.bin" , SIM2T.getReconstruction() , sz );
        //Chrono.toc("saveRecon");
    } catch( std::exception &e ) {
        // error?
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
