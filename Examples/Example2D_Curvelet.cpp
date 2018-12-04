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
        // load lena
        int x, y;
        cxVec lena = sigma::loadAscii2D( "Signals/lena.asc", x, y );
        sigma::point<2> sz( std::array<double,2>{ (double)x , (double)y } ), Fs( sz ),
                        num_chans(std::array<double,2>{7,17});

        // start chronometer
        Chrono.tic();

        // make rectangular window in warped (polar) Fourier domain in (0,1] x (0,pi/16]
        auto rect2D = []( sigma::point<2> const& p ) {
            return (0<p[0] && p[0]<=2) && (0<=p[1] && p[1]<(M_PI/32));
        };

        // make channels
        std::vector<sigma::point<2>> grid = sigma::meshgridN<2>(
						std::array<std::vector<double>,2>{ sigma::linspace( log2(1), log2(64), 4 ) ,
                                                           sigma::linspace( -M_PI/2, M_PI/2,  33 ) });

        // adjust step-width to parabolic scaling
        for( auto& st : grid )
            st[1] *= exp( log(2)/2 * st[0] );

        //construct 2D WaveletTransform transform
        sigma::Curvelet2D    CurveT(
            //NULL ,   // window (point<N>->cmpx (NULL for warped Gaussian))
            rect2D ,   // window (point<N>->cmpx (NULL for warped Gaussian))
            Fs ,       // spatial/temporal sampling rate  ( point<N> )
            sz ,       // signal length ( point<N> )
            grid       // shift-steps in warped Fourier domain ( point<N>-vec )
        );

        // stop and restart chronometer
        Chrono.toc("construct").tic();

        // analyze
        CurveT( lena ); // == CurveT.analyze( lena );
        Chrono.toc("analyze").tic();

        // synthesize
        CurveT.synthesize();
        Chrono.toc("synthesize").tic();

        // save coefficients
        sigma::save2file_bin( "lena_coeff.bin", CurveT.getCoeffs() );
        Chrono.toc("saveCoeffs").tic();

        // save windows
        sigma::save2file_bin( "lena_windows.bin", CurveT.getWindows() );
        Chrono.toc("saveWindows").tic();

        // save reconstruction
        //sigma::save2file( "lena_rec.bin" , CurveT.getReconstruction() , sz );
        //Chrono.toc("saveRecon");
    } catch( std::exception &e ) {
        // error?
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
