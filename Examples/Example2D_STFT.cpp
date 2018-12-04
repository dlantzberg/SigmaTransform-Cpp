// for std::cerr,std::cout
#include <iostream>
// for std::vector
#include <vector>
// for std::complex
#include <complex>
// the class-templace
#include "SigmaTransformN.h"
// specific implementations, like STFT2D, SIM2Transform, ShearletTransform, etc.
#include "SigmaTransform2D.h"

namespace sigma  = SigmaTransform;
using cxVec = std::vector<std::complex<double>>;

sigma::cmpx rect2D( sigma::point<2> const& p ) {
    if( abs( p[0] ) <= 32 && abs( p[1] ) <= 32 ) {
        return 1;
    } else {
        return 0;
    }
}

int main( int argc , char** argv ) {
    try {
        // Chronometer, for benchmarking purposes
        sigma::Chronometer    Chrono;

        // start chronometer
        Chrono.tic();

        // load bat signal
        sigma::point<2> Fs(128), num_chans(5);
        int x, y;
        cxVec lena = sigma::loadAscii2D( "Signals/lena.asc", x, y );
        sigma::point<2>    sz( std::array<double,2>{ (double)x , (double)y } );

        //construct 2D STFT
        sigma::SigmaTransform<2>    sigT2D(
            NULL,       // diffeomorphism ( point<N> -> point<N> ; (NULL for identity))
            rect2D,     // window ( point<N> -> cmplx ; (NULL for warped Gaussian))
            Fs ,        // spatial/temporal sampling rate  ( point<N> )
            sz ,        // signal length ( point<N> )
            {num_chans} // Number of steps in warped Fourier domain ( point<N> )
                        // or: sampling points in warped Fourier domain ( point<N>-vector )
            //sigma::meshgridN<2>( {sigma::linspace( -Fs/2, Fs/2, 400 )} )
        );
        // stop and restart chronometer
        Chrono.toc("construct").tic();

        // analyze
        sigT2D( lena ); // == sigT2D.analyze( lena );
        Chrono.toc("analyze").tic();

        // synthesize
        sigT2D.synthesize();
        Chrono.toc("synthesize").tic();

        // save coefficients
        // sigma::save2file_bin( "bat_coeff.bin", sigT2D.getCoeffs() );
        // Chrono.toc("saveCoeffs").tic();

        // save reconstruction
        // sigma::save2file_bin( "bat_out.bin" , sigT2D.getReconstruction() );
        // Chrono.toc("saveRecon");
    } catch( std::exception &e ) {
        // error?
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
