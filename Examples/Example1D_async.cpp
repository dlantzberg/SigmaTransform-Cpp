// for std::cerr,std::cout
#include <iostream>
// for std::vector
#include <vector>
// for std::complex
#include <complex>
// for std::chrono
#include <chrono>
// for std::thread
#include <thread>
// the class-templace
#include "SigmaTransformN.h"

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
        std::cout   <<"Async STFTransform with " << std::fixed
                    << (int)(numsteps * len / 1000 / 1000) << " million coefficients.\n";

        std::cout<<"\n-------------------- starting to construct transform ------------------\n\n";

        // start chronometer
        Chrono.tic();

        // create STFT
        sigma::SigmaTransform<1> sig(
            NULL,                     // diffeomorphism ( point<N> -> point<N> ; (NULL for identity))
            (sigma::winFunc<1>)NULL,  // window ( point<N> -> cmplx ; (NULL for warped Gaussian))
            Fs ,                      // spatial/temporal sampling rate  ( point<N> )
            len ,                     // signal length ( point<N> )
            {numsteps}                // Number of steps in warped Fourier domain ( point<N> )
                                      // or: sampling points in warped Fourier domain ( point<N>-vector )
            //sigma::meshgridN<1>({sigma::linspace( -Fs/2, Fs/2, 400 )})
        );
        // stop and restart chronometer
        Chrono.toc("Done constructing").tic();

        // analyze async
        std::cout<<"\n------------------ starting to analyze asynchronously -----------------\n\n";
        sig.analyze( bat_signal, [&](sigma::SigmaTransform<1>* obj){
            Chrono.toc("Done analyzing").tic();
        } );

        // synthesize async
        std::cout<<"\n---------------- starting to synthesize asynchronously ----------------\n\n";
        sig.synthesize( [&](sigma::SigmaTransform<1>* obj){Chrono.toc("Done synthesizing").tic();} );

        // do some other stuff
        std::cout << "-- Doing some other stuff" << std::flush;
        for( int k = 0 ; k < 39 ; ++k ) {
            std::this_thread::sleep_for( std::chrono::milliseconds(25) );
            std::cout << "." << std::flush;
        } std::cout << "done --\n";

        // wait for asynchronous jobs to join
        std::cout<<"\n----------------------- waiting for threads to join -------------------\n\n";
        sig.join();

        // save
        std::cout << "\nAll joined. Saving to file." << std::endl;
        sigma::save2file_bin("bat_out_async.bin",sig.getReconstruction() );
    } catch( std::exception &e ) {
        // error?
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
