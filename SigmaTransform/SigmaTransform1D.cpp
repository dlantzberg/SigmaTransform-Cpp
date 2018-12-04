
#include <math.h>

//#include "SigmaTransform_utilities.h"
#include "SigmaTransform1D.h"

//using namespace std;

namespace SigmaTransform {

    SigmaTransform1D::SigmaTransform1D(
        diffFunc<1> sigma, winFunc<1> window, const point<1> &Fs, const point<1> &size, const std::vector<point<1>> &steps, actFunc<1> act )
        : SigmaTransform<1>( sigma, window, Fs , size , steps , act ) { }

    SigmaTransform1D::SigmaTransform1D(
        diffFunc<1> sigma, const point<1> &width, const point<1> &Fs, const point<1> &size, const std::vector<point<1>> &steps, actFunc<1> act )
        : SigmaTransform<1>( sigma, width, Fs , size , steps , act ) { }

        
    STFT1D::STFT1D( winFunc<1> window, const point<1> &Fs, const point<1> &size, const std::vector<point<1>> &steps )
        : SigmaTransform<1>( id<1>, window, Fs , size , steps ) { }
        
    STFT1D::STFT1D( const point<1> &width, const point<1> &Fs, const point<1> &size, const std::vector<point<1>> &steps )
        : SigmaTransform<1>( id<1>, width, Fs , size , steps ) { }

        
    WaveletTransform1D::WaveletTransform1D( winFunc<1> window, const point<1> &Fs, const point<1> &size, const std::vector<point<1>> &steps )
        : SigmaTransform<1>( logpos<1>, window, Fs , size , steps ) { }
 
    WaveletTransform1D::WaveletTransform1D(const point<1> &width, const point<1> &Fs, const point<1> &size, const std::vector<point<1>> &steps )
        : SigmaTransform<1>( logpos<1>, width, Fs , size , steps ) { }
 
 
    CQTransform1D::CQTransform1D( winFunc<1> window, const point<1> &Fs, const point<1> &size, const std::vector<point<1>> &steps )
        : SigmaTransform<1>( cq, window, Fs , size , steps ) { }
 
    CQTransform1D::CQTransform1D( const point<1> &width, const point<1> &Fs, const point<1> &size, const std::vector<point<1>> &steps )
        : SigmaTransform<1>( cq, width, Fs , size , steps ) { }

        
    double CQTransform1D::cq( const point<1> &x ) { return Q * log2( abs( x[0] / f_0 ) + 1E-16 ); };

} // namespace SigmaTransform
