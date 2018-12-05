
//#include <math.h>

#include "SigmaTransform2D.h"

namespace SigmaTransform {

    point<2> polar(const point<2> &p) {
        //return std::move( point<2>({ log2( p.sq().sum()+1E-16 )/2.0, atan2( p[1] , p[0] ) }) );
        return std::move( point<2>(std::array<double,2>{ log2( p.sq().sum()+1E-16 )/2.0, atan( p[1] / p[0] ) }) );
    }

    point<2> shear(const point<2> &p) {
        return std::move( point<2>(std::array<double,2>{ log2( abs(p[0]) + 1E-16 ), p[1] / p[0] }) );
    }

    point<2> parabolicAction(const point<2> &l ,const point<2> &r ) {
        return std::move( point<2>(std::array<double,2>{ l[0] - r[0] , exp( -r[0] / 2.0 * log(2) )  * ( l[1] - r[1] ) }) );
    }

    SigmaTransform2D::SigmaTransform2D(
        diffFunc<2> sigma, winFunc<2> window, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps, actFunc<2> act )
        : SigmaTransform<2>( sigma, window, Fs , size , steps , act ) { }

    SigmaTransform2D::SigmaTransform2D(
        diffFunc<2> sigma, const point<2> &width, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps, actFunc<2> act )
        : SigmaTransform<2>( sigma, width, Fs , size , steps , act ) { }


    STFT2D::STFT2D( winFunc<2> window, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( id<2>, window, Fs , size , steps ) { }

    STFT2D::STFT2D( const point<2> &width, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( id<2>, width, Fs , size , steps ) { }


    WaveletTransform2D::WaveletTransform2D( winFunc<2> window, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( logabs<2>, window, Fs , size , steps ) { }

    WaveletTransform2D::WaveletTransform2D( const point<2> &width, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( logabs<2>, width, Fs , size , steps ) { }


    SIM2D::SIM2D( winFunc<2> window, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( polar, window, Fs , size , steps ) { }

    SIM2D::SIM2D( const point<2> &width, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( polar, width, Fs , size , steps ) { }


    Curvelet2D::Curvelet2D( winFunc<2> window, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( polar, window, Fs , size , steps , parabolicAction ) { }

    Curvelet2D::Curvelet2D( const point<2> &width, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( polar, width, Fs , size , steps , parabolicAction ) { }


    NPShearlet2D::NPShearlet2D( winFunc<2> window, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( shear, window, Fs , size , steps ) { }

    NPShearlet2D::NPShearlet2D( const point<2> &width, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( shear, width, Fs , size , steps ) { }


    Shearlet2D::Shearlet2D( winFunc<2> window, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( shear, window, Fs , size , steps , parabolicAction ) { }

    Shearlet2D::Shearlet2D( const point<2> &width, const point<2> &Fs, const point<2> &size, const std::vector<point<2>> &steps )
        : SigmaTransform<2>( shear, width, Fs , size , steps , parabolicAction ) { }

} // namespace SigmaTransform
