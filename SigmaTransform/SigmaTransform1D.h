#ifndef SIGMATRANSFORM1D_H
#define SIGMATRANSFORM1D_H

#include <math.h>

//#include "SigmaTransform_utilities.h"
#include "SigmaTransformN.h"

//using namespace std;

namespace SigmaTransform {

    class SigmaTransform1D : public SigmaTransform<1> {
        public: 
        SigmaTransform1D(
            diffFunc<1> sigma = NULL,
            winFunc<1> window = NULL,
            const point<1> &Fs = point<1>(0),
            const point<1> &size = point<1>(0),
            const std::vector<point<1>> &steps = std::vector<point<1>>(0),
            actFunc<1> act = minus<1>
        ); 
        SigmaTransform1D(
            diffFunc<1> sigma = NULL,
            const point<1> &width = 0,
            const point<1> &Fs = point<1>(0),
            const point<1> &size = point<1>(0),
            const std::vector<point<1>> &steps = std::vector<point<1>>(0),
            actFunc<1> act = minus<1>
        );
    };
 
    class STFT1D : public SigmaTransform<1> {
        public: STFT1D(
            winFunc<1> window = NULL,
            const point<1> &Fs = point<1>(0),
            const point<1> &size = point<1>(0),
            const std::vector<point<1>> &steps = std::vector<point<1>>(0)
        );
        STFT1D(
            const point<1> &width,
            const point<1> &Fs = point<1>(0),
            const point<1> &size = point<1>(0),
            const std::vector<point<1>> &steps = std::vector<point<1>>(0)
        );
    };

    class WaveletTransform1D : public SigmaTransform<1> {
        public: WaveletTransform1D(
            winFunc<1> window = NULL,
            const point<1> &Fs = point<1>(0),
            const point<1> &size = point<1>(0),
            const std::vector<point<1>> &steps = std::vector<point<1>>(0)
        );
        WaveletTransform1D(
            const point<1> &width,
            const point<1> &Fs = point<1>(0),
            const point<1> &size = point<1>(0),
            const std::vector<point<1>> &steps = std::vector<point<1>>(0)
        );
    };

    class CQTransform1D : public SigmaTransform<1> {
        static const int Q = 8, f_0 = 1;
        static double cq( const point<1> &x );
        
        public: CQTransform1D(
            winFunc<1> window = NULL,
            const point<1> &Fs = point<1>(0),
            const point<1> &size = point<1>(0),
            const std::vector<point<1>> &steps = std::vector<point<1>>(0)
        );
        CQTransform1D(
            const point<1> &width,
            const point<1> &Fs = point<1>(0),
            const point<1> &size = point<1>(0),
            const std::vector<point<1>> &steps = std::vector<point<1>>(0)
        );
    };

} // namespace SigmaTransform

#endif //SIGMATRANSFORM1D_H
