#ifndef SIGMATRANSFORM2D_H
#define SIGMATRANSFORM2D_H

#include <math.h>

//#include "SigmaTransform_utilities.h"
#include "SigmaTransformN.h"

//using namespace std;

namespace SigmaTransform {
    point<2> polar(const point<2> &p);
    point<2> shear(const point<2> &p);
    point<2> parabolicAction(const point<2> &l ,const point<2> &r);
    
    class SigmaTransform2D : public SigmaTransform<2> {
        public: 
        SigmaTransform2D(
            diffFunc<2> sigma = NULL,
            winFunc<2> window = NULL,
            const point<2> &Fs = point<2>(0),
            const point<2> &size = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0),
            actFunc<2> act = minus<2>
        );
        SigmaTransform2D(
            diffFunc<2> sigma = NULL,
            const point<2> &width = 0,
            const point<2> &Fs = point<2>(0),
            const point<2> &size = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0),
            actFunc<2> act = minus<2>
        );
    };

    class STFT2D : public SigmaTransform<2> {
        public: STFT2D(
            winFunc<2> window = NULL,
            const point<2> &Fs = point<2>(0),
            const point<2> &size = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0)
        );
        STFT2D(
            const point<2> &width,
            const point<2> &Fs = point<2>(0),
            const point<2> &size = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0)
        );
    };

    class WaveletTransform2D : public SigmaTransform<2> {
        
        public: 
        WaveletTransform2D(
            winFunc<2> window = NULL,
            const point<2> &Fs = point<2>(0),
            const point<2> &size = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0)
        );
        WaveletTransform2D(
            const point<2> &width,
            const point<2> &Fs = point<2>(0),
            const point<2> &size = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0)
        );
    };

    class SIM2D : public SigmaTransform<2> { 
    
        public: 
        SIM2D(  
            winFunc<2> window = NULL,
            const point<2> &Fs = point<2>(0), 
            const point<2> &size = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0) 
        );
        SIM2D(  
            const point<2> &width,
            const point<2> &Fs = point<2>(0), 
            const point<2> &size = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0) 
        );
    };

    class Curvelet2D : public SigmaTransform<2> {
        
        public: 
        Curvelet2D(
            winFunc<2> window       = NULL,
            const point<2> &Fs      = point<2>(0),
            const point<2> &size    = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0)
        );
        Curvelet2D(
            const point<2> &width,
            const point<2> &Fs      = point<2>(0),
            const point<2> &size    = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0)
        );
    };

    class NPShearlet2D : public SigmaTransform<2> {
        
        public: 
        NPShearlet2D(
            winFunc<2> window       = NULL,
            const point<2> &Fs      = point<2>(0),
            const point<2> &size    = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0)
        );
        NPShearlet2D(
            const point<2> &width,
            const point<2> &Fs      = point<2>(0),
            const point<2> &size    = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0)
        );
    };
    
    class Shearlet2D : public SigmaTransform<2> {
        
        public: 
        Shearlet2D(
            winFunc<2> window       = NULL,
            const point<2> &Fs      = point<2>(0),
            const point<2> &size    = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0)
        );
        Shearlet2D(
            const point<2> &width,
            const point<2> &Fs      = point<2>(0),
            const point<2> &size    = point<2>(0),
            const std::vector<point<2>> &steps = std::vector<point<2>>(0)
        );
    };
} // namespace SigmaTransform

#endif //SIGMATRANSFORM1D_H
