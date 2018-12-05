#ifndef SIGMATRANSFORM2D_H
#define SIGMATRANSFORM2D_H

#include <math.h>

#include "SigmaTransformN.h"

namespace SigmaTransform {
    /** Performs a polar diffeomorphism.
    *   @param      a 2D cartesian point
    *
    *   @return     the 2D point in polar-coordinates
    */
    point<2> polar(const point<2> &p);

    /** Performs a shear diffeomorphism.
    *   @param      a 2D cartesian point
    *
    *   @return     the 2D point in "shearing-coordinates"
    */
    point<2> shear(const point<2> &p);

    /** Performs a parabolic shift.
    *   @param   l  2D point
    *   @param   r  2D point
    *
    *   @return     a 2D point, "l" parabolic shifted by "r"
    */
    point<2> parabolicAction(const point<2> &l ,const point<2> &r);

    /** Class SigmaTransform2D is a two-dimensional instantation of SigmaTransform<N>
    *   and can thus, since no more templates are involved, be used in a precompiled form.
    *
    *   See the documentation of SigmaTransform<N> for more information.
    */
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

    /** Class STFT2D is a two-dimensional STFT-instantation of SigmaTransform<N>,
    *   by setting the spectral diffeomorphism's handle to the (2D-)identity.
    *   Since no more templates are involved, it can be used in a precompiled form.
    *
    *   See the documentation of SigmaTransform<N> for more information.
    */
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

    /** Class WaveletTransform2D is a two-dimensional Wavelet-instantation of SigmaTransform<N>,
    *   by setting the spectral diffeomorphism's handle to the (2D-binary) logarithm.
    *   Since no more templates are involved, it can be used in a precompiled form.
    *
    *   See the documentation of SigmaTransform<N> for more information.
    */
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


    /** Class SIM2D is a two-dimensional SIM(2)-instantation of SigmaTransform<N>,
    *   by setting the spectral diffeomorphism's handle to the 2D polar diffeomorphism.
    *   Since no more templates are involved, it can be used in a precompiled form.
    *
    *   See the documentation of SigmaTransform<N> for more information.
    */
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

    /** Class Curvelet2D is a two-dimensional Curvelet-instantation of SigmaTransform<N>,
    *   by setting the spectral diffeomorphism's handle to the 2D polar diffeomorphism
    *   and the "group action"-handle to a "parabolic shift".
    *   Since no more templates are involved, it can be used in a precompiled form.
    *
    *   See the documentation of SigmaTransform<N> for more information.
    */
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

    /** Class NPShearlet2D is a two-dimensional Shearlet-instantation of SigmaTransform<N>,
    *   by setting the spectral diffeomorphism's handle to the 2D Shearing diffeomorphism.
    *   Since no more templates are involved, it can be used in a precompiled form.
    *
    *   See the documentation of SigmaTransform<N> for more information.
    */
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

    /** Class NPShearlet2D is a two-dimensional Shearlet-instantation of SigmaTransform<N>,
    *   by setting the spectral diffeomorphism's handle to the 2D Shearing diffeomorphism
    *   and the "group action"-handle to a "parabolic shift".
    *   Since no more templates are involved, it can be used in a precompiled form.
    *
    *   See the documentation of SigmaTransform<N> for more information.
    */
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
