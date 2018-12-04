#ifndef SIGMATRANSFORM_H
#define SIGMATRANSFORM_H

#include <vector>
#include <map>
#include <complex>
#include <thread>
#include <array>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <mutex>
#include <functional>
#include <fftw3.h>

#define DEBUG

#include "SigmaTransform_util.h"

namespace SigmaTransform {

    // shorthands
    using   cmpx        = std::complex<double>;
    using   cxVec       = std::vector<cmpx>;

    template<size_t N> using   diffFunc = std::function<point<N>(const point<N>&)>;
    template<size_t N> using   winFunc  = std::function<cmpx(const point<N>&)>;
    template<size_t N> using   actFunc  = std::function<point<N>(const point<N>&,const point<N>&)>;


    template<size_t N>
    class SigmaTransform {

        public:
            SigmaTransform( diffFunc<N> sigma=NULL, const point<N>& winWidth={0}, const point<N> &Fs=point<N>(0), const point<N> &size=point<N>(0),
                            const std::vector<point<N>> &steps=std::vector<point<N>>(0), actFunc<N> action=minus<N> , int const& numThreads = 4 )
            : SigmaTransform( sigma, (winFunc<N>)(NULL), Fs,size,steps,action,numThreads) { m_winWidth = winWidth; }
            SigmaTransform( diffFunc<N> sigma=NULL, winFunc<N> window=NULL, const point<N> &Fs=point<N>(0), const point<N> &size=point<N>(0),
                            const std::vector<point<N>> &steps=std::vector<point<N>>(0), actFunc<N> action=minus<N> , int const& numThreads = 4 )
            : m_window(window),m_sigma(sigma?sigma:id<N>),m_action(action?action:minus<N>),m_windows(0),m_coeff(0),m_reconstructed(0),
              m_size(size),m_fs(Fs) , m_winWidth(0.0) {
                setSteps( steps );
                if( !fftw_init_threads() )
                    std::cerr << "thread error\n";
                setNumThreads( numThreads );
            }

            SigmaTransform& setWindow( winFunc<N> window ) { m_window = window; return *this; }
            SigmaTransform& setSigma( diffFunc<N> sigma ) { m_sigma = sigma?sigma:id<N>; return *this;  }
            SigmaTransform& setAction( actFunc<N> action ) { m_action = action?action:minus<N>; return *this; }
            SigmaTransform& setFs( const point<N> &Fs ) { m_fs = Fs; return *this; }
            SigmaTransform& setSize( const point<N> &size ) { m_size = size; return *this; }
            SigmaTransform& setWinWidth( const double& winWidth ) { m_window = NULL; m_winWidth = winWidth; return *this; }
            SigmaTransform& setNumThreads( const int &numThreads )
                { m_numThreads = numThreads ; fftw_plan_with_nthreads(numThreads); return *this; }
            SigmaTransform& setSteps( const std::vector<point<N>> &steps ) {
                if( steps.empty() ) {
                    throw std::runtime_error("steps must not be empty.");
                }
                // only number of stepps given? Create steps-array
                if( steps.size() == 1 ) {
                    if( m_fs==0.0 || !m_sigma) {
                        throw std::runtime_error("Sampling rate and diffeomorphism must be set for automatic step-array determination");
                    }
                    // make Fourier-domain
                    std::array<std::vector<double>,N> stps;
                    for( int k = 0 ; k < N ; ++k ) {
                        stps[k] = linspace( -m_fs[k] / 2.0 , m_fs[k] / 2.0 , steps[0][k] );
                    }
                    m_steps = meshgridN( stps );
                    // get boundary  points of warped domain
                    point<N> maxi,mini; maxi=mini=m_steps[0];
                    for( auto& step : m_steps ) {
                        step = m_sigma( step );
                        for( int k = 0 ; k < N ; ++k ) {
                            maxi[k]     = (step[k]>maxi[k])?step[k]:maxi[k];
                            mini[k]     = (step[k]<mini[k])?step[k]:mini[k];
                        }
                    }

                    // make linear spacing in warped domain
                    for( int k = 0 ; k < N ; ++k ) {
                        stps[k] = linspace( mini[k] , maxi[k] , steps[0][k] );
                    }
                    m_steps = meshgridN( stps );
                } else {
                    m_steps = steps;
                }
                //for( auto& step : m_steps )
                //    std::cout << step << " " << std::flush;
                return *this;
            }

            cxVec& getCoeffs(){ return m_coeff; }
            cxVec& getWindows(){ return m_windows; }
            cxVec& getReconstruction(){ return m_reconstructed; }
            SigmaTransform& operator()( cxVec const& sig ){ return analyze( sig ); }

            SigmaTransform& analyze( cxVec const& sig , std::function<void(SigmaTransform*)> onFinish = NULL ) {
                if( !m_steps.size() ) {
                    throw std::runtime_error("steps not set");
                }
                if( !m_fs.prod() ) {
                    throw std::runtime_error("Fs not set");
                }
                if( !m_size.prod() ) {
                    throw std::runtime_error("size not set");
                }
                return (onFinish) ? asyncTransform( sig , onFinish ) : applyTransform( sig );
            }

            SigmaTransform& synthesize( std::function<void(SigmaTransform* obj)> onFinish = NULL ){
                return (onFinish) ? asyncInverseTransform( onFinish ) : applyInverseTransform( );
            }

            cxVec& multiplier( cxVec const& sig, cxVec const& mask, std::function<void(SigmaTransform*)> onFinish = NULL  ){
                return (onFinish)?asyncMultiplier( sig, mask, onFinish ) : analyze(sig).applyMask(mask).synthesize().getReconstruction();
            }

            SigmaTransform& applyMask( const cxVec &mask ) {
                // error?
                if( mask.size() != m_coeff.size() ) {
                    throw std::runtime_error("Size of mask does not match size of coefficients.");
                }
                // deallocate threads after leaving scope
                std::unique_ptr<std::thread[]> _threads( new std::thread[m_numThreads] );
                // start threads
                int stepsPerThread = ceil( (double) (m_steps.size()) / m_numThreads );
                for( int k = 0, stepsLeft = m_steps.size() ; k < m_numThreads ; ++k, stepsLeft -= stepsPerThread ) {
                    _threads[k] = std::move( std::thread( [this,&stepsPerThread,&mask,stepsLeft,k]() {
                        // get iterator from correct offset
                        int offset = k*stepsPerThread*m_size.prod(),
                            numval = ((stepsLeft>stepsPerThread)?stepsPerThread:stepsLeft)*m_size.prod();
                        // multiply nth part of coeffs
                        for( int i = 0 ; i < numval ; ++i ) {
                            m_coeff[offset+i] *= mask[offset+i];
                        }
                    } ) );
                }
                // wait for all threads to finish
                for( int k = 0 ; k < m_numThreads ; ++k ) {
                    _threads[k].join();
                }
                // return
                return *this;
            }

            SigmaTransform& applyMask( std::function<cmpx(point<N>const&,point<N>const&)> maskFunc ) {
                // spatial domain
                auto spatialDom = makeSpatialDomain();
                // deallocates itself after leaving scope
                std::unique_ptr<std::thread[]> _threads( new std::thread[m_numThreads] );
                    // start threads
                int stepsPerThread = ceil( (double) (m_steps.size()) / m_numThreads );
                for( int k = 0, stepsLeft = m_steps.size() ; k < m_numThreads ; ++k, stepsLeft -= stepsPerThread ) {
                    _threads[k] = std::move( std::thread( [this,&maskFunc,&spatialDom,&stepsPerThread,stepsLeft,k]() {
                        // get offsets and iterator
                        int  stepsOffset = stepsPerThread*k,
                             numval      = (stepsLeft>stepsPerThread)?stepsPerThread:stepsLeft;
                        auto coeff       = m_coeff.begin() + stepsOffset*m_size.prod();
                        // run thru all steps
                        for( int i=0;i<numval;++i) {
                            // run thru all spatial values
                            for( auto const& x : spatialDom ) {
                                *coeff++ *= maskFunc( x , m_steps[stepsOffset+i] );
                            }
                        }
                    } ) );
                }
                // wait for all threads to finish
                for( int k = 0 ; k < m_numThreads ; ++k ) {
                    _threads[k].join();
                }
                // return
                return *this;
            }

            // make warped domain
            void makeWarpedDomain() {
                // make (fft-shifted) domain...
                std::array<std::vector<double>,N> doms;
                auto itFs = m_fs.begin(),itSz = m_size.begin();
                for(auto& d : doms) {
                    d = FourierAxis( *itFs++ , *itSz++ );
                }
                m_domain = meshgridN( doms );
                // ...and warp the domain
                std::for_each(m_domain.begin(),m_domain.end(),[&](point<N>&x){x=m_sigma(x);});
            }

            // make spatial-domain
            std::vector<point<N>> makeSpatialDomain() {
                std::array<std::vector<double>,N> doms;
                auto itFs = m_fs.begin(), itSz = m_size.begin();
                for(auto& d : doms) {
                    //std::cout << "0 - " << (*itSz-1) << "/"  << *itFs << " - " << *itSz << "\n";
                    d = linspace( 0 , (*itSz-1) / *itFs , *itSz );
                    itFs++; itSz++;
                }
                return std::move( meshgridN( doms ) );
            }

            // makes a gaussian window
            void makeWarpedGaussian() {
                // make adequate standard deviation
                point<N> maxi,mini,num_steps{1};
                maxi=mini=m_steps[0];
                for( auto& step : m_steps ) {
                    for( int k = 0 ; k < N ; ++k ) {
                        if( step[k]>maxi[k] ) {
                            maxi[k] = step[k];
                            num_steps[k]++;
                        }
                        if( step[k]<mini[k] ) {
                            mini[k] = step[k];
                        }
                    }
                }
                // winWidth set?
                if( m_winWidth==0.0 )
                    // set window width of 8
                    m_winWidth = point<N>(8.0);

                auto width = (maxi-mini) / num_steps * m_winWidth;
                #ifdef DEBUG
//                std::cout << "Width = " << width << std::endl;
                #endif
                m_window = [&](const point<N>&x)->cmpx{ return gauss_stddev( x , width ); };
            }

            // create a set of "m_steps.size()" windows in the Fourier domain
            void makeWindows( ) {
                // reserve space for windows..
                m_windows.resize( m_size.prod() * m_steps.size() );
                // make Domain
                makeWarpedDomain();
                // check if window was given, else calculate good width for a warped gaussian window
                if( !m_window )
                    makeWarpedGaussian();
                // deallocates itself after leaving scope
                std::unique_ptr<std::thread[]> _threads( new std::thread[m_numThreads] );
                // start threads
                int stepsPerThread = ceil( (double) (m_steps.size()) / m_numThreads );
                for( int k = 0, stepsLeft = m_steps.size() ; k < m_numThreads ; ++k, stepsLeft -= stepsPerThread ) {
                    _threads[k] = std::move( std::thread( [this,&stepsPerThread,k,stepsLeft]() {
                        // get iterator from correct offset
                        auto win = m_windows.begin() + k*stepsPerThread * m_size.prod();
                        // create n-th part of the windows
                        for( int i = 0 ; i < ((stepsLeft>stepsPerThread)?stepsPerThread:stepsLeft) ; ++i ) {
                            for( auto const& x : m_domain ) {
                                *win++ = m_window( m_action( x , m_steps[k*stepsPerThread + i] ) );
                            }
                        }
                    } ) );
                }
                // wait for all threads to finish
                for( int k = 0 ; k < m_numThreads ; ++k ) {
                    _threads[k].join();
                }
            }

            /* thread-based asynchronous methods */
            SigmaTransform& asyncTransform( cxVec const& sig, std::function<void(SigmaTransform*)> onFinish ) {
                 m_threads.insert( std::make_pair<std::string,std::thread>( "Transform" , std::move( std::thread( [this,sig,onFinish]() {
                    //std::cout << "Starting 'Transform' asynchronously."<<std::endl;
                    // try to acquire mutex
                    std::unique_lock<std::mutex> lk( m_mtx );
                    //std::cout << "Acquired mutex in 'Transform'"<<std::endl;
                    // transform asynchronously
                    this->analyze( sig );
                    // call callback
                    onFinish( this );
                } ) ) ) );
                return *this;
            }

            SigmaTransform& asyncInverseTransform( std::function<void(SigmaTransform*)> onFinish ) {
                 m_threads.insert( std::make_pair<std::string,std::thread>( "Inverse" , std::move( std::thread( [this,onFinish]() {
                    //std::cout << "Starting 'Inverse' asynchronously."<<std::endl;
                    // try to acquire mutex
                    std::unique_lock<std::mutex> lk( m_mtx );
                    //std::cout << "Acquired mutex in 'Inverse'"<<std::endl;
                    // transform asynchronously
                    this->synthesize( );
                    // call callback
                    onFinish( this );
                } ) ) ) );
                return *this;
            }

            SigmaTransform& asyncMultiplier( cxVec const& sig, cxVec const& mask, std::function<void(SigmaTransform*)> onFinish ) {
                 m_threads.insert( std::make_pair<std::string,std::thread>( "Multiplier" , std::move( std::thread( [this,sig,mask,onFinish]() {
                    //std::cout << "Starting 'Multiplier' asynchronously."<<std::endl;
                    // try to acquire mutex
                    std::unique_lock<std::mutex> lk( m_mtx );
                    //std::cout << "Acquired mutex in 'Multiplier'"<<std::endl;
                    // transform asynchronously
                    this->multiplier( sig,mask );
                    // call callback
                    onFinish( this );
                } ) ) ) );
                return *this;
            }

            void join() {
                // get threads to join
                std::vector<std::string> toErase;
                for( auto it = m_threads.begin(); it != m_threads.end() ; ++it ) {
                    if( it->second.joinable() ) {
                        //std::cout << "joining " << it->first << ".\n";
                        it->second.join();
                        toErase.push_back( it->first );
                    }
                }
                // erase joined threads
                for( auto const& name : toErase ) {
                    for( auto it = m_threads.begin(); it != m_threads.end() ; ++it ) {
                        if( name.compare( it->first) == 0 ) {
                            //std::cout << "erasing " << it->first << ".\n";
                            m_threads.erase( it );
                            break;
                        }
                    }
                }
            }

            cxVec fft( cxVec const& in , int const& howmany = 1 ) {
                // get space
                cxVec out( in.size() );
                // fft transform the signal
                fftN( reinterpret_cast<fftw_complex*> (out.data()) ,
                      reinterpret_cast<fftw_complex*>(const_cast<cmpx*> (in.data())) ,
                      m_size , howmany , FFTW_FORWARD );
                // return memory
                return std::move( out );
            }
            void fft_inplace( cxVec& inout , int const& howmany = 1 ) {
                // ifft transform the signal
                fftN( reinterpret_cast<fftw_complex*> (inout.data()) ,
                      reinterpret_cast<fftw_complex*>(const_cast<cmpx*> (inout.data())) ,
                      m_size , howmany , FFTW_FORWARD );
            }

            cxVec ifft( cxVec const& in , int const& howmany = 1 ) {
                // get space
                cxVec out( in.size() );
                // ifft transform the signal
                fftN( reinterpret_cast<fftw_complex*> (out.data()) ,
                      reinterpret_cast<fftw_complex*>(const_cast<cmpx*> (in.data())) ,
                      m_size , howmany , FFTW_BACKWARD );
                // return memory
                return std::move( out );
            }
            void ifft_inplace( cxVec& inout , int const& howmany = 1 ) {
                // ifft transform the signal
                fftN( reinterpret_cast<fftw_complex*> (inout.data()) ,
                      reinterpret_cast<fftw_complex*>(const_cast<cmpx*> (inout.data())) ,
                      m_size , howmany , FFTW_BACKWARD );
            }

            void ifft_inplace( cmpx* inout , int const& howmany = 1 ) {
                // ifft transform the signal
                fftN( reinterpret_cast<fftw_complex*> (inout) ,
                      reinterpret_cast<fftw_complex*> (inout) ,
                      m_size , howmany , FFTW_BACKWARD );
            }

        protected:

            SigmaTransform& applyTransform( const cxVec &in )  {
                // convenient var
                double sigsize  = m_size.prod();
                // make windows
                makeWindows( );
                // fft transform the signal
                cxVec Fsig = fft( in );
                // copy the windows
                m_coeff = m_windows;
                // deallocates itself after leaving scope
                std::unique_ptr<std::thread[]> _threads( new std::thread[m_numThreads] );
                // start threads
                int stepsPerThread = ceil( (double) (m_steps.size()) / m_numThreads );
                for( int k = 0, stepsLeft = m_steps.size() ; k < m_numThreads ; ++k, stepsLeft -= stepsPerThread ) {
                    _threads[k] = std::move( std::thread( [this,&Fsig,&stepsPerThread,stepsLeft,k]() {
                        cxVec::iterator coeff = m_coeff.begin() + (int)( k*stepsPerThread*m_size.prod() );
                        for( int i = 0 ; i < ((stepsLeft>stepsPerThread)?stepsPerThread:stepsLeft) ; ++i ) {
                            for( auto const& val : Fsig ) {
                                *coeff = conj(*coeff++) * val / ((double)Fsig.size());
                            }
                        }
                    } ) );
                }
                // wait for all threads to finish
                for( int k = 0 ; k < m_numThreads ; ++k ) {
                    _threads[k].join();
                }
                // transform back
                ifft_inplace( m_coeff , m_steps.size() );
                // return
                return *this;
            }

            SigmaTransform& applyInverseTransform()  {
                // reserve vectorspace with zeros
                m_reconstructed = cxVec( m_coeff.size() / m_steps.size() , 0 );

                // fft transform the signal
                cxVec temp = fft( m_coeff , m_steps.size() );

                // act on signal
                auto coeff = temp.begin(),
                    window = m_windows.begin();
                for(int k=0;k<m_steps.size();++k) {
                    for( auto &val : m_reconstructed ) {
                        val += (*coeff++) * (*window++);
                    }
                }

                // transform back
                ifft_inplace( m_reconstructed );

                // return
                return *this;
            }

            void fftN( fftw_complex *out, fftw_complex *in, const point<N> &size, const int &howmany = 1, const int& DIR = FFTW_FORWARD ) {
                int sz[N];
                for( int k = 0 ; k < N ; ++k )  sz[k] = (int) size[k];
                // make plan
                fftw_plan p = fftw_plan_many_dft( N , sz , howmany ,  in  , NULL , 1 , (int) size.prod() ,
                                                                      out , NULL , 1 , (int) size.prod() ,
                                                                      DIR , FFTW_ESTIMATE );
                // perform
                fftw_execute( p );
                // destroy plan
                fftw_destroy_plan( p );
            }

            // function handles for the transform
            std::function<cmpx(point<N>const&)>                     m_window;
            std::function<point<N>(point<N>const&)>                 m_sigma;
            std::function<point<N>(point<N>const&,point<N>const&)>  m_action;

            // holds data
            cxVec                                   m_windows;
            cxVec                                   m_coeff;
            cxVec                                   m_reconstructed;

            // holds information about data
            point<N>                                m_size;
            point<N>                                m_fs;
            std::vector<point<N>>                   m_steps;
            std::vector<point<N>>                   m_domain;
            int                                     m_numThreads;
            point<N>                                m_winWidth;

            // for asynchronous computations
            std::map<std::string,std::thread>       m_threads;
            std::mutex                              m_mtx;

    }; // class SigmaTransform

} // namespace SigmaTransform

#endif //SIGMATRANSFORM_H
