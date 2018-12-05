#ifndef SIGMATRANSFORM_UTIL_H
#define SIGMATRANSFORM_UTIL_H

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

#define DEBUG

#ifdef DEBUG

#include <chrono>
#include <iomanip>

#endif

#include <fftw3.h>

namespace SigmaTransform {

	const double PI     	= 3.141592654;
	const double pow2_1_4  	= 1.189207115;

    #ifdef DEBUG

    /** Class to conveniently measure the timing of executed code.
    */
    class Chronometer {
            std::chrono::steady_clock::time_point _tic;
        public:
            Chronometer();

            /** Starts/resets the timer.
             *
             *  @return             void
             */
            void tic();

            /** Stops the timer and prints the elapsed time, along with the message "str".
             *
             *  @param  str         a string, containing a message to print to stdout
             *
             *  @return             void
             */
            Chronometer& toc(std::string const& str );
    };

    #endif

    /** Class template for an N-dimensional point.
    *
    *   This class contains all kinds of overloaded operator for usage as what is
    *   usually perceived as a "point" in the linear algebraic sense, that is,
    *   points may be added, subtracted, multiplied, compared, etc. and implementations for
    *   the modulus, square, sum are available and the standard c++-output streams are overloaded.
    *
    *   Since the names of the methods are self-explanatory, this header shall suffice as a documentation.
    */
    template<size_t N>
    class point {
            // private, internally handled data
            std::array<double,N>    _data;

        public:

            // constructors
            point<N>() : point<N>({0}) {}
            point<N>(double const& x) { for( auto& _x : _data ) _x = x; }
            point<N>( std::array<double,N> dat ) : _data( dat ) {}

            // we may overwrite a point
            void operator=( std::array<int,N> dat ) { for(int k=0;k<N;++k)_data[k]=(double)dat[k]; }

            // overloaded arithmetic operator with const-qualifier
            point<N> operator+( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for( auto& x : p ) x += *R++; return std::move( p ); }
            point<N> operator-( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for( auto& x : p ) x -= *R++; return std::move( p ); }
            point<N> operator*( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for( auto& x : p ) x *= *R++; return std::move( p ); }
            point<N> operator/( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for( auto& x : p ) x /= *R++; return std::move( p ); }
            point<N> operator/( double const& r ) const { point<N> p(_data); for( auto& x : p ) x /= r; return std::move( p ); }
            point<N> operator*( double const& r ) const { point<N> p(_data); for( auto& x : p ) x *= r; return std::move( p ); }
            point<N> operator+( double const& r ) const { point<N> p(_data); for( auto& x : p ) x += r; return std::move( p ); }
            point<N> operator-( double const& r ) const { point<N> p(_data); for( auto& x : p ) x -= r; return std::move( p ); }

            // overloaded arithmetic operator
            void operator+=( point<N> const&  r ) { auto R=r.begin(); for( auto& x : _data ) x += *R++; }
            void operator-=( point<N> const&  r ) { auto R=r.begin(); for( auto& x : _data ) x -= *R++; }
            void operator*=( point<N> const&  r ) { auto R=r.begin(); for( auto& x : _data ) x *= *R++; }
            void operator/=( point<N> const&  r ) { auto R=r.begin(); for( auto& x : _data ) x /= *R++; }

            // sum / prod / abs/ sq(uare)
            double   sum()  const { double s = 0; for( auto const& x : _data ) s+=x; return s; }
            double   prod() const { double p = 1; for( auto const& x : _data ) p*=x; return p; }
            point<N> sq() const { point<N> p( _data ); for( auto& x : p ) x *= x; return std::move( p ); }
            point<N> abs() const { point<N> p( _data ); for( auto& x : p ) x = (x>=0)?x:-x; return std::move( p ); }

            // apply specific functions to each component in the "point-vector"
            point<N> apply( double(*fun)(double const&) ) const { point<N> p; for(int k=0;k<N;++k) p[k] = fun(_data[k]); return std::move( p ); }
            point<N> apply( double(*fun)(double) ) const { point<N> p; for(int k=0;k<N;++k) p[k] = fun(_data[k]); return std::move( p ); }

            // comparison-related operator
            point<N> operator>( double const& d ) const { point<N> p(_data); for( auto& x : p ) x = x>d; return std::move( p ); }
            point<N> operator<( double const& d ) const { point<N> p(_data); for( auto& x : p ) x = x<d; return std::move( p ); }
            point<N> operator<=( double const& d ) const { point<N> p(_data); for( auto& x : p ) x = x<=d; return std::move( p ); }
            point<N> operator>=( double const& d ) const { point<N> p(_data); for( auto& x : p ) x = x>=d; return std::move( p ); }
            point<N> operator==( double const& d ) const { point<N> p(_data); for( auto& x : p ) x = x==d; return std::move( p ); }
            point<N> operator>( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for(auto&x : p) x = x>*R++; return std::move( p ); }
            point<N> operator<( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for(auto&x : p) x = x<*R++; return std::move( p ); }
            point<N> operator>=( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for(auto&x : p) x = x>=*R++; return std::move( p ); }
            point<N> operator<=( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for(auto&x : p) x = x<=*R++; return std::move( p ); }
            point<N> operator==( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for(auto&x : p) x = x==*R++; return std::move( p ); }
            operator bool() { bool b=true; for(auto&x:_data) if(x!=1) b=false; return b; }

            // acquire internal data
            double&         operator[](size_t const& ind)       { return _data[ind]; }
            double const&   operator[](size_t const& ind) const { return _data[ind]; }
            double* get() { return _data; }

            // iterators, for compatibility with the C++ - Std Template Library
            typedef typename std::array<double,N>::iterator                 iterator;
            typedef typename std::array<double,N>::const_iterator           const_iterator;
            typedef typename std::array<double,N>::reverse_iterator         reverse_iterator;
            typedef typename std::array<double,N>::const_reverse_iterator   const_reverse_iterator;

            iterator       begin()                 { return _data.begin(); }
            iterator       end()                   { return _data.end(); }

            const_iterator begin()  const          { return _data.begin(); }
            const_iterator end()    const          { return _data.end(); }

            const_iterator cbegin() const          { return _data.cbegin(); }
            const_iterator cend()   const          { return _data.cend(); }

            reverse_iterator rbegin()              { return _data.rbegin(); }
            reverse_iterator rend()                { return _data.rend(); }

            const_reverse_iterator crbegin() const { return _data.crbegin(); }
            const_reverse_iterator crend()   const { return _data.crend(); }

            // friends for streams
     friend std::ostream& operator<<( std::ostream& os , point<N> const& p ) { os<<"("; for( auto& x : p ) os << x << " "; return os << "\b)"; }
     friend std::istream& operator>>( std::istream& is , point<N>& p )       { for( auto& x : p ) is >> x; return is; }

    };

    // shorthands for cleaner code
    using   cmpx        = std::complex<double>;
    using   cxVec       = std::vector<cmpx>;

    template<size_t N> using   diffFunc = std::function<point<N>(const point<N>&)>;
    template<size_t N> using   winFunc  = std::function<cmpx(const point<N>&)>;
    template<size_t N> using   actFunc  = std::function<point<N>(const point<N>&,const point<N>&)>;

    /** Loads 1D-ASCII file.
     *
     *  @param  filename    the filename
     *
     *  @return             the read file as complex vector
     *
     *  @throws             std::runtime_error
     */
    cxVec loadAscii1D( std::string const& filename );

    /** Loads 2D-ASCII file.
     *
     *  @param  filename    the filename
     *  @param  x           reference to an integer in which to put the x dimension
     *  @param  y           reference to an integer in which to put the y dimension
     *
     *  @return             the read file as complex vector
     *
     *  @throws             std::runtime_error
     */
	cxVec loadAscii2D( std::string const& filename, int &x , int &y );

    /** Makes a linearly spaced vector between L and R in N steps.
     *
     *  @param  L           left boundary
     *  @param  R           right boundary
     *  @param  N           number of points
     *
     *  @return             linearly spaced vector
     */
    std::vector<double> linspace( const double &L , const double &R , const int &N );

    /** Makes a Fourier axis of length "len", with sampling frequency "fs".
     *
     *  @param  fs          sampling frequency
     *  @param  len         number of points
     *
     *  @return             linearly spaced Fourier-axis vector
     */
    std::vector<double> FourierAxis( const double &fs , const unsigned &len );

    /** Starts a recursive loop and applies "toDo" for each point between (0,...,0) and "(max_1,...,max_n)"
    *
    *  @param  max         vector, containing the maximal iteration-indices
    *  @param  toDo        function handle, performing the nested operations, depending on the handed index
    *
    *  @return             void
    */
    void StartRecursiveLoop( const std::vector<int>& max , std::function<void(const std::vector<int>&)> toDo );

    /** Applies a recursive loop, which replaces an arbitrary number of nested loops; internally called by
     *  "StartRecursiveLoop".
     *
     *  @param  index       vector, containing the current iteration-indices
     *  @param  max         vector, containing the maximal iteration-indices
     *  @param  curr_ind    the current nested iteration
     *  @param  toDo        function handle, performing the nested operations, depending on the handed index
     *
     *  @return             void
     */
    void recursiveLoop( std::vector<int> &index, const std::vector<int>& max , const int& curr_ind ,
                        std::function<void(const std::vector<int>&)> toDo );

    /** Makes N-D-meshgrid from N 1D-vectors.
    *
    *  @param  dom         array containg N double-vectors
    *
    *  @return             N-D meshgrid
    */
    template<size_t N>
    std::vector<point<N>> meshgridN( std::array<std::vector<double>,N> const& dom ) {
        // make sizeVector for loop
        std::vector<int> sizeVec(0);
        int p = 1;
        for( auto& d : dom ) {
            p *= d.size();
            sizeVec.push_back( d.size() );
        }
        // reserve space
        std::vector<point<N>> out; out.reserve( p );
        // use recursive loop - very slow, but replaces up to n nested loops
        StartRecursiveLoop( sizeVec , [&](std::vector<int> const& index ) {
            // get current point
            point<N> p;
            for( int k = 0 ; k < index.size() ; ++k )
                p[k] = dom[k][index[k]];
            // insert into vectors
            out.emplace_back( p );
        } );
        // return output vector
        return std::move( out );
    }

    /** Makes symmetric N-D-meshgrid from single 1D-vector.
    *
    *  @param  dom         a double-vectors
    *
    *  @return             symmetric N-D meshgrid
    */
    template<size_t N>
    std::vector<point<N>> meshgridN( std::vector<double> const& dom ) {
        // make sizeVector for loop
        std::vector<int> sizeVec(0);
        int p = 1;
        for(int k = 0 ; k < N;++k) {
            p *= dom.size();
            sizeVec.push_back( dom.size() );
        }
        // reserve space
        std::vector<point<N>> out; out.reserve( p );
        // use recursive loop - very slow, but replaces up to n nested loops
        StartRecursiveLoop( sizeVec , [&](std::vector<int> const& index ) {
            // get current point
            point<N> p;
            for( int k = 0 ; k < index.size() ; ++k )
                p[k] = dom[index[k]];
            // insert into vectors
            out.emplace_back( p );
        } );
        // return output vector
        return std::move( out );
    }

    /** Saves N-dim complex vector to ASCII-file.
    *
    *  @param  filename    name of file to save data to
    *  @param  out         complex vector, containing the data
    *  @param  dim         N-dim point, containing the dimensions of the data
    *
    *  @return             void
    */
    template<size_t N>
    void save2file_asc( std::string const& filename , const cxVec &out , const point<N> & dim = point<N>(0) ) {
        // use a stringstream buffer
        std::stringstream ss;
        // set precision
        ss.precision(6);
		//write to buffer
		for( const auto &c : out )
			ss << c.real() << "," << c.imag() << std::endl;
		// output filestream
		std::ofstream os(filename);
        // save dimensions of vector
        if( dim.prod() == out.size() )
            os 	<< "dim=" << dim << "\n";
        else
            os << "dim=(" <<out.size()<<")\n";
        // write buffer to file
        os<<ss.str();
	}

    /** Saves N-dim complex vector to binary-file.
    *
    *  @param  filename    name of file to save data to
    *  @param  out         complex vector, containing the data
    *
    *  @return             void
    */
    void save2file_bin( std::string const& filename , const cxVec &out );

    /** Simple N-dim addition.
    *
    *  @param  l            left operand
    *  @param  r            right operand
    *
    *  @return              sum
    */
    template<size_t N>
    point<N> plus( const point<N> &l , const point<N> &r )  { return l + r; }

    /** Simple N-dim subtraction.
    *
    *  @param  l            left operand
    *  @param  r            right operand
    *
    *  @return              difference
    */
    template<size_t N>
    point<N> minus( const point<N> &l , const point<N> &r ) { return l - r; }

    /** Logarithm of the modulus.
    *
    *  @param  x            N-dim point
    *
    *  @return              log( abs( x ) )
    */
    template<size_t N>
    point<N> logabs( const point<N> &x ) { return (x.abs()+1E-16).apply( log2 ); }

    /** "Positive" logarithm.
    *
    *  @param  x           N-dim point
    *
    *  @return             logarithm of x if positive, -inf else
    */
    template<size_t N>
    point<N> logpos( const point<N> &x ) { return (x>0).apply( log2 ); }

    /** Identical diffeomorphism.
    *
    *  @param  x            N-dim point
    *
    *  @return              x
    */
    template<size_t N>
    point<N> id( const point<N> &x )  { return x; }

    /** Normalized Gaussian window.
    *
    *  @param  x            N-dim point
    *
    *  @return              normalized N-dim gaussian at x
    */
    template<size_t N>
    cmpx gauss( const point<N> &x  )  {
        return exp( -PI * x.sq().sum() ) * point<N>(pow2_1_4).prod();
    }

    /** Gaussian window of specific std-dev.
    *
    *  @param  x            N-dim point
    *  @param  stddev       N-dim standard deviation
    *
    *  @return              N-dim gaussian with N-dim standarddev "stddev" at x
    */
    template<size_t N>
    cmpx gauss_stddev( const point<N> &x , const point<N> &stddev ) {
        return exp( -PI * (x/stddev).sq().sum() ) * point<N>(pow2_1_4).prod();
    }

} // namespace SigmaTransform

#endif //SIGMATRANSFORM_H
