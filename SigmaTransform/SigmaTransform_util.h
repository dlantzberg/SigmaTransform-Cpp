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

#include <fftw3.h>

#define DEBUG

//using namespace std;

#ifdef DEBUG

#include <chrono>
#include <iomanip>

#endif

namespace SigmaTransform {

	const double PI     	= 3.141592654;
	const double pow2_1_4  	= 1.189207115;

    #ifdef DEBUG

    class Chronometer {
            std::chrono::steady_clock::time_point _tic;
        public:
            Chronometer();
            void tic();
            Chronometer& toc(std::string const& str );
    };

    #endif

    template<size_t N>
    class point {
            // data
            std::array<double,N>    _data;

        public:

            // getData
            double* get() { return _data; }

            // constructor
            point<N>() : point<N>({0}) {}
            point<N>(double const& x) { for( auto& _x : _data ) _x = x; }
            point<N>( std::array<double,N> dat ) : _data( dat ) {}
            //point<N>( std::array<int,N> dat ) { for(int k=0;k<N;++k)_data[k]=(double)dat[k];}

            // overloaded arithmetic operator with const-qualifier
            point<N> operator+( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for( auto& x : p ) x += *R++; return std::move( p ); }
            point<N> operator-( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for( auto& x : p ) x -= *R++; return std::move( p ); }
            point<N> operator*( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for( auto& x : p ) x *= *R++; return std::move( p ); }
            point<N> operator/( point<N> const& r ) const { point<N> p(_data); auto R=r.begin(); for( auto& x : p ) x /= *R++; return std::move( p ); }
            point<N> operator/( double const& r ) const { point<N> p(_data); for( auto& x : p ) x /= r; return std::move( p ); }
            point<N> operator*( double const& r ) const { point<N> p(_data); for( auto& x : p ) x *= r; return std::move( p ); }
            point<N> operator+( double const& r ) const { point<N> p(_data); for( auto& x : p ) x += r; return std::move( p ); }
            point<N> operator-( double const& r ) const { point<N> p(_data); for( auto& x : p ) x -= r; return std::move( p ); }

            void operator=( std::array<int,N> dat ) { for(int k=0;k<N;++k)_data[k]=(double)dat[k]; }

            // overloaded arithmetic operator
            void operator+=( point<N> const&  r ) { auto R=r.begin(); for( auto& x : _data ) x += *R++; }
            void operator-=( point<N> const&  r ) { auto R=r.begin(); for( auto& x : _data ) x -= *R++; }
            void operator*=( point<N> const&  r ) { auto R=r.begin(); for( auto& x : _data ) x *= *R++; }
            void operator/=( point<N> const&  r ) { auto R=r.begin(); for( auto& x : _data ) x /= *R++; }

            double&         operator[](size_t const& ind)       { return _data[ind]; }
            double const&   operator[](size_t const& ind) const { return _data[ind]; }

            // sum / prod / sq(uare)
            double   sum()  const { double s = 0; for( auto const& x : _data ) s+=x; return s; }
            double   prod() const { double p = 1; for( auto const& x : _data ) p*=x; return p; }
            point<N> sq() const { point<N> p( _data ); for( auto& x : p ) x *= x; return std::move( p ); }
            point<N> abs() const { point<N> p( _data ); for( auto& x : p ) x = (x>=0)?x:-x; return std::move( p ); }

            point<N> apply( double(*fun)(double const&) ) const { point<N> p; for(int k=0;k<N;++k) p[k] = fun(_data[k]); return std::move( p ); }
            point<N> apply( double(*fun)(double) ) const { point<N> p; for(int k=0;k<N;++k) p[k] = fun(_data[k]); return std::move( p ); }

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

            // iterators
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

    // shorthands
    using   cmpx        = std::complex<double>;
    using   cxVec       = std::vector<cmpx>;

    template<size_t N> using   diffFunc = std::function<point<N>(const point<N>&)>;
    template<size_t N> using   winFunc  = std::function<cmpx(const point<N>&)>;
    template<size_t N> using   actFunc  = std::function<point<N>(const point<N>&,const point<N>&)>;

    cxVec loadAscii1D( std::string const& filename );
	cxVec loadAscii2D( std::string const& filename, int &x , int &y );

	// make a linearly spaced vector between L and R in N steps
    std::vector<double> linspace( const double &L , const double &R , const int &N );
	// make a Fourier axis of length "len", with sampling frequency "fs"
    std::vector<double> FourierAxis( const double &fs , const unsigned &len );

    // a recursive loop replaces an arbitrary number of nested loops
    void recursiveLoop( std::vector<int> &index, const std::vector<int>& max , const int& curr_ind ,
                        std::function<void(const std::vector<int>&)> toDo );

    // start a recursive loopa and apply "toDo" for each point between (0,...,0) and "(max_1,...,max_n)"
    void StartRecursiveLoop( const std::vector<int>& max , std::function<void(const std::vector<int>&)> toDo );

    // make N-D-meshgrid from N 1D-vectors
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

    // make N-D-meshgrid from 1D-vector
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

    template<size_t N>
    void save2file_asc( std::string const& filename , const cxVec &out , point<N> const& dim = {0} ) {
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

    void save2file_bin( std::string const& filename , const cxVec &out );

    template<size_t N>
    point<N> plus( const point<N> &l , const point<N> &r )  { return l + r; }

    template<size_t N>
    point<N> minus( const point<N> &l , const point<N> &r ) { return l - r; }

    template<size_t N>
    //point<N> logabs( const point<N> &in ) { point<N> out(in); for(auto&x:out)x=log(abs(x)+1E-16); return std::move(out); }
    point<N> logabs( const point<N> &in ) { return (in.abs()+1E-16).apply( log2 ); }
    template<size_t N>
    point<N> logpos( const point<N> &in ) { return (in>0).apply( log2 ); }

    template<size_t N>
    point<N> id( const point<N> &x )  { return x; }

    template<size_t N>
    cmpx gauss( const point<N> &p  )  {
        return exp( -PI * p.sq().sum() ) * point<N>(pow2_1_4).prod();
    }

    template<size_t N>
    cmpx gauss_stddev( const point<N> &p , const point<N> &stddev ) {
        return exp( -PI * (p/stddev).sq().sum() ) * point<N>(pow2_1_4).prod();
    }


} // namespace SigmaTransform

#endif //SIGMATRANSFORM_H
