
#define DEBUG

#ifdef DEBUG

#include <chrono>
#include <iomanip>

#endif

#include "SigmaTransform_util.h"

namespace SigmaTransform {

    #ifdef DEBUG

    Chronometer::Chronometer() : _tic( std::chrono::steady_clock::now() ) { }
    void Chronometer::tic() { _tic = std::chrono::steady_clock::now(); }
    Chronometer& Chronometer::toc(std::string const& str ) {
        // get time
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        double diff = std::chrono::duration_cast<std::chrono::microseconds>(end-_tic).count();
        std::cout << std::fixed << std::showpoint << std::setw(60) << std::setfill('.') << std::setprecision(2)
                  << (diff?diff:1000.0)/1000.0 << " ms elapsed \r["<<str<<"]\n"<<std::flush;

        return *this;
    }

    #endif

    // loads an ascii file, containing 1d data
	cxVec loadAscii1D( std::string const& filename ) {
		// tmp file
		double tmp;
        cxVec out;
		// input filestream
		std::ifstream is(filename);
		if( !is )
            throw std::runtime_error("Error opening File.");
		// read
		while( is >> tmp )
			out.push_back( tmp );
        // return the read data
        return std::move( out );
	}

    // loads an ascii file, containing 2d data
	cxVec loadAscii2D( std::string const& filename , int &x , int &y ) {
		// tmp
		double tmp;
        cxVec out;
		// string for a line
		std::string line;
		// input filestream
		std::ifstream is(filename);
		if( !is )
            throw std::runtime_error("Error opening File.");
		// read file
		y=0;
		while( std::getline(is, line,'\n') ) {
			// read line
			x = 0;
			std::stringstream ss(line);
			while( ss >> tmp ) {
				out.push_back( tmp );
				++x;
			}
			++y;
		}
        // return the read data
        return std::move( out );
	}

    // saves "out" as a binary file "filename"
    void save2file_bin( std::string const& filename , const cxVec &out ) {
		std::ofstream( filename , std::ios::binary ).write((char*)out.data(),sizeof(double)*2*out.size());
	}

	// returns a linearly spaced vector between L and R in N steps
    std::vector<double> linspace( const double &L , const double &R , const int &N ) {
        std::vector<double> out(N);
        double delta = (R-L) / (N-1);
        double Lr = L;
        for( auto &val : out ) {
            val = Lr;
            Lr += delta;
        }
        return out;
    }

	// returns a linearly spaced vector between L and R in N steps
    std::vector<double> FourierAxis( const double &fs , const unsigned &len ) {
        std::vector<double> domain = linspace( 0 , fs , len );
        double offset = fs + domain[1];
        std::for_each( domain.begin() + ceil( len / 2 ) , domain.end() , [&](double& x){ x -= offset; } );
        return domain;
    }

    // this recursive loop replaces nested loops
    void recursiveLoop( std::vector<int> &index, const std::vector<int>& max , const int& curr_ind ,
                        std::function<void(const std::vector<int>&)> toDo ) {
        for( index[curr_ind] = 0 ; index[curr_ind] < max[curr_ind] ; ++index[curr_ind] ) {
            if( curr_ind )
                recursiveLoop( index , max , curr_ind-1 , toDo );
            else
                toDo( index );
        }
    }

    //void StartRecursiveLoop( const std::vector<int>& max , void (*toDo)(const std::vector<int>&) ) {
    void StartRecursiveLoop( const std::vector<int>& max , std::function<void(const std::vector<int>&)> toDo ) {
        std::vector<int> index; index.resize( max.size() , 0 );
        recursiveLoop( index , max , max.size()-1 , toDo );
    }

} // namespace SigmaTransform
