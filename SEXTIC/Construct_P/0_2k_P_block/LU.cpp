/*This code computes the full matrix P(n,m)
 n is range(0,d+1) and m is in range(0,d+1) */

#include <string>
#include <boost/lexical_cast.hpp>

//Boost.Multiprecision mpfr_float
#include <boost/multiprecision/mpfr.hpp>

//for parallel implementation
#include <boost/mpi.hpp>
#include <boost/thread.hpp>

//Boost.Math headers
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>


//Boost.Serialization headers
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/version.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/array.hpp>
#include <stdio.h> // used in creating File object

//Boost.Ublas headers
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

//for serialization and miscellaneous functionalities
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include<stdlib.h>
#include<time.h>

#define digits 3000

using namespace boost::archive;
using namespace boost::serialization;
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;
using namespace std::chrono;

//serialization code for mpfr_float

using realtype = number<mpfr_float_backend<digits, allocate_stack>>;
#define MPFR_BUFFER_SIZE 3010

    namespace boost {
        namespace serialization {
            template<class Archive>
            void save(Archive& ar, const realtype& x, const boost::serialization::version_type&) {
                static char buffer[MPFR_BUFFER_SIZE];
                FILE* fid = fmemopen(buffer, MPFR_BUFFER_SIZE, "wb+");
                mpfr_fpif_export(fid, const_cast<mpfr_ptr>(x.backend().data()));
                fseek(fid, 0L, SEEK_END);
                long length = ftell(fid);
                ar& length;
                ar& boost::serialization::make_array(buffer, length);
                fclose(fid);
            }

            template<class Archive>
            void load(Archive& ar, realtype& x, const boost::serialization::version_type&) {
                static char buffer[MPFR_BUFFER_SIZE];
                long length = 0;

                ar& length;
                ar& boost::serialization::make_array(buffer, length);

                FILE* fid = fmemopen(buffer, length, "r");
                mpfr_fpif_import(x.backend().data(), fid);
                fclose(fid);
            }

            template<class Archive>
            inline void serialize(Archive& ar, realtype& t, const unsigned int file_version) {
                split_free(ar, t, file_version);
            }
        }
    }

namespace mpi = boost::mpi;

int i,j,k,m,n;
const int l = 200, b = 9; // l*(b+1)  = d+1
const int d = l*(b+1) - 1;
const realtype Pi = boost::math::constants::pi<realtype>();
const realtype nu = -realtype(1)/realtype(2);

vector<realtype> term(d+1);
vector<realtype> factor( d+1 );
vector<realtype> factor1( 3*(d+1) + 1 );
vector<realtype> factor2( d+1 );
vector<realtype> factor3( d+1 );

int main()
{

    mpi::environment env;
    mpi::communicator world;
    int num_process = 3;

    int share = (int) ( floor((d + 1)/ num_process) ); /* The rows are divided amongst the processes equally. The last process
                                                            gets the left over rows. */

    // ---------------------------factors--------------------------------

    for (i=0; i < factor.size (); ++i)
    {
        factor(i) = boost::math::factorial <realtype> ( i );// factorial
    }

    for (i=0; i < factor1.size (); ++i)
    {
        factor1(i) = pow( realtype(2), realtype(i+1)-nu) * boost::math::tgamma(realtype(i+1)-nu);
    }

    for (i=0; i < factor2.size (); ++i)
    {
        factor2(i) = pow( realtype(-1), realtype( i )) / pow( factor(i), realtype( 2 ));
        factor3(i) = realtype(1) / factor(i);
    }

    //-------------------------------------------------------------------

    if (world.rank() == 0)
    {
        auto start = high_resolution_clock::now();

        std::vector < matrix<realtype> > ms;
        for ( i = 0; i < num_process - 1; ++i )
        {
            matrix<realtype> P1;
            ms.push_back(P1);
        }

        // create this process's share of the work
        matrix<realtype> P(share, d+1);
        for (n = 0; n < share; ++n)
        {
            for (k=0; k < d+1; ++k)
            {
                term(k) = factor1( 2*n + k )*factor2( k );
            }

            for (m = 0; m < d+1; ++m)
            {
                P(n,m) = factor(m) * inner_prod( project( term, range(0, m+1) ), project( factor3, slice( m, -1, m+1 ) ) );
            }
        }
        // receive the computations of other processes and store them in elements of ms[i]
        for (i = 0; i < num_process -1 ; ++i)
        {
            world.recv( i+1, 17, ms[i] );
        }

        // all parts of the matrix have been received, we can print them

        std::ofstream outfile1;
        outfile1.open("matrix_p.txt", std::ios_base::out);
        outfile1.precision(digits);
        for ( n = 0;n < share; n++ )
        {
            for( m = 0; m < P.size2(); m++)
            {
                outfile1 << std::scientific << P(n,m) << std::endl;
            }
        }

        for ( i = 0; i < num_process-1; ++i )
        {
            for ( n=0; n < ms[i].size1(); ++n )
            {
                for ( m = 0; m < ms[i].size2(); ++m )
                {
                    outfile1 << std::scientific << ms[i](n,m) << std::endl;
                }
            }
        }

        outfile1.close();

        std::cout << "share " << share << std::endl;
        std::cout << "num_process " << num_process << std::endl;
        std::cout << "d+1 is " << d+1 << std::endl;

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<minutes>(stop - start);
        std::cout << duration.count() << "minutes" << std::endl;

    }

    if (world.rank() < num_process - 1  && world.rank() != 0)
    {
        matrix<realtype> P1(share, d+1);

        for (n = 0; n < P1.size1() ; ++n)
        {
            for ( k = 0; k < d+1; ++k )
            {
                term(k) = factor1( k + 2*(n + (world.rank() * share) ) )* factor2( k );
            }

            for (m = 0; m < P1.size2(); ++m)
            {
                P1(n,m) = factor(m) * inner_prod( project( term, range(0, m+1) ), project( factor3, slice( m, -1, m+1 ) ) );
            }
        }

        world.send(0,17, P1);

    }

    if (world.rank() == num_process - 1)
    {

        matrix<realtype> P1( d+1 - (num_process - 1) * share, d+1);

        for (n = 0; n < P1.size1(); ++n)
        {
            for (k=0; k < d+1; ++k)
            {
                term(k) = factor1( k + 2*(n + ( (num_process - 1)  * share) ) )* factor2( k );
            }

            for (m = 0; m < P1.size2(); ++m)
            {
                P1(n,m) = factor(m) * inner_prod( project( term, range(0, m+1) ), project( factor3, slice( m, -1, m+1 ) ) );
            }
        }

        world.send(0,17, P1);
    }

    return 0;


}
