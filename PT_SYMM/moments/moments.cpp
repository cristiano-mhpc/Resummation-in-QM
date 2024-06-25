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
#include <boost/math/special_functions/factorials.hpp>

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
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

//for serialization and miscellaneous functionalities
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <time.h>


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
namespace math = boost::math;

int i,j,k,l,m,n;

realtype waneyt = realtype(1)/sqrt(realtype(8));

realtype Mat( int row , int colm )
{

    realtype r = realtype(row);

    if ( row - colm == -3 )
    {
        return waneyt * sqrt( ( r + realtype(3) )*( r + realtype(2) ) *
                    ( r + realtype(1) ) );
    }
    else if ( row - colm == -1 )
    {
        return  realtype(3) * waneyt * ( r + realtype(1) ) * sqrt( r + realtype(1) ) ;
    }
    else if ( row - colm == 1 )
    {
        return realtype(3) * waneyt * r * sqrt( r );
    }
    else if ( row - colm == 3 )
    {
        return waneyt * sqrt( ( r - realtype(2) ) * ( r - realtype(1) ) * r );
    }

    else
    {
        return realtype(0);
    }

}

int Floor( int x ) //returns  an int floor since return of std::floor() is double.
{
    return (int) std::floor(x);
}

int main()
{

    auto start = high_resolution_clock::now();

    //Declare some variables
    const int d = 10020; // code below works for even d;
    matrix <realtype> a( d+1, 3*(d/2) + 1 );

    vector<realtype> moments( d/2 + 1 );

    //initial values
    moments(0) = realtype(1);

    a (1, 1) = -( realtype(3) * sqrt( realtype(8) ) ) / realtype(16);

    a (1, 3) = -( sqrt( realtype(48) ) / realtype(48) );

    a (2, 2) = -( realtype(1) / ( realtype(2)*realtype(2) ) ) * ( Mat(2 , 2 - 1) * a(2-1 , 2-1) + Mat(2, 2+1) * a(2-1, 2 + 1) );

    a (2, 4) = -( realtype(1) / ( realtype(2)*realtype(4) ) ) * ( Mat(4, 4 - 3) * a(2-1 , 4-3 ) + Mat(4, 4 - 1) * a(2-1 , 4 - 1 ) );

    a (2, 6) = -( realtype(1) / ( realtype(2)*realtype(6) ) ) * ( Mat(6, 6 - 3) * a(2-1, 6 - 3) );

    moments(1) = ( ( realtype(3) / sqrt( realtype(8) ) ) * a ( 2 - 1, 1 ) )
                        + ( (  sqrt( realtype(6) ) / sqrt( realtype(8) )  ) * a ( 2 - 1, 3 ) );

    //upper triangle
    for ( int r = 3; r < d/2 + 1; ++r )
    {
        if ( r%2 == 0 ) //r is even
        {
            //moments(r/2) = E(r)
            moments(r/2) = ( (  realtype(3) / sqrt( realtype(8) ) ) * a (r - 1, 1) )
                        + ( ( sqrt( realtype(6) ) / sqrt( realtype(8) ) ) * a ( r - 1, 3 ) );

            //first three. Requires the all elements a(r, m) above
            //m = 2
            a(r, 2) = -( realtype(1) / ( realtype(2)*realtype(2) ) ) * ( Mat(2 , 2 - 1 ) * a(r-1 , 2 - 1) + Mat(2, 2 + 1 ) * a(r-1, 2 + 1) +  Mat(2, 2 + 3) * a(r-1, 2 + 3)
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 2, 2 , r/2 - 1 ), slice( 2, 0, r/2 - 1 ) ),
                        project( moments , slice ( r/2 - 1 , -1 , r/2 - 1 ) ) ) );
            //m = 4
            a(r, 4) = -( realtype(1) / ( realtype(2)*realtype(4) ) ) * ( Mat(4, 4 - 3) * a(r-1 , 4 - 3 ) + Mat(4 , 4 - 1 )*a(r-1 , 4 - 1) + Mat(4, 4 + 1)*a(r-1, 4 + 1) +  Mat(4, 4 + 3)*a(r-1, 4 + 3)
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 2, 2 , r/2 - 1 ), slice( 4, 0, r/2 - 1 ) ),
                        project( moments , slice ( r/2 - 1 , -1 , r/2 - 1 ) ) ) );
            //m = 8
            a(r, 6) = -( realtype(1) / ( realtype(2)*realtype(6) ) ) * ( Mat( 6, 6 - 3 )*a(r-1 , 6 - 3 ) + Mat(6 , 6 - 1 )*a(r-1 , 6 - 1) + Mat(6, 6 + 1)*a(r-1, 6 + 1) +  Mat(6, 6 + 3)*a(r-1, 6 + 3)
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 2, 2 , r/2 - 1 ), slice( 6, 0, r/2 - 1 ) ),
                        project( moments , slice ( r/2 - 1 , -1 , r/2 - 1 ) ) ) );

            //in-between
            for (m = 8; m < 3*(r - 2) + 1; m+=2) // requires 2 alements to the left and 2 to the right from row above it
            {
                a(r , m) = -( realtype(1) / ( realtype(2)*realtype(m) ) ) * ( Mat(m, m - 3)*a(r - 1, m - 3) + Mat(m , m-1) * a(r - 1, m - 1) +
                        Mat(m , m + 1) * a(r - 1, m + 1) + Mat(m , m + 3) * a(r - 1, m + 3)
                        -  inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 4 + 2*Floor(  ( m - 8 ) / 6 ), 2 , ( r/2 - 2 ) - Floor(  ( m - 8 ) / 6 ) ), slice( m, 0, ( r/2 - 2 ) - Floor(  ( m - 8 ) / 6 )  ) ),
                        project( moments , slice ( r/2 - 2 - Floor(  ( m - 8 ) / 6 )  , -1 , r/2 - 2 - Floor(  ( m - 8 ) / 6 )  ) ) ) );
            }

            //last three. Dont require elements directly above.
            // m = 3r - 4
            a (r , 3*(r - 2) + 2 ) = -( realtype(1) / ( realtype(2)*realtype(3*r - 4) ) ) * ( Mat(3*r - 4, 3*r - 4 - 3)*a(r - 1, 3*r - 4 - 3) + Mat(3*r - 4 , 3*r - 4 - 1 ) * a(r - 1, 3*r - 4 - 1) +  Mat(3*r - 4 , 3*r - 4 + 1) * a(r - 1, 3*r - 4 + 1) );

            // m = 3r - 2
            a (r , 3*r - 2) = -( realtype(1) / ( realtype(2)*realtype(3*r - 2) ) ) * ( Mat(3*r - 2 ,3*r - 2 - 3)*a(r - 1, 3*r - 2 - 3) + Mat(3*r - 2 , 3*r - 2 - 1) * a(r - 1, 3*r - 2 - 1)  );

            // m = 3r
            a (r , 3*r ) = -( realtype(1) / ( realtype(2)*realtype(3*r ) ) ) * ( Mat(3*r ,3*r - 3)*a(r - 1, 3*r - 3) );
        }

        else //r is odd
        {
            //first two. Requires the all elements a(r, m) above
            a(r , 1) = -( realtype(1) / ( realtype(2)*realtype(1) ) ) * ( Mat(1, 1 + 1) * a(r-1, 1 + 1) +  Mat(1, 1 + 3) * a(r-1, 1 + 3)
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 1, 2 , (r+1)/2 - 1 ), slice( 1, 0, (r+1)/2 - 1 ) ),
                        project( moments , slice ( (r+1)/2  - 1  , -1 , (r+1)/2 - 1 ) ) ) );


            a(r , 3) = -( realtype(1) / ( realtype(2)*realtype(3) ) ) * ( Mat(3, 3 - 1) * a(r-1, 3 - 1) + Mat(3, 3 + 1) * a(r - 1, 3 + 1) +  Mat(3, 3 + 3) * a(r - 1, 3 + 3)
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 1, 2 , (r+1)/2 - 1 ), slice( 3, 0, (r+1)/2 - 1 ) ),
                        project( moments , slice ( (r+1)/2  - 1  , -1 , (r+1)/2 - 1 ) ) ) );

            //in-between. requires 2 alements to the left and 2 to the right from row above it
            for (m = 5; m < 3*r - 6 + 1; m+= 2)
            {
                 a(r , m) = -( realtype(1) / ( realtype(2)*realtype(m) ) ) * ( Mat(m ,m - 3)*a(r - 1, m - 3) + Mat(m , m-1) * a(r - 1, m - 1) +
                        Mat(m , m + 1) * a(r - 1, m + 1) + Mat(m , m + 3) * a(r - 1, m + 3)
                        -  inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 3 + 2 * Floor(  ( m - 5 ) / 6 ), 2 , ( (r + 1)/2 - 2 ) - Floor(  ( m - 5 ) / 6 ) ),
                         slice( m, 0, ( (r + 1)/2 - 2 ) - Floor(  ( m - 5 ) / 6 ) ) ),
                        project( moments , slice ( (r + 1)/2 - 2 - Floor(  ( m - 5 ) / 6 ), -1 , (r + 1)/2 - 2 - Floor(  ( m - 5 ) / 6 )  ) ) ) );

            }

            //last three. Dont require elements directly above.

            // m = 3*r - 4
            a(r, 3*r - 4)  = -( realtype(1) / ( realtype(2)*realtype(3*r - 4) ) ) * ( Mat(3*r - 4 ,3*r - 4 - 3)*a(r - 1, 3*r - 4 - 3) + Mat(3*r - 4 , 3*r - 4 - 1) * a(r - 1, 3*r - 4 - 1) + Mat(3*r - 4 , 3*r - 4 + 1) * a(r - 1, 3*r - 4 + 1) );

            // m = 3*r - 2
            a(r, 3*r - 2) = -( realtype(1) / ( realtype(2)*realtype(3*r - 2) ) ) * ( Mat(3*r - 2 ,3*r - 2 - 3)*a(r - 1, 3*r - 2 - 3) + Mat(3*r - 2 , 3*r - 2 -1) * a(r - 1, 3*r - 2 - 1)  );

            // m = 3r
            a (r , 3*r) = -( realtype(1) / ( realtype(2)*realtype(3*r ) ) ) * ( Mat(3*r ,3*r - 3) * a(r - 1, 3*r - 3) );

        }


    }

    // the other half of the triangle. The columns for each row r, is up to m = 3*(d-r)

    for ( int r = (d/2) + 1; r < d; ++r )
    {

        if (r%2 == 0 ) //r is even
        {
            //moments(r/2) = E(r)
            moments(r/2) = (  realtype(3) / sqrt( realtype(8) ) * a (r - 1, 1) )
                        + ( (sqrt( realtype(6) ) / sqrt( realtype(8) ) ) * a ( r - 1, 3 ) );

            //first three. Requires the all elements a(r, m) above
            //m = 2
            a(r, 2) = -( realtype(1) / ( realtype(2)*realtype(2) ) ) * ( Mat(2 , 2 - 1 )*a(r-1 , 2 - 1) + Mat(2, 2+1)*a(r-1, 2 + 1) +  Mat(2, 2+3)*a(r-1, 2 + 3)
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 2, 2 , r/2 - 1 ), slice( 2, 0, r/2 - 1 ) ),
                        project( moments , slice ( r/2 - 1 , -1 , r/2 - 1 ) ) ) );
            //m = 4
            a(r, 4) = -( realtype(1) / ( realtype(2)*realtype(4) ) ) * ( Mat(4, 4 - 3) * a(r-1 , 4-3 ) + Mat(4 , 4 - 1 )*a(r-1 , 4 - 1) + Mat(4, 4 + 1)*a(r-1, 4 + 1) +  Mat(4, 4 + 3)*a(r-1, 4 + 3)
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 2, 2 , r/2 - 1 ), slice( 4, 0, r/2 - 1 ) ),
                        project( moments , slice ( r/2 - 1 , -1 , r/2 - 1 ) ) ) );
            //m = 8
            a(r, 6) = -( realtype(1) / ( realtype(2)*realtype(6) ) ) * ( Mat( 6, 6 - 3 )*a(r-1 , 6-3 ) + Mat(6 , 6 - 1 )*a(r-1 , 6 - 1) + Mat(6, 6 + 1)*a(r-1, 6 + 1) +  Mat(6, 6 + 3)*a(r-1, 6 + 3)
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 2, 2 , r/2 - 1 ), slice( 6, 0, r/2 - 1 ) ),
                        project( moments , slice ( r/2 - 1 , -1 , r/2 - 1 ) ) ) );

            //in-between
            for (m = 8; m < 3*( d - r ) + 1; m+=2) // requires 2 alements to the left and 2 to the right from row above it
            {
                a(r , m) = -( realtype(1) / ( realtype(2)*realtype(m) ) ) * ( Mat(m, m - 3)*a(r - 1, m - 3) + Mat(m , m-1) * a(r - 1, m - 1) +
                        Mat(m , m + 1) * a(r - 1, m + 1) + Mat(m , m + 3) * a(r - 1, m + 3)
                        -  inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 4 + 2*Floor(  ( m - 8 ) / 6 ), 2 , ( r/2 - 2 ) - Floor(  ( m - 8 ) / 6 ) ), slice( m, 0, ( r/2 - 2 ) - Floor(  ( m - 8 ) / 6 )  ) ),
                        project( moments , slice ( r/2 - 2 - Floor(  ( m - 8 ) / 6 ), -1 , r/2 - 2 - Floor(  ( m - 8 ) / 6 )  ) ) ) );
            }

        }


        else //r is odd
        {

            //first two. Requires the all elements a(r, m) above
            a(r , 1) = -( realtype(1) / ( realtype(2)*realtype(1) ) ) * ( Mat(1, 1+1) * a(r-1, 1 + 1) +  Mat(1, 1+3) * a(r-1, 1 + 3)
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 1, 2 , (r+1)/2 - 1 ), slice( 1, 0, (r+1)/2 - 1 ) ),
                        project( moments , slice ( (r+1)/2  - 1  , -1 , (r+1)/2 - 1 ) ) ) );


            a(r , 3) = -( realtype(1) / ( realtype(2)*realtype(3) ) ) * ( Mat(3, 3 - 1) * a(r-1, 3 - 1) + Mat(3, 3 + 1) * a(r - 1, 3 + 1) +  Mat(3, 3 + 3) * a(r - 1, 3 + 3)
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 1, 2 , (r+1)/2 - 1 ), slice( 3, 0, (r+1)/2 - 1 ) ),
                        project( moments , slice ( (r+1)/2  - 1  , -1 , (r+1)/2 - 1 ) ) ) );

            //in-between. requires 2 alements to the left and 2 to the right from row above it
            for (m = 5; m < 3*( d - r ) + 1; m+= 2)
            {
                 a(r , m) = -( realtype(1) / ( realtype(2)*realtype(m) ) ) * ( Mat(m ,m - 3)*a(r - 1, m - 3) + Mat(m , m - 1) * a(r - 1, m - 1) +
                        Mat(m , m + 1) * a(r - 1, m + 1) + Mat(m , m + 3) * a(r - 1, m + 3)
                        -  inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 3 + 2 * Floor(  ( m - 5 ) / 6 ), 2 , ( (r + 1)/2 - 2 ) - Floor(  ( m - 5 ) / 6 ) ),
                         slice( m, 0, ( (r + 1)/2 - 2 ) - Floor(  ( m - 5 ) / 6 ) ) ),
                        project( moments , slice ( (r + 1)/2 - 2 - Floor(  ( m - 5 ) / 6 ) , -1 , (r + 1)/2 - 2 - Floor(  ( m - 5 ) / 6 )  ) ) ) );

            }

        }


    }

    moments(d/2) = ( realtype(3) / sqrt( realtype(8) ) * a (d - 1, 1) )
                        + ( (sqrt( realtype(6) ) / sqrt( realtype(8) ) ) * a ( d - 1, 3 ) );

    // print moments()
    std::ofstream outfile1;
    outfile1.open( "/home/chris.tica/PTSYM/moments_3k/moments.txt",std::ios_base::out );
    outfile1.precision( digits - 5 );


    for ( i = 1; i < moments.size(); ++i )
    {
        outfile1 << std::scientific << realtype(-1) * moments(i) << std::endl;
    }

    outfile1.close();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<hours>(stop - start);
    std::cout << duration.count() << "hours" << std::endl;

    std::cout << "d = " << d << std::endl;

    return 0;

}
