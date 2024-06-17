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

realtype waneyt = realtype(1)/ realtype(8);


realtype Mat( int row , int colm )
{

    realtype r = realtype(row);

    if ( row - colm == -6 )
    {
        return  waneyt * sqrt( ( r + realtype(1) ) * ( r + realtype(2) ) *( r + realtype(3) ) * ( r + realtype(4) ) * ( r + realtype(5) ) * ( r + realtype(6) ) );
    }
    if ( row - colm == -4 )
    {
        return realtype(3) * waneyt * sqrt( ( r + realtype(1) )* ( r + realtype(2) ) *( r + realtype(3) ) * ( r + realtype(4) ) ) * ( realtype(2)*r + realtype(5) );
    }
    else if (row - colm == -2)
    {
        return realtype(15) * waneyt * sqrt( (r + realtype(1) ) * ( r + realtype(2) ) ) * ( pow(r, realtype(2) ) + realtype(3)*r + realtype(3) );
    }
    else if ( row - colm == 0)
    {
        return realtype(5) * waneyt * ( realtype(4) * pow(r , realtype(3) ) + realtype(6)* pow(r , realtype(2) ) + realtype(8)*r  + realtype(3) );
    }
    else if ( row - colm == 2 )
    {
        return realtype(15) * waneyt * sqrt( r * ( r - realtype(1) ) ) * ( pow( r, realtype(2)) - r + realtype(1) );
    }
    else if ( row - colm == 4 )
    {
        return realtype(3) * ( waneyt ) * sqrt( r * ( r - realtype(1) ) * ( r - realtype(2) ) * ( r - realtype(3) ) ) * ( realtype(2)*r - realtype(3) );
    }

    if ( row - colm == 6 )
    {
        return waneyt * sqrt( r * ( r - realtype(5) ) * ( r - realtype(4) ) * ( r - realtype(3) ) * ( r - realtype(2) ) * ( r - realtype(1) ) );
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
    const int d = 5000; // code below works for even d;
    matrix <realtype> a(d + 1, 3*(d/2) + 1);
    matrix <realtype> mat( 3*(d/2) + 1, 3*(d/2) + 1 );// zero index included

    vector<realtype> moments( d+1 );

    /*careful here to use the correct upper limit in populating the matrix.
    I used m<d . Took a while before I caught it.*/

    for ( m = 1; m < 3*(d/2) + 1; ++m )
    {
        for ( n = std::max(1, m - 3); n < std::min(3*(d/2) + 1, m + 4) ; ++n )
        {
            mat(m, n) = Mat( 2*m, 2*n );

        }
    }


    //initial values
    moments(0) = realtype(1);
    moments(1) = realtype(15)/realtype(8);

    a (1, 1) = -( realtype(45) * sqrt( realtype(2) ) ) / realtype(32);
    a (1, 2) = -( realtype(15) * sqrt( math::factorial < realtype > ( 4 ) ) / realtype(64) );
    a (1, 3) = -( sqrt( math::factorial < realtype > ( 6 ) ) / realtype(96) );
    //along each row, the elements can be computed independent of each other. This is an opportunity for parallelization

    ///--------------------------------first half of the triangle-----------------------------///

    for ( int r = 2; r < d/2 + 1; ++r )// assume d is even for now
    {
        moments(r) = ( realtype(45) * sqrt( realtype(2) )  * waneyt * a (r - 1, 1) )
                        + ( realtype(15) * sqrt( math::factorial < realtype > ( 4 ) ) * waneyt * a (r - 1, 2) )
                        + ( sqrt( math::factorial < realtype > ( 6 ) ) * waneyt * a (r - 1, 3) );

        // first column m = 1. Computation requires only 4 terms: on top and 3 to the right.
        a ( r, 1 ) =  -( realtype(1) /  realtype(4) ) * ( mat( 1, 1 )* a(r - 1, 1 ) + mat( 1, 1 + 1 ) * a( r - 1, 1 + 1 ) +  mat( 1, 1 + 2 ) * a( r - 1, 1 + 2 ) + mat( 1, 1 + 3 ) * a( r - 1, 1 + 3 )
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 1, 1 , r - 1 ), slice( 1, 0, r-1 ) ),
                        project( moments , slice ( r - 1 , -1 , r - 1 ) ) ) );



        // second column m = 2. Computation requires 5 terms: on top, one to its left, and 3 elements to the right.
        a ( r, 2 ) =  -( realtype(1) / ( realtype(4)*realtype(2) ) ) * ( mat( 2, 2-1 ) * a( r - 1, 2 - 1 ) + mat( 2, 2 )* a(r - 1, 2 ) + mat( 2, 2 + 1 ) * a( r - 1, 2 + 1 )
                        +  mat( 2, 2 + 2 ) * a( r - 1, 2 + 2 ) + mat( 2, 2 + 3 ) * a( r - 1, 2 + 3 )
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 1, 1 , r - 1 ), slice( 2, 0, r-1 ) )
                        , project( moments , slice ( r - 1 , -1 , r - 1 ) ) ) );



        // third column m = 3. Computation requires 6 terms: on top, 2 to its left, and 3 elements to the right.
        a(r, 3) = -( realtype(1) / ( realtype(4)*realtype(3) ) ) * ( mat( 3, 3-2 ) * a( r - 1, 3 - 2 ) + mat( 3, 3-1 ) * a( r - 1, 3 - 1 ) + mat( 3, 3 )* a(r - 1, 3 ) + mat( 3, 3 + 1 ) * a( r - 1, 3 + 1 )
                        +  mat( 3, 3 + 2 ) * a( r - 1, 3 + 2 ) + mat( 3, 3 + 3 ) * a( r - 1, 3 + 3 )
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 1, 1 , r - 1 ), slice( 3, 0, r-1 ) )
                        , project( moments , slice ( r - 1 , -1 , r - 1 ) ) ) );



        for ( m = 4 ; m < 3*( r - 1 ) + 1 ; ++m )
        {
            a ( r, m ) = -( realtype(1) / ( realtype(4)*realtype(m) ) ) * ( mat( m, m-3 ) * a( r - 1, m - 3 ) + mat( m, m-2 ) * a( r - 1, m - 2 )  +  mat( m, m-1 ) * a( r - 1, m - 1 )  + mat( m, m )* a(r - 1, m )
                        + mat( m, m + 1 ) * a( r - 1, m + 1 ) +  mat( m, m + 2 ) * a( r - 1, m + 2 ) + mat( m, m + 3 ) * a( r - 1, m + 3 )
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 2 + Floor(  ( m - 4 ) / 3 ), 1 , ( r - 2 ) - Floor(  ( m - 4 ) / 3 ) ),
                        slice( m, 0, (r - 2) - Floor(  ( m - 4 ) / 3 ) ) ),
                        project( moments , slice ( r - 2 - Floor(  ( m - 4 ) / 3 ) , -1 , ( r - 2 ) - Floor(  ( m - 4 ) / 3 ) ) ) ) );


        }

        // third to the last column m = 3*(r - 1) + 1. Computation requires only 3 terms to the left
        a ( r, 3*(r - 1) + 1 ) =  -( realtype(1) / ( realtype(4)* ( realtype(3)*realtype(r) - realtype(2) ) ) )  * ( mat( 3*(r - 1) + 1, 3*(r - 1) + 1 - 3 ) * a( r - 1, 3*(r - 1) + 1 - 3 ) +
                                    mat( 3*(r - 1) + 1, 3*(r - 1) + 1 - 2 ) * a( r - 1, 3*(r - 1) + 1 - 2 )+  mat( 3*(r - 1) + 1, 3*(r - 1) + 1 - 1 ) * a( r - 1, 3*(r - 1) + 1 - 1 ) );



        // second to the last column m = 3*(r - 1) + 2. Computation requires only the left-most 2 terms
        a ( r, 3*(r - 1) + 2 ) =  -( realtype(1) / ( realtype(4)* ( realtype(3)*realtype(r) - realtype(1) ) ) )  * ( mat( 3*(r - 1) + 2, 3*(r - 1) + 2 - 3 ) * a( r - 1, 3*(r - 1) + 2 - 3 ) +
                                    mat( 3*(r - 1) + 2, 3*(r - 1) + 2 - 2 ) * a( r - 1, 3*(r - 1) + 2 - 2 ) );



         // last column m = 3*r. Computation requires only 1 term farthest the left.
        a ( r, 3*r ) = - ( realtype(1) / ( realtype(4) * realtype(3) * realtype(r) ) ) * ( mat( 3*r, 3*r - 3 ) * a( r - 1, 3*r - 3 ) );


    }

    ///--------------------------------the other half of the triangle-------------------

    for ( int r = (d/2) + 1; r < d; ++r )// assume d is even for now
    {
        moments(r) = ( realtype(45) * sqrt( realtype(2) )  * waneyt * a (r - 1, 1) )
                        + ( realtype(15)*sqrt( math::factorial < realtype > ( 4 ) ) * waneyt * a (r - 1, 2) )
                        + ( sqrt( math::factorial < realtype > ( 6 ) ) * waneyt * a (r - 1, 3) );

        //first column m = 1. Computation requires only 4 terms: on top and 3 to the right.
        a ( r, 1 ) =  -( realtype(1) /  realtype(4) ) * ( mat( 1, 1 )* a(r - 1, 1 ) + mat( 1, 1 + 1 ) * a( r - 1, 1 + 1 ) +  mat( 1, 1 + 2 ) * a( r - 1, 1 + 2 ) + mat( 1, 1 + 3 ) * a( r - 1, 1 + 3 )
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 1, 1 , r - 1 ), slice( 1, 0, r-1 ) ),
                        project( moments , slice ( r - 1 , -1 , r - 1 ) ) ) );



        // second column m = 2. Computation requires 5 terms: on top, one to its left, and 3 elements to the right.
        a ( r, 2 ) =  -( realtype(1) / ( realtype(4)*realtype(2) ) ) * ( mat( 2, 2-1 ) * a( r - 1, 2 - 1 ) + mat( 2, 2 )* a(r - 1, 2 ) + mat( 2, 2 + 1 ) * a( r - 1, 2 + 1 )
                        +  mat( 2, 2 + 2 ) * a( r - 1, 2 + 2 ) + mat( 2, 2 + 3 ) * a( r - 1, 2 + 3 )
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 1, 1 , r - 1 ), slice( 2, 0, r-1 ) )
                        , project( moments , slice ( r - 1 , -1 , r - 1 ) ) ) );



        // third column m = 3. Computation requires 6 terms: on top, 2 to its left, and 3 elements to the right.
        a(r, 3) = -( realtype(1) / ( realtype(4)*realtype(3) ) ) * ( mat( 3, 3-2 ) * a( r - 1, 3 - 2 ) + mat( 3, 3-1 ) * a( r - 1, 3 - 1 ) + mat( 3, 3 )* a(r - 1, 3 ) + mat( 3, 3 + 1 ) * a( r - 1, 3 + 1 )
                        +  mat( 3, 3 + 2 ) * a( r - 1, 3 + 2 ) + mat( 3, 3 + 3 ) * a( r - 1, 3 + 3 )
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 1, 1 , r - 1 ), slice( 3, 0, r-1 ) )
                        , project( moments , slice ( r - 1 , -1 , r - 1 ) ) ) );



        //other columns beyond 3: 4 <= m <= 3(d -r) .Requires all 7 relevant elements.
        for ( m = 4 ; m <  3*(d - r)  + 1; ++m )
        {
            a ( r, m ) = -( realtype(1) / ( realtype(4)*realtype(m) ) ) * ( mat( m, m-3 ) * a( r - 1, m - 3 ) + mat( m, m-2 ) * a( r - 1, m - 2 )  +  mat( m, m - 1 ) * a( r - 1, m - 1 )  + mat( m, m )* a(r - 1, m )
                        + mat( m, m + 1 ) * a( r - 1, m + 1 ) +  mat( m, m + 2 ) * a( r - 1, m + 2 ) + mat( m, m + 3 ) * a( r - 1, m + 3 )
                        - inner_prod( matrix_vector_slice< matrix<realtype> > ( a, slice ( 2 + Floor(  ( m - 4 ) / 3 ), 1 , ( r - 2 ) - Floor(  ( m - 4 ) / 3 ) ),
                        slice( m, 0, (r - 2) - Floor(  ( m - 4 ) / 3 ) ) ),
                        project( moments , slice ( r - 2 - Floor(  ( m - 4 ) / 3 ) , -1 , ( r - 2 ) - Floor(  ( m - 4 ) / 3 ) ) ) ) );


        }

    }


    //since the entire matrix has now been constructed, the computation of the terms below can now be done in parallel. In the future, this can be done using multiple processes.
    moments(d) = ( realtype(45) * sqrt( realtype(2) )  * waneyt * a (d - 1, 1) )
                        + ( realtype(15)*sqrt( math::factorial < realtype > ( 4 ) ) * waneyt * a (d - 1, 2) )
                        + ( sqrt( math::factorial < realtype > ( 6 ) ) * waneyt * a (d - 1, 3) );

    // print moments()
    std::ofstream outfile1;
    outfile1.open( "moments.txt",std::ios_base::out );
    outfile1.precision( digits - 5 );

    for (i = 1; i < moments.size(); ++i)
    {
        outfile1 << std::scientific << pow( realtype(-1), realtype(i+1) ) * moments(i) << std::endl;
    }

    outfile1.close();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<hours>(stop - start);
    std::cout << duration.count() << "hours" << std::endl;

    return 0;

}
