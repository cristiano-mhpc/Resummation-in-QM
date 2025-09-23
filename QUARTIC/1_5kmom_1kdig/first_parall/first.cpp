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

#define digits 1000

using namespace boost::archive;
using namespace boost::serialization;
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;
using namespace std::chrono;

//serialization code for mpfr_float

using realtype = number<mpfr_float_backend<digits, allocate_stack>>;
#define MPFR_BUFFER_SIZE 1010

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

int i,j,k,m,n,l;
int omegas = 30; // number of omegas to test
const int s = 150, b = 9; // s*(b+1)  = d+1
const int d = s*(b+1) - 1;
const realtype Pi = boost::math::constants::pi<realtype>();
const realtype nu = realtype(1)/realtype(3);
realtype sine = Pi/( pow (realtype(2), nu) * sin(Pi * nu) );

vector<realtype> constants(d+1), term(d+1), term1(d+1);
vector<realtype> factor (d+1);
vector<realtype> factor1(d+1);
vector<realtype> factor2(d+1);

int adder(int number) //used later to determine if LO == 0 or otherwise.
{
    return (number > 0 ? 1 : 0);
}

int main()
{
    mpi::environment env;
    mpi::communicator world;

    auto start = high_resolution_clock::now();

    int num_process = 5;// specify the number of process to use


    int share = (int) ( floor( ( d + 1 ) / num_process ) );
    int LO = (d + 1) - (num_process * share); // number of left over k's beyond what can be covered by (num_process*share)


    vector<realtype> omega ( omegas );//stores the omegas
    for ( i = 0; i < omegas-2; ++i )
    {
        omega(i) = realtype(1) / pow( realtype(10), realtype( i - 1 ) );
    }
    omega(omegas-1) = realtype(5);
    omega(omegas-1) = realtype(1)/realtype(4);
    //--------------------factor--------------------------------------------
    /* here we compute the all the relevant factors appearing repeatedly in the
     each iteration of the multiple sum. This factors are computed here once and
     for all and stored later for use. */

    for ( i=0; i < factor.size (); ++i )
    {
        factor(i) = realtype(1)/ math::factorial < realtype > ( d - i ); //factorial
        factor1(i) = pow( realtype(1)/realtype(2), realtype( d - i ) ) /
                        ( boost::math::tgamma( realtype( d + 1 - i ) + nu ) );
        factor2(i) =  realtype(1) /  pow ( math::factorial < realtype > ( i ), realtype(2) );
    }

    //---------------read-in the constants------------------------------------

    std::string st1;
    std::ifstream infile;
    infile.open("../Constants/Constant.txt");
    for ( i=0; i < constants.size(); ++i )
    {
        std::getline(infile, st1 );
        constants(i) = realtype( st1 );

    }

    infile.close();

    //---------------------------------------------------------------------------
    if ( world.rank() == 0 )
    {
        std::vector < vector<realtype> > ms1;
        for ( i = 0; i < num_process - 1; ++i )
        {
            vector<realtype> first1; //for storing the first1's from other processes
            ms1.push_back( first1 );
        }


        std::vector < vector<int> > ms2;
        for ( i = 0; i < num_process - 1; ++i )
        {
            vector<int> index1; //for storing the index1's from other processes
            ms2.push_back( index1 );
        }

        vector<realtype>  firstter( omegas );
        vector<realtype>  total(d+1);


        vector<int> index( share + adder(LO) ); // will contain the indices assigned to this process
        vector<realtype>  first( index.size() );

        for (k =0; k < index.size() ; ++k)
        {
            index(k) = world.rank() + ( num_process )*k;
        }

        for ( int k = 0; k < index.size(); ++k )
        {
            for(int l = 0; l < index(k)+1; ++l)//I dont think this should be up to d+1 for all k
            {
                term(l) = factor2(l) * factor1( d - index(k) + l );
            }

            for( int m = 0; m < index(k)+1; ++m )
            {
                term1(m) = inner_prod( project( term, range( 0, m+1 ) ), project( factor, range (d-m, d + 1) ) );
            }

            first( k ) = inner_prod( project ( constants, range (0, index(k)+1) ), project (term1 , range (0, index(k)+1) ) );
        }

        for ( i = 0; i < num_process - 1 ; ++i ) // receive all the first1's and indices and store them in elements of ms1[i] and ms2[i], respectively.
        {
            world.recv( i + 1, 15, ms1[i] );
            world.recv( i + 1, 16, ms2[i] );
        }

        std::ofstream outfile1;
        outfile1.open("../results/FIRST.txt", std::ios_base::out);
        outfile1.precision(digits);
        for (i = 0; i < omegas ; ++i )
        {
            for (j = 0; j < num_process - 1 ; ++j)
            {
                for (k = 0; k < ms1[j].size(); ++k)
                {
                    total( ms2[j][k] ) = pow( omega(i), realtype(  ms2[j][k]  ) ) * ms1[j][k];
                }
            }

            for (k = 0; k < index.size(); ++k)
            {
                total( index(k) ) = pow( omega(i), realtype( index(k) )) * first(k);
            }

            firstter(i) = sum(total);
            outfile1 << std::scientific << ( -sine )*firstter(i) << std::endl;
        }

        //print the total run-time
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<minutes>(stop - start);
        std::cout << duration.count() << "minutes" << std::endl;
    }


    //processes assigned with share + 1 number of indices. The extra 1 comes from the left over indices.
    if ( world.rank() < LO  && world.rank() > 0  )
    {
        //this will not execute for the case L0 == 0

        vector<int> index(share + 1); // will contain the indices assigned to this process

        for ( k =0; k < index.size(); ++k )
        {
            index(k) = world.rank() + ( num_process )*k;
        }

        vector<realtype> first1( index.size() );

        for (int k = 0; k < index.size(); ++k)
        {
            for(int l = 0; l < index(k)+1; ++l)//I dont think this should be up to d+1 for all k
            {
                term(l) = factor2(l)*factor1(d - index(k) + l);
            }

            for( int m = 0; m < index(k)+1; ++m )
            {
                term1(m) = inner_prod( project( term, range( 0, m+1 ) ), project( factor, range (d-m, d + 1) ) );
            }

            first1(k) = inner_prod( project ( constants, range (0, index(k)+1) ), project (term1 , range (0, index(k)+1) ) );
        }

        world.send( 0 ,15 , first1 ); //send the computed values for the assigned k's
        world.send( 0 ,16 , index ); // send the indices
    }


    //processes assigned with share number of indices.

    if ( world.rank() > LO - 1 && world.rank() != 0 )
    {
        vector<int> index(share); // will contain the indices assigned to this process

        for (k =0; k < index.size()  ; ++k)
        {
            index(k) = world.rank() + (num_process)*k;
        }

        vector<realtype> first1( index.size() );
        for ( int k = 0; k < index.size(); ++k )
        {
            for(int l = 0; l < index(k)+1; ++l)//I dont think this should be up to d+1 for all k
            {
                term(l) = factor2(l)*factor1(d - index(k) + l);
            }

            for( int m = 0; m < index(k)+1; ++m )
            {
                term1(m) = inner_prod( project( term, range( 0, m+1 ) ), project( factor, range (d-m, d + 1) ) );
            }

            first1(k) = inner_prod( project ( constants, range (0, index(k)+1 ) ), project (term1 , range ( 0, index(k)+1 ) ) );
        }

        world.send( 0 ,15 , first1 ); //send the computed values for the assigned k's
        world.send( 0 ,16 , index ); // send the indices

    }

    return 0;

}
