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
#include<stdlib.h>
#include<time.h>

#define digits 1200

using namespace boost::archive;
using namespace boost::serialization;
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;
using namespace std::chrono;

//serialization code for mpfr_float

using realtype = number<mpfr_float_backend<digits, allocate_stack>>;
#define MPFR_BUFFER_SIZE 1210

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
int omegas = 41; // number of omegas to test
int num_terms = 5000;
const int s = 125, b = 9; // l*(b+1)  = d+1
const int d = s*(b+1) - 1;
const realtype Pi = boost::math::constants::pi<realtype>();
const realtype nu = realtype(2)/realtype(3);
realtype sine = Pi/( pow (realtype(2), nu) * sin(Pi * nu) );

vector<realtype> term (d+1), term1(d+1);
vector<realtype> constants( d+1 );

vector<realtype> factor( d+1 );
vector<realtype> factor1( num_terms );
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

    int num_process = 20;// specify the number of process to use
    int num_threads = boost::thread::hardware_concurrency() ;//use all available threads

    int share = (int) ( floor( ( num_terms - d ) / num_process ) );
    int LO = ( num_terms - d ) - (num_process * share); // number of left over k's beyond what can be covered by (num_process*share)

    vector<realtype> omega(omegas);
    for ( i = 0; i < omegas - 10 ; ++i )
    {
        omega(i) = realtype(1) / pow( realtype(10), realtype( i - 1 ) );
    }
    omega(omegas-10) = realtype(2);
    omega(omegas-9) = pow(realtype(2)/realtype(10),realtype(3))/realtype(4);
    omega(omegas-8) = pow(realtype(4)/realtype(10),realtype(3))/realtype(4);
    omega(omegas-7) = pow(realtype(6)/realtype(10),realtype(3))/realtype(4);
    omega(omegas-6) = pow(realtype(8)/realtype(10),realtype(3))/realtype(4);
    omega(omegas-5) = realtype(1)/realtype(4);
    omega(omegas-4) = pow(realtype(12)/realtype(10),realtype(3))/realtype(4);
    omega(omegas-3) = pow(realtype(14)/realtype(10),realtype(3))/realtype(4);
    omega(omegas-2) = pow(realtype(16)/realtype(10),realtype(3))/realtype(4);
    omega(omegas-1) = pow(realtype(18)/realtype(10),realtype(3))/realtype(4);

    //--------------------factors--------------------------------------------
    /* here we compute the all the relevant factors appearing repeatedly in the
     each iteration of the multiple sum. This factors are computed here once and
     for all and stored later for use. */

    for (i=0; i < factor.size (); ++i)
    {
        factor(i) = realtype(1)/ math::factorial < realtype > ( d - i ); //factorial

        factor2(i) =  realtype(1) /  pow ( math::factorial < realtype > ( i ), realtype(2) );
    }

    for (i = 0; i < num_terms ; ++i)
    {
        factor1(i) = pow( realtype(1)/realtype(2), realtype( i + 1 ) ) /
                        ( boost::math::tgamma( realtype(  i + 2  ) + nu ) );
    }

    //-------------read-in the constants------------------------------------

    std::string st1;
    std::ifstream infile;
    infile.open( "/scratch3/chris.tica/FUNNEL/1_25k_moments_2kdig/Constants/Constant.txt" );
        for (i=0; i < constants.size(); ++i)
    {
        std::getline(infile, st1);
        constants(i) = realtype(st1);
    }

    infile.close();

    //----------------------------------------------------------------------------------

    if (world.rank() == 0)
    {
        std::vector < vector<realtype> > ms1;
        for (i = 0; i < num_process - 1; ++i)
        {
            vector<realtype> fourth1; //for storing the first1's from other processes
            ms1.push_back( fourth1 );
        }

        std::vector < vector<int> > ms2;
        for (i = 0; i < num_process - 1; ++i)
        {
            vector<int> index1; //for storing the index1's from other processes
            ms2.push_back( index1 );
        }

        vector<realtype> fourthter(omegas);
        vector<realtype> total(num_terms - d);

        vector<int> index(share + adder(LO) ); // will contain the indices assigned to this process
        vector<realtype> fourth( index.size() );

        for (k =0; k < index.size() ; ++k)
        {
            index(k) = world.rank() + (num_process)*k;
        }

        for (k = 0; k < index.size(); ++k)
        {
            for( l = 0; l < term.size(); ++l)
            {
                term(l) =  factor2(l)*factor1( index(k) + d  - l );
            }

            for( m = 0; m < d+1; ++m)
            {
                term1(m) = inner_prod( project( term, range(0, m+1) ), project( factor, range (d-m , d + 1) ) );
            }

            fourth( k ) = inner_prod( constants, term1 );
        }

        /* receive all the first1's and indices and store them
        in elements of ms1[i] and ms2[i], respectively.*/
        for ( i = 0; i < num_process - 1 ; ++i )
        {
            world.recv( i + 1, 15, ms1[i] );
            world.recv( i + 1, 16, ms2[i] );
        }

        std::ofstream outfile1;
        outfile1.open("/scratch3/chris.tica/FUNNEL/1_25k_moments_2kdig/results/FOURTH.txt",std::ios_base::out);
        outfile1.precision(digits);

        for (i = 0; i < omegas ; ++i )
        {
            for (j = 0; j < num_process - 1 ; ++j)
            {
                for (k = 0; k < ms1[j].size(); ++k)
                {
                    total( ms2[j][k] ) = pow( omega(i), realtype( d + 1 + ms2[j][k]  ) ) * ms1[j][k];
                }
            }

            for (k = 0; k < index.size(); ++k)
            {
                total( index(k) ) = pow( omega(i), realtype( d + 1 + index(k) ) ) * fourth(k);
            }

            fourthter(i) = sum(total);
            outfile1 << std::scientific << (-sine)*fourthter(i) << std::endl;
        }

    }

    //this will not execute for the case L0 == 0
    if ( world.rank() < LO  && world.rank() > 0 )
    {

        vector<int> index(share + 1); // will contain the indices assigned to this process

        for (k =0; k < index.size()  ; ++k)
        {
            index(k) = world.rank() + (num_process)*k;
        }

        vector<realtype> fourth1( index.size() );

        for ( k = 0; k < fourth1.size(); ++k )
        {
            for( l = 0; l < term.size(); ++l)
            {
                term(l) =  factor2(l)*factor1( index(k) + d  - l );
            }

            for( m = 0; m < d+1; ++m)
            {
                term1(m) = inner_prod( project( term, range(0, m+1) ), project( factor, range (d-m , d + 1) ) );
            }

            fourth1( k ) = inner_prod( constants, term1 );
        }

        world.send( 0 ,15 , fourth1 ); //send the computed values for the assigned k's
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


        vector<realtype> fourth1( index.size() );

        for ( k = 0; k < fourth1.size(); ++k )
        {
            for( l = 0; l < term.size(); ++l)
            {
                term(l) =  factor2(l)*factor1( index(k) + d  - l );
            }

            for( m = 0; m < d+1; ++m)
            {
                term1(m) = inner_prod( project( term, range(0, m+1) ), project( factor, range (d-m , d + 1) ) );
            }

            fourth1( k ) = inner_prod( constants, term1 );
        }

        world.send( 0 ,15 , fourth1 ); //send the computed values for the assigned k's
        world.send( 0 ,16 , index ); // send the indices

    }

    return 0;

}
