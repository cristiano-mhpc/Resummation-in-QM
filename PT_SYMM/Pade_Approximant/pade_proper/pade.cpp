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

int i,j,k,m,n,l;

int betas = 30; // number of omegas to test
const realtype e = boost::math::constants::e<realtype>();
const int N = 1000;
const int M = N+1;
int main()
{

    vector<realtype>  pade(betas), num(betas), total(N + 1), denom(betas), total2(M + 1), A(N + 1), coeff1( N + 1 ), B(M + 1);

    vector<realtype> beta(betas);
    for (i = 0; i < betas-2; ++i)
    {
        beta(i) =  pow( realtype(10), realtype( i - 1 ) );
    }
    beta(betas-2) = realtype(1)/realtype(5);
    beta(betas-1) = realtype(4);

    //------------------read-in the first half of coefficients--------------------------------------
    std::string st;
    std::ifstream infile;
    infile.open( "moments.txt", std::ios_base::in );
    for ( int i = 0; i < coeff1.size(); ++i )
    {
        std::getline(infile, st);
        coeff1(i) = pow(realtype(-1), realtype(i))*realtype(st);
    }

    infile.close();
    infile.clear();

    //------------------read-in the B(i)-----------------------------------------
    infile.open("/scratch3/chris.tica/PTSYM/Pade/Constants/Constant.txt", std::ios_base::in);
    B(0) = realtype(1);
    std::string st1;
    for (i = 1; i < B.size(); ++i)
    {
        std::getline(infile, st1);
        B(i) = realtype( st1 );
    }

    infile.close();
    infile.clear();

    //------------------compute A(i)--------------------------------
    for (n = 0; n < A.size(); ++n )
    {
        A(n) = inner_prod( project(coeff1, slice(n, -1, n + 1) ), project( B, range(0, n+1) ) );
    }

    //----compute the Pade Approximants evaluated at various beta and print them----

    std::ofstream outfile1;
    outfile1.open("pade.txt",std::ios_base::out );
    outfile1.precision( digits - 5 );

    for ( i=0; i < betas; ++i )
    {
        for ( k = 0; k < A.size(); ++k )
        {
            total(k) = pow( beta(i), realtype(k) ) * A(k);//numerator

        }

        for ( k = 0; k < B.size(); ++k )
        {
            total2(k) = pow( beta(i), realtype(k) ) * B(k);//denominator
        }

        pade(i) = sum(total)/sum(total2);

        outfile1 << std::scientific << realtype(1) + (beta(i) * pade(i)) << std::endl;
    }

    return 0;

}
