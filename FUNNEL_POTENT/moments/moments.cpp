/*We compute here the perturbation coefficients for the gound-state energy

    E = Sum_{k=0}^{\infty} E_k (-\mu)^k

of the funnel like potential V(r) = 1/r + gr, \mu = g^{1/2} here is the coupling parameter. The funnel
potential is an example of a screened-Coulumb potential.  The Hamiltonian is
H(r) = (p^2)/2 + V(r). The recurrence relation defining the coefficients E(k) are

E(k) = -3/2 a(k,1) = -3 a(k,2)

The elements of the matrix a(k,j) is

a(k,k-1) = f_k , where f_1 = 1, f_2 = -1 f_k=0 for k>=3.

a(k,j)= (j+3)a(k, j+1) + Sum_{l=2}^{k-2}Sum_{p+q=j} a(l,p)a(k-l,q), for 1<j<=k-2

*/


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


int main()
{

    auto start = high_resolution_clock::now();

    //Declare some variables
    const int d = 1250;
    matrix <realtype> a(d + 1, d + 1);
    a(0,2) = realtype(1)/realtype(3);
    a(2,1) = realtype(1);
    a(4,2) = realtype(1)/realtype(2);
    a(4,1) = realtype(2)*a(4,2);
    ///--------------------------------We compute the matrix a(k,j) here-----------------------------///
    vector<realtype> prod(d+1);

    for (i = 3; i < (d/2)+1 ; ++i) // k = 2i we only compute along even rows
    {
        for (n = i; n < 2*i-1 ; ++n) //from a(2i,i) to a(2i,2). The comlumn in a(k,j) here is j = 2i-n. So j is specified by n.
        {

            //here we separate the cases l = 2 and l = k-2

            prod(2) = inner_prod( matrix_vector_slice< matrix<realtype> > (a, slice(2, 0, std::max( 0, 1 +  std::min(2*i-n-1, 2-1 ) - std::max(1, 2*i-n-(2*i-2-1) ) ) ),
                                slice( std::max(1, 2*i - n -( 2*i-2-1) ), 1 , std::max( 0, 1 +  std::min(2*i-n-1, 2-1 ) - std::max(1, 2*i-n-(2*i-2-1) ) )  ) ),
                                    matrix_vector_slice< matrix<realtype> > ( a, slice( 2*i-2, 0, std::max( 0, 1 +  std::min(2*i-n-1, 2-1 ) - std::max(1, 2*i-n-(2*i-2-1) ) )  ),
                                        slice( 2*i - n - std::max(1, 2*i-n-( 2*i-2-1) ), -1, std::max( 0, 1 +  std::min(2*i-n-1, 2-1 ) - std::max(1, 2*i-n-(2*i-2-1) ) )  ) ) );
            for (l = 3; l < 2*i-2; ++l)
            {
                prod(l) = inner_prod( matrix_vector_slice< matrix<realtype> > (a, slice(l, 0, std::max( 0, 1 +  std::min(2*i-n-1, l-2) - std::max(1, 2*i-n-(2*i-l-2) ) ) ),
                                slice( std::max(1, 2*i - n -( 2*i-l-2) ), 1 , std::max( 0, 1 +  std::min(2*i-n-1, l-2 ) - std::max(1, 2*i-n-(2*i-l-2) ) )  ) ),
                                    matrix_vector_slice< matrix<realtype> > ( a, slice( 2*i-l, 0, std::max( 0, 1 +  std::min(2*i-n-1, l-2 ) - std::max(1, 2*i-n-(2*i-l-2) ) )  ),
                                        slice( 2*i - n - std::max(1, 2*i-n-( 2*i-l-2) ), -1, std::max( 0, 1 +  std::min(2*i-n-1, l-2 ) - std::max(1, 2*i-n-(2*i-l-2) ) )  ) ) );
            }

            prod(2*i-2) = inner_prod( matrix_vector_slice< matrix<realtype> > (a, slice(2*i-2, 0, std::max( 0, 1 +  std::min(2*i-n-1, 2*i-2-1 ) - std::max(1, 2*i-n-(2-1) ) ) ),
                                slice( std::max(1, 2*i - n -(2-1) ), 1 , std::max( 0, 1 +  std::min(2*i-n-1, 2*i-2-1 ) - std::max(1, 2*i-n-(2-1) ) )  ) ),
                                    matrix_vector_slice< matrix<realtype> > ( a, slice( 2, 0, std::max( 0, 1 +  std::min(2*i-n-1, 2*i-2-1 ) - std::max(1, 2*i-n-(2-1) ) )  ),
                                        slice( 2*i - n - std::max(1, 2*i-n-( 2-1) ), -1, std::max( 0, 1 +  std::min(2*i-n-1, 2*i-2-1 ) - std::max(1, 2*i-n-(2-1) ) )  ) ) );

            a(2*i, 2*i - n) = (realtype(1)/realtype(2)) * (realtype(2*i - n + 3) * a(2*i, 2*i - n + 1) + sum( project(prod, range(2, 2*i-1) ) ) ); //Along each rows, we computing the elements near the diagonal first.
        }

        a(2*i,1) = realtype(2)*a(2*i,2);
    }


    // print moments()
    std::ofstream outfile1;
    std::ofstream outfile2;

    outfile1.open("moments.txt",std::ios_base::out);
    outfile1.precision( digits - 5 );


    //outfile1 << std::scientific << -realtype(1) << std::endl;
    outfile1 << std::scientific << -realtype(3) << std::endl;

    for (k = 2; k < floor(d/2) + 1; ++k)
    {
        outfile1 << std::scientific << -realtype(6)*a(2*k, 2) << std::endl;
    }

    outfile1.close();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<hours>(stop - start);
    std::cout << duration.count() << "hours" << std::endl;

    return 0;

}
