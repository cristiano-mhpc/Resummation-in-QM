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
std::string st1;
int omegas = 41; // number of omegas to test

int main()
{

    vector<realtype> firstter(omegas), secondter(omegas), thirdter(omegas), fifthter(omegas), fourthter(omegas), rezult(omegas);

    //------------------read-in the terms--------------------------------------

    std::ifstream infile;
    infile.open("FIRST.txt");
    for (i=0; i < omegas; ++i)
    {
        infile >> st1;
        firstter(i) = realtype(st1);
    }

    infile.close();
    infile.clear();

    infile.open("SECOND.txt");
    for (i=0; i < omegas; ++i)
    {
        infile >> st1;
        secondter(i) = realtype(st1);
    }
    infile.close();
    infile.clear();

    infile.open("THIRD.txt");
    for (i=0; i < omegas; ++i)
    {
        infile >> st1;
        thirdter(i) = realtype(st1);
    }
    infile.close();
    infile.clear();

    infile.open("FOURTH.txt");
    for (i=0; i < omegas; ++i)
    {
        infile >> st1;
        fourthter(i) = realtype(st1);
    }
    infile.close();
    infile.clear();


    infile.open("FIFTH.txt");
    for (i=0; i < omegas; ++i)
    {
        infile >> st1;
        fifthter(i) = realtype(st1);
    }
    infile.close();
    infile.clear();

    std::ofstream outfile1;
    outfile1.open("disectres.txt",std::ios_base::out);
    outfile1.precision(digits);
    for (i=0; i < omegas; ++i)
    {
        rezult(i) = -(realtype(1)/realtype(2)) * (realtype(1) -  (firstter(i) + secondter(i) + fourthter(i)  + thirdter(i) + fifthter(i)) );
        outfile1 << std::scientific << rezult(i) << std::endl;
    }

    return 0;

}
