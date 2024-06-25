#include <fstream>
#include <iostream>
#include "unsupported/Eigen/MPRealSupport"
#include "Eigen/Dense"
#include "mpreal.h"
#include <chrono>

using namespace std::chrono;
using namespace mpfr;
using namespace Eigen;
using namespace std;


int main()
{
    auto start = high_resolution_clock::now();
    using mpfr::mpreal;


    int i,k,m,n;
    const int d = 1249; //number of moments used
    const int digits = 1200; //for setting the precision


    /*Setup default precision for all subsequent computations
    MPFR accepts precision in bits - so we do the conversion*/
    const int bits = mpfr::digits2bits(digits);
    mpreal::set_default_prec(bits);

    //some mp literals
    const mpreal pi = mpfr::const_pi();
    const mpreal nu = mpfr::mpreal(2)/mpfr::mpreal(3);

    //Declare matrix and vector types with multiprecision scalar types
    typedef Matrix<mpreal,Dynamic,Dynamic> MatrixXmp;
    typedef Matrix<mpreal,Dynamic,1>      VectorXmp;
    typedef Matrix<mpreal,1,Dynamic>   RowVectorXmp;


    //Allocate matrices

    MatrixXmp P(d+1,d+1);

    std::string st1;
    std::ifstream infile2;
    infile2.open("/scratch3/chris.tica/FUNNEL/construct_P_2kdig/0_1_25k_P_block/matrix_p.txt", std::ios_base::in );
    for ( n = 0; n < d+1; n++ )
    {
        for ( m = 0; m < d+1; m++ )
        {
            std::getline( infile2, st1 );
            P(n,m) = mpreal(st1);
        }
    }

    infile2.close();

    // allocate and read-in the moments
    VectorXmp moments(d+1);

    std::string st;
    std::ifstream infile;
    infile.open( "moments.txt", std::ios_base::in );
    for ( int i = 0; i < d+1; ++i )
    {
        std::getline(infile, st);
        moments(i) = -mpreal(st);
    }

    infile.close();
    //solve the system P constant = moments
    VectorXmp constant = P.partialPivLu().solve(moments);

    //print out the constant
    ofstream outfile1;
    outfile1.open("Constant.txt",ios_base::out);
    outfile1.precision(digits);
    for (i=0 ;i < d+1;i++)
    {
        outfile1 << scientific << constant(i)*fac_ui(i)<< endl;
    }
    outfile1.close();

    P.resize(0,0);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<hours>(stop - start);
    std::cout << duration.count() << "hours" << std::endl;

    return 0;
}
