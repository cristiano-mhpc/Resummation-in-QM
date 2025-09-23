/*This code implements LU factorization of matrix using the BLOCK ALgorithm
in order to implement the factorization in parallel. Two processes are used and
as many logical threads as possible. To create and manage threads, I used the
Boost.Thread Library. To implement application level parallelization,
I used Boost.MPI.
1. A_00 = L_00 U_00
2. A_10 = L_10 U_00
3. A_01 = L_00 U_01
4. A_11 = L_10 U_01 + L_11 U_11

Steps 2 and 3 can be done in parallel by two processes. While the matrix
update in Step 4 can be done row by row by each of the threads.

*/

#include <boost/thread.hpp>
#include <string>
#include <boost/lexical_cast.hpp>

//Boost.Multiprecision mpfr_float
#include <boost/multiprecision/mpfr.hpp>

//for parallel implementation
#include <boost/mpi.hpp>

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

#include <algorithm>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <time.h>

#define digits 800

using namespace boost::archive;
using namespace boost::serialization;
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;

//serialization code for mpfr_float

using realtype = number<mpfr_float_backend<digits, allocate_stack>>;
#define MPFR_BUFFER_SIZE 810

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
const int l = 100, b = 9; // l*(b+1)  = d+1
const int d = l*(b+1) - 1;
const realtype Pi = boost::math::constants::pi<realtype>();
const realtype nu = realtype(1)/realtype(3);

matrix<realtype> Pr(b+1, b+1);
triangular_matrix<realtype, unit_lower> tul(b+1, b+1);
triangular_matrix<realtype, lower> tup(b+1, b+1);
matrix<realtype> P(d+1, d+1);
vector<realtype> moments(d+1);


void u_solver(boost::numeric::ublas::vector<realtype> mvs, int row , int colm) // solves for U_01 column per column
{
    mvs = matrix_vector_slice<matrix<realtype> >  (P, slice(row, 1, b+1), slice( colm,  0, b + 1));

    inplace_solve(tul, mvs, unit_lower_tag());

    matrix_vector_slice<matrix<realtype> > (P, slice(row, 1, b+1), slice(colm, 0, b+1)) = mvs; //updates P
}

void l_solver(boost::numeric::ublas::vector<realtype> mvs ,int roww, int colm) // solves for L_10 row per row
{
    mvs = matrix_vector_slice<matrix<realtype> > (P, slice(roww, 0, b+1), slice(colm, 1, b+1));

    inplace_solve (tup, mvs, lower_tag ());

    matrix_vector_slice<matrix<realtype> > (P, slice(roww, 0, b+1), slice(colm, 1, b+1)) = mvs; //updates P
}

void update (boost::numeric::ublas::vector<realtype> mvs1, boost::numeric::ublas::vector<realtype> mvs2, const matrix<realtype> *U01 ,int row, int start)
{
    /*We perform the update A_11 - L10 U01  = L_11 U_11 = A'_11, a chunk of rows at a time. Each row of A'_11 is updated by a thread. */

    /* start specify starting points for defining the various matrices used in this function. The starting points
     depend on the stage of the iteration. start = n*(b+1), start = 0 at first iteration: n = 0*/

    mvs1 = matrix_vector_slice<matrix<realtype> > ( P, slice(row, 0, b+1), slice(start, 1, b+1) );

    mvs2 = matrix_vector_slice<matrix<realtype> >  ( P, slice(row, 0, P.size2() - ( start + b+1 )), slice( start + b + 1 , 1, P.size2() - ( start + b+1 ) ) );

    matrix_vector_slice<matrix<realtype> > ( P, slice(row, 0, P.size2() - ( start + b+1 )), slice( start + b + 1 , 1, P.size2() - ( start + b+1 ) ) )

    =  mvs2 - prod(mvs1, *U01);
}

void up_tri() // for obtaining the upper triangular stored in Pr
{
    for(int i = 0; i < tup.size1(); ++i)
    {
        for(int j = 0; j <  i+1 ; ++j)
        {
            tup(i,j) = Pr(j,i);
        }
    }

}

void low_tri() // for obtaining the unit lower triangular stored in Pr
{
    for(int i = 0; i < tul.size1(); ++i)
    {
        for(int j = 0; j < i ; ++j)
        {
            tul(i,j) = Pr(i,j);
        }
    }
}

int main()
{
    using namespace std::chrono;
    mpi::environment env;
    mpi::communicator world;

    int num_process = 2; //specify the number of processes to be used

    int num_threads = boost::thread::hardware_concurrency(); /* gets the maximum number of physically available
                                                                logical threads. */

    auto start = high_resolution_clock::now();

    std::string st;
    std::ifstream infile;
    infile.open( "moments.txt", std::ios_base::in );
    for ( int i = 0; i < d+1; ++i )
    {
        std::getline(infile, st);
        moments(i) = realtype(st);
    }

    infile.close();

//---------------------------------read-in P(n,m)----------------------------------//
    std::string st2;
    std::ifstream infile2;
    infile2.open("../0_1k_P_block/matrix_p.txt", std::ios_base::in );
    for ( n = 0; n < d+1 ; ++n )
    {
        for ( m = 0; m < d+1; ++m )
        {
            std::getline( infile2, st2 );
            P(n,m) = realtype(st2);
        }
    }

    infile2.close();

//--------------------------------------LU Factorization-------------------------------------//

    if (world.rank() == 0) // l rounds of Block LU factorization to perform
    {

        for ( n = 0; n < l ; ++n )
        {
            permutation_matrix<std::size_t> PI( b+1 );//Records the row indices placed on top in the partial pivoting during LU factorization of Pr.

            Pr = matrix_range<matrix<realtype> > ( P, range( n*(b+1), ( n+1 )*( b+1 ) ), range( n*( b+1 ), ( n+1 )*( b + 1 )  ) );

            lu_factorize(Pr,PI); //unit lower and upper triangular factors of Pr are now stored in Pr
            world.send(1, 16, Pr); //send Pr to process 1

            /*row swapping was made in the submatrix Pr, we do the same for the
            corresponding elements of the  moments(i) and rows of P(i,j) adjacent with Pr. Row swaps made in Pr
            are encoded in PI */

            for (i = 0; i < PI.size(); ++i )
            {
                std::swap(moments(i + n*(b+1)), moments( PI(i) + n*(b+1) )); //swap corresponding rows of moments

                for (j = (n+1)*(b+1); j < P.size2(); ++j)
                {
                    std::swap(row(P, i  + n*(b+1))(j), row(P, PI(i) + n*(b+1))(j)); //swap row elements (not rows) of P to the right of Pr, j specifies the column
                }
                for (j = 0; j < n*(b+1); ++j)
                {
                    std::swap(row(P, i + n*(b+1))(j), row(P, PI(i) + n*(b+1))(j));//swap row elements (not rows) of P to the left of Pr
                }
            }

            //--------------------------------------------------------------------------------------------------//

            //update the relevant submatrix of P with the factorizations of Pr

            matrix_range<matrix<realtype> > (P, range(n*(b+1), (n+1)*(b+1)), range(n*(b+1), (n+1)*(b+1))) = Pr;

            //--------------------------------------------------------------------------------------------------//

            boost::thread t7(up_tri); // obtain the upper triangular factor stored in Pr
            t7.join();

            //matrix_range<matrix<realtype> > U_01( P, range( Start , Start + ( b + 1 )  ), range( Start + ( b+1 ) , P.size2() ) ); // if not allowed

            matrix<realtype> A_01( b+1, P.size2() - (n+1)*(b+1) );

            //matrix<realtype> U_01;

            A_01 = matrix_range<matrix<realtype> > ( P, range( n*(b+1), ( n+1 ) * ( b + 1 )), range( ( n+1) * ( b+1 ), P.size2() ) );

            world.send(1, 17, A_01); // send A_01 to process 1

            //-----------process 1 will proceed to compute U_01-----------------//

            //----------------Here we compute L_10 -----------------------------//

            int batch = (int) floor( ( P.size2() - (n+1)*(b+1))  / num_threads );

            int Start = n*(b+1); /*defines starting points for column updates in L_10
                                and row updates in L_10. */

            boost::numeric::ublas::vector<realtype> mhs(b+1);

            for (k = 0; k < batch; ++k)
            {
                boost::thread t[num_threads];
                for (i = 0 ; i < num_threads; ++i)
                {
                    int row = (n+1)*(b+1) + i + (k*num_threads);
                    t[i] = boost::thread{l_solver, mhs, row, Start};
                }

                for (int i = 0; i < num_threads; ++i)
                {
                    t[i].join();
                }

            }

            //last batch
            boost::thread t4[P.size2() - (n+1)*(b + 1) - (batch*num_threads)];
            for (i = 0; i < P.size2() - (n+1)*(b + 1) - (batch*num_threads); ++i)
            {
                int row = i + (n+1)*(b + 1) + (batch*num_threads);

                t4[i] =boost::thread{l_solver, mhs, row, Start};
            }

            for ( int i = 0; i < P.size2() - (n+1)*(b + 1) - (batch*num_threads); ++i )
            {
                t4[i].join();
            }

            //-----------------solve for A'_11 and replace the current trailing sub matrix-----------------------------//

            matrix<realtype> U_01;

            //U_01 = matrix_range<matrix<realtype> > ( P, range(n*(b+1), ( n+1 )*( b+1 )), range(( n+1) * ( b+1 ), P.size2() ) );

            world.recv(1, 18 , U_01);// receive U_01 from process 1

            //then assign to the appropriate matrix block
            matrix_range<matrix<realtype> > ( P, range( n*(b+1), ( n+1 )*( b+1 ) ), range( ( n+1 ) * ( b+1 ), P.size2() ) ) = U_01;

            boost::numeric::ublas::vector<realtype> mvh1( P.size2()-(n+1)*(b+1) );

            for (k = 0; k < batch; ++k) // divide the rows/columns into batches
            {
                boost::thread t[num_threads];

                for (i = 0 ; i < num_threads; ++i)
                {
                    int row = (n+1)*(b+1) + i + (k*num_threads);

                    t[i] = boost::thread{ update, mhs, mvh1, &U_01, row, Start };

                }

                for ( int i = 0; i < num_threads; ++i )
                {
                    t[i].join();
                }
            }

            // remaining rows to be taken care of by the last batch of threads
            boost::thread t5[P.size2() - (n+1)*(b + 1) - (batch*num_threads)];

            for ( i = 0; i < P.size2() - (n+1)*(b + 1) - (batch*num_threads); ++i )
            {
                int row = i + (n+1)*(b + 1) + (batch*num_threads);

                t5[i] =boost::thread{ update, mhs, mvh1, &U_01, row , Start };
            }

            for ( int i = 0; i < P.size2() - (n+1)*(b + 1) - (batch*num_threads); ++i )
            {
                t5[i].join();
            }

            PI.clear();
            U_01.clear();
            A_01.clear();
            Pr.clear();
            tup.clear();

        }

        //----------------------------------------------------------------------------------------------//
        //P has now been LU factorized. We can now solve the system of Linear equation.

        /*Constants is used as both input and output container. The system is
        Ax = (LU)x = b. Let Ux = x' so Lx' = b. Solve this for x'
        by forward substitution. Then Solve Ux = x' by backward substitution.
        */
        vector<realtype> Constants(d+1);
        Constants = moments;
        permutation_matrix<size_t> PI3(d+1);
        lu_substitute(P, PI3, Constants);

        //print out the constant
        std::ofstream outfile1;
        outfile1.open("Constant.txt",std::ios_base::out);
        outfile1.precision(digits);

        for (i = 0; i < Constants.size(); ++i)
        {
            outfile1 << std::scientific << Constants(i) * boost::math::factorial <realtype> ( i ) << std::endl;
        }

        outfile1.close();

        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<minutes>(stop - start);

        std::cout << duration.count() << "minutes" << std::endl;

        return 0;

    }

    if (world.rank() == 1)
    {
        for ( n = 0; n < l ; ++n )
        {
            //matrix<realtype> U_01 ( b+1, P.size2() - (n+1)*(b+1) );

            matrix<realtype> U_01;

            int batch = (int) floor( ( P.size2() - (n+1)*(b+1))  / num_threads);// the number of batches of threads to use

            int Start = n*(b+1); /*defines starting points for column updates in U_01
                                and row updates in L_10. */

            //matrix_range<matrix<realtype> > U_01( P, range( Start , Start + ( b + 1 )  ), range( Start + ( b+1 ) , P.size2() ) );

            world.recv(0, 16, Pr); // receive Pr from process 0
            world.recv( 0, 17, U_01 ); //receive A_01  from process 0 and store it in U_01
            boost::thread t6(low_tri); // obtain the unit lower triangular factor stored in Pr
            t6.join();

            matrix_range<matrix<realtype> > ( P, range( n*(b+1), ( n+1 )*( b+1 ) ), range( ( n+1) * ( b+1 ), P.size2()) ) = U_01;

            boost::numeric::ublas::vector<realtype> mhs(b+1);
            //Solve for U_01 block column by column. Each thread gets ones column.
            for (k = 0; k < batch; ++k)
            {
                boost::thread t[num_threads];
                for (i = 0 ; i < num_threads; ++i)
                {
                    int colm = (n+1)*(b+1) + i + (k*num_threads);

                    t[i] =boost::thread{u_solver, mhs, Start, colm }; // pass the m
                }

                for (int i = 0; i < num_threads; ++i)
                {
                    t[i].join();
                }
            }

            //last batch
            boost::thread t[ P.size2() - ( n+1 )*(b + 1) - (batch*num_threads) ];

            for (i = 0; i < P.size2() - (n+1)*(b + 1) - (batch*num_threads); ++i)
            {
                int colm = i + (n+1)*(b + 1) + (batch*num_threads);

                t[i] =boost::thread{u_solver, mhs, Start, colm};
            }

            for (int i = 0; i < P.size2() - (n+1)*(b + 1) - (batch*num_threads); ++i)
            {
                t[i].join();
            }

            U_01 = matrix_range<matrix<realtype> > ( P, range( n*(b+1), ( n+1 )*( b+1 ) ), range( ( n+1) * ( b+1 ), P.size2()) );

            //send back updated U_10 to process 0

            world.send( 0, 18, U_01 );

            U_01.clear();
            Pr.clear();
            tul.clear();

        }

    }

    return 0;

}
