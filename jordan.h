//Gvozdev 17.02.2019
#include "mpi.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring> //for memcpy


#define ERR_CANNOT_OPEN      -1
#define ERR_READ             -2
#define ERR_LITTLE_MATRIX    -3
#define BLOCK_DEGENERATE     -4
#define MATRIX_DEGENERATE    -5

#define EPS 1e-16
#define PRINT_SIZE 15

#define LOG(...) std::cerr<< #__VA_ARGS__<< " = " <<(__VA_ARGS__)<< "\n"




void
data_init ( double *&a , double *&x , double *&b , double *&save,
            int n , int k , int m , int l , int p , int my_rank,
            int local_lines);
int
read_matrix (double *a, int n, int m, int k, int l, int p, int my_rank,
             char * name, int local_lines, double *save);
void
init_matrix( double *a , int n , int m , int k , int l , int p , int my_rank,
             int local_lines);

void
print_matrix ( double *a , int n , int m , int k , int l , int p , int my_rank ,
                   int size , double * save, int local_lines);
double
get_full_time();
double
get_time();
void
matrix_mult_vector (double *a, double *b, double *x, int n, int k, int m, int l, int p,
                    int my_rank, double *save, int local_lines);



int solve( double *a , double *b , double *x , double *save , int n , int k , int m ,
              int l , int p , int my_rank, int local_lines);
void make_0( int hight , int width , int real_size_c , double * c );
void make_E( int real_size_c , double * c );
double norm_matrix ( double * a , int n , int k, int m , int l,
                     int p , int my_rank, int local_lines);
void get_block(int i, int j, int width, int hight, int real_size_c,
               int n, double *c, double *a, int k, int m,
               int local_lines, int my_rank, int kol_proc);
void put_block(int i, int j, int width, int hight, int real_size_c ,
               int n, double *c, double *a, int k, int m,
               int local_lines, int my_rank, int kol_proc);
int num_min_inverse_norm_in_column( double * a , int num_column, double * c0 ,
                                    double * c1 , int n , int m , int k , int p,
                                    double norm , int * index_in_block,
                                    int * index, int local_lines, int my_rank);
void multiply_line_on_inverse( int i , int column , double * a , double * b , double * c0 ,
                               double * c1 , double * c2 , int n , int k , int m , int l ,
                               double norm , int * index_in_block , int my_rank, int p,
                               int local_lines);
void substract_lines( int i , int column , double * a , double * b , double * c0 ,
                      double * c1 , double * c2 , int n , int m , int k , int l,
                      double * save, int local_lines, int my_rank, int p);
void substract_block( double * c0 , double * c1 , double * c2 , int hight ,
                      int width , int real_size_c);
void block_mult_block(double *c, double *a, double *b, int n);
int search_inverse_matrix ( double * c0 , double * c1 , int m ,
                            double norm , int * index_in_block );
void recovery_solution( double * b , double * x , int k , int m , int l ,
                        int p , int my_rank, int local_lines, int * index);
double search_error( double * x , double * b , int k, int m,
                     int l, int p, int local_lines, int my_rank);
double search_residual( double * a , double * x , double * save , int n , int m, int k, int l,
                        int p, int my_rank, int local_lines);
template< typename T >
void print_vector( T *b, int n);
void print_block( double *a, int size);
void print_br( double *a, int hight, int width);

using namespace std;


class Double_int
{
public:
    double elem = 0;
    int num = 0;
    Double_int()
    {
        elem=0; num=0;
    }
};
