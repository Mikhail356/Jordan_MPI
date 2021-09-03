//Gvozdev 29.02.20
#include "jordan.h"
int
solve (double *a, double *b, double */*x*/, double *save, int n, int k, int m,
       int l, int p, int my_rank, int local_lines)
{
    int *index_in_block;
    double *c0 , *c1 , *c2 ;
    double norm, y = 0;
    int * index, error = 0, fatal_error = 0;
//    int local_lines = 0;

//    local_lines = (k-my_rank)/p;
    index_in_block = new int [m];
    c0 = new double [m*m];
    c1 = new double [m*m];
    c2 = new double [m*m];
    index = new int [k];
    memset(index, -1, k*sizeof (int));
//LOG(k%p);
//LOG(local_lines);

    //print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
    //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);

    norm = EPS*norm_matrix (a, n, k, m, l, p, my_rank, local_lines);
    //LOG(norm);

    for( int i = 0 ; i < m ; i ++ ) { index_in_block[i] = i ; }

    make_0(0,0,m,c2);
    make_0(0,0,m,c1);
    make_0(0,0,m,c0);
//LOG("1");
    for( int i = 0 ; i < k ; i ++ ) // n == k * m + l // loop through columns
    {
        //print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
        //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);
        error = num_min_inverse_norm_in_column(a, i, c0, c1, n, m, k, p, norm,
                                               index_in_block, index,
                                               local_lines, my_rank);
        //LOG("1.1");
        if( error == MATRIX_DEGENERATE )
        {
            delete [] index_in_block ;
            delete [] c0 ;
            delete [] c1 ;
            delete [] c2 ;
            delete [] index ;
            return MATRIX_DEGENERATE ;
        }
        else {
            index[i] = error;
        }
        //LOG(index[i]);
        multiply_line_on_inverse (index[i], i, a, b, c0, c1, c2, n, k, m, l,
                                  norm, index_in_block, my_rank, p,
                                  local_lines);
        //LOG("1.2");

        //print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
        //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);

        substract_lines (index[i], i, a, b, c0, c1, c2, n, m, k, l, save,
                         local_lines, my_rank, p);
        //LOG("1.3");
        //print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
        //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);
    }
//LOG("2");
    //print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
    //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);

    error = 0;
    local_lines*=m;
    //LOG(norm);
    for( int i = my_rank/*0*/ ; i < k/*local_lines*/ ; i+=p /*++*/ )
    {
        for(int t = 0; t < m; t ++)
        {
            y = 0;
            for( int j = (index[i])*m/*0*/ ; j < n ; j ++ )/*need start from
                                                                index[t+i]*m*/
            {
                y += fabs( a[((i/p)*m+t)*n+j] );
//                LOG(a[((i/p)*m+t)*n+j]);
//                LOG(fabs( a[((i/p)*m+t)*n+j] ));
            }

            if(y <= norm)
            {
                //LOG(y);
                error = MATRIX_DEGENERATE;
            }
        }
    }
    if(my_rank==k%p/*0*/)
    {
        for(int i = 0; i < l; i++)
        {
            y=0;
//            LOG(y);
            for(int j=k*m ; j<n; j++)
            {
//                LOG(local_lines);
//                LOG((local_lines+i)*n+j);
//                LOG(a[(local_lines+i)*n+j]);
//                LOG(fabs( a[(local_lines+i)*n+j] ));
                y += fabs( a[(local_lines+i)*n+j] );
//                LOG(y);
//                LOG(norm);
            }
            if(y <= norm)
            {
                error = MATRIX_DEGENERATE;
//                fatal_error = MATRIX_DEGENERATE;
//                LOG(error);
//                LOG(fatal_error);
            }
        }
    }
    MPI_Allreduce(&error, &fatal_error, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
//    LOG(error);
//    LOG(fatal_error);
//LOG("3");
    if(fatal_error == MATRIX_DEGENERATE)
    {
//        LOG("yo");
//        delete [] save ;
        delete [] index_in_block ;
        delete [] c0 ;
        delete [] c1 ;
        delete [] c2 ;
        delete [] index ;
        return MATRIX_DEGENERATE;
    }
    local_lines/=m;

    /* start last step of algoritm
     * for down right (block l*l) side matrix a and vector b
     * search inverse matrix and multiply on it block l*l(for what?)*/
    if( l != 0 )
    {
        //LOG("4");

        //print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
        //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);
        if(my_rank == k%p/*0*/)
        {
            get_block(local_lines*n*m , k*m , l , l , m , n , c0 , a,
                       k, m, local_lines, my_rank, 0);

            int schetchik = 0;
            for ( int i = 0 ; i < m ; i++ )
            {
                for( int j =0 ; j < m ; j++)
                {
                    if( i < l && j < l)
                    {
                        c0[ schetchik ]=c0[ i*m+j ];
                        schetchik ++ ;
                    }
                }
            }
//LOG("5");
//print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
//print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);
            error = search_inverse_matrix(c0, c1, l, norm, index_in_block);
            if( error == BLOCK_DEGENERATE) error = MATRIX_DEGENERATE;
            else
            {
                schetchik=0;
                for ( int i = 0 ; i < m ; i++ )
                {
                    for( int j =0 ; j < m ; j++)
                    {
                        c1[ i*m+j ]=c0[ schetchik ];
                        if( i < l && j < l)
                        {
                            schetchik ++ ;
                        }
                    }
                }
                make_0(l, l, m, c1);
                get_block(local_lines*m, 0, 1, l, m, 1, c0, b, k, m,
                          local_lines, my_rank, 0);
                make_0(l, 1, m, c0);

                block_mult_block(c2, c1, c0, m);//c2=c1*c0
                put_block(local_lines*m, 0, 1, l, m, 1, c2, b, k, m,
                          local_lines, my_rank, 0);

                make_E( m , c2);
                put_block(local_lines*n*m , k*m , l , l , m , n , c2 , a,
                          k, m, local_lines, my_rank, 0);//t k A^(-1)*A==E
            }
//LOG("6");
//print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
//print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);
        }

        MPI_Allreduce(&error, &fatal_error, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
        if(fatal_error==MATRIX_DEGENERATE)
        {
            delete [] index_in_block ;
            delete [] c0 ;
            delete [] c1 ;
            delete [] c2 ;
            delete [] index ;
            return MATRIX_DEGENERATE ;
        }

//LOG("7");
//print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
//print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);
//print_block(c0,m);
//print_block(c1,m);
//print_block(c2,m);
//LOG("go");
//LOG(k);
    }

    if(my_rank == k%p/*0*/)
    {
        //LOG("hello");
        memcpy(save+n*m, b+local_lines*m, l*sizeof(double));
        //print_vector(save, l);
        //LOG(save[n*m]);
    }
    //print_vector(save, l);
    MPI_Bcast(save+n*m, l, MPI_DOUBLE, k%p/*0*/, MPI_COMM_WORLD);
//LOG("7.1");
    /*  continue last step of algoritm
     *  substract block l*l  */
    if(l!=0)
    {
        for( int i = 0 ; i < local_lines ; i ++ )
        {
            get_block(i*n*m, k*m, l, m, m, n, c0, a, k, m, local_lines,
                      my_rank, 0);
//        print_block(c0,m);
//            get_block( k*m  ,  -1  , 1 , l , m , 1 , c1 , b,
//                       k, m, local_lines, my_rank, p);
            get_block(0, n*m, 1, l, m, 1, c1, save, k, m, local_lines,
                      my_rank, 0);
//        print_block(c1,m);
            make_0(m, l, m, c0);
            make_0(l, 1, m, c1);
            //LOG("c0");
//        print_block(c0,m);
        //LOG("c1");
//        print_block(c1,m);
            block_mult_block(c2, c0, c1, m);
//print_block(c2,m);
            get_block(i*m, 0, 1, m, m, 1, c1, b, k, m, local_lines, my_rank,
                      0);
//print_block(c1,m);
            substract_block(c0 ,c1 ,c2 ,m ,1 ,m);//c0=c1-c2
//print_block(c0,m);
            put_block(i*m, 0, 1, m, m, 1, c0, b, k, m, local_lines, my_rank,
                      0);
        }

    }
//LOG("8");

    //print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
    //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);

//    if(my_rank == 0)print_vector( index, k);
    recovery_solution( b , save , k , m , l , p ,my_rank ,local_lines ,index);
    //LOG("helol");
    //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);
//LOG("9");
    delete [] index_in_block ;
    delete [] c0 ;
    delete [] c1 ;
    delete [] c2 ;
    delete [] index ;

    return 0;
}

void get_block(int i, int j, int width, int hight, int real_size_c,
               int n, double *c, double *a, int /*k*/, int /*m*/,
               int /*local_lines*/, int /*my_rank*/, int /*p*/)
/*i,j koordinaty verhnego levogo ugla bloka c v matrice a
width, hight shirina i vysota bloka c; n size of matrix a */
{
    int count=i+j;
    for(int t = 0 ; t < hight ; t ++ )
    {
        for( int r = 0 ; r < width ; r ++ )
        {
            c[ t * real_size_c + r ] = a[ count + r ];
        }
        count+=n;
    }
}

void put_block(int i, int j, int width, int hight, int real_size_c ,
               int n, double *c, double *a, int /*k*/, int /*m*/,
               int /*local_lines*/, int /*my_rank*/, int /*p*/)
{
    int count=0;
    count=i+j;
    for( int t = 0 ; t < hight ; t ++ )
    {
        for( int r = 0 ; r < width ; r ++ )
        {
            a[ count + r ] = c[ t * real_size_c + r ];
        }
        count+=n;
    }
}

void make_0( int hight , int width , int real_size_c , double * c )
{
    for( int i = 0 ; i < real_size_c ; i ++ )
    {
        for( int j = 0 ; j < real_size_c ; j ++ )
        {
            if( i >= hight || j >= width )
            {
                c[ i * real_size_c + j ] = 0 ;
            }
        }
    }
}

void make_E( int real_size_c , double * c )
{
    for( int i = 0 ; i < real_size_c ; i ++ )
    {
        for( int j = 0 ; j < real_size_c ; j ++ )
        {
            c[ i * real_size_c + j ] = ( i != j ? 0 : 1 ) ;
        }
    }
}

double norm_matrix ( double * a , int n , int k, int m , int l,
                     int p , int my_rank, int local_lines)
// a matrix   n * n , norm by line
{
    double norm_of_matrix = 0;
    double max1 = 0 , max2 = 0 ;
    int count = 0;

    //LOG("I in norm matrix");
    for( int i = 0 ; i < local_lines ; i ++ )
    {
        count = i * m ;
        for( int temp = 0 ; temp < m ; temp ++ )
        {
            max1 = 0;
            for( int j = 0 ; j < n ; j ++ )
            {
                max1 += fabs(a[ (count + temp)*n + j ]);
            }

            if ( max1 > max2 ) {max2 = max1 ;}
        }
    }
    if(my_rank == k%p && l != 0)
    {
        //LOG(my_rank);
        //LOG(local_lines);
        count = local_lines * m ;
        for( int temp = 0 ; temp < l ; temp ++ )
        {
            max1 = 0;
            for( int j = 0 ; j < n ; j ++ )
            {
                max1 += fabs(a[ (count + temp)*n + j ]);
                //LOG(fabs(a[ (count + temp)*n + j ]));
                //LOG( (count + temp)*n + j );
            }

            if ( max1 > max2 ) {max2 = max1 ;}
        }
    }

    MPI_Allreduce(&max2, &norm_of_matrix, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    return norm_of_matrix;
}
double norm_block( double * a , int n )
{
    double max1 , max2 = 0 ;
    int schet;

    for( int i = 0 ; i < n ; i ++ )
    {
        max1 = 0;
        schet = i * n ;

        for( int j = 0 ; j < n ; j ++ )
        {
            max1 = max1 + fabs(a[ schet + j ]);
        }

        if ( max1 > max2 ) {max2 = max1 ;}
    }
    return max2;
}

int search_inverse_matrix ( double * c0 , double * c1 , int m ,
                            double norm , int * index_in_block )
//write inverse matrix for c0 in c0 . Jordan
{
    double y,y0;
    int www=0 , scheti , schetj ;
    for( int i = 0 ; i < m ; i ++ )
    {
        index_in_block[i]=i;
        for ( int  j = 0 ; j < m ; j ++ )
        {
            c1[ i * m + j ] = ( i != j ? 0 : 1 );
        }
    }

    for ( int i = 0 ; i < m ; i ++ )
    {
        y=0;
        for( int j = i ; j < m ; j ++ )
        {
            y0 = fabs(  c0[ index_in_block[j] * m + i ] );
            if ( y0 > y )
            {
                y = y0;
                www = j;
            }
        }
        if( fabs(y) <= norm ) { return BLOCK_DEGENERATE ; }
        swap(index_in_block [ i ] , index_in_block [ www ]) ;

        scheti = index_in_block [ i ] * m;
        y = 1./c0 [ scheti + i];

        for ( int  j = 0 ; j < m ; j ++ )
        {
            c0 [ scheti + j] *= y ;
            c1 [ scheti + j] *= y ;
        }

        for( int j = 0 ; j < m ; j ++ )
            //vichest is stroki j stroku www umnojuyu na koef j stroki
        {
            schetj = index_in_block[j] * m;
            y=c0[ schetj + i ];
            www=0;

            for ( int t = 0 ; t < m ; t ++ )
            {
                if( scheti != schetj )
                {
                    c0[ schetj + t ] = c0[ schetj + t ]
                            - ( y * c0[ scheti + t ] );
                    c1[ schetj + t ] = c1[ schetj + t ]
                            - ( y * c1[ scheti + t ] );
                }
                if( fabs(c0[ schetj + t ]) < norm )
                {
                    www ++ ;
                }
            }

            if( www == m )
            {
                return BLOCK_DEGENERATE;
            }
        }
    }
    for ( int i = 0 ; i < m ; i ++ )
    {

        for( int j = 0 ; j < m ; j ++ )
        {
            c0 [ i * m + j ] = c1 [ index_in_block[ i ] * m + j ] ;
        }
    }
    return 0;
}

int num_min_inverse_norm_in_column( double * a , int num_column, double * c0 ,
                                    double * c1 , int n , int m , int k ,
                                    int p, double norm , int * index_in_block,
                                    int * index, int local_lines, int my_rank)
{
    int h = 0;
    Double_int pair, res_pair;
    double f1/*,f2 = 0*/;
    int count = 0, ind = 0, error = 0;

    //LOG("in num inverse");
    for(int i = my_rank; i < k; i +=p)
    {
//        if(index[i]==-1)
//        {
            get_block( (i/p)*n*m , num_column*m , m , m , m , n , c0 , a,
                       k, m, local_lines, my_rank, p);
            //print_block(c0,m);
            h = search_inverse_matrix( c0 , c1 , m , norm , index_in_block );

            //print_block(c0,m);
            if( h == BLOCK_DEGENERATE )
            {
                count++;
                //LOG(h); LOG(i);
                //LOG(count);
                h=0;
            }
            else
            {
                f1 = norm_block( c0 , m );
                //LOG(f1);    LOG(i); LOG(pair.elem);
                //print_vector(index, k);
                if(ind == 0)
                {
                    ind = 1;
                    pair.elem = f1;
                    pair.num = i;
                }
                else if(pair.elem > f1 && index [i] == -1)
                {
                    pair.elem = f1;
                    pair.num = i;
                }
            }
//        }
//        else
//        {
//            count++;
//        }
    }
//LOG(local_lines);
    if(count == local_lines)
    {
        pair.elem = 0;
        pair.num = 0;
    }
    else
    {
        pair.elem = 1./pair.elem;
//        pair.num = my_rank + pair.num*p;
    }
//LOG(pair.elem);
//LOG(pair.num);
    MPI_Allreduce( &count, &error, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //LOG(error);
    //LOG(k);

    //LOG("uhozhu iz num_inverse");
    if ( error == k )
    {
        return MATRIX_DEGENERATE ;
    }

    MPI_Allreduce( &pair, &res_pair, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
                   MPI_COMM_WORLD);

    return res_pair.num;
}

void multiply_line_on_inverse (int i, int column, double *a, double *b,
                               double *c0, double *c1, double *c2, int n,
                               int k, int m, int l, double norm,
                               int *index_in_block, int my_rank, int p,
                               int local_lines) // n == k * m + l
{
    //LOG("multiply line");
    if(my_rank == i%p)
    {
        get_block( (i/p)*n*m , column*m , m , m , m , n , c0 , a ,
                   k, m, local_lines, my_rank, 0);
//        print_block(c0, m);
        search_inverse_matrix( c0 , c1 , m , norm , index_in_block );
        //upper inverse matrix in c0
//        print_block(c0, m);

        for( int j = column+1/*my_rank*/ ; j < k ; j ++ )
        {
//            if(j>column)
//            {
                get_block( (i/p) * n * m , j * m , m , m , m , n , c2 , a ,
                           k, m, local_lines, my_rank, 0);
                block_mult_block( c1 , c0 , c2 , m ); // c1 = c0 * c2
                put_block( (i/p) * n * m, j * m , m , m , m , n , c1 , a ,
                           k, m, local_lines, my_rank, 0);
//            }
        }

//        if(my_rank==(i%p))
//        {
            if(l!=0)
            {
                get_block( (i/p) * n * m, k * m , l , m , m , n , c2 , a ,
                           k, m, local_lines, my_rank, 0);
//                print_block(c0, m);
//                print_block(c2, m);
                block_mult_block( c1 , c0 , c2 , m );
//                print_block(c1, m);
                put_block( (i/p) * n * m, k * m , l , m , m , n , c1 , a ,
                           k, m, local_lines, my_rank, 0);
            }
            get_block((i/p)*m,0,1,m,m,1, c2, b,
                      k, m, local_lines, my_rank, 0);
//            print_block(c0, m);
//            print_block(c2, m);
            block_mult_block( c1 , c0 , c2 , m );
//            print_block(c1, m);
            put_block((i/p)*m,0,1,m,m,1, c1, b,
                      k, m, local_lines, my_rank, 0);
//        }
    }
}

void
substract_lines (int i, int column, double *a, double *b, double *c0,
                 double *c1, double *c2, int n, int m, int k, int l,
                 double *save, int local_lines, int my_rank, int p)
//n == k*m + l
// a[t][j] = a[t][j] - a[t][i]*a[i][j]
{
    //LOG("substract_lines()");
    //LOG("before MPI_Bcast(...)");
    //print_br(a, (my_rank==k%p?local_lines*m+l:local_lines*m) ,n);
    //print_br(save, m, n);
    //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);

    //LOG(i); LOG(p);
    if(i%p==my_rank)
    {
        //LOG("befor memcpy(...)"); LOG(my_rank);
        memcpy(save, a+((i/*-my_rank*/)/p)*n*m, n*m*sizeof(double));
        memcpy(save+n*m, b+((i/*-my_rank*/)/p)*m, m*sizeof(double));
        //print_br(a, (my_rank==k%p?local_lines*m+l:local_lines*m) ,n);
        //print_br(save, m, n);
        //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);
        //LOG("after memcpy(...)");
    }

    /*int error = -100;
    error = */MPI_Bcast( save, ((n*m)+m), MPI_DOUBLE, i%p, MPI_COMM_WORLD);
//    MPI_Bcast( save, ((n*m)+m), MPI_DOUBLE, i%p, MPI_COMM_WORLD);
    //LOG(error);
    //LOG("After MPI_Bcast(...)");
//    LOG(n); LOG(m);
    //print_br(a, (my_rank==k%p?local_lines*m+l:local_lines*m) ,n);
    //print_br(save, m, n);
    //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);

//LOG(save[n*m]);
//LOG(save[0]);
    for( int t = my_rank ; t < k ; t += p )
    {
        if(t!=i)
        {
            //LOG("for a");
            for( int j = column+1 ; j < k ; j ++)
            {
                get_block( (t/p)*n*m , column*m , m , m , m , n , c2 , a ,
                           k, m, local_lines, my_rank, 0);
                //print_block(c2,m);

                //for a
                get_block( 0 , j*m , m , m , m , n/*size*/ , c1 , save ,
                           k, m, local_lines, my_rank, 0);
                //print_block(c1,m);

                block_mult_block( c0 , c2 , c1 , m );//c0 = c2 * c1
                //print_block(c0,m);

                get_block( (t/p)*n*m , j*m , m , m , m , n , c1 , a ,
                           k, m, local_lines, my_rank, 0);
                //print_block(c1,m);
                substract_block( c2 , c1 , c0 , m , m , m );//c2=c1-c0

                put_block( (t/p)*n*m , j*m , m , m , m , n , c2 , a ,
                           k, m, local_lines, my_rank, 0);
                //print_block(c2,m);
            }

            //LOG("for b------------------------------");
            //for b
            get_block( (t/p)*n*m , column*m , m , m , m , n , c2 , a ,
                       k, m, local_lines, my_rank, 0);
            //print_block(c2, m);

            get_block( 0 , n*m , 1 , m , m , 1 , c1 , save ,
                       k, m, local_lines, my_rank, 0);
            //print_block(c1, m);

            block_mult_block( c0 , c2 , c1 , m );
            //print_block(c0, m);
            get_block((t/p)*m,0,1,m,m,1, c1, b,
                      k, m, local_lines, my_rank, 0);
            //print_block(c1,m);
            substract_block( c2 , c1 , c0 , m , 1 , m);
            //print_block(c2, m);
            put_block((t/p)*m,0,1,m,m,1, c2, b,
                      k, m, local_lines, my_rank, 0);
            //LOG("end for b------------------------------");
            if( l != 0 && (my_rank == (t%p)) )/* for right side a and block l*l
                                                 except t-line */
            {
                get_block( (t/p)*n*m , column*m , m , m , m , n , c2 , a ,
                           k, m, local_lines, my_rank, 0);//m x m

                //for a

                get_block( 0 , k*m , l , m , m , n/*size*/ , c1 , save ,
                           k, m, local_lines, my_rank, 0);//m x l
                //print_block(c1,m);

                block_mult_block( c0 , c2 , c1 , m );//c0 = c2 * c1
                get_block( (t/p)*n*m , k*m , l , m , m , n , c1 , a ,
                           k, m, local_lines, my_rank, 0);//m x l
                substract_block( c2 , c1 , c0 , m , l , m);//c2=c1-c0
                put_block( (t/p)*n*m , k*m , l , m , m , n , c2 , a ,
                           k, m, local_lines, my_rank, 0);

                //for b already make
            }
            //LOG("end for a");
        }
    }
//LOG("substract i-lines...-------------------------------------------------");
    //print_br(a, (my_rank==k%p?local_lines*m+l:local_lines*m) ,n);
    //print_br(save, m, n);
    //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);

    if(l != 0 && i != k && my_rank == k%p)/*substract i-lines from down-right
                                            block l*l matrix *a and block l*1
                                            vector *b*/
    {
        get_block( local_lines*n*m , column*m , m , l , m , n , c2 , a ,
                   k, m, local_lines, my_rank, 0);//l x m
//print_block(c2,m);
        //for a
//LOG("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
        get_block( 0 , k*m , l , m , m , n/*size*/ , c1 , save ,
                   k, m, local_lines, my_rank, 0);//m x l
//LOG("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
//        print_block(c1,m);
        block_mult_block( c0 , c2 , c1 , m );//c0 = c2 * c1
//print_block(c0,m);
        get_block( local_lines*n*m , k*m , l , l , m , n , c1 , a ,
                   k, m, local_lines, my_rank, 0);//m x l
//print_block(c1,m);
        substract_block( c2 , c1 , c0 , l , l , m);//c2=c1-c0
//print_block(c2,m);
        put_block( local_lines*n*m , k*m , l , l , m , n , c2 , a ,
                   k, m, local_lines, my_rank, 0);
//LOG("for b2... -----------------------------------------------------------");
        //for b
        get_block( local_lines*n*m , column*m , m , l , m , n , c2 , a ,
                   k, m, local_lines, my_rank, 0);
//print_block(c2,m);
        get_block( 0 , n*m , 1 , m , m , 1/*size*/ , c1 , save ,
                   k, m, local_lines, my_rank, 0);
//print_block(c1,m);
        block_mult_block( c0 , c2 , c1 , m );
//print_block(c0,m);
        get_block( local_lines*m , 0 , 1 , l , m , 1 , c1 , b ,
                   k, m, local_lines, my_rank, 0);
//print_block(c1,m);
        substract_block( c2 , c1 , c0 , l , 1 , m);
//print_block(c2,m);
        put_block( local_lines*m , 0 , 1 , l , m , 1 , c2 , b ,
                   k, m, local_lines, my_rank, 0);
        //LOG(" end for b2... ----------------------------------------------");
    }
    //LOG("end substract i-lines...-----------------------------------------");
    //print_br(a, (my_rank==k%p?local_lines*m+l:local_lines*m) ,n);
    //print_br(save, m, n);
    //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);

//LOG("for bottom...--------------------------------------------------------");
    if(l != 0 && my_rank == k%p) /* for bottom line of block except columns
                                    with num < column+1 and also except last
                                    column with blocks size of m*l and l*l */
    {
        for(int j = column+1 ; j < k ; j ++ )
        {
            get_block (local_lines*n*m, column*m, m, l, m, n, c2, a, k, m,
                       local_lines, my_rank, 0);// l x m
//print_block(c2,m);
            //for a
            get_block (0 ,j*m ,m ,m ,m ,n ,c1 ,save , k, m, local_lines,
                       my_rank, p);
//print_block(c1,m);
            block_mult_block (c0 ,c2 ,c1 ,m);//c0 = c2 * c1
//print_block(c0,m);
            get_block (local_lines*n*m, j*m, m, l, m, n, c1, a, k, m,
                       local_lines, my_rank, 0);
//            print_block(c1,m);
            substract_block (c2 ,c1, c0, l, m, m);//c2=c1-c0 l x m
//print_block(c2,m);
            put_block (local_lines*n*m, j*m, m, l, m, n, c2, a, k, m,
                       local_lines, my_rank, 0);
        }
    }
    //LOG("end for bottom...----------------------------------------------");
    //print_br(a, (my_rank==k%p?local_lines*m+l:local_lines*m) ,n);
    //print_br(save, m, n);
    //print_vector(b, my_rank==k%p?local_lines*m+l:local_lines*m);
}

void substract_block( double * c0 , double * c1 , double * c2 , int hight ,
                      int width , int real_size_c)//c0=c1-c2 hight*width
{
    for(int i = 0 ; i < hight ; i ++ )
    {
        for( int j = 0 ; j < width ; j ++ )
        {
            c0[ i*real_size_c+j ] = c1[ i*real_size_c+j ]
                    - c2[ i*real_size_c+j ];
        }
    }
}

void block_mult_block(double *c, double *a, double *b, int n)
// n%3==0  c=a*b
{
        int i,j,k;
        double *pc = c , *pa = a , *pb = b , sum[9] ;

        for(i=0;i<n;i++)
    {
                for(j=0;j<n;j++)
        {
            c[i*n+j]=0.;
        }
    }

        for(i=0;i<n;i+=3)
    {
                for(j=0;j<n;j+=3)
        {
                        sum[0]=sum[1]=sum[2]=sum[3]=sum[4]=sum[5]=sum[6]
                                =sum[7]=sum[8]=0.;

                        for(k=0;k<n;k++)
                        {
                                pa=a+i*n+k;
                                pb=b+k*n+j;

                                sum[0]+=pa[0]*pb[0];
                                sum[1]+=pa[0]*pb[1];
                                sum[2]+=pa[0]*pb[2];
                                sum[3]+=pa[n]*pb[0];
                                sum[4]+=pa[n]*pb[1];
                                sum[5]+=pa[n]*pb[2];
                                sum[6]+=pa[2*n]*pb[0];
                                sum[7]+=pa[2*n]*pb[1];
                                sum[8]+=pa[2*n]*pb[2];
                        }

                        pc=c+i*n+j;

                        pc[0]			+=sum[0];
                        pc[1]			+=sum[1];
                        pc[2]			+=sum[2];
                        pc[n]			+=sum[3];
                        pc[n+1]			+=sum[4];
                        pc[n+2]			+=sum[5];
                        pc[2*n]			+=sum[6];
                        pc[2*n+1]		+=sum[7];
                        pc[2*n+2]		+=sum[8];
                }
    }
}

void
recovery_solution (double *b, double *save, int k, int m, int l, int p,
                   int my_rank, int local_lines, int * index)
{
    //LOG(m);
    if(my_rank==k%p/*0*/)
    {
        MPI_Status status;
        for( int i = 0 ; i < k ; i ++ )
        {
            if(index[i]%p!=my_rank/*0*/)
                MPI_Recv(save+i*m, m, MPI_DOUBLE, index[i]%p, index[i]/*/p*/,
                         MPI_COMM_WORLD, &status);
            else
            {
                memcpy(save+i*m, b+((index[i]/p)*m), m*sizeof(double));
            }
        }
        if( l!= 0 )
        {
            for( int i = 0 ; i < l ; i ++ )
            {
                save[k*m+i]=b[local_lines*m+i];
            }
        }
    }
    else
    {
        for(int i = my_rank; i < k; i += p)
        {
            MPI_Send(b+(i/p)*m, m, MPI_DOUBLE, k%p/*0*/, i, MPI_COMM_WORLD);
        }
    }
    MPI_Bcast(save, k*m+l, MPI_DOUBLE, k%p/*0*/, MPI_COMM_WORLD);
}

double
search_error (double *x, double *save, int k, int m, int l, int p,
              int local_lines, int my_rank)
{
    double y1, y2=0;

    for(int i=0; i<local_lines; i++)
    {
        for(int j=0; j<m; j++)
        {
            y1 = fabs(save[(my_rank+i*p)*m+j]-x[(my_rank+i*p)*m+j]);
            if(y1 > y2)
                y2 = y1;
        }
    }
    if(my_rank == p-1)
    {
        for(int i = 0 ; i < l; i ++)
        {
            y1 = fabs(save[k*m+i] - x[k*m+i]);
            if(y1 > y2)
                y2 = y1;
        }
    }
    MPI_Allreduce(&y2, &y1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return y1;
}

double
search_residual (double *a, double *x, double *save, int n, int m, int k,
                 int l, int p, int my_rank, int local_lines)
{
    double y1=0,y2=0,save1=0,save2=0;
    int schet = 0 ;
    MPI_Bcast(save, n, MPI_DOUBLE, k%p/*0*/, MPI_COMM_WORLD);
    for(int i = 0 ; i < local_lines ; i ++ )
    {
        schet = i*m;
        for( int j = 0 ; j < m ; j ++ )
        {
            save1 = 0; save2 = 0;
            for( int t = 0 ; t < n ; t ++ )
            {
                save1 += a[(schet+j)*n+t] * x[t];
                save2 += a[(schet+j)*n+t] * save[t];
            }
            y1=fabs(save1-save2);
            if( y1 > y2 ){ y2 = y1 ; }
        }
    }
    if(my_rank==k%p/*0*/&&l!=0)
    {
        schet = local_lines*m;
        for( int j = 0 ; j < l ; j ++ )
        {
            save1 = 0; save2 = 0;
            for( int t = 0 ; t < n ; t ++ )
            {
                save1 += a[(schet+j)*n+t] * x[t];
                save2 += a[(schet+j)*n+t] * save[t];
            }
            y1=fabs(save1-save2);
            if( y1 > y2 ){ y2 = y1 ; }
        }
    }
    MPI_Allreduce(&y2, &y1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return y1;
}

template< typename T >
void print_vector( T *b, int n)
{
    for ( int i = 0 ; i < n/*min(n,PRINT_SIZE)*/ ; i++ )
    {
        cout<<b[i]<<" ";
    }
    cout<<"\n\n";
}

void print_block( double *a, int size)
{
    for(int i = 0 ; i < min(size, PRINT_SIZE) ; i ++)
    {
        for (int j = 0; j < min(size, PRINT_SIZE); j++) {
            cout<<a[i*size+j]<<" ";
        }
        cout<<"\n";
    }
    cout<<"\n";
}

void print_br( double *a, int hight, int width)
{
    for(int i = 0; i < hight; i++)
    {
        for(int j = 0; j < width; j++)
        {
            cout<<a[i*width+j]<<" ";
        }
        cout<<"\n";
    }
    cout<<"\n";
}
