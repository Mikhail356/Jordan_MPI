//Gvozdev 17.02.2019
#include "jordan.h"

int
main(int argc, char *argv[])
{
    MPI_Init( &argc , &argv );
    int n , m , k , p , l , my_rank , err = 0, local_lines = 0; //Ax=b
    double *a , *b , *x , *save , total_time = 0, error=0, residual=0;
    char * name = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if ( argc < 3 || argc > 4 || ((n = atoi(argv[1])) <= 0) || ((m = atoi(argv[2])) <= 0) )
    {
        if(my_rank==0)  printf("usage: %s n m [name]\n", argv[0]); //n division matrix a
        MPI_Finalize();
        return 0;
    }
    else if(argc == 4) { name = argv[3]; }
    if( m%3 != 0 )
    {
        if(my_rank==0)  printf("Error: m/3 is not integer\n");
        MPI_Finalize();
        return 0;
    }
    if( n < m )
    {
        if(my_rank==0)  printf("Error: n < m \n");
        MPI_Finalize();
        return 0;
    }

    k = n / m ;
    l = n - m*k ;
    for(int i = my_rank; i < k; i+=p)
    {
        local_lines++;
    }
//LOG(my_rank);
    data_init( a , x , b , save , n , k , m , l , p , my_rank, local_lines);
    //print_br(save, m, n*m);
    //print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
    //print_br(b, my_rank==k%p?local_lines*m+l:local_lines*m, 1);
    //print_br(x, n, 1);

    if(name!=0)
    {
        //LOG(my_rank); LOG("read");
        err = read_matrix(a, n, m, k, l, p, my_rank, name, local_lines, save);
        if( err < 0 )
        {
            if( my_rank == 0 )
            {
                switch (err)
                {
                    case ERR_CANNOT_OPEN:
                        printf("Cannot open file %s \n",name);
                        break;
                    case ERR_READ:
                        printf("Error read in file %s \n", name);
                        break;
                    default:
                        printf("Unknown error\n");
                }
            }

            delete [] a; delete [] b; delete []x; delete []save;
            MPI_Finalize();
            return 0;
        }
    }
    else
    {
        //LOG(my_rank); LOG("init");
        init_matrix( a , n , k , m , l , p , my_rank, local_lines);
    }
    //print_br(save, m, n*m);
    //print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
    //print_br(b, my_rank==k%p?local_lines*m+l:local_lines*m, 1);
    //print_br(x, n, 1);
    for( int i = 0 ; i < n ; i ++ ) {
        x[i] = i%2;
    }

    matrix_mult_vector (a, b, x, n, k, m, l, p, my_rank, save, local_lines);

    //print_br(save, m, n*m);
    //print_br(a, my_rank==k%p?local_lines*m+l:local_lines*m, n);
    //print_br(b, my_rank==k%p?local_lines*m+l:local_lines*m, 1);
    //print_br(x, n, 1);
    print_matrix( a , n , m , k , l , p , my_rank , min(n,PRINT_SIZE) , save , local_lines);
//    if(my_rank==0) print_vector(x, n);
//LOG("Vhoju v solve");
    if(my_rank==0){total_time = get_full_time();}
    err = solve( a , b , x , save , n , k , m , l , p , my_rank, local_lines );
    if(my_rank==0){total_time = get_full_time() - total_time;}

//MPI_Bcast(save, k*m+l, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//LOG("vishel is solva");
//print_vector(save, n);
    if(err == MATRIX_DEGENERATE)
    {
        if(my_rank == 0)
            printf("Matrix degenerate\n");
        delete [] a; delete [] b; delete []x; delete []save;
        MPI_Finalize();
        return 0;
    }

//    print_vector(save, n);
    //LOG("x");print_vector(x,n);
    error = search_error(x, save, k, m, l, p, local_lines, my_rank);
    memcpy (x, save, n*sizeof(double));
//LOG(error);
//    print_vector(save, n);
    //LOG("init");
    if(name!=0)
    {
        //LOG("Goo");
        read_matrix( a , n , m , k , l , p , my_rank , name , local_lines, save);
        //LOG("gle");
    }
    else
    {
        init_matrix( a , n , k , m , l , p , my_rank, local_lines );
    }
//LOG(10);
    memcpy (save, x, n*sizeof(double));
    for(int i = 0; i < n ; i ++)
    {
        x[i] = i % 2;
    }
//    print_vector(save, n);
    //LOG(10.1);
//    print_vector(x,n);
    //LOG(11);
    residual = search_residual(a, x, save, n, m, k, l, p, my_rank, local_lines);
//LOG(12);
    if(my_rank==0)
    {
        printf("Residual = %e  Error = %e  Total time: %.2lf  n = %d  m = %d  p = %d\n",
               residual, error, total_time, n, m, p);
    }

    /*if(local_lines>0){*/delete [] a; delete [] b; delete []x; delete []save;/*}*/
    MPI_Finalize();
    return 0;
}

void
data_init ( double *&a , double *&x , double *&b , double *&save,
            int n , int k , int m , int l , int p , int my_rank,
            int local_lines)
{
    int kol = 0 ;

    kol = local_lines * m;
    if(my_rank==k%p) kol+=l;

    a = new double [(kol) * n];
    b = new double [kol];
    x = new double [n];
    memset(a, 0, kol*n);
    memset(b, 0, kol);
    memset(x, 0, n);

    kol = ((n * m) + m);
    save = new double [ kol ] ;
    memset(save, 0, kol);
}

void
init_matrix ( double *a , int n , int k , int m , int l , int p , int my_rank,
              int local_lines)
{
    for(int i = my_rank ; i < k ; i += p)
    {
        for( int j = 0 ; j < m ; j ++ )
        {
            for( int t = 0 ; t < n ; t ++ )
            {
                a[((i/p)*m+j)*n+t]=abs((i*m)+j-t);
            }
        }
    }
    if(my_rank==k%p&&l!=0)
    {
        for( int j = 0 ; j < l ; j ++ )
        {
            for( int t = 0 ; t < n ; t ++ )
            {
                a[(local_lines*m+j)*n+t]=abs(k*m+j-t);
            }
        }
    }
}

int
read_matrix (double *a, int n, int m, int k, int l, int p, int my_rank,
             char * name, int local_lines, double *save)
{
    int err_read = 0 , nm = 0/* , loc = 0*/;
    MPI_Status status;
    ifstream file;

    nm = n * m;

    if(my_rank == k%p/*0*/)
    {
        file.open(name);
        if (!(file.is_open())) err_read = ERR_CANNOT_OPEN;
    }

    MPI_Bcast( &err_read , 1 , MPI_INT , k%p/*0*/ , MPI_COMM_WORLD );
    if(err_read < 0 ) return err_read;

    if(my_rank == k%p)
    {
        for(int j = 0 ; j < k ; j ++ )
        {
//            if( j % p != my_rank )
//            {
//                loc = (((j-my_rank)/p)+1)*nm;
//            }
//            else
//            {
//                loc = ((j-my_rank)/p)*nm;
//            }

            for(int i = 0 ; i < m && err_read == 0 ; i++ )
            {
                for( int t = 0 ; t < n ; t ++ )
                {
                    if(!(file>>save[i*n+t]/*a[loc + i*n + t]*/))
                    {
                        err_read = ERR_READ;
                        break;
                    }
//                     LOG(loc + i*n + t);
                }
            }

            if(j%p != my_rank)
            {
                MPI_Send( &err_read , 1 , MPI_INT , j%p , 0 , MPI_COMM_WORLD );
                if(err_read == 0)
                {
                    MPI_Send(save , nm , MPI_DOUBLE , j%p , 1+(j/p) , MPI_COMM_WORLD );
                }
            }
            else
            {
                memcpy (a+(j/p)*nm, save, n*m*sizeof(double));
            }
        }

        if( l != 0 && err_read == 0 )
        {
//            loc = ((k-my_rank)/p)*nm;
            for(int i = 0 ; i < l && err_read==0 ; i++ )
            {
                for( int t = 0 ; t < n ; t ++ )
                {
                    if(!(file>>save[i*n+t]/*a[local_lines*n*m + i*n + t]*/))
                    {
                        err_read = ERR_READ;
                        break;
                    }
//                    LOG(local_lines*n*m + i*n + t);
                }
            }
            memcpy(a+local_lines*m*n, save, l*n*sizeof(double));
        }

        file.close();
    }
    else
    {
//        loc=0;
//        while((k-loc)%p!=0)
//            loc++;

        for( int i = my_rank ; i < k ; i += p)
        {
            MPI_Recv(&err_read, 1, MPI_INT, k%p, 0, MPI_COMM_WORLD, &status);
            if(err_read == 0)
            {
                MPI_Recv(a+((i/p)*nm), m*n, MPI_DOUBLE, k%p, 1+(i/p),
                         MPI_COMM_WORLD, &status);
            }
            else break;
        }
    }
    MPI_Bcast(&err_read, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return err_read;
}

void
print_matrix ( double *a , int n , int m , int k , int l , int p , int my_rank,
                   int size , double * save, int local_lines)
{
    MPI_Status status;
    if(my_rank == k%p)
    {
        //LOG(my_rank);
        for( int i = 0 ; i < k ; i ++ )
        {
            if(i%p != my_rank)
            {
                MPI_Recv( save , m*n , MPI_DOUBLE , i%p , i/p , MPI_COMM_WORLD , &status );
            }
            else memcpy( save , a+((i/p)*n*m) , n*m*sizeof(double));
            for( int j = 0 ; j < m ; j ++ )
            {
                if(i*m+j >= size) break;
                for( int t = 0 ; t < size ; t ++ )
                {
                    cout<<save[j*n+t]<<" ";
                }
                cout<<"\n";
            }
        }
        if(l!=0)
        {/*
            int kol = k/p;
            if(k%p!=0) kol++;*/
            for( int j = 0 ; j < l && k*m+j < size ; j ++ )
            {
                for( int t = 0 ; t < size ; t ++ )
                {
                    cout<<a[(local_lines*m+j)*n+t]<<" ";
                }
                cout<<"\n";
            }
        }
        cout<<"\n\n";
    }
    else
    {
        //LOG(my_rank);
        for( int i = my_rank ; i < k ; i += p )
        {
            MPI_Send( a+((i/p)*m*n) , m*n , MPI_DOUBLE , k%p/*0*/ , i/p , MPI_COMM_WORLD );
        }
    }
}

void
init_x (double *x, int /*n*/, int k, int m, int l, int p, int my_rank)
{
    int count = 0;
    for(int i = my_rank ; i < k ; i+=p )
    {
        count = (i/p)*m;
        for( int j = 0 ; j < m ; j++ )
        {
            x[count+j] = (k*m+j)%2;
        }
    }

    if(my_rank == 0)
    {
        count = (k/p)*m;
        for( int j = 0 ; j < l ; j++ )
        {
            x[count+j] = (k*m+j)%2;
        }
    }
}

void
matrix_mult_vector (double *a, double *b, double *x, int n, int k, int m, int l, int p,
                    int my_rank, double* /*save*/, int local_lines)
{
    double count = 0 ;
//    if(my_rank==0)
//        memcpy(save, x, n*sizeof(double));
//    MPI_Bcast(save, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for(int i = 0 ; i < local_lines ; i ++ )
    {
        for( int j = 0 ; j < m ; j ++ )
        {
            count = 0;
            for( int t = 0 ; t < n ; t ++ )
            {
                count += a[(i*m+j)*n+t] * /*save*/x[t];
            }
            b[i*m+j] = count;
        }
    }
    if( my_rank==k%p && l!=0 )
    {
        for( int j = 0 ; j < l ; j ++ )
        {
            count = 0 ;
            //LOG(j);
            //LOG(((local_lines)*m+j)*n);
            //LOG(a[((local_lines)*m+j)*n]);
            //LOG(save[0]);
            for( int t = 0 ; t < n ; t ++ )
            {
                count += a[((local_lines)*m+j)*n+t] * /*save*/x[t];
                //LOG(t);
                //LOG(save[t]);
            }
            b[(local_lines)*m+j] = count;
        }
    }
}

template< typename T >
void print_vector( T *b, int n)
{
    for ( int i = 0 ; i < min(n,PRINT_SIZE) ; i++ )
    {
        cout<<b[i]<<" ";
    }
    cout<<"\n\n";
}
