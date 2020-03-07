#include "mpi.h"
#include <iostream>
#include <cstdlib>

inline int fact(const int N)
{
	if(N == 1)
		return 1;
	else
		return N * fact(N - 1);
}

inline void print(int *P, const int N)
{
	for(int i = 0; i < N; i++) std::cout << P[i];
}

void swap(int *P2, int i1, int j1)
{
	int temp = P2[i1];
	P2[i1] = P2[j1];
	P2[j1] = temp;
}


int distance(int *P2, const int N, int **D2)
{
	int d = D2[0][P2[1]] + D2[P2[N-1]][0];
	for(int j = 2; j < N; j++)
		d += D2[P2[j-1]][P2[j]];
	return d;
}

void perm(int *P1, const int N, int **D1, int rk, int sz, MPI_Status stat)
{
	int *min_dist = new int[sz];
	for(int i = 0; i < sz; i++)
        min_dist[i] = distance(P1, N, D1);
	int **res_P = new int*[sz];
	for(int i = 0; i < sz; i++)
        res_P[i] = new int[N];
    double t1, t2;
    t1 = MPI_Wtime();

	for(int k = rk + 1; k < N; k += sz){
        if(k != 1){
            for(int p = 2; p <= k; p++)
                P1[p] = p - 1;
        }
        if(k != N - 1){
            for(int p = k + 1; p < N; p++)
                P1[p] = p;
        }
        P1[1] = k;
        print(P1, N);
        std::cout << "\tDistance = " << distance(P1, N, D1) << " rank = " << rk << std::endl;

        for(int q = 1; q < fact(N-2); q++){
            int f = 2;
            int r = N - 1;
            int i = r - 1;
            while(i > f && P1[i] >= P1[i+1]) i--;

            int j = r;
            while(P1[i] >= P1[j]) j--;
            swap(P1, i, j);

            f = i + 1;
            while(f < r) swap(P1, f++, r--);
            print(P1, N);
            std::cout << "\tDistance = " << distance(P1, N, D1) << " rank = " << rk << std::endl;

            if(distance(P1, N, D1) < min_dist[rk]){
                min_dist[rk] = distance(P1, N, D1);
                for(int v = 0; v < N; v++)
                    res_P[rk][v] = P1[v];
            }
        }
	}

    int *min_dist_z = new int[sz];
    MPI_Gather(&min_dist[rk], 1, MPI_INT, &min_dist_z[rk], 1, MPI_INT, 0, MPI_COMM_WORLD);

    int x;
    if(rk == 0){
        int t = 100;
        for(int i = 0; i < sz; i++){
            if(min_dist_z[i] < t){
                x = i;
                t = min_dist_z[i];
            }
        }
        std::cout << "MINIMAL DISTANCE = " << t;

    }
    MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(rk == x){
        std::cout << " PERMUTATION WITH MINIMAL DISTANCE = ";
        for(int i = 0; i < N; i++)
            std::cout << res_P[rk][i];
        std::cout << std::endl;
    }

    if(rk == 0){
        t2 = MPI_Wtime();
        printf("TIME OF EXECUTION\t%f", t2 - t1);
    }



}

int main(int argc, char *argv[])
{
	const int N = 6;

    int rank, size;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	int **D = new int *[N];
	for(int i = 0; i < N; i++)
		D[i] = new int[N];

	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			if(i == j)
				D[i][j] = -1;
			else if(i > 0 && j < i)
				D[i][j] = 1 + rand() % 20;

	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			if(i > 0 && j < i)
				D[j][i] = D[i][j];

	if(rank == 0){
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                    std::cout << D[i][j] << "\t";
            }
            std::cout << std::endl;
        }
	}
	int *P = new int [N];
	for(int i = 0; i < N; i++)
		P[i] = i;

	perm(P, N, D, rank, size, status);
	delete P;
	delete [] D;
	MPI_Finalize();
	return 0;
}
