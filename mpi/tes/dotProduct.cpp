#include <iostream>
#include <vector>
#include <math.h>
#include <mpi.h>

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	int id, nproc;

	uint64_t total_size = 1 << 24;
	std::vector<int> A;
	std::vector<int> B;
	uint64_t result = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	const unsigned int begin = floor((float) id * total_size / nproc);
	const unsigned int end = floor((float) (id + 1) * total_size / nproc);
	const unsigned int size = end - begin;

	if(!id)
	{
		for(int p = nproc - 1; p > 0; p--)
		{
			unsigned int begin_pos = floor((float) p * total_size / nproc);
			unsigned int end_pos = floor((float) (p + 1) * total_size / nproc);
			unsigned int n_elem = end_pos - begin_pos;

			A.clear();
			B.clear();
			A.resize(n_elem, 2);
			B.resize(n_elem, 5);

			MPI_Send(A.data(), n_elem, MPI_INT, p, p, MPI_COMM_WORLD);
			MPI_Send(B.data(), n_elem, MPI_INT, p, p, MPI_COMM_WORLD);
		}

		A.clear();
		B.clear();
		A.resize(size, 2);
		B.resize(size, 5);

	}else
	{
        MPI_Status status;
        
		A.resize(size);
		B.resize(size);

		MPI_Recv(A.data(), size, MPI_INT, 0, id, MPI_COMM_WORLD, &status);
		MPI_Recv(B.data(), size, MPI_INT, 0, id, MPI_COMM_WORLD, &status);
	}

	// Dot product
	for(unsigned long int i = 0; i < size; i++)
		result += A[i] * B[i];

    if(!id) MPI_Reduce(MPI_IN_PLACE, &result, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    else MPI_Reduce(&result, &result, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
        
	if(!id) std::cout << "Error: " << total_size * 10 - result << std::endl;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}

