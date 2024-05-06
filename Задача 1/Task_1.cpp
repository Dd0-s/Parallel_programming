#include <iostream>
#include <vector>
#include <math.h>
#include <mpi.h>

float func(double x){
    return 4 / (1 + x * x);
}

double integral(double begin, double end, double delta){
    double result = 0;
    for(double i = 0; begin + i < end; i+= delta){
        result += (func(begin + i) + func(begin + i + delta)) * delta / 2;
    }
    return result;
}

int main(int argc, char *argv[]){
    int rank, size;
    double N, consistently;
    std::cin >> N;
    double begin, end, delta;
    double time_b1, time_b2, time_e1, time_e2;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0){
        time_b1 = MPI_Wtime();
        consistently = integral(0, 1, 1. / N);
        time_e1 = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    time_b2 = MPI_Wtime();
    if(rank == 0){
        int number = trunc(N / size);
        int difference = N - number * size;
        end = 0;
        delta = 1. / N;
        for(int i = 1; i < size; i++){
            MPI_Send(&delta, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
            begin = end;
            MPI_Send(&begin, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            if(difference != 0){
                end += (number + 1) * delta;
                difference -= 1;
            }else{
                end += number * delta;
            }
            MPI_Bsend(&end, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
        begin = end;
        end = 1;
    }else{
        MPI_Recv(&begin, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&end, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&delta, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    double result = integral(begin, end, delta);
    MPI_Barrier(MPI_COMM_WORLD);
    time_e2 = MPI_Wtime();

    if (rank == 0) {
        double answer;
        std::cout << "Process: " << 0 << "  I = " << result << std::endl;
        for (int i = 1; i < size; i++){
            MPI_Recv(&answer, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
            std::cout << "Process: " << i << "  I = " << answer << std::endl;
            result += answer;
        }
        std::cout << "Parallel: " << result << " Time: " << time_e2 - time_b2 << std::endl;
        std::cout << "Consistently: " << consistently << " Time: " << time_e1 - time_b1 << std::endl;
    }else{
        MPI_Bsend(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    return 0;
}
