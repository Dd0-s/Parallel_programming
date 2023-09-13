#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>

double exact_solution_dot(double x, double t){
    double k1 = 4 / M_PI;
    double k2 = - pow(M_PI, 2) * t;
    double k3 = M_PI * x;
    double temperature = 0;
    for(int i = 0; i <= 1000; i++){
        temperature += k1 * exp(k2 * pow(2 * i + 1, 2)) / (2 * i  + 1) * sin(k3 * (2 * i + 1));
    }
    return temperature;
}

std::vector<double> exact_solution(int N, double T){
    std::vector<double> temperatures(11);
    int number = static_cast<int>(trunc(N / 10));
    int difference = N - number * 10;
    int step = 0;

    temperatures[0] = 0;
    for(int i = 1; i < 10; i++){
        if(difference != 0){
            step += number + 1;
            difference -= 1;
        }else{
            step += number;
        }
        temperatures[i] = exact_solution_dot(1. / N * step, T);
    }
    temperatures[10] = 0;

    return temperatures;
}

double u_next(double left, double bottom, double right, double c){
    return bottom + c * (left - 2 * bottom + right);
}

void calculate(double data[2], std::vector<double>& steps, int size, double c){
    double help_1;
    double help_2 = steps[0];
    steps[0] = u_next(data[0], steps[0], steps[1], c);
    for(int i = 1; i < size - 1; i++){
        help_1 = steps[i];
        steps[i] = u_next(help_2, steps[i], steps[i + 1], c);
        help_2 = help_1;
    }
    steps[size - 1] = u_next(help_2, steps[size - 1], data[1], c);
}

int main(int argc, char *argv[]){
    int rank, size, N = 50;
    double T = 0.1;
    double data[2];
    double time_b, time_e;
    //std::cin >> N;

    MPI_Status Status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    time_b = MPI_Wtime();

    if(rank == 0){
        for(int i = 1; i < size; i++){
            MPI_Send(&N, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }else{
        MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
    }

    int process_size;
    std::vector<int> process_sizes(size);
    if(rank == 0){
        int number_per_process = static_cast<int>(trunc((N - 1) / size));
        int difference = (N - 1) - number_per_process * size;

        for(int i = size - 1; i > 0; i--){
            if(difference != 0){
                process_size = number_per_process + 1;
                difference -= 1;
            }else{
                process_size = number_per_process;
            }
            process_sizes[i] = process_size;
            MPI_Send(&process_size, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
        }

        process_size = number_per_process;
        process_sizes[0] = process_size;
        data[0] = 0;
        data[1] = 1;
        if(size == 1){
            data[1] = 0;
        }
    }else{
        MPI_Recv(&process_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &Status);
        data[0] = 1;
        data[1] = 1;
        if(rank == size - 1){
            data[1] = 0;
        }
    }

    std::vector<double> steps(process_size);
    for(int i = 0; i < process_size; i++){
        steps[i] = 1;
    }
    double h = 1. / N;
    double dt = h * h / 2;
    int iterations = floor(T / dt);
    double c = dt / (h * h);

    for(int i = 0; i < iterations; i++){
        calculate(data, steps, process_size, c);

        if(rank > 0){
            MPI_Send(&steps[0], 1, MPI_DOUBLE, rank - 1, 2, MPI_COMM_WORLD);
        }

        if(rank < size - 1){
            MPI_Recv(&data[1], 1, MPI_DOUBLE, rank + 1, 2, MPI_COMM_WORLD, &Status);
        }

        if(rank < size - 1){
            MPI_Send(&steps[process_size - 1], 1, MPI_DOUBLE, rank + 1, 3, MPI_COMM_WORLD);
        }

        if(rank > 0){
            MPI_Recv(&data[0], 1, MPI_DOUBLE, rank - 1, 3, MPI_COMM_WORLD, &Status);
        }

    }

    std::vector<double> full_temperatures;
    if(rank == 0){
        for(int j = 0; j < process_sizes[0]; j++){
            full_temperatures.push_back(steps[j]);
        }
        for(int i = 1; i < size; i++){
            std::vector<double> process_answer(process_sizes[i]);
            MPI_Recv(&process_answer[0], process_sizes[i], MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &Status);
            for(int j = 0; j < process_sizes[i]; j++){
                full_temperatures.push_back(process_answer[j]);
            }
        }

        std::vector<double> temperatures(11);
        temperatures[0] = 0;
        int number = static_cast<int>(trunc(N / 10));
        int difference = N - number * 10;
        int step = -1;

        for(int i = 1; i < 10; i++){
            if(difference != 0){
                step += number + 1;
                difference -= 1;
            }else{
                step += number;
            }
            temperatures[i] = full_temperatures[step];
        }
        temperatures[10] = 0;

        time_e = MPI_Wtime();

        std::cout << std::endl;
        std::vector<double> ex_solution = exact_solution(N, T);
        std::cout << "The exact solution:" << std::endl;
        for(int i = 0; i < 11; i++){
            std::cout << ex_solution[i] << ' ';
        }
        std::cout << std::endl << std::endl;

        std::cout << "The approximation solution:" << std::endl;
        for(int i = 0; i < 11; i++){
            std::cout << temperatures[i] << ' ';
        }
        std::cout << std::endl << std::endl;

        std::cout << "Mean squared error: ";
        double error = 0.0;
        for(int i = 0; i < 11; i++){
            error +=  pow(temperatures[i] - ex_solution[i], 2);
        }
        std::cout << pow(error, 0.5) << std::endl << std::endl;

        std::cout << "Time: " << time_e - time_b << std::endl;

    }else{
        MPI_Send(&steps[0], process_size, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
