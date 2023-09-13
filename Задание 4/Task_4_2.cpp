#include <pthread.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <semaphore.h>

#define NUM_THREADS 8
sem_t sem_read[NUM_THREADS];
sem_t sem_write[NUM_THREADS];

const double c = 1.0;
const double T = 0.0001;
const int N = 1000000;
const double dt = 0.0000001;
const double h = 1.0 / N;
const long iterations = static_cast<long>(T / dt);
const double nu = dt * c / h;
double result[10 * N + 1];

typedef struct{
    double left;
    long count;
    int number;
    int index;
} input;

double round_2(double x){
    if (std::round(x * 100) / 100 == 0){
        return 0;
    }else{
        return std::round(x * 100) / 100;
    }
}

double g(double x){
    if((0 < x) and (x < 2)){
        return x * (x - 2);
    }
    return 0;
}

double exact_solution_dot(double x, double t){
    return g(x - c * t);
}

std::vector<double> exact_solution(){
    std::vector<double> solution(10 * N + 1);

    for(int i = 0; i < 10 * N + 1; i++){
        solution[i] = exact_solution_dot(h * i, c * T);
    }

    return solution;
}

double u_next(double now, double left){
    return now - nu * (now - left);
}

void* compute(void* param){
    input params = *static_cast<input *>(param);

    for(int j = 0; j < iterations; j++) {
        sem_wait(&sem_write[params.index]);
        double help_1 = result[params.number + 1];
        double help_2;
        result[params.number + 1] = u_next(result[params.number + 1], params.left);
        for (int i = params.number + 2; i <= params.number + params.count; i++) {
            help_2 = result[i];
            result[i] = u_next(result[i], help_1);
            help_1 = help_2;
        }

        if(params.index < NUM_THREADS - 1) {
            sem_post(&sem_read[params.index + 1]);
        }else{
            sem_post(&sem_write[params.index]);
        }

        if(params.index > 0) {
            sem_wait(&sem_read[params.index]);
            params.left = result[params.number];
            sem_post(&sem_write[params.index - 1]);
        }

    }
    pthread_exit(nullptr);
}

int main(){

    struct timespec start{}, end{};
    double time;
    clock_gettime(CLOCK_REALTIME, &start);

    pthread_t threads[NUM_THREADS];
    input information[NUM_THREADS];

    for(int i = 0; i < 10 * N + 1; i++){
        result[i] = g(h * i);
    }

    int number_per_thread = static_cast<int>(trunc(10 * N / NUM_THREADS));
    int difference = 10 * N - number_per_thread * NUM_THREADS;
    int process_size;
    int length = 0;

    for(auto& info : information) {

        if(difference != 0){
            process_size = number_per_thread + 1;
            difference -= 1;
        }else{
            process_size = number_per_thread;
        }

        info.count = process_size;
        info.left = g(length * h);
        info.number = length;
        length += process_size;
    }

    for(int i = 0; i < NUM_THREADS; i++) {
        sem_init(&sem_read[i], 0, 0);
        sem_init(&sem_write[i], 0, 1);
    }

    for(int i = 0; i < NUM_THREADS; i++) {
        information[i].index = i;
        pthread_create(&threads[i], nullptr, compute, &information[i]);
    }

    for(auto& thread : threads){
        pthread_join(thread, nullptr);
    }

    clock_gettime(CLOCK_REALTIME, &end);
    time = static_cast<double>(end.tv_sec - start.tv_sec);
    time += static_cast<double>(end.tv_nsec - start.tv_nsec) / 1000000000;

    std::vector<double> ex_solution = exact_solution();
    std::cout << "The exact solution:" << std::endl;
    for(int i = 0; i < 10 * N + 1; i += N / 2){
        std::cout << round_2(ex_solution[i]) << ' ';
    }
    std::cout << std::endl << std::endl;

    std::cout << "The approximation solution:" << std::endl;
    for(int i = 0; i < 10 * N + 1; i += N / 2){
        std::cout << round_2(result[i]) << ' ';
    }
    std::cout << std::endl << std::endl;

    std::cout << "Mean squared error: ";
    double error = 0.0;
    for(int i = 0; i < 10 * N + 1; i++){
        error +=  pow(result[i] - ex_solution[i], 2);
    }
    std::cout << pow(error, 0.5) << std::endl << std::endl;

    std::cout << "Time: " << time << std::endl;

    for(int i = 0; i < NUM_THREADS; i++) {
        sem_destroy(&sem_write[i]);
        sem_destroy(&sem_read[i]);
    }

    return 0;
}
