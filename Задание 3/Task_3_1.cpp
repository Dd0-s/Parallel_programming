#include <pthread.h>
#include <iostream>
#include <cmath>

#define NUM_THREADS 1
pthread_mutex_t mutex;

long part = 0;

typedef struct{
    double x, y, z;
} point;

void get_random(point& result, unsigned int& seed){
    result.x = (double(rand_r(&seed)) / double(RAND_MAX)) * M_PI;
    result.y = (double(rand_r(&seed)) / double(RAND_MAX));
    result.z = (double(rand_r(&seed)) / double(RAND_MAX)) * M_PI;
}

typedef struct{
    long count;
    unsigned int seed;
} input;

void* compute(void* param){
    input params = *(input*)param;
    point tmp;
    long result = 0;

    for(int i = 0; i < params.count; i++){
        get_random(tmp, params.seed);
        if((tmp.y <= sin(tmp.x)) and (tmp.z <= tmp.x * tmp.y)){
            result += 1;
        }
    }

    pthread_mutex_lock(&mutex);
    part += result;
    pthread_mutex_unlock(&mutex);

    pthread_exit(nullptr);
}

int main(){

    pthread_mutex_init(&mutex, nullptr);

    struct timespec start{}, end{};
    double time;
    clock_gettime(CLOCK_REALTIME, &start);

    long N = 1000000000;

    pthread_t threads[NUM_THREADS];
    input inits[NUM_THREADS];

    int number_per_thread = static_cast<int>(trunc((N - 1) / NUM_THREADS));
    int difference = static_cast<int>(N - 1) - number_per_thread * NUM_THREADS;
    long process_size;

    for(int i = 0; i < NUM_THREADS; i++) {

        if(difference != 0){
            process_size = number_per_thread + 1;
            difference -= 1;
        }else{
            process_size = number_per_thread;
        }

        inits[i].count = process_size;
        inits[i].seed = random();
        pthread_create(&threads[i], nullptr, compute, &inits[i]);
    }

    for(auto& thread : threads) {
        pthread_join(thread, nullptr);
    }

    clock_gettime(CLOCK_REALTIME, &end);
    time = static_cast<double>(end.tv_sec - start.tv_sec);
    time += static_cast<double>(end.tv_nsec - start.tv_nsec) / 1000000000;

    std::cout << "Result: " << static_cast<double>(part) / static_cast<double>(N) * M_PI * M_PI << std::endl;
    std::cout << "Time: " << time << std::endl;

    pthread_mutex_destroy(&mutex);

    return 0;
}
