#include <iostream>
#include <pthread.h>

using namespace std;

// Define a mutex to synchronize output
pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;

void* foo(void* arg) {
    pthread_detach(pthread_self());

    // Lock the mutex before printing
    pthread_mutex_lock(&output_mutex);
    cout << "Inside the thread" << endl;
    // Unlock the mutex after printing
    pthread_mutex_unlock(&output_mutex);

    // Exit the current thread
    pthread_exit(NULL);
}

int main(int argc, char* argv[]) {
    pthread_t ptid;
    pthread_create(&ptid, NULL, &foo, NULL);

    // Lock the mutex before printing
    pthread_mutex_lock(&output_mutex);
    cout << "Main, Foo, and bar executed" << endl;
    // Unlock the mutex after printing
    pthread_mutex_unlock(&output_mutex);

    pthread_join(ptid, NULL);
    pthread_exit(NULL);
    return 0;
}
