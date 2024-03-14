#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

// This function will be executed by the first thread
void *thread_function1(void *arg) {
    for (int i = 0; i < 5; i++) {
        printf("Thread 1: %d\n", i);
        // Sleep for a short time to simulate some work
        usleep(100000); // 100 millisecondsgcc
    }
    return NULL;
}

// This function will be executed by the second thread
void *thread_function2(void *arg) {
    for (int i = 0; i < 5; i++) {
        printf("Thread 2: %d\n", i);
        // Sleep for a short time to simulate some work
        usleep(100000); // 100 milliseconds
    }
    return NULL;
}

int main() {
    pthread_t thread1, thread2;

    // Create two threads
    if (pthread_create(&thread1, NULL, thread_function1, NULL) != 0) {
        perror("pthread_create");
        return 1;
    }
    if (pthread_create(&thread2, NULL, thread_function2, NULL) != 0) {
        perror("pthread_create");
        return 1;
    }

    // Wait for both threads to finish
    pthread_join(thread1, NULL);
    pthread_join(thread2, NULL);

    printf("Both threads have finished.\n");

    return 0;
}
