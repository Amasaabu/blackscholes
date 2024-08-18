#include "Thread_Pool.h"
#include <iostream>
#include <functional>


/**
 * Worker method that runs in each thread of the thread pool.
 * Continuously pops tasks from the work queue and executes them.
 * If the shutdown flag is set and the queue is empty, the worker thread exits.
 */
void Thread_Pool::worker() {
    while (true) {
        Func task;
        if (!wrk_queue.wait_and_pop(task, done)) {
            if (done.load() && wrk_queue.isEmpty()) {
                break;  // Exit if the shutdown flag is set and the queue is empty
            }
            continue;  // Retry if no task was popped and not shutting down
        }
        // Calling supplied lambda function
        task();
    }
}

/**
 * Constructor for the Thread_Pool class.
 * Initializes the thread pool with the specified number of items.
 * Spins up worker threads based on the number of CPU cores or the number of items, whichever is smaller.
 * @param number_of_items The number of items to process in the thread pool.
 */
Thread_Pool::Thread_Pool(int number_of_items) {
    done.store(false);
    this->num_of_threads = std::thread::hardware_concurrency();
    this->num_of_itmes = number_of_items;
    num_of_threads = std::min(this->num_of_threads, number_of_items);
    // Spin up workers based on CPU cores or number of items if less than CPU cores
    for (int i = 0; i < this->num_of_threads; i++) {
        threads.push_back(std::thread(&Thread_Pool::worker, this));
    }
}

/**
 * Submits a task to the thread pool for execution.
 * Calculates the workload per thread to ensure even distribution.
 * Divides the number of items by the number of threads and assigns any remaining work to additional threads.
 * @param F The function to be executed by the thread pool.
 */
void Thread_Pool::submit(Func2 F) {
    // Calculate workload per thread to ensure even distribution
    int work_load_per_thread = num_of_itmes / num_of_threads;
    int remaining_work = num_of_itmes % num_of_threads;

    int start_row = 0;
    for (int i = 0; i < num_of_threads; i++) {
        int end_row = start_row + work_load_per_thread + (remaining_work > 0 ? 1 : 0);
        // Adjust remaining work after allocating an extra item to this thread
        if (remaining_work > 0) {
            remaining_work--;
        }

        // Capture the current start_row and end_row for this task
        int task_start = start_row;
        int task_end = end_row;

        // Submit a new task to the work queue
        wrk_queue.push([F, task_start, task_end]() {
            F(task_start, task_end);
        });

        // Update start_row for the next iteration
        start_row = end_row;
    }
}

/**
 * Destructor for the Thread_Pool class.
 * Ensures that all threads are joined before the thread pool is destroyed.
 */
Thread_Pool::~Thread_Pool() {
    shutdown();
}

/**
 * Shuts down the thread pool.
 * Sets the shutdown flag to true and waits for all worker threads to join.
 */
void Thread_Pool::shutdown() {
    {
        done.store(true);
    }
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }
}
