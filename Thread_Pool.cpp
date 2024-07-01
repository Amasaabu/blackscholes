#include "Thread_Pool.h"
#include <iostream>
#include <functional>

void Thread_Pool::worker() {
	while (true) {

		Func task;
		if (!wrk_queue.wait_and_pop(task, done)) {
			if (done.load() && wrk_queue.isEmpty()) {
				break;  // Exit if the shutdown flag is set and the queue is empty
			}
			continue;  // Retry if no task was popped and not shutting down
		}
		//calling supplied lambd
		task();

	}
}


Thread_Pool::Thread_Pool(int number_of_items) {
	done.store(false);
	//this->start_point = start_point;
	//this->end_point = end_point;
	this->num_of_threads = std::thread::hardware_concurrency();
	//this->work_load_per_thread = (end_point-start_point) / this->num_of_threads;
	//this->total_zize_of_data = end_point - start_point;
	//dont use more threads than data#
	this->num_of_itmes = number_of_items;
	num_of_threads = std::min(this->num_of_threads,number_of_items);
	//spin up workers based on CPU cores
	for (int i = 0; i < this->num_of_threads; i++) {
		threads.push_back(std::thread(&Thread_Pool::worker, this));
	}
}

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

//Thread_Pool::Thread_Pool(int total_zize_of_data) {
//	done.store(false);
//	this->num_of_threads = std::thread::hardware_concurrency();
//	this->work_load_per_thread = total_zize_of_data / this->num_of_threads;
//	this->total_zize_of_data = total_zize_of_data;
//	//spin up workers based on CPU cores
//	for (int i = 0; i < this->num_of_threads; i++) {
//		int start_row = i * this->work_load_per_thread;
//		int end_row = (i == this->num_of_threads - 1) ? total_zize_of_data : start_row + this->work_load_per_thread;
//		threads.push_back(std::thread(&Thread_Pool::worker, this));
//	}
//}

//void Thread_Pool::submit(Func2 F) {
//	for (int i = 0; i < this->num_of_threads; i++) {
//		int start_row = i * this->work_load_per_thread;
//		int end_row = (i == this->num_of_threads - 1) ? this->total_zize_of_data : start_row + this->work_load_per_thread;
//		wrk_queue.push([start_row, end_row, F]() {F(start_row, end_row); });
//
//	}
//	//wrk_queue.push(F);
//}

//ensure all threads are joined
Thread_Pool::~Thread_Pool() {
	shutdown();
}

//shut down the pool
void Thread_Pool::shutdown() {

	{
		done.store(true);
	}
	for (auto& t : threads) {
		if (t.joinable()) t.join();
	}
}
