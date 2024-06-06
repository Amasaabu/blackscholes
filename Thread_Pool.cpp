#include "Thread_Pool.h"
#include <iostream>
#include <functional>




void Thread_Pool::worker() {
	while (true) {
		Func task;
		wrk_queue.pop(task);
		task();
	}
}

Thread_Pool::Thread_Pool() {
	//get number of cores
	int num_cores = std::thread::hardware_concurrency();
	for (int i = 0; i < num_cores; i++) {
		threads.push_back(std::thread(&Thread_Pool::worker, this));
	}
}

void Thread_Pool::submit(Func F) {
	wrk_queue.push(F);
}

//ensure all threads are joined
Thread_Pool::~Thread_Pool() {
	for (auto& t : threads) {
		t.join();
	}
}