#include "Thread_Pool.h"
#include <iostream>
#include <functional>




void Thread_Pool::worker() {
	while (true) {
		/*std::unique_lock<std::mutex> lck(mtx);
		cv.wait(lck, [this] { return done.load() || !wrk_queue.isEmpty();});*/

	/*	if (done.load() && wrk_queue.isEmpty()) {
			break;
		}
		Func task;
		wrk_queue.pop(task);
		task();*/

		Func task;
		if (!wrk_queue.wait_and_pop(task, done)) {
			if (done.load() && wrk_queue.isEmpty()) {
				break;  // Exit if the shutdown flag is set and the queue is empty
			}
			continue;  // Retry if no task was popped and not shutting down
		}
		task();
		
	}
}

Thread_Pool::Thread_Pool() {
	done.store(false);
	//get number of cores
	int num_cores = std::thread::hardware_concurrency();
	//spin up workers based on CPU cores
	for (int i = 0; i < num_cores; i++) {
		threads.push_back(std::thread(&Thread_Pool::worker, this));
	}
}

void Thread_Pool::submit(Func F) {
	wrk_queue.push(F);
	//cv.notify_one();
}

//ensure all threads are joined
Thread_Pool::~Thread_Pool() {
	shutdown();
}

//shut down the pool
void Thread_Pool::shutdown() {

	{
		done.store(true);
		cv.notify_all();
	}
	for (auto& t : threads) {
		t.join();
		std::cout << "After joinning" << std::endl;
	}
}