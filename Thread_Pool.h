#pragma once
#include <vector>
#include <thread>
#include <functional>
#include "Queue.h"
using Func = std::function<void()>;
using Func2 = std::function<void(int, int)>;

class Thread_Pool
{
	std::vector<std::thread> threads;
	Queue<Func> wrk_queue;
	void worker();
	std::atomic<bool> done;
	std::mutex mtx;
	int num_of_threads;
	int work_load_per_thread;
	int total_zize_of_data;
	int start_point;
	int end_point;

	int num_of_itmes;
public:
	Thread_Pool(int);
	~Thread_Pool();
	void submit(Func2 f);
	void shutdown();

};