#pragma once
#include <vector>
#include <thread>
#include <functional>
#include "Queue.h"
using Func = std::function<void()>;
class Thread_Pool
{
	std::vector<std::thread> threads;
	Queue<Func> wrk_queue;
	void worker();
	std::atomic<bool> done;
	std::mutex mtx;
public:
	Thread_Pool();
	~Thread_Pool();
	void submit(Func f);
	void shutdown();
};

