#pragma once
#include <queue>
#include <mutex>
#include <condition_variable>

template <class T>
class Queue
{
	std::mutex m;
	std::queue<T> q;
	std::condition_variable cv;
public:
	void push(T value) {
		std::lock_guard<std::mutex> lg(m);
		q.push(value);
		//notify thread to re-evaluate predicate
		cv.notify_all();
	}
	void pop(T& value) {
		std::unique_lock<std::mutex> lg(m);
		//release lg mutex since queue is empty and wait for task to be added to the queue
		cv.wait(lg, [this] {return !q.empty();});
		value = q.front();
		q.pop();
	}

};

