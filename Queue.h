#pragma once
#include <queue>
#include <mutex>
#include <condition_variable>
#include <iostream>
#include <atomic>

template <class T>
class Queue
{
	std::mutex m;
	std::queue<T> q;
	std::condition_variable cv;
	//std::atomic<bool> done;
public:
	//Queue() : done(false) {}
	void push(T value) {
		std::lock_guard<std::mutex> lg(m);
		q.push(value);
		//notify thread to re-evaluate predicate
		cv.notify_all();
	}

	bool wait_and_pop(T& value, std::atomic<bool>& done) {
		std::unique_lock<std::mutex> lg(m);
		cv.wait_for(lg, std::chrono::milliseconds(100), [this, &done] { return !q.empty() || done.load(); });

		if (q.empty() && done.load()) {
			return false;  // Indicate shutdown if no tasks and done is true
		}

		if (!q.empty()) {
			value = std::move(q.front());
			q.pop();
			return true;
		}

		return false;
	}
	bool pop(T& value) {
		std::unique_lock<std::mutex> lg(m);
		//release lg mutex since queue is empty and wait for task to be added to the queue
		cv.wait(lg, [this] {return !q.empty();});
		if ( q.empty()) {
			return false;
		}
		value = q.front();
		q.pop();
		return true;
	}

	bool try_pop(T& value) {
		std::lock_guard<std::mutex> lg(m);
		if (q.empty()) {
			return false;  // Return false if no task is available
		}
		value = q.front();
		q.pop();
		return true;
	}
	bool isEmpty() {
		std::lock_guard<std::mutex> lg(m);
		return q.empty();
	}
	/*void shutdown() {
			std::lock_guard<std::mutex> lock(m);
			done = true;
			cv.notify_all();
	}*/

};

