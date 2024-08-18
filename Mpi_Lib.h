#pragma once
#include "Mpi_Lib.h"
#include "mpi.h"
#include <vector>
#include "Thread_Pool.h"
using namespace std;
using Funca = std::function<void(int, int)>;

template <class T, class U, class R>
class Mpi_Lib
{
public:
	Mpi_Lib(int&, char**);
	void init(int&, int elements_per_unit = 1);
	void broadcast(U*, int, int, MPI_Datatype);
	int get_world_rank();
	vector<int> get_sendcounts();
	void scatterV(T*, T*, MPI_Datatype, Funca f);
	void gather_v(R* local_result, R* result, MPI_Datatype, bool use_element_per_unit_scaling = false);
	int* get_displs();
	void barrier();
	~Mpi_Lib();

private:
	int world_rank;
	int world_size;
	int remainder;
	int size_of_data;
	int elements_per_unit;
	vector<int> sendcounts;
	vector<int> displs;
	//int* local_data;

	int elem_per_process;
	int extra_rows;
};


/**
 * @brief Constructor for the Mpi_Lib class.
 * 
 * @param argc The number of command line arguments.
 * @param argv The command line arguments.
 */
template <class T, class U, class R>
Mpi_Lib<T, U, R>::Mpi_Lib(int& argc, char** argv)
{
	// Initialize the MPI environment
	MPI_Init(&argc, &argv);
	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &this->world_rank);
	// Get the number of processors
	MPI_Comm_size(MPI_COMM_WORLD, &this->world_size);
}


/**
 * @brief Initializes the Mpi_Lib object.
 * 
 * @param size_v The size of the data.
 * @param elements_per_unit The number of elements per unit.
 */
template <class T, class U, class R>
void Mpi_Lib<T, U, R>::init(int& size_v, int elements_per_unit) {


	//initialize the sendcounts and displs
	this->sendcounts.resize(this->world_size);
	this->displs.resize(this->world_size);
	this->size_of_data = size_v;
	this->elements_per_unit = elements_per_unit;
	this->elem_per_process = this->size_of_data / world_size;
	this->extra_rows = this->size_of_data % world_size;
	int current_displ = 0;
	int sum = 0;
	for (int i = 0; i < world_size; ++i) {

		// Calculate the number of rows for this process
		int rows_for_this_process = elem_per_process + (i < extra_rows ? 1 : 0);

		// Calculate the number of elements (rows * columns) for this process
		this->sendcounts[i] = rows_for_this_process * this->elements_per_unit;
		// Set the displacement for this process
		displs[i] = current_displ;

		// Update the displacement for the next process
		current_displ += this->sendcounts[i];
	}
}

/**
 * @brief Broadcasts data from the root process to all other processes.
 *
 * @param data_to_brodcast The data to be broadcasted.
 * @param count The number of elements in the data.
 * @param root The rank of the root process.
 * @param type_of_data_to_brodcast The MPI datatype of the data.
 */
template <class T, class U, class R>
void Mpi_Lib<T, U, R>::broadcast(U* data_to_brodcast, int count, int root, MPI_Datatype type_of_data_to_brodcast)
{
	MPI_Bcast(data_to_brodcast, count, type_of_data_to_brodcast, root, MPI_COMM_WORLD);
}

	/**
	 * @brief Returns the rank of the current process in the MPI_COMM_WORLD communicator.
	 * 
	 * @return The rank of the current process.
	 */
template <class T, class U, class R>
int Mpi_Lib<T, U, R>::get_world_rank() {
	return this->world_rank;
}



template <class T, class U, class R>
	/**
	 * @brief Returns the vector of sendcounts used for scatterv operation.
	 * 
	 * @return The vector of sendcounts.
	 */
vector<int> Mpi_Lib<T, U, R>::get_sendcounts() {
	return this->sendcounts;
}




/**
 * @brief Scatters data from the root process to all other processes using MPI_Scatterv.
 *
 * This method scatters the data from the root process to all other processes using the MPI_Scatterv function.
 * It divides the data_to_scatter array into smaller chunks and distributes them to the local_data array of each process.
 * After scattering the data, it splits the work across multiple threads using a thread pool and submits the provided function f for execution.
 * Finally, it ensures that all threads are joined before returning.
 *
 * @param data_to_scatter The data to be scattered.
 * @param local_data The local data array to store the scattered data.
 * @param type_to_scatter The MPI datatype of the data.
 * @param f The function to be executed by the thread pool.
 */
template <class T, class U, class R>
void Mpi_Lib<T, U, R>::scatterV(T* data_to_scatter, T* local_data, MPI_Datatype type_to_scatter, Funca f) {
	// Scatter the data from the root process to all other processes
	MPI_Scatterv(data_to_scatter, this->sendcounts.data(), displs.data(), type_to_scatter, local_data, sendcounts[world_rank], type_to_scatter, 0, MPI_COMM_WORLD);

	// Split the work across threads
	int threadPoolSize = sendcounts[world_rank] / this->elements_per_unit;
	
	// Ensure every processor has work to do
	if (sendcounts[world_rank] == 0) {
		cout << "-------------Process " << world_rank << " received zero data." << endl;
		return;
	}
	
	// Create a thread pool and submit the function for execution
	Thread_Pool thread_pool(threadPoolSize);
	thread_pool.submit(f);
	
	// Ensure all threads are joined
	thread_pool.shutdown();
}

//Custom Gathrev
/**
 * @brief Gathers the results from all processors using MPI_Gatherv.
 *
 * This method gathers the results from all processors using the MPI_Gatherv function.
 * It calculates the sendcounts and displs arrays based on the number of rows each process is responsible for.
 * The gathered results are stored in the result array.
 *
 * @param local_result The local result array containing the results from each process.
 * @param result The array to store the gathered results.
 * @param type_of_result The MPI datatype of the result.
 * @param use_element_per_unit_scaling Flag indicating whether to scale the sendcounts based on the number of elements per unit.
 */
template <class T, class U, class R>
void Mpi_Lib<T, U, R>::gather_v(R* local_result, R* result, MPI_Datatype type_of_result, bool use_element_per_unit_scaling) {
	//gather the results from all processors
	vector<int> sendcounts_gather(this->world_size, 0);
	vector<int> displs_gather(this->world_size, 0);

	int curr_displ = 0;
	for (int i = 0; i < this->world_size; ++i) {
		// Calculate the number of rows this process is responsible for
		int rows_for_this_proc = elem_per_process + (i < extra_rows ? 1 : 0);
		if (use_element_per_unit_scaling) {
			// Each process sends back rows_for_this_proc * size_v elements
			sendcounts_gather[i] = rows_for_this_proc * this->elements_per_unit;
		}
		else {
			sendcounts_gather[i] = rows_for_this_proc;
		}

		// Displacement for this process's data in the gathered array
		displs_gather[i] = curr_displ;

		// Update displacement for the next process
		curr_displ += sendcounts_gather[i];
	}
	MPI_Gatherv(local_result, sendcounts_gather[world_rank], type_of_result, result, sendcounts_gather.data(), displs_gather.data(), type_of_result, 0, MPI_COMM_WORLD);
}

/**
 * @brief Returns a pointer to the displacements array used for scatterv and gatherv operations.
 * 
 * @return A pointer to the displacements array.
 */
template <class T, class U, class R>
int* Mpi_Lib<T, U, R>::get_displs() {
	return displs.data();
}

/**
 * @brief Blocks until all processes in the MPI_COMM_WORLD communicator have reached this point.
 */
template <class T, class U, class R>
void Mpi_Lib<T, U, R>::barrier() {
	MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * @brief Destructor for the Mpi_Lib class.
 * 
 * This method finalizes the MPI environment and cleans up any resources used by the Mpi_Lib object.
 */
template <class T, class U, class R>
Mpi_Lib<T, U, R>::~Mpi_Lib()
{
	// Finalize the MPI environment
	MPI_Finalize();
	//delete local_data
}