#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <mpi.h>
#include <cmath>

using namespace std;

//Input
vector<int> LoadFromFile();

//Serial QuickSort
void quicksort(vector<int>& values);
int compare(const void * a, const void * b);

//Hquicksort
// vector< vector<int> > GetSubArraysPos(int valuesSize, int qty);
void HyperQSort();		//Temp signature

//outputs
void Print(vector<int>& values);
void PrintChrono(int &nodes, size_t qty, double &duration);

int d;
int p;
int idle;
vector<int> values;
int arraySize;



//TODO:
//Use ScatterV
//Broadcast pivot
//Sort and echange
//Create a group with only a multiple of 8 processes
//Create smaller groups for dividing among the hypercube.

int main(int argc, char* argv[])
{		
	int mpiWorldRank;
	int mpiWorldSize;
	 
	int mpiRank; 
	int mpiSize;

	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiWorldRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);
	MPI_Comm MPI_COMM_HYPERCUBE;
	
	
	double startTime;
	
	d = log2(mpiWorldSize);		
	p = pow(2, d);					//Number of sorting processes
	idle = mpiWorldSize - p;		//number of idle processes 
	int toExclude[idle];			//idle processes to exclude
	int sendCounts[p];				//array of values sizes to be sent to processes
	
	int i = 0;
	while (i < idle)
	{
		toExclude[i] = mpiWorldSize - 1 - i;
		++i;
	}
	
	
	vector < vector<int> > pos;		//Input is stored in this vector.
	
	int transferDataSize;		
	
	//CREATING HYPERCUBE GROUP: Group of size of power of 2 -----------------
	// Obtain the group of processes in the world communicator
	MPI_Group world_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	// Remove all unnecessary ranks
	if (idle > 0)
	{
		MPI_Group newGroup;		
		MPI_Group_excl(world_group, idle, toExclude, &newGroup);

		// Create a new communicator
		MPI_Comm_create(MPI_COMM_WORLD, newGroup, &MPI_COMM_HYPERCUBE);

		//Abort any processor not part of the hypercube.	
		if (mpiWorldRank >= p)
		{
			cout << "aborting: " << mpiWorldRank <<endl;
			MPI_Finalize();
			return 0;
		}	
	}	
	else 
	{
		MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_HYPERCUBE);
	}
	MPI_Comm_rank(MPI_COMM_HYPERCUBE, &mpiRank);
	MPI_Comm_size(MPI_COMM_HYPERCUBE, &mpiSize);
	//REFACTOR: Organize code in different methods --------------------------

	
	if (mpiRank == 0)
	{
		//STEP0: Start Chrono
		startTime = MPI_Wtime();
	
		//STEP1: Read input
		values = LoadFromFile();
		arraySize = values.size();

		//STEP3: Select and distribute pivot
		int pivot0 = values[0];		//Arbitrary pivot
		//TODO: Broadcast Pivot! SORT BEFORE??
		
	}
	
	//STEP2: Scatter values
	MPI_Bcast(&arraySize, 1, MPI_INT, 0, MPI_COMM_HYPERCUBE);
	cout << "rank " << mpiRank << " - array size broadcasted: " << arraySize << endl;

    int nmin = arraySize / p;
	int remainingData = arraySize % p;
	int displs[p];
	int recvCount;

    int k = 0;
    for (i=0; i<p; i++)
    {
		sendCounts[i] = i < remainingData
			? nmin+1
			: nmin;
        displs[i] = k;
		k += sendCounts[i];
    }

	recvCount = sendCounts[mpiRank];
	int recvValues[recvCount];
	MPI_Scatterv(&values[0], sendCounts, displs, MPI_INT, recvValues , recvCount, MPI_INT, 0, MPI_COMM_HYPERCUBE);
	
	//TOREMOVE: checking values. That's all.
	int j = 0;
	while (j < recvCount)
	{
		cout << "rank " << mpiRank << " received: " << recvValues[j] << endl;
		++j;
	}
		
	//void HyperQSort();	

	MPI_Barrier(MPI_COMM_HYPERCUBE);
	if (mpiRank == 0)
	{
		//STEPX: Stop chrono and print results
		double endTime = MPI_Wtime();
		double duration = endTime - startTime;

		cout << values.size() << " numbers (quick)sorted" << endl;
		cout << "Duration: " << duration << " seconds" << endl;
		Print(values);
		PrintChrono(mpiSize, values.size(), duration);
	}

	MPI_Finalize();
	return 0;
}

void HyperQSort(int array[], int arraySize)
{
	//shall I create groups per dimensions? (1024 proc. is 10 groups only)
	//pivot should be sent by lowest id of the group. p0 to all, then to 3, then to 1. (don't send to himself)
	//first time to all, 
	//then to ones with same most significant bit, 
	//then to one with same most significant 2 bits, etc....
	
	
	int i = 0;
	// while (i < d)
	// {
	// 	if (mpiRank == 0)
	// 	{
	// 		//Arbitrary pivot value taken in the middle in case values are already ordered.
	// 		int pivot = array[arraySize / 2];
	//
	// 		//broadcast pivot
	// 	}
	//
	// 	vector<int> arrayL;
	// 	vector<int> arrayU;
	//
	// 	int j = 0;
	// 	while (j < arraySize)
	// 	{
	// 		if(array[j] <= pivot)
	// 		{
	// 			arrayL.push_back(array[j]);
	// 		}
	// 		else
	// 		{
	// 			arrayU.push_back(array[j]);
	// 		}
	// 		++j;
	// 	}
	//
	// 	int destId = mpiRank ^ pow(2, i);
	//
	// 	if (rank & (1<<i))	//checks if nth bit is set
	// 	{
	// 		//send lower
	// 		//receive from upper
	// 		//array = arrayU U receivedArray
	// 	}
	// 	else
	// 	{
	// 		//send upper
	// 		//etc...
	// 	}
	// 	++i;
	// }
	
}

//TOREMOVE: Not usefull anymore!!!!!!!!!!!!!!!!!!!!
// vector< vector<int> > GetSubArraysPos(int valuesSize, int processes)
// {
// 	int i = 0;
// 	arraySize = valuesSize / p;
// 	vector< vector<int> > pos;
// 	pos.resize( p , vector<int>(2) );
//
// 	while (i < p)
// 		{
// 			//each value of pos is the position of the beginning and end of each subarray
// 			pos[i][0] =  i * arraySize;
// 			pos[i][1] = i == (p - 1)	//Ending postion. Last p contains remaining values.
// 				? valuesSize - 1
// 				: pos[i][0] + arraySize - 1;
// 			cout << "pos[" << i << "]: " <<  pos[i][0] << " | " << pos[i][1] << endl;
// 			++i;
// 		}
//
// 	return pos;
// }


//Helpers, output and serial operations --------------------------

void quicksort(vector<int>& values)
{
	qsort(&values[0], values.size(), sizeof(int), compare);
}

vector<int> LoadFromFile()
{
	ifstream input("./input.txt", ios::in);

	int number;
	vector<int> numbers;
	
	if (input.is_open())
	{
		while(input >> number)
		{
			numbers.push_back(number);
			input.get();
		}
		input.close();
	}
	return numbers;
}

int compare(const void * a, const void * b)
{
	return (*(int*)a - *(int*)b);
}

void Print(vector<int>& values)
{
	std::ofstream myfile;
	myfile.open("./output.txt");

	int i = 0;
	while (i < values.size())
	{
		myfile << values[i] << endl;
		++i;
	}

	myfile.close();
}

void PrintChrono(int &nodes, size_t qty, double &duration)
{
    ofstream myfile;
	myfile.open("./HquicksortChrono.csv", ios::app);
	myfile << nodes << "\t" << qty << "\t" << duration << "\n";
	myfile.close();
}
