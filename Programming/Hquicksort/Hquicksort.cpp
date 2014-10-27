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
int subArraySize;


//TODO:
//Use ScatterV
//Broadcast pivot
//Sort and echange
//Create a group with only a multiple of 8 processes
//Create smaller groups for dividing among the hypercube.

int main(int argc, char* argv[])
{		
	int mpiRank; 
	int mpiSize;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	
	double startTime;
	
	d = log2(mpiSize);		
	p = pow(2, d);					//Number of sorting processes
	idle = mpiSize - p;				//number of idle processes 
	
	vector < vector<int> > pos;		//Input is stored in this vector.
	
	int transferDataSize;			
	
	
	//EXAMPLE ON HOW TO CREATE GROUPS: -----------------
	// Obtain the group of processes in the world communicator
	// MPI_Group world_group;
// 	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
//
// 	// Remove all unnecessary ranks
// 	MPI_Group new_group;
// 	int ranges[3] = { process_limit, size-1, 1 };
// 	MPI_Group_range_excl(world_group, 1, ranges, &new_group);
//
// 	// Create a new communicator
// 	MPI_Comm newworld;
// 	MPI_Comm_create(MPI_COMM_WORLD, new_group, &newworld);

		
		
		
	//REFACTOR: Organize code in different methods --------------------------

	
	if (mpiRank == 0)
	{
		//STEP0: Start Chrono
		startTime = MPI_Wtime();
	
		//STEP1: Read input
		values = LoadFromFile();
		arraySize = values.size();
		
		
		// //SEND remaining values to last process.
// 		int remainingVSize = values.size() % p;
// 		int remainingValues[remainingVSize];
// 		copy(values.end() - remainingVSize, values.end(), remainingValues);
// 		MPI_Send(&remainingValues, remainingVSize, MPI_INT, p-1, 0, MPI_COMM_WORLD);

		//STEP3: Select and distribute pivot
		int pivot0 = values[0];		//Arbitrary pivot
		//TODO: Broadcast Pivot! SORT BEFORE??
		
	}
	
	//PERF: work with an array from beginning instead of transforming vector into array
	int valuesArray[values.size()];
	if(mpiRank == 0)
	{
		copy(values.begin(), values.end(), valuesArray);
	}
	
	//STEP2: Scatter values
	//TODO: ScatterV in order to distribute all values!!
	MPI_Bcast(&arraySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	subArraySize = arraySize / p;
	int subValues[subArraySize];
	
	if (mpiRank !=0)
	{
		cout << "arraySize: " << arraySize << endl;
		cout << "subArraySize: " << subArraySize << endl;
	}


	//TODO: Need to define groups. At the moment idle processors receive values!!!!
	MPI_Scatter(&valuesArray, subArraySize, MPI_INT, subValues , subArraySize, MPI_INT, 0, MPI_COMM_WORLD);
	cout << "toSend 0 of rank " << mpiRank << " : " << subValues[0] << endl;	
	
	// //last process receives remaining values.
// 	if (mpiRank == p-1 && arraySize % p > 0)
// 	{
// 		int remainingValues[arraySize % p];
// 		MPI_Recv(&remainingValues, arraySize % p, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
// 		//need to merge arrays! or use scatterV
//
// 		int i = 0;
// 		while (i < arraySize % p)
// 		{
// 			cout << "remaining: " << remainingValues[i] << endl;
// 			++i;
// 		}
// 	}
	
	//TOREMOVE: checking values. That's all.
	int j = 0;
	while (j < 10)
	{
		cout << "toSend all values of rank " << mpiRank << " : " << subValues[j] << endl;		
		++j;			
	}
		
	//void HyperQSort();	

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
