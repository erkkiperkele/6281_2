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
void GetSubArraysPos(int pos[], int valuesSize, int qty);
void HyperQSort();		//Temp signature

//outputs
void Print(vector<int>& values);
void PrintChrono(int &nodes, size_t qty, double &duration);

int d;
int p;
int idle;

int main(int argc, char* argv[])
{
	
	vector<int> values = LoadFromFile();
		
	int mpiRank; 
	int mpiSize;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	
	d = log2(mpiSize - 1);		
	p = pow(2, d);				//Number of sorting processes
	idle = mpiSize - 1 - p;		//number of idle processes 

	//positions is an array of array[2] containing beginning and end of sub array.
	int* positions[p];
	GetSubArraysPos(positions, values.size(), p);
	
	double startTime = MPI_Wtime();
	
	if (mpiRank == 0)
	{
		//TODO
		int rank = 0;
		while (rank < p)		//PERF: Use Boradcast?
		{
			int size = positions[rank-1][1] - positions[rank-1][0];
			MPI_Send(&values[rank-1][0], size, MPI_INT, id, 0, MPI_COMM_WORLD);
		}
		
		HyperQSort(&values[0], values.size());
	}

	if (mpiRank != 0)
	{
		//TODO
	}
		
	//void HyperQSort();	

	double endTime = MPI_Wtime();
	double duration = endTime - startTime;

	cout << values.size() << " numbers (quick)sorted" << endl;
	cout << "Duration: " << duration << " seconds" << endl;
	Print(values);
	PrintChrono(mpiSize, values.size(), duration);

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
	while (i < d)
	{
		if (mpiRank == 0)
		{
			//Arbitrary pivot value taken in the middle in case values are already ordered.
			int pivot = array[arraySize / 2];
			
			//broadcast pivot
		}
		
		vector<int> arrayL;
		vector<int> arrayU;
		
		int j = 0;
		while (j < arraySize)
		{
			if(array[j] <= pivot)
			{
				arrayL.push_back(array[j]);
			}
			else
			{
				arrayU.push_back(array[j]);
			}
			++j;
		}
		
		int destId = mpiRank ^ pow(2, i);
		
		if (rank & (1<<i))	//checks if nth bit is set
		{
			//send lower
			//receive from upper
			//array = arrayU U receivedArray
		}
		else 
		{
			//send upper
			//etc...
		}
		++i;
	}
	
}

void GetSubArraysPos(int* pos[], int valuesSize, int qty)
{
	int i = 0;
	int arraySize = valuesSize / (qty);
	
	while (i < qty)
		{
			//each value of pos is the position of the beginning and end of each subarray
			int arrayLimits[2];
			arrayLimits[0] = i * arraySize;			//Starting position
			arrayLimits[1] = i == (qty - 1)			//Ending postion. Last array contains any remaining values.
				? valuesSize - 1
				: arrayLimits[0] + arraySize
			positions[i] = arrayLimits;
			++i;
		}
}

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
