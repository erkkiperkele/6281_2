#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <cmath>

using namespace std;

//Input
vector<int> LoadFromFile();

//Serial QuickSort
void quicksort(vector<int>& values);
int compare(const void * a, const void * b);

//Hquicksort
int* HyperQSort(int currentd, int *currentValues);
void SplitList(vector<int> &lower, vector<int> &upper, int* valueSet, int valueSetSize, int pivot);

//outputs
void Print(int* values, int valuesSize);
void PrintChrono(int &nodes, size_t qty, double &duration);

//communicator(s)
void HypercubeInit(int toExclude[]);

int mpiWorldRank;
int mpiWorldSize;
MPI_Comm MPI_COMM_HYPERCUBE;
 
int mpiRank; 
int mpiSize;

int d;
int p;
int idle;
vector<int> values;
int arraySize;

int pivot;
int* currentValues;
int currentValuesSize;

int* toSend;
vector<int> lower;
vector<int> upper;

//TODO:
//remove global variable and pass it as parameters instead

int main(int argc, char* argv[])
{		
	//STEP0: initialize MPI, global variables and start chrono
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiWorldRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);
	double startTime = MPI_Wtime();
	
	//Processes infos
	d = log2(mpiWorldSize);			//Hypercube dimension
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
	
	HypercubeInit(toExclude);
	
	//Abort any processor not part of the hypercube (better way to do it).	
	if (mpiWorldRank >= p)
	{
		cout << "aborting: " << mpiWorldRank <<endl;
		MPI_Finalize();
		return 0;
	}	
	MPI_Comm_rank(MPI_COMM_HYPERCUBE, &mpiRank);
	MPI_Comm_size(MPI_COMM_HYPERCUBE, &mpiSize);
	
	//STEP1: Master read input
	if (mpiRank == 0)
	{
		values = LoadFromFile();
		arraySize = values.size();
	}
	
	//STEP2: Scatter values
	MPI_Bcast(&arraySize, 1, MPI_INT, 0, MPI_COMM_HYPERCUBE);

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
			
	//STEP3: Compare and exchange
	int currentd = 0;
	currentValues = &recvValues[0];
	currentValuesSize = recvCount;
	
	toSend = currentValues;
	while (currentd < d)
	{
		toSend = HyperQSort(currentd, &toSend[0]);
		++currentd;
	}
	
	//STEP4: Gather results
	int resultSizes[p];
	int finalResults[arraySize];
	int receive_displacements[p];
	
	//Get each process result size
	MPI_Gather(&currentValuesSize, 1, MPI_INT, resultSizes , 1, MPI_INT, 0, MPI_COMM_HYPERCUBE);
	
	i = 0;
	k = 0;
	while (i < p)
	{
		k += i == 0
			? 0
			: resultSizes[i-1];
		receive_displacements[i] = k;
		++i;
	}

	//Concatenate data in proper order since each process result size is known.
	MPI_Gatherv(&toSend[0], currentValuesSize, MPI_INT, &finalResults[0], &resultSizes[0], receive_displacements, MPI_INT, 0, MPI_COMM_HYPERCUBE);

	if (mpiRank == 0)
	{	
		//STEP5: Stop chrono and print results
		double endTime = MPI_Wtime();
		double duration = endTime - startTime;

		cout << values.size() << " numbers (quick)sorted" << endl;
		cout << "Duration: " << duration << " seconds" << endl;
		Print(finalResults, currentValuesSize);
		PrintChrono(mpiSize, values.size(), duration);
	}

	MPI_Finalize();
	return 0;
}

//TODO: to implement properly...
vector<int> GetGroup(int currentd)
{

	if (currentd == 0)
	{
		static const int arr[] = {0,1,2,3,4,5,6,7};
		 vector<int> group (arr, arr + sizeof(arr) / sizeof(arr[0]) );
		 return group;
	}
	if (currentd == 1)
	{
		static const int arr1[] = {0,2,4,6};
		static const int arr2[] = {1,3,5,7};
		
		if (mpiRank % 2 == 0)
		{
			vector<int> group (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
			return group;
		}
		else 
		{
			vector<int> group (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
			return group;
		}
	}
	if (currentd == 2)
	{
		static const int arr1[] = {0,4};
		static const int arr2[] = {1,5};
		static const int arr3[] = {2,6};
		static const int arr4[] = {3,7};
		
		if (mpiRank == 0 || mpiRank == 4)
		{
			vector<int> group (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
			return group;
		}
		if (mpiRank == 1 || mpiRank == 5)
		{
			vector<int> group (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
			return group;
		}
		if (mpiRank == 2 || mpiRank == 6)
		{
			vector<int> group (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]) );
			return group;
		}
		if (mpiRank == 3 || mpiRank == 7)
		{
			vector<int> group (arr4, arr4 + sizeof(arr4) / sizeof(arr4[0]) );
			return group;
		}
	}
    vector<int> group;
	 return group;
}

void PushArrayToVector(vector<int> &toExtend, int* received, int receivedSize)
{
	int i = 0;
	while (i < receivedSize)
	{
		toExtend.push_back(received[i]);
		++i;
	}
}

int* HyperQSort(int currentd, int* toSort)
{
	//STEP A: Select pivot and Broadcast	
	vector<int> myGroup = GetGroup(currentd);
	int broadcaster = myGroup[0];
	int sendSize = (p / (pow(2,currentd)));		//Send pivot to n processes
	
	//TODO: fallback broadcaster!!
	if (broadcaster == mpiRank)
	{
		pivot = toSort[currentValuesSize-1];

		int i = 0;
		while(i < sendSize)
		{
			if (myGroup[i] != mpiRank)
			{
				int receiver = myGroup[i];
				MPI_Send(&pivot, 1, MPI_INT, receiver, 0, MPI_COMM_HYPERCUBE);
			}
			++i;	
		}	
	}
	
	else
	{
		MPI_Recv(&pivot, 1, MPI_INT, broadcaster, 0, MPI_COMM_HYPERCUBE, MPI_STATUS_IGNORE);	
	}

	SplitList(lower, upper, toSort, currentValuesSize, pivot);		
	
	//STEP B: Echange data with neighbour
	int power = pow(2, currentd);
	int destId = mpiRank ^ power;
	int maxSize = arraySize; // / p +1;
	int received[maxSize];
	
	
	MPI_Request sendRequest;
	MPI_Request recvRequest;
	MPI_Status status;

	bool isUpperKeeper = mpiRank & (1<<currentd);
	int* toSend = isUpperKeeper
		? &lower[0]
		: &upper[0];
	
	int toSendSize = isUpperKeeper
		? lower.size()
		: upper.size();
	
	MPI_Isend(toSend, toSendSize, MPI_INT, destId, 0, MPI_COMM_HYPERCUBE, &sendRequest);
	MPI_Irecv(&received[0], maxSize, MPI_INT, destId, 0, MPI_COMM_HYPERCUBE, &recvRequest);
	
	int receivedSize = 0;
	MPI_Wait(&recvRequest, &status);
	MPI_Get_count( &status, MPI_INT, &receivedSize );
	
	//STEPC: merge data kept with data received
	if (isUpperKeeper)
	{
		PushArrayToVector(upper, received, receivedSize);
	}
	else
	{
		PushArrayToVector(lower, received, receivedSize);
	}
	
	//STEPD: Serial quickSort on last dimension
	if (currentd == d-1)
	{
		quicksort(isUpperKeeper ? upper	: lower);
	}

	toSort = isUpperKeeper
			? &upper[0]
			: &lower[0];
	currentValuesSize = isUpperKeeper
			? upper.size()
			: lower.size();
	
	//Make sure every process is done before starting on next dimension
	MPI_Barrier(MPI_COMM_HYPERCUBE);
	return toSort;
}

void HypercubeInit(int toExclude[])
{
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
	}	
	else 
	{
		MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_HYPERCUBE);
	}
	
	//PERF: Should create all hypercube's subgroups
}

void SplitList(vector<int> &lower, vector<int> &upper, int* valueSet, int valueSetSize, int pivot)
{
	int i = 0;
	lower.clear();
	upper.clear();
	while (i < valueSetSize)
	{
		bool isLower = valueSet[i] <= pivot;
		if (isLower)
		{
			lower.push_back(valueSet[i]);
		}
		else
		{
			upper.push_back(valueSet[i]);			
		}
		++i;
	}
}

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

void Print(int* values, int valuesSize)
{
	std::ofstream myfile;
	myfile.open("./output.txt");

	int i = 0;
	while (i < valuesSize)
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