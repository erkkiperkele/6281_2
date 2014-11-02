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
int* HyperQSort(int currentd, int *currentValues);
void SplitList(vector<int> &lower, vector<int> &upper, int* valueSet, int valueSetSize, int pivot);

//outputs
void Print(vector<int>& values);
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
	
	//Abort any processor not part of the hypercube.	
	if (mpiWorldRank >= p)
	{
		cout << "aborting: " << mpiWorldRank <<endl;
		MPI_Finalize();
		return 0;
	}	
	MPI_Comm_rank(MPI_COMM_HYPERCUBE, &mpiRank);
	MPI_Comm_size(MPI_COMM_HYPERCUBE, &mpiSize);
	
	//STEP1: Read input
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
	
	// //TOREMOVE: checking values. That's all.
	// int j = 0;
	// while (j < recvCount)
	// {
	// 	cout << "rank " << mpiRank << " received: " << recvValues[j] << endl;
	// 	++j;
	// }
			
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

//TODO: Implement properly: far too ugly!
vector<int> GetGroup(int currentd)
{
	// vector<int> group;
	// int step = p / pow(2,currentd);
	//
	// int size = mpiRank / step;
	// int start = step * size;
	// int end = start + step;
	// while (start < end)
	// {
	// 	group.push_back(start);
	// 	++start;
	// }
	
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
	
	// return group;
}

void PushArrayToVector(vector<int> &toExtend, int* received, int receivedSize)
{
	int i = 0;
	// if (mpiRank == 0)
	// {
		// cout << "before size: " <<toExtend.size() << endl;
		// cout << "before receivedSize: " <<receivedSize << endl;
	// }
		
	while (i < receivedSize)
	{
		toExtend.push_back(received[i]);
		++i;
	}
	// if (mpiRank == 0)
	// {
		// cout << "before size: " <<toExtend.size() << endl;
		// cout << "before receivedSize: " <<receivedSize << endl;
	// }
}

int* HyperQSort(int currentd, int* toSort)
{
	// if (mpiRank == 0 || mpiRank == 4 || mpiRank == 2 || mpiRank == 6)
// 	{
// 		cout << currentd << " -- zzz receiving toSort address (" << mpiRank << "): " << &toSort[0] << endl;
// 		cout << currentd << "-- zzz receiving toSort[0] (" << mpiRank << "): " << toSort[0] << endl;
// 	}

	cout << currentd << " -- aaa rank: " << mpiRank << " toSortSize: " << currentValuesSize << endl;	
	//STEP3: Select pivot and Broadcast	
	vector<int> myGroup = GetGroup(currentd);
	int broadcaster = myGroup[0];
	if (mpiRank == 5)
		cout << currentd << " -- ooo rank5 - broadcaster: " << broadcaster << endl;

	int sendSize = (p / (pow(2,currentd)));		//Send pivot to n processes

	//Check every bit to determine the broadcaster. 
	//(Ex: 000 has none of its 3 bits set, and is broadcaster of dimension 0)
	//(Ex: 100 has none of its 2 last bits set and is broadcaster of dimension 3-1 = 2) etc...
	bool isBroadcaster = true;
	int i = 0;
	while ((i < d - currentd) && isBroadcaster)
	{
		if ((mpiRank & (1<<i)))
		{
			isBroadcaster = false;
		}
		++i;
	}
	
	//TODO: fallback broadcaster!!
	// if(isBroadcaster)
	if (broadcaster == mpiRank)
	{
		pivot = toSort[currentValuesSize-1];
		// pivot = toSort[1];
		// pivot = (toSort[0] + toSort[1]) / 2;
		// if (mpiRank == 0 || mpiRank == 4)
		// {
			// cout << currentd << " -- zzz toSort[0] pivot (" << mpiRank << ")= " << toSort[0] << endl;
			// cout << "currentValuesSize (" << mpiRank << ")= " << currentValuesSize << endl;
			// cout << "toSort address (" << mpiRank << ")= " << &toSort[0] << endl;
			
		// }

		int i = 0;
		while(i < sendSize)
		{
			if (myGroup[i] != mpiRank)
			{
				int receiver = myGroup[i];
				cout << currentd << " -- zzz rank: " << mpiRank << " send pivot " << pivot << " to: " << receiver << endl;
				cout << currentd << " -- zzz rank: " << mpiRank 
					<< " pivot - toSort[0] " << toSort[0] << " toSort[currentValuesSize-1] " << toSort[currentValuesSize-1]  
					<< " pivot " << pivot << endl;
				MPI_Send(&pivot, 1, MPI_INT, receiver, 0, MPI_COMM_HYPERCUBE);
			}
			++i;	
		}	
	}
	
	else
	{
		MPI_Recv(&pivot, 1, MPI_INT, broadcaster, 0, MPI_COMM_HYPERCUBE, MPI_STATUS_IGNORE);	
	}
	// cout << currentd << " -- zzz rank: " << mpiRank << " received pivot: " << pivot << endl;

	
	//STEP4: Split values in 2
	// vector<int> lower;
	// vector<int> upper;
	//TEST
	// if (currentd == 0)
	// {
	SplitList(lower, upper, toSort, currentValuesSize, pivot);		
	cout << currentd << " -- yyy rank: " << mpiRank << " lowerSize + upperSize: " << lower.size() + upper.size() << endl;	
	// }

	//VERIF ONLY:
	// int l = 0;
	// if (mpiRank == 0 && currentd == 0)
	// {
	// 	cout << "SEPARATED INTO: --------" << endl;
	// 	while (l < lower.size())
	// 	{
	// 		cout << lower[l] << endl;
	// 		++l;
	// 	}
	// 	cout << "upper than: " << pivot << endl;
	// 	int u = 0;
	// 	while (u < upper.size())
	// 	{
	// 		cout << upper[u] << endl;
	// 		++u;
	// 	}
	// }
	
	//STEP5: Echange data with neighbour
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
	
	int keepSize = !isUpperKeeper
		? lower.size()
		: upper.size();
	
	cout << currentd << " -- bbb rank: " << mpiRank << " send to: " << destId << " toSendSize " << toSendSize << endl;
	MPI_Isend(toSend, toSendSize, MPI_INT, destId, 0, MPI_COMM_HYPERCUBE, &sendRequest);
	MPI_Irecv(&received[0], maxSize, MPI_INT, destId, 0, MPI_COMM_HYPERCUBE, &recvRequest);
	cout << currentd << " -- ccc rank: " << mpiRank << " receive from destId: " << destId <<  endl;
	
	int receivedSize = 0;
	MPI_Wait(&recvRequest, &status);
	MPI_Get_count( &status, MPI_INT, &receivedSize );
	cout << currentd << " -- ddd rank: " << mpiRank << " got Count: " << receivedSize << endl;	
	cout << currentd << " -- ddd rank: " << mpiRank << " totalSize " << keepSize + receivedSize << endl;
	
	// //VERIF ONLY:
	// int m = 0;
	// if (mpiRank == 0 && currentd == 0)
	// {
	// 	cout << "RECEIVED values: --------" << endl;
	// 	while (m < receivedSize)
	// 	{
	// 		cout << received[m] << endl;
	// 		++m;
	// 	}
	// 	cout << "END OF RECEIVED values: --------" << endl;
	// }
	
	//STEP6: merge data kept with data received
	if (isUpperKeeper)
	{
		// if (mpiRank == 0) 
			// cout << "push uppersize" << upper.size() << endl;
		PushArrayToVector(upper, received, receivedSize);
	}
	else
	{
		// if (mpiRank == 0)
			// cout << "push lowersize" << lower.size() << endl;
		PushArrayToVector(lower, received, receivedSize);
		// if (mpiRank == 0)
			// cout << "--------------- 0 IS LOWER KEEPER!!! --------------" << endl;
	}
	
	//STEP7: Serial quickSort on last dimension
	if (currentd == d-1)
	{
		quicksort(isUpperKeeper ? upper	: lower);
	}
	
	// if (mpiRank < 3 && currentd == 0)
	if (currentd == 2)
	{
		cout << endl <<"MPI_get_count: " << receivedSize << endl;
		cout << "destId :" << destId << endl;
		cout << "MERGED values (" << mpiRank << ")--------" << endl;
		int m = 0;
		int size = isUpperKeeper
			? upper.size()
			: lower.size();

		while (m < size)
		{
			int result = isUpperKeeper
			? upper[m]
			: lower[m];
			cout << "yyy " << result << endl;
			++m;
		}
		cout << "END OF MERGED values: --------" << endl;
	}

	toSort = isUpperKeeper
			? &upper[0]
			: &lower[0];
	currentValuesSize = isUpperKeeper
			? upper.size()
			: lower.size();
	
	
	
	// if (mpiRank == 0 || mpiRank == 4 || mpiRank ==2 || mpiRank == 6)
	// {
	// 	cout << "DIMENSION " << currentd << " toSort currentValuesSize (" << mpiRank << ") " << currentValuesSize << endl;
	// 	cout << "DIMENSION " << currentd << " toSort sent (" << mpiRank << ") " << upper[0] << endl;
	// 	cout << "DIMENSION " << currentd << " toSort sent size (" << mpiRank << ") " << upper.size() << endl;
	// 	cout << "DIMENSION " << currentd << " toSort[0] (" << mpiRank << ") " << toSort[0] << endl;
	// 	cout << "DIMENSION " << currentd << " toSort[0] address (" << mpiRank << ") " << &toSort[0] << endl;
	// }
		

	//ALMOST THERE! Just need my groups to be right! At the moment, they don't share properly
	//exit and restart
	
	
	
	
	//shall I create groups per dimensions? (1024 proc. is 10 groups only)
	//pivot should be sent by lowest id of the group. p0 to all, then to 3, then to 1. (don't send to himself)
	//first time to all, 
	//then to ones with same most significant bit, 
	//then to one with same most significant 2 bits, etc....
	MPI_Barrier(MPI_COMM_HYPERCUBE);
	cout << currentd << " -- eee rank: " << mpiRank << " valueSize (end of loop): " << currentValuesSize << endl;
	
	
	
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
		// cout << "value compared to " << valueSet[i] << endl;
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