#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <mpi.h>

using namespace std;

void quicksort(vector<int>& values);
vector<int> LoadFromFile();
int compare(const void * a, const void * b);
void Print(vector<int>& values);
void PrintChrono(int &nodes, size_t qty, double &duration);

int main(int argc, char* argv[])
{
	
	vector<int> values = LoadFromFile();
		
	int mpiRank; 
	int mpiSize;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	
	double startTime = MPI_Wtime();
	
	if (mpiRank == 0)
	{
		//TODO
	}

	if (mpiRank != 0)
	{
		//TODO
	}
	

	double endTime = MPI_Wtime();
	double duration = endTime - startTime;

	cout << values.size() << " numbers (quick)sorted" << endl;
	cout << "Duration: " << duration << " seconds" << endl;
	Print(values);
	PrintChrono(mpiSize, values.size(), duration);

	return 0;
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
