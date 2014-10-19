#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>

using namespace std;
using namespace std::chrono;

void quicksort(vector<int>& values);
vector<int> LoadFromFile();
int compare(const void * a, const void * b);
void Print(vector<int>& values);
void PrintChrono(double &duration, int qty);

int main(int argc, char* argv[])
{
	high_resolution_clock::time_point startTime = high_resolution_clock::now();

	vector<int> values = LoadFromFile();

	quicksort(values);

	high_resolution_clock::time_point endTime = high_resolution_clock::now();
	double time_span = duration_cast<duration<double> >(endTime - startTime).count();

	cout << values.size() << " numbers (quick)sorted" << endl;
	cout << "Duration: " << time_span << " seconds" << endl;
	Print(values);
	PrintChrono(time_span, values.size());

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

void PrintChrono(double &duration, int qty)
{
	std::ofstream myfile;
	myfile.open("./SquicksortChrono.csv", ios::app);
	myfile << "1," << qty << "," << duration << "\n";
	myfile.close();
}
