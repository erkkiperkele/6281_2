#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <climits>

using namespace std;

void Print(int qty, int* item);

int main(int argc, char* argv[])
{
	int i;
	int qty = atoi(argv[1]);
	int item[qty];
	cout << "number of elements generated: " << qty << endl;

	srandom(time(NULL));
	
	for (i=0; i<qty; i++)
	{
		item[i] = rand() % INT_MAX + 1;
	}
	
	Print(qty, item);
	
	return 0;
}

void Print(int qty, int* item)
{
	std::ofstream myfile;
	myfile.open("./input.txt");

	int i = 0;
	while (i < qty)
	{
		myfile << item[i] << endl;
		i++;
	}

	myfile.close();
}

