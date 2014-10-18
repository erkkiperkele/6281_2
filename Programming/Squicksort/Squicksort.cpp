#include <stdio.h>
#include <stdlib.h>

int compare(const void * a, const void * b)
{
	return (*(int*)a - *(int*)b);
}

int main()
{
	quicksort(values, valuesSize);
	return 0;
}

void quicksort(int values[], valuesSize)
{
	int n;
	qsort(values, valuesSize, sizeof(int), compare);
	for (n=0; n<valuesSize; n++)
		printf("%d", values[n]);
}
