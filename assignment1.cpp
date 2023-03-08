#include <iostream>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <cmath>
#include "mt19937-64.c"

using namespace std;

#define ull unsigned long long

void printArray(ull* arr, int size)
{
    for (int i = 0; i < size; i++)
    {
        cout << arr[i] << " ";
    }
    cout << endl;
}

void merge(ull* arr, int start, int middle, int end)
{
    int arr1_size = middle - start + 1;
    int arr2_size = end - middle;

    ull* arr1 = new ull[arr1_size];
    ull* arr2 = new ull[arr2_size];

    for (int i = 0; i < arr1_size; i++)
    {
        arr1[i] = arr[start + i];
    }

    for (int i = 0; i < arr2_size; i++)
    {
        arr2[i] = arr[middle + 1 + i];
    }

    int i = 0, j = 0, k = start;

    while (i < arr1_size && j < arr2_size)
    {
        if (arr1[i] <= arr2[j])
        {
            arr[k] = arr1[i];
            i++;
        }
        else
        {
            arr[k] = arr2[j];
            j++;
        }
        k++;
    }

    while (i < arr1_size)
    {
        arr[k] = arr1[i];
        i++;
        k++;
    }

    while (j < arr2_size)
    {
        arr[k] = arr2[j];
        j++;
        k++;
    }

    delete[] arr1;
    delete[] arr2;
}

void mergeSort(ull* arr, int start, int end)
{
    if (start < end)
    {
        int middle = start + (end - start) / 2;

        mergeSort(arr, start, middle);
        mergeSort(arr, middle + 1, end);

        merge(arr, start, middle, end);
    }
}

ull sum(ull* arr, int size)
{
    ull sum = 0;
    for (int i = 0; i < size; i++)
    {
        sum += arr[i];
    }
    return sum;
}

ull** getPartitions(ull* arr, int size, int p)
{
    ull** partitions = new ull * [p];

    for (int i = 0; i < p; i++)
    {
        int start = i * size;
        int end = start + size - 1;

        // allocate memory for partition
        partitions[i] = new ull[end - start + 1];

        // copy partition from original array to partition array
        for (int j = start, k = 0; j <= end; j++, k++)
        {
            partitions[i][k] = arr[j];
        }
    }

    return partitions;
}

ull* getDividers(ull** partitions, int div, int dividersSize, int size, int p)
{
    int nonDiv = size - div;
    ull* dividers = new ull[dividersSize];

    int step = nonDiv / p + 1;

    int divIndex = 0;

    for (int i = 0; i < p; i++) // for each bucket
    {
        int r = nonDiv % p;
        int q = -1;

        for (int j = 1; j < p; j++) // for each div
        {
            q += step;
            dividers[divIndex++] = partitions[i][q];
            if (r > 0)
            {
                q++;
                r--;
            }
        }
    }

    return dividers;
}

ull* getBucketDel(ull* dividers, int div, int dividersSize, int p)
{
    ull* bucketDel = new ull[div];

    int notDel = dividersSize - div;
    int delStep = notDel / p + 1;
    int rem = notDel % p;
    int k = -1;

    int bucketDelIndex = 0;

    for (int i = 1; i < p; i++)
    {
        k += delStep;
        bucketDel[bucketDelIndex++] = dividers[k];
        if (rem > 0)
        {
            k++;
            rem--;
        }
    }

    return bucketDel;
}

ull** getSizeMatrix(ull** partitions, ull* bucketDel, int partitionSize, int div, int p)
{
    ull** sizeMat = new ull * [p];

    for (int i = 0; i < p; i++)
    {
        sizeMat[i] = new ull[p];
        for (int k = 0; k < div; k++)
        {
            sizeMat[i][k] = 0;
            for (int j = 0; j < partitionSize; j++)
            {
                if (partitions[i][j] < bucketDel[k])
                {
                    sizeMat[i][k] += 1;
                }
            }
        }
    }

    for (int i = 0; i < p; i++)
    {
        for (int j = 1; j < div; j++)
        {
            sizeMat[i][j] -= sum(sizeMat[i], j);
        }
        sizeMat[i][div] = partitionSize - sum(sizeMat[i], div);
    }

    return sizeMat;
}

ull* sampleSort(ull* arr, int n, int p) {
    // Divide the array into p partitions
    int partitionSize = n / p;
    int div = p - 1;
    int dividersSize = p * div;

    ull** partitions = getPartitions(arr, partitionSize, p);
    for (int i = 0; i < p; i++)
    {
        mergeSort(partitions[i], 0, partitionSize - 1);
        printArray(partitions[i], partitionSize);
    }
    cout << endl;

    ull* dividers = getDividers(partitions, div, dividersSize, partitionSize, p);
    mergeSort(dividers, 0, dividersSize - 1);
    for (int i = 0; i < dividersSize; i++)
    {
        cout << dividers[i] << " ";
    }
    cout << endl << endl;

    ull* bucketDel = getBucketDel(dividers, div, dividersSize, p);
    for (int i = 0; i < div; i++)
    {
        cout << bucketDel[i] << " ";
    }
    cout << endl << endl;

    ull** sizeMat = getSizeMatrix(partitions, bucketDel, partitionSize, div, p);
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < p; j++)
        {
            cout << sizeMat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    ull* bucketSize = new ull[p];
    for (int i = 0; i < p; i++)
    {
        int colSum = 0;
        for (int j = 0; j < p; j++)
        {
            colSum += sizeMat[j][i];
        }
        bucketSize[i] = colSum;
    }
    for (int i = 0; i < p; i++)
    {
        cout << bucketSize[i] << " ";
    }
    cout << endl << endl;

    ull** bucket = new ull * [p];
    for (int i = 0; i < p; i++)
    {
        bucket[i] = new ull[bucketSize[i]];
    }

    bool** flags = new bool* [p];
    for (int i = 0; i < p; i++)
    {
        flags[i] = new bool[partitionSize];
        for (int j = 0; j < partitionSize; j++)
        {
            flags[i][j] = false;
        }
    }

    ull* bucketIndices = new ull[p];
    for (int i = 0; i < p; i++)
    {
        bucketIndices[i] = 0;
    }

    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < div; j++)
        {
            for (int k = 0; k < partitionSize; k++)
            {
                if (partitions[i][k] < bucketDel[j] && !flags[i][k])
                {
                    int bucketElemIndex = bucketIndices[j];
                    bucket[j][bucketElemIndex] = partitions[i][k];
                    flags[i][k] = true;
                    bucketIndices[j] += 1;
                }
            }
        }
    }

    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < partitionSize; j++)
        {
            if (partitions[i][j] >= bucketDel[div - 1])
            {
                int bucketElemIndex = bucketIndices[p - 1];
                bucket[div][bucketElemIndex] = partitions[i][j];
                bucketIndices[div] += 1;
            }
        }
    }

    for (int i = 0; i < p; i++)
    {
        cout << "Bucket [" << i + 1 << "] : ";
        for (int j = 0; j < bucketSize[i]; j++)
        {
            cout << bucket[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    for (int i = 0; i < p; i++)
    {
        mergeSort(bucket[i], 0, bucketSize[i]);
    }

    ull* sortedArray = new ull[n];
    int ind = 0;
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < bucketSize[i]; j++)
        {
            sortedArray[ind++] = bucket[i][j];
        }
    }

    for (int i = 0; i < p; i++)
    {
        delete[] partitions[i];
        delete[] sizeMat[i];
        delete[] bucket[i];
        delete[] flags[i];
    }

    delete[] partitions;
    delete[] dividers;
    delete[] bucketDel;
    delete[] sizeMat;
    delete[] bucketSize;
    delete[] bucket;
    delete[] flags;
    delete[] bucketIndices;

    return sortedArray;
}

int main()
{
    int n, p;

    /*cout << "Enter n: ";
    cin >> n;
    cout << "Enter p: ";
    cin >> p;*/

    n = 24;
    p = 3;

    int seed = 123;
    init_genrand64(seed);

    ull* arr = new ull[n];
    for (int i = 0; i < n; i++)
    {
        arr[i] = genrand64_int64() % 100;
    }

    printArray(arr, n);
    cout << endl;

    double startTime = omp_get_wtime();
    ull* sortedArr = sampleSort(arr, n, p);
    double endTime = omp_get_wtime();

    cout << "Sorted array: ";
    printArray(sortedArr, n);
    cout << endl;

    // Verify that the array is sorted
    for (int i = 0; i < n - 1; i++)
    {
        if (sortedArr[i] > sortedArr[i + 1])
        {
            cout << "Error: Array not sorted" << endl;
            break;
        }
    }

    double time_used = endTime - startTime;
    cout << "Time used: " << time_used << " seconds" << endl;

    delete[] arr;
    return 0;
}