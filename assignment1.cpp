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

        partitions[i] = new ull[end - start + 1];

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

    for (int i = 0; i < p; i++)
    {
        int r = nonDiv % p;
        int q = -1;

        for (int j = 1; j < p; j++)
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

    double p1t1Start = omp_get_wtime();
    ull** partitions = getPartitions(arr, partitionSize, p);
    for (int i = 0; i < p; i++)
    {
        mergeSort(partitions[i], 0, partitionSize - 1);
        // printArray(partitions[i], partitionSize);
    }
    double p1t1End = omp_get_wtime();

    double p1t2Start = omp_get_wtime();
    ull* dividers = getDividers(partitions, div, dividersSize, partitionSize, p);
    mergeSort(dividers, 0, dividersSize - 1);
    // printArray(dividers, dividersSize);
    double p1t2End = omp_get_wtime();

    double p1t3Start = omp_get_wtime();
    ull* bucketDel = getBucketDel(dividers, div, dividersSize, p);
    ull** sizeMat = getSizeMatrix(partitions, bucketDel, partitionSize, div, p);

    // printArray(bucketDel, div);
    /* for (int i = 0; i < p; i++)
    {
        printArray(sizeMat[i], p);
    }*/

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
    // printArray(bucketSize, p);

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

    /* for (int i = 0; i < p; i++)
    {
        printArray(bucket[i], bucketSize[i]);
    } */

    double p1t3End = omp_get_wtime();

    double p1t4Start = omp_get_wtime();
    for (int i = 0; i < p; i++)
    {
        mergeSort(bucket[i], 0, bucketSize[i]);
        // printArray(bucket[i], bucketSize[i]);
    }
    double p1t4End = omp_get_wtime();

    ull* sortedArray = new ull[n];
    int ind = 0;
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < bucketSize[i]; j++)
        {
            sortedArray[ind++] = bucket[i][j];
        }
    }

    cout << "1. " << p1t1End - p1t1Start << "s" << endl;
    cout << "2. " << p1t2End - p1t2Start << "s" << endl;
    cout << "3. " << p1t3End - p1t3Start << "s" << endl;
    cout << "4. " << p1t4End - p1t4Start << "s" << endl;

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

    cout << "Enter n: ";
    cin >> n;
    cout << "Enter p: ";
    cin >> p;
    cout << endl;

    if (n % p == 0)
    {
        // omp_set_num_threads(p);

        int seed = 123;
        init_genrand64(seed);

        ull* arr = new ull[n];
        for (int i = 0; i < n; i++)
        {
            arr[i] = genrand64_int64() % 100;
        }

        // printArray(arr, n);

        double startTime = omp_get_wtime();
        ull* sortedArr = sampleSort(arr, n, p);
        double endTime = omp_get_wtime();

        // printArray(sortedArr, n);

        // verify that the array is sorted
        for (int i = 0; i < n - 1; i++)
        {
            if (sortedArr[i] > sortedArr[i + 1])
            {
                cout << "Error: Array not sorted" << endl;
                break;
            }
        }

        double time_used = endTime - startTime;
        cout << "\nTime used in total: " << time_used << "s" << endl;

        delete[] arr;
    }
    else
    {
        cout << "Please choose p and n such that p is divisible by n" << endl;
    }

    return 0;
}