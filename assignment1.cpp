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

//int binarySearch(int arr[], int l, int r, int x)
//{
//    if (r >= l) 
//    {
//        int mid = l + (r - l) / 2;
//        if (arr[mid] == x)
//            return mid;
//
//        if (arr[mid] > x)
//            return binarySearch(arr, l, mid - 1, x);
//
//        return binarySearch(arr, mid + 1, r, x);
//    }
//    return -1;
//}

//void bucketSort(float* array, int size) {
//    vector<float> bucket[size];
//    for (int i = 0; i < size; i++) 
//    {
//        bucket[int(size * array[i])].push_back(array[i]);
//    }
//    for (int i = 0; i < size; i++)
//    {
//        sort(bucket[i].begin(), bucket[i].end());
//    }
//    int index = 0;
//    for (int i = 0; i < size; i++)
//    {
//        while (!bucket[i].empty()) {
//            array[index++] = *(bucket[i].begin());
//            bucket[i].erase(bucket[i].begin());
//        }
//    }
//}

void sampleSort(ull* arr, int n, int p) {
    // Divide the array into p partitions
    int partitionSize = n / p;
    ull** partitions = new ull * [p];

    for (int i = 0; i < p; i++)
    {
        int start = i * partitionSize;
        int end = start + partitionSize - 1;

        // allocate memory for partition
        partitions[i] = new ull[end - start + 1];

        // copy partition from original array to partition array
        for (int j = start, k = 0; j <= end; j++, k++)
        {
            partitions[i][k] = arr[j];
        }
    }

    cout << "Sorted partitions:" << endl;
    for (int i = 0; i < p; i++)
    {
        mergeSort(partitions[i], 0, partitionSize - 1);

        for (int j = 0; j < partitionSize; j++)
        {
            cout << partitions[i][j] << " ";
        }
        cout << endl;
    }

    int div = p - 1;
    int nonDiv = partitionSize - div;
    int dividersSize = p * div;
    ull* dividers = new ull[dividersSize];

    // cout << "\ndiv / non-div : " << div << " / " << nonDiv << endl;

    double step = (double)nonDiv / (double)(div + 1);
    int startIndex = (step < 1) ? 0 : (int)round(step);
    step = (step < 1) ? 1 : (int)round(step);

    int divIndex = 0;

    // cout << "start index, step : " << startIndex << ", " << step << endl << endl;

    for (int i = 0; i < p; i++) // for each bucket
    {
        int bucketDivCounter = 0;

        for (int j = startIndex; j < partitionSize; j += step + 1) // for each div
        {
            if (bucketDivCounter == div) break;
            dividers[divIndex++] = partitions[i][j];
            bucketDivCounter++;
        }

        if (bucketDivCounter != div)
        {
            int j = startIndex + 1;
            while (bucketDivCounter != div)
            {
                dividers[divIndex++] = partitions[i][j];

                j += step + 1;
                bucketDivCounter++;
            }
        }
    }

    cout << "\nSorted dividers: ";
    mergeSort(dividers, 0, dividersSize - 1);
    printArray(dividers, dividersSize);
    cout << endl;

    int bucketDelSize = p - 1;
    ull* bucketDel = new ull[bucketDelSize];

    int notDel = dividersSize - bucketDelSize;
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

    cout << "Pick every p'th divider: ";
    printArray(bucketDel, bucketDelSize);
    for (int i = 0; i < p; i++)
    {
        delete[] partitions[i];
    }

    delete[] partitions;
    delete[] dividers;
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

    // mergeSort(arr, 0, n - 1);

    double startTime = omp_get_wtime();
    sampleSort(arr, n, p);
    cout << endl;
    double endTime = omp_get_wtime();

    //printArray(arr, n);
    //cout << endl;

    // Verify that the array is sorted
    for (int i = 0; i < n - 1; i++)
    {
        if (arr[i] > arr[i + 1])
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