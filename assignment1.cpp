#include <iostream>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <cmath>
#include "mt19937-64.c"

using namespace std;

#define ull unsigned long long

void printArray(vector<ull> arr)
{
    for (int i = 0; i < arr.size(); i++)
    {
        cout << arr.at(i) << " ";
    }
    cout << endl;
}

ull sum(vector<ull> arr, int end)
{
    ull sum = 0;
    for (int i = 0; i < end; i++)
    {
        sum += arr[i];
    }
    return sum;
}

vector<vector<ull> > getPartitions(vector<ull> arr, int size, int p)
{
    vector<vector<ull> > partitions(p);

    for (int i = 0; i < p; i++)
    {
        int start = i * size;
        int end = start + size - 1;

        partitions[i].resize(end - start + 1);

        for (int j = start, k = 0; j <= end; j++, k++)
        {
            partitions[i][k] = arr[j];
        }
    }

    return partitions;
}

vector<ull> getDividers(vector<vector<ull> > partitions, int div, int divSize, int size, int p)
{
    int nonDiv = size - div;
    vector<ull> dividers(divSize);

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

vector<ull> getBucketDel(vector<ull> dividers, int div, int divSize, int p)
{
    vector<ull> bucketDel(div);

    int notDel = divSize - div;
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

vector<vector<ull> > getSizeMatrix(vector<vector<ull> > partitions, vector<ull> bucketDel, int partitionSize, int div, int p)
{
    vector<vector<ull> > sizeMat(p);

    for (int i = 0; i < p; i++)
    {
        sizeMat[i].resize(p);
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

vector<ull> sampleSort(const vector<ull>& arr, int n, int p)
{
    int partitionSize = n / p;
    int div = p - 1;
    int divSize = p * div;

    double p1t1Start = omp_get_wtime();
    vector<vector<ull> > partitions = getPartitions(arr, partitionSize, p);
    for (int i = 0; i < p; i++)
    {
        sort(partitions[i].begin(), partitions[i].begin() + partitionSize);
        printArray(partitions[i]);
    }
    double p1t1End = omp_get_wtime();

    double p1t2Start = omp_get_wtime();
    vector<ull> dividers = getDividers(partitions, div, divSize, partitionSize, p);
    sort(dividers.begin(), dividers.begin() + divSize);
    printArray(dividers);
    double p1t2End = omp_get_wtime();

    double p1t3Start = omp_get_wtime();
    vector<ull> bucketDel = getBucketDel(dividers, div, divSize, p);
    vector<vector<ull> > sizeMat = getSizeMatrix(partitions, bucketDel, partitionSize, div, p);

    printArray(bucketDel);
    for (int i = 0; i < p; i++)
    {
        printArray(sizeMat[i]);
    }
    
    vector<ull> bucketSize(p, 0);
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < p; j++)
        {
            bucketSize[i] += sizeMat[j][i];
        }
    }
    printArray(bucketSize);

    vector<vector<ull> > bucket(p);
    for (int i = 0; i < p; i++)
    {
        bucket[i].resize(bucketSize[i]);
    }
    
    vector<vector<bool> > flags(p);
    for (int i = 0; i < p; i++)
    {
        flags[i].resize(partitionSize);
        for (int j = 0; j < partitionSize; j++)
        {
            flags[i][j] = false;
        }
    }

    
    vector<ull> bucketIndices(p, 0);
    printArray(bucketIndices);
    
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
                flags[i][j] = true;
                bucketIndices[div] += 1;
            }
        }
    }
    
    for (int i = 0; i < p; i++)
    {
        printArray(bucket[i]);
    }
    
    double p1t3End = omp_get_wtime();

    double p1t4Start = omp_get_wtime();
    for (int i = 0; i < p; i++)
    {
        sort(bucket[i].begin(), bucket[i].begin() + bucketSize[i]);
        printArray(bucket[i]);
    }
    double p1t4End = omp_get_wtime();

    vector<ull> sorted;
    for (int i = 0; i < p; i++) 
    {
        for (int j = 0; j < bucketSize[i]; j++)
        {
            sorted.push_back(bucket[i][j]);
        }
    }

    cout << "1. " << p1t1End - p1t1Start << "s" << endl;
    cout << "2. " << p1t2End - p1t2Start << "s" << endl;
    cout << "3. " << p1t3End - p1t3Start << "s" << endl;
    cout << "4. " << p1t4End - p1t4Start << "s" << endl;

    return sorted;
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
        omp_set_num_threads(p);

        int seed = 123;
        init_genrand64(seed);

        vector<ull> arr;
        for (int i = 0; i < n; i++)
        {
            arr.push_back(genrand64_int64() % 100);
        }

        printArray(arr);
    
        double startTime = omp_get_wtime();
        if (p == 1)
        {
            sort(arr.begin(), arr.begin() + n);
        }
        else
        {
            arr = sampleSort(arr, n, p);
        }
        double endTime = omp_get_wtime();

        printArray(arr);

        // verify that the array is sorted
        for (int i = 0; i < n - 1; i++)
        {
           if (arr[i] > arr[i + 1])
           {
               cout << "Error: Array not sorted" << endl;
               cout << "index [" << i << "] : " << arr[i] << " " << arr[i + 1] << endl;
               break;
           }
        }

        double time_used = endTime - startTime;
        cout << "\nTime used in total: " << time_used << "s" << endl;
    }
    else
    {
        cout << "Please choose p and n such that p is divisible by n" << endl;
    }

    return 0;
}