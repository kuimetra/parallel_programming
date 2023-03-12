#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <omp.h>
#include <cmath>
#include "mt19937-64.c"

using namespace std;

#define ull unsigned long long
ofstream file;

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
    ull s = 0;
    #pragma omp parallel for reduction(+:s)
    for (int i = 0; i < end; i++)
    {
        s += arr[i];
    }
    return s;
}

void fillArrayWithDividers(vector<ull> &del, vector<ull> &arr, int p, int nonDivAmount, int step)
{
    int rem = nonDivAmount % p;
    int k = -1;

    for (int j = 1; j < p; j++)
    {
        k += step;
        del.push_back(arr[k]);
        if (rem > 0)
        {
            k++;
            rem--;
        }
    }
}

vector<ull> sampleSort(vector<ull> &arr, int n, int p)
{
    int subseqSize = n / p;
    int div = p - 1;
    int divTotal = p * div;

    double p1t1Start = omp_get_wtime();

    vector<vector<ull> > subseqVector(p);

    #pragma omp parallel num_threads(p)
    {
        int tid = omp_get_thread_num();
        int start = tid * subseqSize;
        int end = start + subseqSize - 1;

        for (int j = start; j <= end; j++)
        {
            subseqVector[tid].push_back(arr[j]);
        }

        sort(subseqVector[tid].begin(), subseqVector[tid].begin() + subseqSize);
    }

    double p1t1End = omp_get_wtime();

    double p1t2Start = omp_get_wtime();

    vector<ull> divVector(divTotal);

    int nonDivOfBucket = subseqSize - div;
    int stepOfBucket = nonDivOfBucket / p + 1;

    #pragma omp parallel for
    for (int i = 0; i < p; i++)
    {
        vector<ull> divs;
        fillArrayWithDividers(divs, subseqVector[i], p, nonDivOfBucket, stepOfBucket);

        #pragma omp critical
        {
            copy(divs.begin(), divs.end(), divVector.begin() + i * divs.size());
        }
    }

    sort(divVector.begin(), divVector.begin() + divTotal);

    double p1t2End = omp_get_wtime();

    double p1t3Start = omp_get_wtime();

    vector<ull> bucketDel;

    int nonDiv = divTotal - div;
    int step = nonDiv / p + 1;
    fillArrayWithDividers(bucketDel, divVector, p, nonDiv, step);

    vector<vector<ull> > sizeMat(p);

    #pragma omp parallel for shared(sizeMat, subseqVector, bucketDel) num_threads(p)
    for (int i = 0; i < p; i++)
    {
        sizeMat[i].resize(p, 0);
        for (int k = 0; k < div; k++)
        {
            for (int j = 0; j < subseqSize; j++)
            {
                if (subseqVector[i][j] < bucketDel[k])
                {
                    #pragma omp atomic update
                    sizeMat[i][k] += 1;
                }
            }
        }
    }

    #pragma omp parallel for shared(sizeMat) num_threads(p)
    for (int i = 0; i < p; i++)
    {
        for (int j = 1; j < div; j++)
        {
            sizeMat[i][j] -= sum(sizeMat[i], j);
        }
        sizeMat[i][div] = subseqSize - sum(sizeMat[i], div);
    }
    
    vector<ull> bucketSize(p, 0);

    #pragma omp parallel for shared(sizeMat, bucketSize) num_threads(p)
    for (int i = 0; i < p; i++)
    {
        ull s = 0;
        #pragma omp parallel for reduction(+:s)
        for (int j = 0; j < p; j++)
        {
            s += sizeMat[j][i];
        }
        bucketSize[i] = s;
    }

    vector<vector<ull> > bucket(p);
    #pragma omp parallel for shared(bucketSize, bucket) num_threads(p)
    for (int i = 0; i < p; i++)
    {
        bucket[i].resize(bucketSize[i]);
    }

    vector<vector<bool> > flags(p);
    
    for (int i = 0; i < p; i++)
    {
        flags[i].resize(subseqSize);
        for (int j = 0; j < subseqSize; j++)
        {
            flags[i][j] = false;
        }
    }

    vector<ull> bucketIndices(p, 0);

    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < div; j++)
        {
            for (int k = 0; k < subseqSize; k++)
            {
                if (subseqVector[i][k] < bucketDel[j] && !flags[i][k])
                {
                    int bucketElemIndex = bucketIndices[j];
                    bucket[j][bucketElemIndex] = subseqVector[i][k];
                    flags[i][k] = true;
                    bucketIndices[j] += 1;
                }
            }
        }
    }
    

    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < subseqSize; j++)
        {
            if (subseqVector[i][j] >= bucketDel[div - 1])
            {
                int bucketElemIndex = bucketIndices[p - 1];
                bucket[div][bucketElemIndex] = subseqVector[i][j];
                flags[i][j] = true;
                bucketIndices[div] += 1;
            }
        }
    }

    double p1t3End = omp_get_wtime();

    double p1t4Start = omp_get_wtime();

    #pragma omp parallel num_threads(p)
    {
        int tid = omp_get_thread_num();
        sort(bucket[tid].begin(), bucket[tid].begin() + bucketSize[tid]);
    }

    double p1t4End = omp_get_wtime();

    vector<ull> sorted;

    for (int i = 0; i < p; i++)
    {
        sorted.insert(sorted.end(), bucket[i].begin(), bucket[i].end());
    }

    cout << p1t1End - p1t1Start << "\t";
    cout << p1t2End - p1t2Start << "\t";
    cout << p1t3End - p1t3Start << "\t";
    cout << p1t4End - p1t4Start << endl;

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

    if (p > 0 && n % p == 0)
    {
        omp_set_num_threads(p);

        int seed = 123;
        init_genrand64(seed);

        vector<ull> arr;
        for (int i = 0; i < n; i++)
        {
            arr.push_back(genrand64_int64() % 100);
        }

        // printArray(arr);

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

        for (int i = 0; i < n - 1; i++)
        {
            if (arr[i] > arr[i + 1])
            {
                cout << "Error: Array not sorted (arr[" << i << "] > arr[" << i + 1 << "] : " << arr[i] << " > " << arr[i + 1] << ")" << endl;
                break;
            }
        }

        // printArray(arr);

        double time_used = endTime - startTime;
        cout << time_used << endl;
    }
    else
    {
        cout << "Please choose p and n such that p is divisible by n" << endl;
    }
    
    return 0;
}