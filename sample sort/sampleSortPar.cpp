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
    // allocation of memory and initialization

    int subseqSize = n / p;
    int div = p - 1;
    int divTotal = p * div;

    int nonDivOfBucket = subseqSize - div;
    int stepOfBucket = nonDivOfBucket / p + 1;
    int nonDiv = divTotal - div;
    int step = nonDiv / p + 1;

    vector<vector<ull> > subseqVector(p);
    vector<ull> divVector;
    vector<ull> bucketDel;
    vector<vector<ull> > sizeMat(p);
    vector<ull> bucketSize(p, 0);
    vector<vector<ull> > bucket(p);
    vector<vector<bool> > flags(p);
    vector<ull> bucketIndices(p, 0);
    vector<ull> sorted;

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int start = tid * subseqSize;
        int end = start + subseqSize - 1;

        for (int j = start; j <= end; j++)
        {
            subseqVector[tid].push_back(arr[j]);
        }
    }

    double startTime = omp_get_wtime(), t1Start = omp_get_wtime();

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        sort(subseqVector[tid].begin(), subseqVector[tid].end());
    }

    double t1End = omp_get_wtime();

    double t2Start = omp_get_wtime();

    for (int i = 0; i < p; i++)
    {
        fillArrayWithDividers(divVector, subseqVector[i], p, nonDivOfBucket, stepOfBucket);
    }

    sort(divVector.begin(), divVector.begin() + divTotal);

    double t2End = omp_get_wtime();

    double t3Start = omp_get_wtime();

    fillArrayWithDividers(bucketDel, divVector, p, nonDiv, step);

    #pragma omp parallel for shared(sizeMat, subseqVector, bucketDel)
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

    #pragma omp parallel for shared(sizeMat)
    for (int i = 0; i < p; i++)
    {
        for (int j = 1; j < div; j++)
        {
            sizeMat[i][j] -= sum(sizeMat[i], j);
        }
        sizeMat[i][div] = subseqSize - sum(sizeMat[i], div);
    }

    #pragma omp parallel for shared(sizeMat, bucketSize)
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

    for (int i = 0; i < p; i++)
    {
        bucket[i].resize(bucketSize[i]);
    }

    for (int i = 0; i < p; i++)
    {
        flags[i].resize(subseqSize);
        for (int j = 0; j < subseqSize; j++)
        {
            flags[i][j] = false;
        }
    }
    
    #pragma omp parallel for
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < div; j++)
        {
            for (int k = 0; k < subseqSize; k++)
            {
                if (subseqVector[i][k] < bucketDel[j] && !flags[i][k])
                {
                    int bucketElemIndex;
                    #pragma omp atomic capture
                    {
                        bucketElemIndex = bucketIndices[j];
                        bucketIndices[j] += 1;
                    }
                    
                    bucket[j][bucketElemIndex] = subseqVector[i][k];
                    flags[i][k] = true;
                   
                }
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < subseqSize; j++)
        {
            if (subseqVector[i][j] >= bucketDel[div - 1])
            {
                int bucketElemIndex;
                #pragma omp atomic capture
                {
                    bucketElemIndex = bucketIndices[p - 1];
                    bucketIndices[div] += 1;
                }
                
                bucket[div][bucketElemIndex] = subseqVector[i][j];
                flags[i][j] = true;
            }
        }
    }

    double t3End = omp_get_wtime();

    double t4Start = omp_get_wtime();

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        sort(bucket[tid].begin(), bucket[tid].end());
    }

    double t4End = omp_get_wtime();

    for (int i = 0; i < p; i++)
    {
        sorted.insert(sorted.end(), bucket[i].begin(), bucket[i].end());
    }

    double endTime = omp_get_wtime();

    file << t1End - t1Start << ",";
    file << t2End - t2Start << ",";
    file << t3End - t3Start << ",";
    file << t4End - t4Start << ",";
    file << endTime - startTime << endl;

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

    file.open("log_p.csv");
    
    if (p > 0 && n % p == 0)
    {
        file << n << "," << p << ",";
        
        omp_set_num_threads(p);

        int seed = 123;
        init_genrand64(seed);

        vector<ull> arr;
        for (int i = 0; i < n; i++)
        {
            arr.push_back(genrand64_int64() % 100000);
        }

        // printArray(arr);

        if (p == 1)
        {
            double startTime = omp_get_wtime();
            sort(arr.begin(), arr.end());
            double endTime = omp_get_wtime();
            file << ",,,," << endTime - startTime << endl;
        }
        else
        {
            arr = sampleSort(arr, n, p);
        }

        for (int i = 0; i < n - 1; i++)
        {
            if (arr[i] > arr[i + 1])
            {
                cout << "Error: Array not sorted (arr[" << i << "] > arr[" << i + 1 << "] : " << arr[i] << " > " << arr[i + 1] << ")" << endl;
                break;
            }
        }

        // printArray(arr);
    }
    else
    {
        cout << "Please choose p and n such that p is divisible by n" << endl;
    }

    file.close();

    return 0;
}