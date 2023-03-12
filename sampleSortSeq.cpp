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

    for (int i = 0; i < p; i++)
    {
        int start = i * subseqSize;
        int end = start + subseqSize - 1;

        for (int j = start; j <= end; j++)
        {
            subseqVector[i].push_back(arr[j]);
        }

        sort(subseqVector[i].begin(), subseqVector[i].begin() + subseqSize);
    }

    double p1t1End = omp_get_wtime();

    double p1t2Start = omp_get_wtime();

    vector<ull> divVector;

    int nonDivOfBucket = subseqSize - div;
    int stepOfBucket = nonDivOfBucket / p + 1;

    for (int i = 0; i < p; i++)
    {
        fillArrayWithDividers(divVector, subseqVector[i], p, nonDivOfBucket, stepOfBucket);
    }

    sort(divVector.begin(), divVector.begin() + divTotal);

    double p1t2End = omp_get_wtime();

    double p1t3Start = omp_get_wtime();

    vector<ull> bucketDel;

    int nonDiv = divTotal - div;
    int step = nonDiv / p + 1;
    fillArrayWithDividers(bucketDel, divVector, p, nonDiv, step);

    vector<vector<ull> > sizeMat(p);

    for (int i = 0; i < p; i++)
    {
        sizeMat[i].resize(p, 0);
        for (int k = 0; k < div; k++)
        {
            for (int j = 0; j < subseqSize; j++)
            {
                if (subseqVector[i][j] < bucketDel[k])
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
        sizeMat[i][div] = subseqSize - sum(sizeMat[i], div);
    }
    
    vector<ull> bucketSize(p, 0);

    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < p; j++)
        {
            bucketSize[i] += sizeMat[j][i];
        }
    }

    vector<vector<ull> > bucket(p);
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

    for (int i = 0; i < p; i++)
    {
        sort(bucket[i].begin(), bucket[i].begin() + bucketSize[i]);
    }

    double p1t4End = omp_get_wtime();

    vector<ull> sorted;

    for (int i = 0; i < p; i++)
    {
        sorted.insert(sorted.end(), bucket[i].begin(), bucket[i].end());
    }

    file << p1t1End - p1t1Start << ",";
    file << p1t2End - p1t2Start << ",";
    file << p1t3End - p1t3Start << ",";
    file << p1t4End - p1t4Start << ",";

    return sorted;
}

int main()
{
    int n, p;

    // cout << "Enter n: ";
    // cin >> n;
    // cout << "Enter p: ";
    // cin >> p;
    // cout << endl;

    file.open("log.csv");
    
    int nArr[3] = {10000000, 20000000, 50000000};
    int pArr[8] = {1, 2, 5, 10, 20, 40, 50, 80};

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0 ; j < 8; j++)
        {
            int n = nArr[i], p = pArr[j];

            if (p > 0 && n % p == 0)
            {
                file << n << "," << p << ",";

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
                file << time_used << endl;
            }
            else
            {
                cout << "Please choose p and n such that p is divisible by n" << endl;
            }
        }
    }

    file.close();
    
    return 0;
}