#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <fstream>
#include <omp.h>
#include "matrix_generator.h"
#include "gale_shapley.h"
using namespace std;

void write_matrix_to_file(int n, int **mat)
{
    ofstream f("f.txt", ios::app);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            f << mat[i][j] << " ";
        }
        f << endl;
    }
    f << endl;
    f.close();
}

int main()
{
    int n = 8;
    int **man_mat = new int *[n];
    int **wom_mat = new int *[n];
    int **ranking_mat = new int *[n];

    for (int i = 0; i < n; i++)
    {
        man_mat[i] = new int[n];
        wom_mat[i] = new int[n];
        ranking_mat[i] = new int[n];
    }

    double startTime1 = omp_get_wtime();
    generate_pref_matrix(n, man_mat, wom_mat);
    double endTime1 = omp_get_wtime();
    double time1 = endTime1 - startTime1;

    double startTime2 = omp_get_wtime();
    srand(time(NULL));
    shuffle_rows(n, man_mat);
    shuffle_rows(n, wom_mat);
    double endTime2 = omp_get_wtime();
    double time2 = endTime2 - startTime2;

    double startTime3 = omp_get_wtime();
    fill_ranking_matrix(n, ranking_mat, wom_mat);
    double endTime3 = omp_get_wtime();
    double time3 = endTime3 - startTime3;

    double startTime4 = omp_get_wtime();
    vector<int> res = gsa_seq(n, man_mat, wom_mat, ranking_mat);
    double endTime4 = omp_get_wtime();
    double time4 = endTime4 - startTime4;

    cout << n << "," << time1 << "," << time2 << "," << time3 << "," << time4 << endl;

    write_matrix_to_file(n, man_mat);
    write_matrix_to_file(n, wom_mat);
    write_matrix_to_file(n, ranking_mat);

    cout << "Pairs:" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << i << " - " << res[i] << endl;
    }
    return 0;
}
