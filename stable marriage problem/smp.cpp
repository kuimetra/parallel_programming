#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <fstream>
#include <omp.h>
#include "matrix_generator.h"
#include "gale_shapley.h"
using namespace std;

int main()
{
    int p;
    cout << "Enter amount of threads: ";
    cin >> p;

    omp_set_num_threads(p);

    int option;
    cout << "Enter:\n1 - to generate matrix\n2 - to read from file\n--> ";
    cin >> option;

    if (option != 1 && option != 2)
    {
        cerr << "Error: invalid input!" << endl;
        exit(1);
    }

    ofstream f_w;
    ifstream f_r;
    string filename;

    int n;

    if (option == 1)
    {
        cout << "Enter n: ";
        cin >> n;

        cout << "Enter filename: ";
        cin >> filename;
    }
    else if (option == 2)
    {
        cout << "Enter filename: ";
        cin >> filename;

        f_r.open(filename.c_str());

        if (!f_r)
        {
            cerr << "Error opening file!" << endl;
            exit(1);
        }

        f_r >> n; // read size of matrices
        cout << "n = " << n << endl;
    }

    int **man_mat = new int *[n];
    int **wom_mat = new int *[n];
    int **ranking_mat = new int *[n];

    for (int i = 0; i < n; i++)
    {
        man_mat[i] = new int[n];
        wom_mat[i] = new int[n];
        ranking_mat[i] = new int[n];
    }

    if (option == 1)
    {
        double startTime_gen = omp_get_wtime();

        generate_pref_matrix(n, man_mat, wom_mat);
        srand(time(NULL));
        shuffle_rows(n, man_mat);
        shuffle_rows(n, wom_mat);
        fill_ranking_matrix(n, ranking_mat, wom_mat);

        double endTime_gen = omp_get_wtime();
        double time_gen = endTime_gen - startTime_gen;

        f_w.open(filename.c_str());

        f_w << n << endl;
        write_matrix_to_file(f_w, n, man_mat);
        cout << "✓ matrix with men's preferences" << endl;
        write_matrix_to_file(f_w, n, wom_mat);
        cout << "✓ matrix with women's preferences" << endl;
        write_matrix_to_file(f_w, n, ranking_mat);
        cout << "✓ ranking matrix" << endl;

        cout << "Matrix generation time: " << time_gen << "s" << endl;
        f_w.close();
    }
    else if (option == 2)
    {
        // read matrices from file
        for (int k = 0; k < 3; k++)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    int val;
                    f_r >> val;
                    if (k == 0)
                    {
                        man_mat[i][j] = val;
                    }
                    else if (k == 1)
                    {
                        wom_mat[i][j] = val;
                    }
                    else if (k == 2)
                    {
                        ranking_mat[i][j] = val;
                    }
                }
                string dummy;
                getline(f_r, dummy); // consume extra newline character
            }
        }

        cout << "✓ file reading is completed" << endl;

        // print matrices to console
        // print_matrix(n, man_mat);
        // print_matrix(n, wom_mat);
        // print_matrix(n, ranking_mat);
    }

    double startTime_gsa_seq = omp_get_wtime();
    vector<int> pairs_gsa_seq = gsa_seq(n, man_mat, ranking_mat);
    double endTime_gsa_seq = omp_get_wtime();
    double time_gsa_seq = endTime_gsa_seq - startTime_gsa_seq;

    double startTime_gsa_par = omp_get_wtime();
    vector<int> pairs_gsa_par = gsa_par(n, man_mat, ranking_mat);
    double endTime_gsa_par = omp_get_wtime();
    double time_gsa_par = endTime_gsa_par - startTime_gsa_par;

    cout << "[SEQ] Gale-Shapley algorithm time: " << time_gsa_seq << "s" << endl;
    cout << "[PAR] Gale-Shapley algorithm time: " << time_gsa_par << "s" << endl;

    // cout << "[SEQ] Pairs:" << endl;
    // for (int i = 0; i < n; i++)
    // {
    //     cout << i << " - " << pairs_gsa_seq[i] << endl;
    // }

    // cout << "[PAR] Pairs:" << endl;
    // for (int i = 0; i < n; i++)
    // {
    //     cout << i << " - " << pairs_gsa_par[i] << endl;
    // }
    return 0;
}
