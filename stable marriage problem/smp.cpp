#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <omp.h>
using namespace std;

vector<vector<int> > generate_pref_matrix(int n)
{
    vector<vector<int> > pref_mat(n, vector<int>(n));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            pref_mat[i][j] = j;
        }
    }

    return pref_mat;
}

void shuffle_rows(int n, vector<vector<int> >& mat) 
{
    for (int i = 0; i < n; i++) 
    {
        random_shuffle(mat[i].begin(), mat[i].end());
    }
}

void print_matrix(int n, vector<vector<int> > mat)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
}

vector<vector<int> > get_ranking_matrix(int n, vector<vector<int> > wom_mat)
{
    vector<vector<int> > rank(n, vector<int>(n));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            rank[i][wom_mat[i][j]] = j;
        }
    }
    return rank;
}

vector<int> gale_shapley(int n, vector<vector<int> > man_mat, vector<vector<int> > wom_mat, vector<vector<int> > rank)
{
    vector<int> curr_partner_W(n, -1), curr_rank_M(n, 0);
    queue<int> free_M;

    // Initialize all men as free
    for (int i = 0; i < n; i++)
    {
        free_M.push(i);
    }
    // Start the algorithm
    while (!free_M.empty())
    {
        int m = free_M.front();
        free_M.pop();
        // Get the man's preferred woman who hasn't rejected him
        int w = man_mat[m][curr_rank_M[m]++];
        // Check if the woman is free
        if (curr_partner_W[w] == -1)
        {
            curr_partner_W[w] = m;
        }
        else
        {
            // The woman is already engaged
            int m1 = curr_partner_W[w];
            if (rank[w][m] < rank[w][m1])
            {
                // The woman prefers this man over her current partner
                curr_partner_W[w] = m;
                free_M.push(m1); // Add the rejected man back to the queue
            }
            else
            {
                // The woman prefers her current partner over this man
                free_M.push(m); // Add this man back to the queue
            }
        }
    }

    return curr_partner_W;
}

int main()
{
    int n;
    int mat_size[] = {125, 250, 500, 1000, 2000, 4000, 8000, 16000, 32000, 64000};
    for (int i = 0; i < 10; i++)
    {
        n = mat_size[i];

        double startTime1 = omp_get_wtime();
        vector<vector<int> > mat = generate_pref_matrix(n);
        double endTime1 = omp_get_wtime();
        double time1 = endTime1 - startTime1;

        double startTime2 = omp_get_wtime();
        shuffle_rows(n, mat);
        vector<vector<int> > man_pref = mat;
        double endTime2 = omp_get_wtime();
        double time2 = endTime2 - startTime2;

        double startTime3 = omp_get_wtime();
        shuffle_rows(n, mat);
        vector<vector<int> > woman_pref = mat;
        double endTime3 = omp_get_wtime();
        double time3 = endTime3 - startTime3;

        double startTime4 = omp_get_wtime();
        vector<vector<int> > ranking_mat = get_ranking_matrix(n, woman_pref);
        double endTime4 = omp_get_wtime();
        double time4 = endTime4 - startTime4;
        
        double startTime5 = omp_get_wtime();
        vector<int> res = gale_shapley(n, man_pref, woman_pref, ranking_mat);
        double endTime5 = omp_get_wtime();
        double time5 = endTime5 - startTime5;

        cout << n << "," << time1 << "," << time2 << "," << time3 << "," << time4 << "," << time5 << endl;
    }

    // cout << "Man preferences:" << endl;
    // print_matrix(n, man_pref);

    // cout << "Woman preferences:" << endl;
    // print_matrix(n, woman_pref);

    // cout << "Ranking:" << endl;
    // print_matrix(n, ranking_mat);

    // cout << "Pairs:" << endl;
    // for (int i = 0; i < n; i++)
    // {
    //     cout << i << " - " << res[i] << endl;
    // }

    return 0;
}
