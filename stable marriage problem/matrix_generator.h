using namespace std;

void print_matrix(int n, int **mat)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void write_matrix_to_file(ofstream &f, int n, int **mat)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            f << mat[i][j] << " ";
        }
        f << endl;
    }
    f << endl;
}

void generate_pref_matrix(int n, int **man_mat, int **wom_mat)
{
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            man_mat[i][j] = j;
            wom_mat[i][j] = j;
        }
    }
}

void shuffle_rows(int n, int **mat)
{
    for (int i = 0; i < n; i++)
    {
        random_shuffle(mat[i], mat[i] + n);
    }
}

void fill_ranking_matrix(int n, int **ranking_mat, int **wom_mat)
{
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            ranking_mat[i][wom_mat[i][j]] = j;
        }
    }
}