// Parallel Breadth First Search
// -----------------------------
// Performs a BFS starting from vertex 1
// The parent of each vertex in the BFS tree along with its distance from the starting
// vertex is computed.
//
// The algorithm should gather all discovered vertices from round i, so that they can be
// distributed among the threads before the search in round i+1.
//
// Parameters:
// n     : number of vertices
// ver   : array of length n. ver[i] points to the start of the neighbor list of vertex i in edges
// edges : array containing lists of neighbors for each vertex, each edge is listed in both direction
// p     : array of length n used for parent pointers
// dist  : array of length n used for distance from starting vertex
// S     : array of length n used for maintaining queue of vertices to be processed
// T     : array of length n where n >> number of threads.
//
// Note that the vertices are numbered from 1 to n (inclusive). Thus there is
// no vertex 0.

void pbfs(int n, int *ver, int *edges, int *p, int *dist, int *S, int *T)
{
    int i, j;
    int v, w;

    const int n_t = omp_get_num_threads();
    const int t = omp_get_thread_num();

    int *loc_T = (int *)malloc(n * sizeof(int));
    int loc_num_w;

#pragma omp single
    {
        for (i = 1; i <= n; i++)
        {
            p[i] = -1;
            dist[i] = -1;
        }

        S[1] = 1;
        p[1] = 1;
        dist[1] = 0;

        S[0] = 1;
    }

    while (S[0] != 0)
    {
        loc_num_w = 0;

#pragma omp for
        for (i = 1; i <= S[0]; i++)
        {
            v = S[i];
            for (j = ver[v]; j < ver[v + 1]; j++)
            {
                w = edges[j];
                if (p[w] == -1)
                {
                    p[w] = v;
                    dist[w] = dist[v] + 1;
                    loc_T[loc_num_w++] = w;
                }
            }
        }
        T[t] = loc_num_w;

#pragma omp barrier

#pragma omp single
        {
            S[0] = 0;
            for (i = 0; i < n_t; i++)
            {
                S[0] += T[i];
            }

            for (i = 1; i < n_t; i++)
            {
                T[i] += T[i - 1];
            }
        }

        for (i = 0; i < loc_num_w; i++)
        {
            j = T[t] - loc_num_w + 1 + i; // +1 because S[0] is sum
            S[j] = loc_T[i];
        }

        T[t] = 0;

#pragma omp barrier
    }
}