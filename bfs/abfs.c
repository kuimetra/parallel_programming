// Parallel Breadth First Search
// -----------------------------
// Berforms a BFS starting from vertex 1
// The parent of each vertex in the BFS tree along with its distance from the starting
// vertex is computed.
//
// The algorithm should first perform some rounds of sequential BFS before starting a parallel
// execution. In the parallel part each thread should be allocated a part of the vertices from the
// last round of the sequential algorithm. Any discovered vertices in the parallel part should
// remain with the thread that discovered them. This continues until the entire graph has been
// explored.
//
// Parameters:
// n     : number of vertices
// ver   : ver[i] points to the start of the neighbor list of vertex i in edges
// edges : lists of neighbors of each vertex, each edge is listed in both direction
// p     : array of length n used for parent pointers
// dist  : array of length n used for distance from starting vertex
// S     : array of length n used for maintaining queue of vertices to be processed, only used in the
//         sequential part.
// T     : array of length n where n >> number of threads.
//
// Note that the vertices are numbered from 1 to n (inclusive). Thus there is
// no vertex 0.

void abfs(int n, int *ver, int *edges, int *p, int *dist, int *S, int *T)
{
    int i, j;
    int v, w;
    int num_r, num_w;
    int *temp;

    int seq_depth = 200;
    int depth = 1;

#pragma omp single
    {
        int *seq_S = S;
        int *seq_T = T;

        for (i = 1; i <= n; i++)
        {
            p[i] = -1;
            dist[i] = -1;
        }

        p[1] = 1;
        dist[1] = 0;
        seq_S[1] = 1; // add the starting vertex to S 

        seq_S[0] = 1; // number of vertices in seq_S
        seq_T[0] = 0; // number of vertices in single_T

        while (seq_S[0] != 0)
        {
            if (depth >= seq_depth)
                break;

            for (i = 0; i < seq_S[0]; i++)
            {
                v = seq_S[i + 1];
                for (j = ver[v]; j < ver[v + 1]; j++)
                {
                    w = edges[j];
                    if (p[w] == -1)
                    {
                        p[w] = v;
                        dist[w] = dist[v] + 1;
                        seq_T[seq_T[0]++ + 1] = w; // +1 because single_T[0] is sum
                    }
                }
            }
            temp = seq_S;
            seq_S = seq_T;
            seq_T = temp;
            seq_T[0] = 0;
            depth++;
        }

        for (i = 0; i <= seq_S[0]; i++)
        {
            S[i] = seq_S[i];
        }
    }

    const int n_t = omp_get_num_threads();
    const int t = omp_get_thread_num();
    int k = 10;

    int *loc_S = (int *)malloc(n * sizeof(int));
    int *loc_T = (int *)malloc(n * sizeof(int));
    int loc_num_w;

// #pragma omp single
//     for (i = 0; i <= loc_S[0]; i++)
//     {
//         loc_S[i] = S[i];
//     }

    int iter_count = 0;
    while (S[0] != 0)
    {
        loc_num_w = 0;

#pragma omp for
        for (i = 1; i <= S[0]; i++)
        {
            v = S[i];
            // v = loc_S[i];
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
        iter_count++;

#pragma omp barrier

#pragma omp single
        {
            S[0] = 0;
            // loc_S[0] = 0;
            for (i = 0; i < n_t; i++)
            {
                S[0] += T[i];
                // loc_S[0] += T[i];
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

            // if (iter_count % k == 0)
            // {
            //     S[j] == loc_S[j];
            // }
        }

        T[t] = 0;

#pragma omp barrier
    }
}