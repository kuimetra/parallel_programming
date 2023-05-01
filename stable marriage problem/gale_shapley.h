vector<int> gsa_seq(int n, int **man_mat, int **wom_mat, int **rank)
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