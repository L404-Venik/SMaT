#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

#include "Domain.h"


bool domain::FindOptimalPartitionRC(int P, int M, int N, int& R, int& C)
{
    // Checking for all pairs R and C, so R * C = P.
    for (int r_trial = 1; r_trial * r_trial <= P; ++r_trial) 
    {
        if (P % r_trial == 0) 
        {
            int R1 = r_trial;
            int C1 = P / r_trial;

            int R2 = C1;
            int C2 = R1;

            std::vector<std::pair<int, int>> rc_pairs;
            rc_pairs.push_back({ R1, C1 });
            if (R1 != R2) 
            {
                rc_pairs.push_back({ R2, C2 });
            }

            for (const auto& pair : rc_pairs) 
            {
                int R_cur = pair.first;
                int C_cur = pair.second;

                int min_nx = M / R_cur;
                int max_nx = (M + R_cur - 1) / R_cur;

                int min_ny = N / C_cur;
                int max_ny = (N + C_cur - 1) / C_cur;

                if (min_nx == 0 || min_ny == 0) 
                    continue;

                // [1/2, 2] sides ratio constrain check:
                double min_ratio = (double)min_nx / max_ny;

                double max_ratio = (double)max_nx / min_ny;

                if (min_ratio >= 0.5 && max_ratio <= 2.0) 
                {
                    R = R_cur;
                    C = C_cur;
                    return true; // found
                }
            }
        }
    }

    // not found
    R = 0; C = 0;
    return false; 
}

std::vector<Domain> domain::SplitDomain2D(int P, const Domain& InitialDomain)
{
    int M = InitialDomain.Nx_local;
    int N = InitialDomain.Ny_local;

    int R = 0, C = 0; // R - доменов по X, C - доменов по Y
    const int OVERLAP_NODES = 2; // Желаемый нахлёст в узлах
    // Нахлёст в узлах (O) = Нахлёст в сегментах (O_seg) + 1
    const int OVERLAP_SEGMENTS = OVERLAP_NODES - 1; // 1 сегмент

    if (!FindOptimalPartitionRC(P, M, N, R, C)) 
    {
        std::cerr << "error finding partition" << std::endl;
        return {};
    }

    // Base segment distribution for load balancing
    int M_rem = M % R;
    int M_base = M / R;
    int N_rem = N % C;
    int N_base = N / C;

    // Grid steps
    double dx = (InitialDomain.x_max - InitialDomain.x_min) / M;
    double dy = (InitialDomain.y_max - InitialDomain.y_min) / N;

    // Unique segment sizes for each domain
    std::vector<int> M_seg_sizes(R);
    std::vector<int> N_seg_sizes(C);

    // The number of segments differs by at most 1 (M_base or M_base + 1)
    for (int i = 0; i < R; ++i) 
    {
        M_seg_sizes[i] = M_base + (i < M_rem ? 1 : 0); // M_segments
    }
    for (int j = 0; j < C; ++j) 
    {
        N_seg_sizes[j] = N_base + (j < N_rem ? 1 : 0); // N_segments
    }

    // Unique segment boundaries (global segment indices [0, M])
    // These are the ideal split points that ensure load balancing.
    std::vector<int> x_unique_boundaries(R + 1);
    std::vector<int> y_unique_boundaries(C + 1);
    x_unique_boundaries[0] = 0;
    y_unique_boundaries[0] = 0;

    for (int i = 0; i < R; ++i) 
    {
        x_unique_boundaries[i + 1] = x_unique_boundaries[i] + M_seg_sizes[i];
    }
    for (int j = 0; j < C; ++j) 
    {
        y_unique_boundaries[j + 1] = y_unique_boundaries[j] + N_seg_sizes[j];
    }

    // 4. Create Domain structures, applying overlap/halo
    std::vector<Domain> domains;
    domains.reserve(P);

    for (int j = 0; j < C; ++j) 
    {
        for (int i = 0; i < R; ++i) 
        {
            Domain d;

            // Start index of the UNIQUE segments (i.e., where the previous domain ended)
            int x_unique_start = x_unique_boundaries[i];
            int x_unique_end = x_unique_boundaries[i + 1];

            // 1. Calculate total index range (x_start_idx, x_end_idx) including overlap
            int x_start_idx = x_unique_start + 1;
            int x_end_idx = x_unique_end;

            if (i > 0) 
            {
                // Add overlap *before* the start for all domains except the first
                x_start_idx -= OVERLAP_SEGMENTS;
            }
            if (i < R - 1) 
            {
                // Add overlap *after* the end for all domains except the last
                x_end_idx += OVERLAP_SEGMENTS;
            }

            // Apply global constraints (must not exceed [0, M])
            if (i == 0) x_start_idx = 0;
            if (i == R - 1) x_end_idx = M;

            // 2. Set Domain physical coordinates and total size
            d.x_min = InitialDomain.x_min + x_start_idx * dx;
            d.x_max = InitialDomain.x_min + x_end_idx * dx;
            d.Nx_total = x_end_idx - x_start_idx + 1; // Total nodes = segments + 1

            // 3. Set Nx_local (Unique nodes that belong to this domain)
            // The number of unique nodes is (M_segments + 1).
            d.Nx_local = M_seg_sizes[i] + 1;
            if (i > 0) 
            {
                // The first node of a non-first domain is an overlap node 
                // counted in the previous domain, so we subtract it from Nx_local.
                d.Nx_local -= (OVERLAP_NODES - 1); // Subtract 1 node
            }

            // --- Y-Axis Calculations (Analogous) ---
            int y_unique_start = y_unique_boundaries[j];
            int y_unique_end = y_unique_boundaries[j + 1];

            int y_start_idx = y_unique_start + 1;
            int y_end_idx = y_unique_end;

            if (j > 0) 
            {
                y_start_idx -= OVERLAP_SEGMENTS;
            }
            if (j < C - 1) 
            {
                y_end_idx += OVERLAP_SEGMENTS;
            }

            if (j == 0) y_start_idx = 0;
            if (j == C - 1) y_end_idx = N;

            d.y_min = InitialDomain.y_min + y_start_idx * dy;
            d.y_max = InitialDomain.y_min + y_end_idx * dy;
            d.Ny_total = y_end_idx - y_start_idx + 1;

            d.Ny_local = N_seg_sizes[j] + 1;
            if (j > 0) 
            {
                d.Ny_local -= (OVERLAP_NODES - 1);
            }

            domains.push_back(d);
        }
    }

    return domains;
}