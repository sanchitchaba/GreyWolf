#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <chrono>
#include <thread>

const int INF = std::numeric_limits<int>::max();

std::pair<int, std::vector<int> > tsp(int n, const std::vector<std::vector<int> >& dist) {
    std::vector<std::vector<int> > dp(1 << n, std::vector<int>(n, INF));
    std::vector<std::vector<int> > parent(1 << n, std::vector<int>(n, -1));

    dp[1][0] = 0;

    for (int visited = 1; visited < (1 << n); ++visited) {
        for (int current = 0; current < n; ++current) {
            if (!(visited & (1 << current))) {
                continue;
            }

            int prev_visited = visited & ~(1 << current);
            for (int prev = 0; prev < n; ++prev) {
                if (prev == current || !(visited & (1 << prev))) {
                    continue;
                }

                if (dp[prev_visited][prev] != INF) {
                    int new_cost = dp[prev_visited][prev] + dist[prev][current];
                    if (new_cost < dp[visited][current]) {
                        dp[visited][current] = new_cost;
                        parent[visited][current] = prev;
                    }
                }
            }
        }
    }

    int min_cost = INF;
    int last_city = -1;
    for (int i = 1; i < n; ++i) {
        if (dp[(1 << n) - 1][i] != INF) {
            int new_cost = dp[(1 << n) - 1][i] + dist[i][0];
            if (new_cost < min_cost) {
                min_cost = new_cost;
                last_city = i;
            }
        }
    }

    std::vector<int> optimal_path;
    int visited = (1 << n) - 1;
    int current = last_city;

    while (current != -1) {
        optimal_path.push_back(current);
        int prev = parent[visited][current];
        visited &= ~(1 << current);
        current = prev;
    }

    std::reverse(optimal_path.begin(), optimal_path.end());
    optimal_path.push_back(0);

    return std::make_pair(min_cost, optimal_path);
}

int main() {
    std::ifstream input_file("input2.txt");
    if (!input_file) {
        std::cerr << "Error: Unable to open input file." << std::endl;
        return 1;
    }

    int n;
    input_file >> n;

    std::vector<std::vector<int> > dist(n, std::vector<int>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            input_file >> dist[i][j];
        }
    }

    auto start = std::chrono::steady_clock::now();
    std::pair<int, std::vector<int> > result = tsp(n, dist);
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;
    double desired_time = elapsed_seconds.count() + 0.01;

    std::cout << "Execution time: " << elapsed_seconds.count() << "s" << std::endl;

    std::cout << "The shortest tour length is: " << result.first << std::endl;
    std::cout << "The optimal tour path is: ";
    for (int i = 0; i < result.second.size(); ++i) {
        std::cout << result.second[i] + 1 << " ";
    }
    std::cout << std::endl;

    return 0;
}
