#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include "/usr/local/opt/libomp/include/omp.h"

const int INF = std::numeric_limits<int>::max();
const int NUM_WOLVES = 1;
const int MAX_ITERATIONS = 10000;
const double a = 2.0; // Alpha constant
const double b = 2.0; // Beta constant

double distance(int x1, int y1, int x2, int y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

int tour_distance(const std::vector<int>& tour, const std::vector<std::vector<int>>& dist) {
    int total_distance = 0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        total_distance += dist[tour[i]][tour[i + 1]];
    }
    total_distance += dist[tour.back()][tour[0]];
    return total_distance;
}

void initialize_wolves(std::vector<std::vector<int>>& wolves, int num_cities) {
    for (int i = 0; i < NUM_WOLVES; ++i) {
        std::vector<int> tour(num_cities);
        for (int j = 0; j < num_cities; ++j) {
            tour[j] = j;
        }
        std::random_shuffle(tour.begin() + 1, tour.end());
        wolves.push_back(tour);
    }
}

void update_wolves(std::vector<std::vector<int>>& wolves, const std::vector<std::vector<int>>& dist) {
    #pragma omp parallel for
    for (size_t i = 0; i < wolves.size(); ++i) {
        int j_rand = rand() % wolves[i].size();
        int k_rand = rand() % wolves[i].size();
        while (j_rand == k_rand) {
            k_rand = rand() % wolves[i].size();
        }
        std::swap(wolves[i][j_rand], wolves[i][k_rand]);
    }
}

void two_opt(std::vector<int>& tour, const std::vector<std::vector<int>>& dist) {
    bool improvement = true;
    while (improvement) {
        improvement = false;
        for (size_t i = 1; i < tour.size() - 2; ++i) {
            for (size_t j = i + 1; j < tour.size() - 1; ++j) {
                int delta = dist[tour[i]][tour[j]] + dist[tour[i + 1]][tour[j + 1]] - dist[tour[i]][tour[i + 1]] - dist[tour[j]][tour[j + 1]];
                if (delta < 0) {
                    std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                    improvement = true;
                }
            }
        }
    }
}

std::pair<int, std::vector<int>> tsp_gwo(int num_cities, const std::vector<std::vector<int>>& dist) {
    std::vector<std::vector<int>> wolves;
    initialize_wolves(wolves, num_cities);

    int best_distance = INF;
    std::vector<int> best_tour;

    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        update_wolves(wolves, dist);

        int local_best_distance = INF;
        std::vector<int> local_best_tour;

        #pragma omp parallel for
        for (size_t i = 0; i < wolves.size(); ++i) {
            int current_distance = tour_distance(wolves[i], dist);
            two_opt(wolves[i], dist);
            if (current_distance < local_best_distance) {
                local_best_distance = current_distance;
                local_best_tour = wolves[i];
            }
        }

        #pragma omp critical
        {
            if (local_best_distance < best_distance) {
                best_distance = local_best_distance;
                best_tour = local_best_tour;
            }
        }
    }

    return std::make_pair(best_distance, best_tour);
}

int main() {
    srand(time(nullptr));

    std::ifstream input_file("input2.txt");
    if (!input_file) {
        std::cerr << "Error: Unable to open input file." << std::endl;
        return 1;
    }

    int num_cities;
    input_file >> num_cities;

    std::vector<std::vector<int>> dist(num_cities, std::vector<int>(num_cities));

    for (int i = 0; i < num_cities; ++i) {
        for (int j = 0; j < num_cities; ++j) {
            input_file >> dist[i][j];
        }
    }

    auto start = std::chrono::steady_clock::now();
    std::pair<int, std::vector<int>> result = tsp_gwo(num_cities, dist);
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Execution time: " << elapsed_seconds.count() << "s" << std::endl;

    std::cout << "The shortest tour length is: " << result.first << std::endl;
    std::cout << "The optimal tour path is: ";
    for (size_t i = 0; i < result.second.size(); ++i) {
        std::cout << result.second[i] + 1 << " ";
    }
    std::cout << std::endl;

    return 0;
}
