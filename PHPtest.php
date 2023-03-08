#include <mpi.h>

#include <iostream>
#include <random>
#include <vector>

// dot: ans = (x, y)
// See https://www.mpich.org/static/docs/v3.3/www3/MPI_Allreduce.html

// utils
void print_vec(std::vector<double>& vec, int rank);
void random_array(std::vector<double>& vec);

double dot(size_t size, double* x, double* y, MPI_Comm comm) {
    double ans = 0;
    for (size_t i = 0; i < size; i++) {
        ans += x[i] * y[i];
    }
    MPI_Allreduce(&ans, &ans, 1, MPI_DOUBLE, MPI_SUM, comm);
    return ans;
}

// main and utils are written in C++
int main(int argc, char* argv[]) {
    size_t M = 3;
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double ans = 0;
    std::vector<double> x(M);
    std::vector<double> y(M);

    random_array(x);
    random_array(y);

    ans = dot(y.size(), x.data(), y.data(), MPI_COMM_WORLD);
    std::cout << ans << std::endl;

    MPI_Finalize();
}

void print_vec(std::vector<double>& vec, int rank) {
    for (auto i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << "," << rank << std::endl;
    }
}

void random_array(std::vector<double>& array) {
    std::random_device random;
    std::mt19937 mt(random());  // FIXME
    std::uniform_real_distribution<> rand(0.0, 1.0);

    for (size_t i = 0; i < array.size(); i++) {
        array[i] = rand(mt);
    }
}
