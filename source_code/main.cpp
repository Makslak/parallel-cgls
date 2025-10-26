#include <mpi.h>
#include <vector>
#include <random>
#include <iostream>
#include <chrono>
#include <string>
#include <stdexcept>

#include "ParallelBlockNormalCG.h"

bool ParseArgs(int argc, char** argv, int& M, int& N, double& alpha, int& iters, int& procs_expect, bool& want_help) 
{
    want_help = false;
    procs_expect = -1;
    for (int i = 1; i < argc; ++i) 
    {
        std::string a = argv[i];
        auto need = [&](const char* name) {
            if (i + 1 >= argc) 
                throw std::runtime_error(std::string("Missing value for ") + name);
            return std::string(argv[++i]);
        };

        if (a == "-h" || a == "--help")
            want_help = true;
        else if (a == "-M" || a == "--M")
            M = std::stoi(need("-M"));
        else if (a == "-N" || a == "--N")
            N = std::stoi(need("-N"));
        else if (a == "-alpha" || a == "--alpha")
            alpha = std::stod(need("-alpha"));
        else if (a == "-iters" || a == "--iters")
            iters = std::stoi(need("-iters"));
        else if (a == "-procs" || a == "--procs" || a == "-np")
            procs_expect = std::stoi(need("-procs"));
        else
            throw std::runtime_error("Unknown arg: " + a);
    }
    return true;
}

int main(int argc, char** argv)
{
    int M = 1000, N = 1000;
    double alpha = 1e-3;
    int iters = 1000;

    MPI_Init(&argc, &argv);
    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int procs_expect = -1;
    bool help = false;
    try 
    {
        ParseArgs(argc, argv, M, N, alpha, iters, procs_expect, help);
    }
    catch (const std::exception& e) 
    {
        if (world_rank == 0)
            std::cerr << e.what() << "\nUse -h for help\n";
        MPI_Finalize();
        return 1;
    }

    if (help) 
    {
        if (world_rank == 0) 
        {
            std::cout <<
                "Usage: mpiexec -n <P> ./Solver [options]\n"
                "  -M <rows>    (default: 1000)\n"
                "  -N <cols>    (default: 1000)\n"
                "  -alpha <val> (default: 1e-3)\n"
                "  -iters <k>   (default: 1000)\n"
                "  -procs <P>   (optional check: compare with actual MPI -n)\n";
        }
        MPI_Finalize();
        return 0;
    }

    if (procs_expect > 0 && procs_expect != world_size && world_rank == 0) 
    {
        std::cerr << "[warn] -procs " << procs_expect
            << " != actual world size " << world_size
            << " (MPI процессы задаются mpiexec -n ...)\n";
    }

    {
        ParallelBlockNormalCG Solver(MPI_COMM_WORLD, M, N, alpha, iters);

        std::vector<double> A, b, x0;
        if (Solver.rank() == 0) {
            std::mt19937_64 rng(42);
            std::uniform_real_distribution<double> d(-1.0, 1.0);
            A.resize((size_t)M * N);
            b.resize(M);
            x0.resize(N);
            for (double& v : A)  v = d(rng);
            for (double& v : b)  v = d(rng);
            for (double& v : x0) v = d(rng);
        }

        Solver.SetGlobalData(A, b, x0);

        MPI_Barrier(MPI_COMM_WORLD);
        auto t0 = std::chrono::high_resolution_clock::now();

        Solver.Solve();

        auto t1 = std::chrono::high_resolution_clock::now();
        double local_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        double wall_ms = 0.0;
        MPI_Reduce(&local_ms, &wall_ms, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        auto x = Solver.GatherXToRoot();
        if (world_rank == 0) {
            std::cout   << " M=" << M
                        << " N=" << N
                        << " alpha=" << alpha
                        << " iters=" << iters
                        << " procs=" << world_size
                        << " Solve_ms=" << wall_ms
                        << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
