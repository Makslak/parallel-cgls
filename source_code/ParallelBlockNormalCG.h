#pragma once

#include <mpi.h>
#include <vector>


class ParallelBlockNormalCG
{
public:

	ParallelBlockNormalCG(MPI_Comm world, int M, int N, double alpha = 0.0, int iters = -1);
	~ParallelBlockNormalCG();

	int rank() const { return rank_; }
	int size() const { return P_; }
	int R() const { return R_; }
	int row() const { return row_; }
	int col() const { return col_; }
	int M() const { return M_; }
	int N() const { return N_; }
	int m_loc() const { return m_loc_; }
	int n_loc() const { return n_loc_; }

	void SetGlobalData(const std::vector<double>& A_full,
		const std::vector<double>& b_full,
		const std::vector<double>& x0_full);

	void SetLocalDataWithLayout(const std::vector<double>& A_part,
		const std::vector<double>& b_part,
		const std::vector<double>& x_part,
		const std::vector<int>& rcounts_M,
		const std::vector<int>& displs_M,
		const std::vector<int>& rcounts_N,
		const std::vector<int>& displs_N);

	void Solve();

	std::vector<double> GatherXToRoot() const;

	const std::vector<double>& APart() const { return A_part_; }
	const std::vector<double>& BPart() const { return b_part_; }
	const std::vector<double>& XPart() const { return x_part_; }

	ParallelBlockNormalCG(const ParallelBlockNormalCG&) = delete;
	ParallelBlockNormalCG& operator=(const ParallelBlockNormalCG&) = delete;


private:
	MPI_Comm world_{ MPI_COMM_NULL };
	MPI_Comm comm_row_{ MPI_COMM_NULL };
	MPI_Comm comm_col_{ MPI_COMM_NULL };
	int P_{ 0 }, R_{ 0 }, rank_{ 0 }, row_{ 0 }, col_{ 0 };

	int M_{ 0 }, N_{ 0 };
	double alpha_{ 0.0 };
	int iters_{ 0 };

	std::vector<int> rcounts_M_, displs_M_, rcounts_N_, displs_N_;
	int m_loc_{ 0 }, n_loc_{ 0 };

	std::vector<double> A_part_, b_part_, x_part_;

	bool has_data_{ false };
	static constexpr int tag_A_ = 777;
};