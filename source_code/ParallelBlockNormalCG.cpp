#include "ParallelBlockNormalCG.h"

#include <algorithm>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <iostream>


namespace 
{
	void counts_offsets(int total, int R, std::vector<int>& counts, std::vector<int>& offs) 
	{
		counts.assign(R, total / R);
		int rem = total % R;
		for (int i = 0; i < rem; ++i) 
			counts[i]++;

		offs.resize(R); 
		offs[0] = 0; 
		for (int i = 1; i < R; ++i) 
			offs[i] = offs[i - 1] + counts[i - 1];
	}

	double dot_local(const std::vector<double>& a, const std::vector<double>& b) 
	{
		double s = 0.0; size_t n = a.size(); 
		for (size_t i = 0; i < n; ++i) 
			s += a[i] * b[i];
		return s;
	}

	void axpy_inplace(double a, const std::vector<double>& x, std::vector<double>& y) 
	{
		size_t n = x.size(); 
		for (size_t i = 0; i < n; ++i) 
			y[i] += a * x[i];
	}
	void matvec(const std::vector<double>& A, int m, int n, const std::vector<double>& x, std::vector<double>& y) 
	{
		std::fill(y.begin(), y.end(), 0.0);
		for (int i = 0; i < m; ++i) 
		{ 
			const double* Ai = &A[(size_t)i * n];
			double sum = 0.0; 
			for (int j = 0; j < n; ++j) 
			{
				sum += Ai[j] * x[j];
				y[i] = sum;
			}
		}
	}

	void tmatvec(const std::vector<double>& A, int m, int n, const std::vector<double>& y, std::vector<double>& z)
	{
		std::fill(z.begin(), z.end(), 0.0);
		for (int i = 0; i < m; ++i) 
		{ 
			const double* Ai = &A[(size_t)i * n];
			double yi = y[i]; 
			for (int j = 0; j < n; ++j) 
				z[j] += Ai[j] * yi; 
		}
	}

	void pack_block(const std::vector<double>& A, int M, int N, int r0, int mr, int c0, int nc, std::vector<double>& out) 
	{
		out.resize((size_t)mr * nc);
		for (int i = 0; i < mr; ++i) 
		{ 
			const double* src = &A[(size_t)(r0 + i) * N + c0];
			double* dst = &out[(size_t)i * nc];
			std::memcpy(dst, src, sizeof(double) * nc); 
		}
	}
}


ParallelBlockNormalCG::ParallelBlockNormalCG(MPI_Comm world, int M, int N, double alpha, int iters)
	: world_(world), M_(M), N_(N), alpha_(alpha), iters_(iters) 
{
	if (M_ <= 0 || N_ <= 0) 
		throw std::runtime_error("M and N must be > 0");
	MPI_Comm_size(world_, &P_);
	MPI_Comm_rank(world_, &rank_);
	R_ = (int)std::lround(std::sqrt((double)P_));

	if (R_ * R_ != P_) 
	{
		if (rank_ == 0)
		{
			std::cerr << "[FATAL] Number of processes must be a perfect square (got " << P_ << ")";
			MPI_Abort(world_, 1);
		}
	}

	row_ = rank_ / R_;
	col_ = rank_ % R_;
	MPI_Comm_split(world_, row_, col_, &comm_row_);
	MPI_Comm_split(world_, col_, row_, &comm_col_);
	if (iters_ < 0) 
		iters_ = N_;


	counts_offsets(M_, R_, rcounts_M_, displs_M_);
	counts_offsets(N_, R_, rcounts_N_, displs_N_);
	m_loc_ = rcounts_M_[row_];
	n_loc_ = rcounts_N_[col_];


	A_part_.resize((size_t)m_loc_ * n_loc_);
	b_part_.resize(m_loc_);
	x_part_.resize(n_loc_);
}


ParallelBlockNormalCG::~ParallelBlockNormalCG() 
{
	if (comm_row_ != MPI_COMM_NULL) 
		MPI_Comm_free(&comm_row_);
	if (comm_col_ != MPI_COMM_NULL) 
		MPI_Comm_free(&comm_col_);
}


void ParallelBlockNormalCG::SetGlobalData(const std::vector<double>& A_full,
	const std::vector<double>& b_full,
	const std::vector<double>& x0_full)
{
	if (rank_ == 0) 
	{
		if (!A_full.empty() && (int)A_full.size() != M_ * N_) 
			throw std::runtime_error("A_full wrong size");
		if (!b_full.empty() && (int)b_full.size() != M_)
			throw std::runtime_error("b_full wrong size");
		if (!x0_full.empty() && (int)x0_full.size() != N_) 
			throw std::runtime_error("x0_full wrong size");
	}


	if (rank_ == 0) 
	{
		for (int r = 0; r < R_; ++r) 
		{
			int mr = rcounts_M_[r]; 
			int i0 = displs_M_[r];
			for (int c = 0; c < R_; ++c) 
			{
				int nc = rcounts_N_[c]; 
				int j0 = displs_N_[c]; 
				int dst = r * R_ + c;

				std::vector<double> buf; 
				pack_block(A_full, M_, N_, i0, mr, j0, nc, buf);
				if (dst == 0)
					std::copy(buf.begin(), buf.end(), A_part_.begin());
				else 
					MPI_Send(buf.data(), (int)buf.size(), MPI_DOUBLE, dst, tag_A_, world_);
			}
		}
	}
	else 
	{
		MPI_Status st{};
		MPI_Recv(A_part_.data(), m_loc_ * n_loc_, MPI_DOUBLE, 0, tag_A_, world_, &st);
	}

	if (col_ == 0) 
	{
		if (rank_ == 0) 
		{
			MPI_Scatterv(b_full.data(), rcounts_M_.data(), displs_M_.data(), MPI_DOUBLE,
				b_part_.data(), m_loc_, MPI_DOUBLE, 0, comm_col_);
		}
		else 
		{
			MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DOUBLE,
				b_part_.data(), m_loc_, MPI_DOUBLE, 0, comm_col_);
		}
	}
	MPI_Bcast(b_part_.data(), m_loc_, MPI_DOUBLE, 0, comm_row_);

	if (row_ == 0) 
	{
		if (rank_ == 0) 
		{
			MPI_Scatterv(x0_full.data(), rcounts_N_.data(), displs_N_.data(), MPI_DOUBLE,
				x_part_.data(), n_loc_, MPI_DOUBLE, 0, comm_row_);
		}
		else 
		{
			MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DOUBLE,
				x_part_.data(), n_loc_, MPI_DOUBLE, 0, comm_row_);
		}
	}
	MPI_Bcast(x_part_.data(), n_loc_, MPI_DOUBLE, 0, comm_col_);

	has_data_ = true;
}

void ParallelBlockNormalCG::SetLocalDataWithLayout(const std::vector<double>& A_part,
	const std::vector<double>& b_part,
	const std::vector<double>& x_part,
	const std::vector<int>& rcounts_M,
	const std::vector<int>& displs_M,
	const std::vector<int>& rcounts_N,
	const std::vector<int>& displs_N)
{
	if ((int)rcounts_M.size() != R_ || (int)displs_M.size() != R_ || (int)rcounts_N.size() != R_ || (int)displs_N.size() != R_)
		throw std::runtime_error("counts/displs must have size R");

	rcounts_M_ = rcounts_M; displs_M_ = displs_M;
	rcounts_N_ = rcounts_N; displs_N_ = displs_N;
	m_loc_ = rcounts_M_[row_];
	n_loc_ = rcounts_N_[col_];


	if ((int)A_part.size() != m_loc_ * n_loc_) 
		throw std::runtime_error("A_part wrong size");
	if ((int)b_part.size() != m_loc_) 
		throw std::runtime_error("b_part wrong size");
	if ((int)x_part.size() != n_loc_) 
		throw std::runtime_error("x_part wrong size");

	A_part_ = A_part; b_part_ = b_part; x_part_ = x_part;
	has_data_ = true;
}

void ParallelBlockNormalCG::Solve() 
{
	if (!has_data_) 
		throw std::runtime_error("Call SetGlobalData or SetLocalDataWithLayout before Solve()");

	std::vector<double> r_part(n_loc_, 0.0);
	std::vector<double> p_part(n_loc_, 0.0);
	std::vector<double> q_part(n_loc_, 0.0);
	std::vector<double> tmp_m(m_loc_, 0.0);
	std::vector<double> Ap_part(m_loc_, 0.0);

	double dot_pq_prev = 0.0;
	for (int s = 0; s < iters_; ++s) 
	{
		if (s == 0) 
		{
			matvec(A_part_, m_loc_, n_loc_, x_part_, Ap_part);
			MPI_Allreduce(MPI_IN_PLACE, Ap_part.data(), m_loc_, MPI_DOUBLE, MPI_SUM, comm_row_);
			for (int i = 0; i < m_loc_; ++i) 
				tmp_m[i] = Ap_part[i] - b_part_[i];
			tmatvec(A_part_, m_loc_, n_loc_, tmp_m, r_part);
			MPI_Allreduce(MPI_IN_PLACE, r_part.data(), n_loc_, MPI_DOUBLE, MPI_SUM, comm_col_);
			if (alpha_ != 0.0) axpy_inplace(alpha_, x_part_, r_part);
		}
		else 
		{
			if (dot_pq_prev == 0.0) break;
			axpy_inplace(-1.0 / dot_pq_prev, q_part, r_part);
		}

		double rr_local = dot_local(r_part, r_part);
		double rr = 0.0; MPI_Allreduce(&rr_local, &rr, 1, MPI_DOUBLE, MPI_SUM, comm_row_);
		if (rr == 0.0) break;
		for (int j = 0; j < n_loc_; ++j)
			p_part[j] += r_part[j] / rr;

		matvec(A_part_, m_loc_, n_loc_, p_part, Ap_part);
		MPI_Allreduce(MPI_IN_PLACE, Ap_part.data(), m_loc_, MPI_DOUBLE, MPI_SUM, comm_row_);
		tmatvec(A_part_, m_loc_, n_loc_, Ap_part, q_part);
		MPI_Allreduce(MPI_IN_PLACE, q_part.data(), n_loc_, MPI_DOUBLE, MPI_SUM, comm_col_);
		if (alpha_ != 0.0)
			axpy_inplace(alpha_, p_part, q_part);

		double dpq_local = dot_local(p_part, q_part);
		double dpq = 0.0;
		MPI_Allreduce(&dpq_local, &dpq, 1, MPI_DOUBLE, MPI_SUM, comm_row_);
		if (dpq == 0.0) 
			break;
		for (int j = 0; j < n_loc_; ++j) 
			x_part_[j] -= p_part[j] / dpq;
		dot_pq_prev = dpq;
	}
}

std::vector<double> ParallelBlockNormalCG::GatherXToRoot() const 
{
	std::vector<double> x_full;
	if (row_ == 0) 
	{
		if (rank_ == 0)
		{
			x_full.assign(N_, 0.0);
			MPI_Gatherv(const_cast<double*>(x_part_.data()), n_loc_, MPI_DOUBLE,
				x_full.data(), const_cast<int*>(rcounts_N_.data()), const_cast<int*>(displs_N_.data()), MPI_DOUBLE,
				0, comm_row_);
		}
		else 
		{
			MPI_Gatherv(const_cast<double*>(x_part_.data()), n_loc_, MPI_DOUBLE,
				nullptr, nullptr, nullptr, MPI_DOUBLE,
				0, comm_row_);
		}
	}
	return x_full;
}