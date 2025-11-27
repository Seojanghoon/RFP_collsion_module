#pragma once

class Gyrokinetic : public Indexer
{
	public:
		Gyrokinetic(const Sokuri &, const MPI_Comm &comm = MPI_COMM_WORLD);
		Gyrokinetic() = delete;
		~Gyrokinetic() = default;

		//main loop module
		void procs_mult(void);
		vector<vector<value_type>> sp_dg_coeff;

		protected:
		MPI_Comm comm;
		int rank, nproc;

		Collision_dg *collision_dg = nullptr;
		int gk_col_implicit_op;

		//setup module for M matrix
		void procs_mult_setup(int);
		int sp_tot_num;
		int nstep;
		value_type dt;
		int nx, nv, num_elm, dof;
};

