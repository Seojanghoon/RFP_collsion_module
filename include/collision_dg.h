#pragma once

#include <mpi.h>
#include <iostream>
#include <memory>
#include <random>
#include <fstream>
#include <ctime>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <deque>

#include <Eigen/PardisoSupport>
#include<Eigen/IterativeLinearSolvers>

#include "config.h"
#include "eq_reader.h"
#include "species.h"
#include "functor.h"
#include "type.h"
#include "mesh.h"
#include "basis_dg.h"
#include "transform_dg.h"
#include "integration_dg.h"
#include "flux.h"
#include "quadrature_dg.h"
#include "diagnosis.h"
#include "postprocess.h"
#include "elliptic_integral.hpp"
#include "svd.h"   

using namespace std;


class Collision_dg
{
	public:
		Collision_dg()=delete;
		Collision_dg(Mesh &, const EQ_reader &, const Config &, const Species_data &, const MPI_Comm &comm = MPI_COMM_WORLD);
		~Collision_dg();

		//collision module setup part
		void collision_dg_general_setup(void); //collision general setup
		void collision_dg_coeff_setup(Mesh &, Quadrature &, Basis &, Flux &, int); //collision species setup
		void col_consv_mat_setup(SparseMatrix &, int, int); //setup for col_diag_cell_consv_mat[sp_id][ix] : needed for conservation 
		void col_time_out(); //collision module speed diagnostic

		//lowest order part of f update module
		void initial_h0g0_setup(vector<vector<value_type>> &); //called when Gyrokinetic class is generated
		void h0g0fM_update(vector<vector<value_type>> &, const int &, const int &); //called periodically

		//implicit time scheme collision module
		void RK_implicit_col(vector<vector<value_type>> &, const int &);

		//additional modules
		//Coulomb logarithm functions
		value_type coulomb_log_ii(const value_type &, const value_type &, const value_type &, const value_type &);
		value_type coulomb_log_ee(const value_type &ne_N, const value_type &Te_N);
		value_type coulomb_log_Ii(const value_type &ZI, const value_type &MI, const value_type &nI_N, const value_type &TI_N, const value_type &Zi, const value_type &Mi, const value_type &ni_N, const value_type &Ti_N);

		//collision frequency functions
		value_type taui_coll_time_cgs(const value_type &, const value_type &, const value_type &, const value_type &);
		value_type taui_ab_coll_time_cgs(const int sp_a, const value_type &na_N, const value_type &Ta_N, const int sp_b, const value_type &nb_N, const value_type &Tb_N);
			value_type tau_ii_coll_time_norm(const value_type &, const value_type &, const value_type &, const value_type &);
		value_type taui_ab_coll_time_norm(const int sp_a, const value_type &na_N, const value_type &Ta_N, const int sp_b, const value_type &nb_N, const value_type &Tb_N);

		//Gamma coefficients for RFP operator
		value_type gamma_ab_cgs(const value_type &, const value_type &, const value_type &, const value_type &);
		value_type gamma_ab_norm(const value_type &, const value_type &, const value_type &, const value_type &);
		value_type gamma_ab_gen_cgs(int sp_a, const value_type &na_N, const value_type &Ta_N, int sp_b, const value_type &nb_N, const value_type &Tb_N);
		value_type gamma_ab_gen_norm(int sp_a, const value_type &na_N, const value_type &Ta_N, int sp_b, const value_type &nb_N, const value_type &Tb_N);

		//Electon <-> ion momentum, energy compensation routines
		value_type ie_mom_en_transfer_fn(int, int, int, int);
		void ie_mom_en_transfer_write_fn(int, int, int, int, value_type); //momentum, energy conservation write for restart of the simulation

	protected:
		MPI_Comm comm;
		int rank, nproc, nstep;

		Mesh *mesh_sp1 = nullptr;
		vector<Mesh*> mesh_arr;
		vector<Quadrature*> quadrature_arr;
		vector<Basis*> basis_arr;
		vector<Flux*> flux_arr;

		const EQ_reader *eq_reader = nullptr;
		const Config *config = nullptr;
		const Species_data *species_data = nullptr;
		const Functor *bfield = nullptr, *bsfield = nullptr;

		vector<Integration *> integration_col_arr;

		void sp_h0g0_SE_update(const int, const int); //core function for the lowest order part of f update module
		void sp_hMbgMb_with_UaTa_update(const int, const int); //Correction term for non-zero equilibrium collision effect


		void RK_implicit_col_single(const vector<Vector> &, const int &, const int &, const int &, const int &, const value_type &, vector<Vector> &); //core function for the collision operation
		void M_C_mat_fa_to_hg_ab(const int &sp_a, const int &ix, const Vector &h_ab, const Vector &g_ab, SparseMatrix &M_C_mat, int flag); //Big matrix construction for the collision operation

		value_type get_w_loc_mod(const value_type &);

		//general variable
		int col_on_op, col_method, col_ion_period, col_elec_period, col_fM_update_period, col_f_2D_single_cell_diag_period, col_diag_1d_period, col_adjustable_quad_op;
		value_type col_adjustable_quad_tor;
		int col_consv_can_ang_mom, col3_hMgM_analytic_op, col3_hMbgMb_with_UaTa_op, col3_error_dump_method, col_consv_onoff_op, col3_hMbgMb_with_UaTa_wo_self_op;
		int col_ei_pitch_angle_op, col_ei_pitch_no_v_in_nu_ei_op, col_ei_pitch_hMbgMb_with_UaTa_op;
		value_type col_ei_pitch_no_v_in_nu_ei_v_e_min;
		value_type system_dt;
		value_type max_del_t_ov_tau_i, ratio_del_col_t_ov_delt_max_theo;
		value_type col_fhat_vol_wid;
		value_type col_mult_fac;
		value_type col3_T_eV_min_lim, col3_norm_T0_eV;

		int tot_species_num, col_diag_quantity_num;
		int nx, tor_wedge_n;
		vector<value_type> vol_cell, vol_B_cell, R_cen_cell, Z_cen_cell;

		int col3_vdim_factor, col3_diff_sp_nl_col_op;
		int col3_h0_g0_op;

		value_type ni00_norm_mks;


		//species variable
		vector<int> fhat_dim_per_edge_arr;
		vector<vector<SparseMatrix>> col_fhat_cell_mat;

		vector<value_type> Ms_arr, Zs_arr, Ms_ov_Mp_arr, Zs_ov_e_arr;
		vector<vector<value_type>> den_arr;

		vector<int> dof_arr, dof_x_arr, dof_v_arr, nve_arr, nv, f_data_size_per_ix, sp_kinetic_op_arr;
		vector<vector<SparseMatrix>> col_diag_cell_mat, col_diag_cell_consv_mat;

		vector<vector<SparseMatrix>> col_f_dg_to_valid_fhat_mat_arr;

		vector<vector<vector<value_type>>> col3_gamma_ab_arr;
		vector<vector<vector<int>>> col3_smaller_vspace_ab_arr;
		vector<vector<vector<value_type>>> col3_vspace_domain_ratio_ab_arr;

		int vspace_edge_num_of_type = 4;
		vector<vector<vector<int>>> vspace_edge_list_arr;
		vector<vector<int>> vspace_edge_ele_arr;

		vector<int> col3_vp_n_arr, col3_u_n_arr, col3_hg_inner_n_arr, col3_hg_tot_n_arr;
		vector<value_type> col3_vp_min_arr, col3_vp_max_arr, col3_u_min_arr, col3_u_max_arr;
		vector<value_type> col3_dvp_arr, col3_du_arr, col3_inv_dvp_arr, col3_inv_du_arr, col3_vp_tot_min_arr, col3_u_tot_min_arr;

		vector<vector<vector<int>>> vp_tot_min_index_arr, vp_tot_max_index_arr, u_tot_min_index_arr, u_tot_max_index_arr;

		vector<int> col3_h_bc_n_arr, col3_g_bc_n_arr;
		vector<vector<SparseMatrix>> col3_f_to_h_bc_mat, col3_f_to_g_bc_mat;
		vector<vector<SparseMatrix>> col3_f_to_h_source_mat, col3_h_to_g_source_mat;

		vector<vector<SparseMatrix>> col3_h_stiffness_mat, col3_g_stiffness_mat;
		vector<vector<Eigen::PardisoLU<SparseMatrix>*>> col3_h_solver_stiffness, col3_g_solver_stiffness;
		vector<vector<Eigen::PardisoLU<SparseMatrix>*>> col3_hg_formula_to_hg0_solver;

		vector<vector<Eigen::PardisoLLT<SparseMatrix>*>> col3_M_solver_arr;

		vector<vector<vector<vector<int>>>> col3_adjustable_quad_n_arr;
		int ilower_ele;

		int col_implicit_solve_op;
		value_type col_implicit_tor;
		int col_fast_to_slow_moment_op;
		value_type col_fast_to_slow_moment_ratio_max;
		vector<vector<vector<Vector>>> col_fast_to_slow_moment_hg_Vec;

		vector<vector<vector<Matrix>>> col_f_dg_to_f_hat_mat_arr; //used in M_C_mat_fa_to_hg_ab

		vector<vector<vector<SparseMatrix>>> consv_fa_fb0_mat_arr;
		vector<vector<SparseMatrix>> consv_fa_del_hg_b_mat_arr;
		vector<vector<vector<Vector>>> consv_fa_fb0_UaTa_mat_arr;

		int col_implicit_linear_ion_op, col_implicit_linear_electron_op;
		vector<int> col_implicit_linear_sp_op;

		vector<vector<vector<Eigen::PardisoLU<SparseMatrix>*>>> col_implicit_linear_solver_arr;
		vector<vector<Vector>> col3_C_fMa_to_fMb_wUaTa_summed_arr, col3_C_fMa_to_fMb_diff_UbTb_UaTa_summed_arr;
		vector<vector<SparseMatrix>> col_sp_Minv_arr, col_sp_M_arr, consv_fa_fb0_mat_sum_arr;

		int col3_fM_coeff_for_hg_init_flag = 0;
		vector<Vector> col3_fM_coeff_for_hg_arr, col3_h0_g0_avged_Q_arr, max_t_step_ratio;
		vector<vector<Vector>> col3_fM_coeff_for_hg_nx_arr, col3_nUT_arr, col3_M_C_fMa_to_fMb_wUaTa_summed_arr, col3_nUT_before_arr;
		vector<int> col3_fM_coeff_update_before_time_arr;
		vector<vector<vector<Vector>>> col3_f0_h0g0_col_arr, col3_ei_mom_en_transfer_arr, col3_ie_mom_en_transfer_arr;
		vector<vector<int>> col3_f0_h0g0_col_flag_arr;
		vector<vector<vector<int>>> col3_sp_ab_ix_col_flag_arr;
		vector<vector<vector<int>>> col3_fa_to_fMb_col_flag_arr;
		vector<int> col3_valid_ix_flag_arr;
		vector<vector<vector<SparseMatrix>>> col3_hg_b_to_hg_a_source_arr, col3_hg_b_to_hg_a_bc_arr;

		vector<vector<value_type>> tau_aa_init;
		vector<vector<vector<value_type>>> col3_hg_bc_points_vp_u_arr;


		vector<vector<vector<Vector>>> col3_hg_Xmat_arr, col3_h_LF_Vmat_arr;
		vector<vector<SparseMatrix>> col3_h_to_Qv_mat, col3_g_to_Qv_mat;
		vector<vector<SparseMatrix>> col3_h_to_Qv_mat_whole;

		vector<vector<vector<vector<Vector>>>> col3_Qv_h0_ab_arr, col3_Qv_g0_ab_arr;
		vector<vector<Vector>> col3_Qv_h0_ab_summed_arr, col3_Qv_g0_ab_summed_arr;



		vector<int> col_edge_flux_qd_num;
		vector<Vector> col_edge_flux_qd_points;

		vector<vector<vector<Vector>>> col3_stored_del_quantity_arr;

					
		//col3_f_to_mom_consv_mat_arr is a parallel shift in V_parallel direction : change momentum mainly
		//col3_f_to_en_consv_mat_arr is a isotropic diffusion in v space : change energy mainly
		vector<vector<vector<SparseMatrix>>> col3_f_to_mom_consv_mat_arr;
		vector<vector<SparseMatrix>> col3_f_to_en_consv_mat_arr;

		vector<int> S_lp_f_index_arr;
		vector<vector<vector<int>>> Q_S_op_flags_arr;

		//Collision module speed diagnostic variables
		clock_t tot_dc1=0.0, tot_dc2=0.0, tot_dc3=0.0, tot_dc4=0.0, tot_dc5=0.0, tot_dc6=0.0;
		clock_t col_dc1=0.0, col_dc2=0.0, col_dc2_1=0.0, col_dc2_2=0.0, col_dc2_3_1_1=0.0, col_dc2_3_1_2=0.0, col_dc2_3_2=0.0, col_dc2_4=0.0, col_dc2_5=0.0, col_dc2_6;  

		clock_t col_dc3 = 0.0, col_dc4 = 0.0, col_dc5 = 0.0, col_dc6 = 0.0, col_dc7 = 0.0, col_dc7_1 = 0.0;
		clock_t col_dc7_2 = 0.0, col_dc7_3 = 0.0, col_dc7_4 = 0.0;



};
