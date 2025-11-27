#pragma once

#include <cmath>

#include "integration.h"
#include "basis_dg.h"
#include "transform_dg.h"

using namespace std;

class Integration_MSE_mat_cal : public Integration
{
	public:
		Integration_MSE_mat_cal() = default;
		~Integration_MSE_mat_cal() = default;

		void init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const EQ_reader &, const Functor &, const int);
		void ME_mat_cal_f_init(const ElementX &, const ElementV &, const int &, vector<value_type> &, vector<value_type> &, vector<value_type> &);

	protected:
		const Basis *basis_cg = nullptr;
		const Quadrature *quadrature_cg = nullptr;
};

class Integration_Col : public Integration
{
	public:
		Integration_Col() = default;
		~Integration_Col() = default;

		void init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const EQ_reader &, const int, const int &, const value_type &);

		void col3_init(const int &, const int &, const vector<int> &, const vector<int> &, const vector<value_type> &, const vector<value_type> &, const vector<value_type> &, const vector<value_type> &, const int &, const int &, const value_type &);
		void vol_cal_weight(const ElementX &, const int, const int, vector<value_type> &);
		void col_cell_mat_cal(const ElementX &, const ElementV &, vector<value_type> &);
		void col_cell_consv_mat_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &, const int &);
		void col_fhat_mat_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, vector<value_type> &, vector<int> &);
		void col_fhat_mat_cal_bd(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, vector<value_type> &, vector<int> &);
		void col_f_dg_to_valid_fhat_mat_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &);
		void col3_bc_mat_cal(const ElementX &, const ElementV &, const int &, const int &, const vector<value_type> &, const vector<value_type> &, const int &, vector<value_type> &);
		void col3_hg_vol_mat_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &, vector<value_type> &, vector<value_type> &);
		void col3_hg_b_to_hg_a_vol_mat_cal(const ElementX &, const ElementV &, const int &, const int &, vector<value_type> &, vector<value_type> &);
		void col3_small_hg_b_to_hg_a_vol_mat_cal(const ElementX &, const ElementV &, const int &, const int &, vector<value_type> &);
		void col3_slow_moment_to_fast_hg_cal(const ElementX &, const ElementV &, const int &, const int &, const value_type &, const int &, Vector &);
		void col3_slow_moment_to_fast_hg_adj_n_cal(const ElementX &, const ElementV &, const int &, const int &, const value_type &, const int &, Vector &, int &);
		void col3_small_fb_to_big_h_a_source_mat_cal(const ElementX &, const ElementV &, const int &, const int &, vector<value_type> &);
		void col3_edge_bd_hg_stiffness(const ElementX &, const EdgeV &, const int &, const int &, vector<value_type> &);
		Vector col_fM_coeff_vol_cal(const ElementX &, const ElementV &, const int &, const Vector &);

		void col3_hg_to_Qx_cal(const ElementX &, Vector &);
		void col3_hg_to_Qv_vol_cal(const ElementX &, const ElementV &, Vector &);
		void col3_hMgM_to_Qv_vol_cal(const int &, const ElementX &, const ElementV &, const Vector &, const int &, Vector &);
		void col3_hMgM_to_Qv_vol_adj_n_cal(const int &, const ElementX &, const ElementV &, const Vector &, const int &, Vector &, int &, value_type);
		void col3_hg_S_LF_Vmat_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &, Vector &, Vector &);
		void col3_hMgM_S_Vmat_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &, Vector &);
		void col3_hMgM_S_Vmat_adj_n_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &, Vector &, int &);

		void col3_g_S_Vmat_bd_cal(const ElementX &,  const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &);
		void col3_gM_S_Vmat_bd_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &);
		void col3_gM_S_Vmat_bd_adj_n_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &, int &);
		void col3_h_flux_qd_cal(const ElementX &,  const EdgeV &, const int &, const Vector &, Vector &);
		void col3_hM_flux_qd_cal(const int &, const ElementX &, const EdgeV &, const int &, const Vector &, const Vector &, const int &, Vector &);
		void col3_f_to_consv_vol_cal(const ElementX &, const ElementV &, Vector &);
		void col3_f_to_mom_consv_surf_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, Vector &);
		void col3_f_to_en_consv_inner_surf_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &);
		void col3_f_to_en_consv_bd_surf_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &);

	protected:
		int tot_species_num, sp_own_id;
		int ei_pitch_angle_col_op, ei_pitch_no_v_in_nu_ei_op;
		value_type ei_pitch_no_v_in_nu_ei_v_e_min;

		int diag_num;
		int data_size_0d, data_size_1d;
		value_type fhat_vol_wid;

		int col_edge_type_num, dof_vp_1d, dof_u_1d, dof_v_1d, fhat_dim_per_edge;
		value_type dvp, du;
		vector<int> col_fhat_mat_num;
		vector<vector<int>> col_fhat_mat_half_size, col_fhat_org_mat_group, col_fhat_org_mat_index, col_fhat_new_mat_group, col_fhat_new_mat_index;
		vector<vector<vector<int>>> col_fhat_org_mat_inv_info, col_fhat_new_mat_inv_info;

		vector<Vector> basisX_qd, basisV_qd;
		vector<vector<vector<Vector>>> basisXV_qd;
		vector<vector<vector<vector<vector<Vector>>>>> basisXVE_qd;
		vector<vector<vector<vector<Vector>>>> col_vol_terms_qd;
		vector<vector<vector<vector<Vector>>>> col_edge_terms_qd, col_edge_terms_fhat_qd;

		int col3_int_vp_n, col3_int_u_n;
		int col3_int_org_f_vp_n, col3_int_org_f_u_n;
		int col3_int_hg_tot_n, col3_int_hg_inner_n;
		value_type col3_int_vp_min, col3_int_vp_max, col3_int_u_min, col3_int_u_max;
		value_type col3_int_vp_tot_min, col3_int_u_tot_min;
		value_type col3_int_dvp, col3_int_du;
		value_type col3_int_inv_dvp, col3_int_inv_du;
		vector<value_type> Mb_ov_Ma_arr;

		value_type col3_int_spline_half_wid;

		vector<int> vp_n_arr, u_n_arr, hg_inner_n, hg_total_n;
		vector<value_type> vp_min_arr, vp_max_arr, u_min_arr, u_max_arr, dvp_arr, du_arr, inv_dvp_arr, inv_du_arr;

		int col3_int_h_interpol_num, col3_int_g_interpol_num;
		vector<vector<SparseMatrix>> col3_int_vol_hg_interpol_mat;
		vector<vector<SparseMatrix>> col3_int_edge_hg_interpol_mat;

};

