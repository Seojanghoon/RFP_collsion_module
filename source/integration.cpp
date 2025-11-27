#include "integration.h"

#define __false_assert__ { assert(false && "A abstract function was called."); }

void Integration::set_eq_reader(const EQ_reader &eq)
{
	eq_reader = &eq;
}


void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &) __false_assert__
void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Functor &) __false_assert__
void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Functor &, const Functor &) __false_assert__
void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Functor &) __false_assert__
void Integration::init(const Mesh &, const Basis &, const Quadrature &, const function<value_type(const Point2 &)> &) __false_assert__
void Integration::init(const Mesh &, const Basis &, const Quadrature &) __false_assert__
void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, int) __false_assert__
void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const Efld_data &, int) __false_assert__


vector<vector<value_type>> Integration::operator()(const EdgeX &, const ElementV &, const vector<value_type> &, const vector<value_type> &) __false_assert__
vector<vector<value_type>> Integration::operator()(const EdgeX &, const ElementV &, const vector<value_type> &) __false_assert__
vector<vector<value_type>> Integration::operator()(const ElementX &, const EdgeV &, const vector<value_type> &) __false_assert__
vector<vector<value_type>> Integration::operator()(const EdgeX &, const ElementV &, const Vector &) __false_assert__
vector<vector<value_type>> Integration::operator()(const ElementX &, const EdgeV &, const Vector &) __false_assert__
vector<vector<value_type>> Integration::operator()(const EdgeX &, const ElementV &, const value_type *) __false_assert__
vector<vector<value_type>> Integration::operator()(const ElementX &, const EdgeV &, const value_type *) __false_assert__
vector<Matrix> Integration::operator()(const ElementX &, const EdgeV &) const __false_assert__
vector<Matrix> Integration::operator()(const EdgeX &, const ElementV &) const __false_assert__
value_type Integration::operator()(const ElementX &ti, const ElementV &tj, int i, int j) const __false_assert__
value_type Integration::operator()(const ElementX &ti, const ElementV &tj, int i) const __false_assert__
value_type Integration::operator()(const EdgeX &, const ElementV &, int i) const __false_assert__
value_type Integration::operator()(const ElementX &, const EdgeV &, int i) const __false_assert__
value_type Integration::operator()(const EdgeX &, const ElementV &, int i, int j) const __false_assert__
value_type Integration::operator()(const ElementX &, const EdgeV &, int i, int j) const __false_assert__
value_type Integration::operator()(const ElementX &, const ElementX &, const ElementV &, int i, int j) const __false_assert__
value_type Integration::operator()(const ElementX&, const ElementV &, const ElementV &, int i, int j) const __false_assert__
value_type Integration::operator()(const ElementX&, int i, int j) const __false_assert__
value_type Integration::operator()(const ElementX&, int i) const __false_assert__
Matrix Integration::operator()(const ElementX &) const __false_assert__

void Integration::set_sp(int _sp_id)
{
	sp_id = _sp_id;
	flag_f0 = species_data->f0_type(sp_id);
	flag_fcm = species_data->fcm_type(sp_id);
}
#ifdef SEO_MOD3
void Integration::set_timer_zero(int _dummy_num)
{
	tot_dc1=0, tot_dc2=0, tot_dc3=0, tot_dc4=0;
	tot_dc5=0, tot_dc6=0, tot_dc7=0;
}

void Integration::timer_value_out(int _dummy_num)
{
	cout << "CPU time used : volume element cal1 (E field cal) : " << double(tot_dc1)/CLOCKS_PER_SEC << " s" << endl;
	cout << "CPU time used : volume element cal2 (U cal) : " << double(tot_dc2)/CLOCKS_PER_SEC << " s" << endl;
	cout << "CPU time used : volume element cal3 (U.dot cal) : " << double(tot_dc3)/CLOCKS_PER_SEC << " s" << endl;
	cout << "CPU time used : volume element cal4 (ccf load) : " << double(tot_dc4)/CLOCKS_PER_SEC << " s" << endl;
	cout << "CPU time used : volume element cal5 (basis val load) : " << double(tot_dc5)/CLOCKS_PER_SEC << " s" << endl;
	cout << "CPU time used : volume element cal6 (mat ele construct) : " << double(tot_dc6)/CLOCKS_PER_SEC << " s" << endl << endl;
	tot_dc1=0, tot_dc2=0, tot_dc3=0, tot_dc4=0;
	tot_dc5=0, tot_dc6=0, tot_dc7=0;
}
#endif


//#ifdef SEO_MOD
#ifdef SEO_MOD_BASIS
void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const EQ_reader &, const Functor &, const int)
#else
void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const Functor &, const int)
#endif
{
	assert(false && "An abstract function was called.");
}


#ifdef SEO_SE_MOD
void Integration::init_2D_efld(const Basis &, const Quadrature &)
{
	assert(false && "An abstract function was called.");
}
#endif


void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const Functor &, const EQ_reader &, const int, const int &)
{
	assert(false && "An abstract function was called.");
}

void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const EQ_reader &, const int, const int, const int, const value_type, const value_type, const int &, const vector<int> &, const int &)
{
	assert(false && "An abstract function was called.");
}
void Integration:: S1_mat_cal(const EdgeX &, const ElementV &, const int &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
void Integration:: S2_mat_cal(const ElementX &, const EdgeV &, const int &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::ME_mat_cal_f_init(const ElementX &, const ElementV &, const int &, vector<value_type> &, vector<value_type> &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::vol_cal_weight(const ElementX &, const int, const int, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::diag_cell_mat_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
#ifdef SEO_DIAG_TEST
void Integration::diag_cell_test_diag_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
#endif

#ifdef SEO_MULT
void Integration::sh_cell_diag_cal(const ElementX &, const ElementV &, const ElementV &, const int &, vector<value_type> &, vector<value_type> &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
#ifdef SEO_DEBUG_SH3
void Integration::sh_cell_diag_moments_cal(const ElementX &, const ElementV &, const ElementV &, const int &, vector<value_type> &, int )
{
	assert(false && "An abstract function was called.");
}
#endif
#endif

value_type Integration::diag_th_fn(const value_type, const value_type)
{
	assert(false && "An abstract function was called.");
}
value_type Integration::diag_mod_th(const value_type)
{
	assert(false && "An abstract function was called.");
}

//#endif
#ifdef SEO_1D_EFLD
void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const EQ_reader &, const Efld_1d_data &, const int, const int, const value_type, const value_type, const int, const vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
/*
   void Integration::integration_S_efld_1d_setup(const int, const value_type, const value_type, const int, const EQ_reader &)
   {
   assert(false && "An abstract function was called.");
   }
   */
void Integration::cell_f_to_1d_den_spline_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &, vector<int> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::cell_f_to_1d_den_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &, vector<int> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::del_E_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
Vector Integration::pot_1d_eval_test(const Point2 &)
{
	assert(false && "An abstract function was called.");
}
SparseMatrix Integration::efld_spline_M_cal(const int)
{
	assert(false && "An abstract function was called.");
}
#endif
#ifdef SEO_2D_EFLD
void Integration::init(const Mesh &, const Basis &, const Basis &, const Quadrature &, const Quadrature &, const Flux &, const Species_data &, const EQ_reader &, const int &, const int &)
{
	assert(false && "An abstract function was called.");
}
void Integration::f_to_source(const ElementX &, const ElementV &, const int &, const int &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::ij_val(const ElementX &, int , int , vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
value_type Integration::ele_to_node_vol(const ElementX &, const int &)
{
	assert(false && "An abstract function was called.");
}
#endif

#if defined(SEO_1D_EFLD) || defined(SEO_2D_EFLD)
void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const EQ_reader &, const Efld_data &, const int, const int, const int)
{
	assert(false && "An abstract function was called.");
}
#endif

#ifdef SEO_COL
void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const EQ_reader &, const int, const int &, const value_type &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_init(const int &, const int &, const vector<int> &, const vector<int> &, const vector<value_type> &, const vector<value_type> &, const vector<value_type> &, const vector<value_type> &, const int &, const int &, const value_type &)
{
	assert(false && "An abstract function was called.");
}
void Integration::col_cell_mat_cal(const ElementX &, const ElementV &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::col_cell_consv_mat_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &, const int &)
{
	assert(false && "An abstract function was called.");
}
#ifdef SEO_MOD_COL_CORR
void Integration::col_cell_consv_mat_mod_corr_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
#endif
value_type Integration::col_L2_norm_f_dg_fM_coeff_cal(const ElementX &, const ElementV &, const int &, const Vector &, const Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col_fhat_mat_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, vector<value_type> &, vector<int> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::col_fhat_mat_cal_bd(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, vector<value_type> &, vector<int> &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col_f_dg_to_valid_fhat_mat_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col_ions_vol_nl(const ElementX &ti, const ElementV &tj, const int &sp_a_id, const vector<vector<int>> &sp_ab_ion_col_flag_arr, const Vector &coeff, const vector<vector<Vector>> &h_arr, const vector<vector<Vector>> &g_arr, const vector<Vector> &col3_nUvthsq_fM, const vector<Vector> &col3_ccf_h0_fM, const int &hMgM_analytic_op, vector<Vector> &evol_coeff_arr, vector<Vector> &evol_correction, vector<vector<Vector>> &sum_consv)
{
	assert(false && "An abstract function was called.");
}

void Integration::col_electron_vol_nl(const ElementX &ti, const ElementV &tj, const int &sp_a_id, const vector<vector<int>> &sp_ab_ion_col_flag_arr, const Vector &coeff, const vector<vector<Vector>> &h_arr, const vector<vector<Vector>> &g_arr, const vector<Vector> &col3_nUvthsq_fM, const vector<Vector> &col3_ccf_h0_fM, const int &hMgM_analytic_op, vector<Vector> &evol_coeff_arr, vector<Vector> &evol_correction, vector<vector<Vector>> &sum_consv)
{
	assert(false && "An abstract function was called.");
}
void Integration::col_ions_vol_fMa_nl(const ElementX &ti, const ElementV &tj, const int &sp_a_id, const vector<vector<int>> &sp_ab_ion_col_flag_arr, const Vector &fMa_coeff, const Vector &del_fa_coeff, const vector<vector<Vector>> &h_arr, const vector<vector<Vector>> &g_arr, const vector<Vector> &col3_nUvthsq_fM, const vector<Vector> &col3_ccf_h0_fM, const int &hMgM_analytic_op, vector<Vector> &evol_coeff_arr, vector<Vector> &evol_correction, vector<vector<Vector>> &sum_consv)
{
	assert(false && "An abstract function was called.");
}

void Integration::col_ions_inner_edge_nl(const ElementX &ti, const EdgeV &te, const ElementV &tj0, const ElementV &tj1, const int &i_op, const vector<vector<int>> &sp_ab_ion_col_flag_arr, const Vector &coeff_part2, const Vector &fhat_part, const vector<vector<Vector>> &gamma_aa_h_arr, const vector<vector<Vector>> &gamma_aa_g_arr, const vector<Vector> &col3_nUvthsq_fM, const vector<Vector> &col3_ccf_h0_fM, const int &hMgM_analytic_op, vector<Vector> &tmp_evol_coeff_arr, vector<Vector> &tmp_evol_correction)
{
	assert(false && "An abstract function was called.");
}
		
void Integration::col_electron_inner_edge_nl(const ElementX &ti, const EdgeV &te, const ElementV &tj0, const ElementV &tj1, const int &i_op, const int &sp_a_id, const vector<vector<int>> &sp_ab_ion_col_flag_arr, const Vector &coeff_part2, const Vector &fhat_part, const vector<vector<Vector>> &gamma_aa_h_arr, const vector<vector<Vector>> &gamma_aa_g_arr, const vector<Vector> &col3_nUvthsq_fM, const vector<Vector> &col3_ccf_h0_fM, const int &hMgM_analytic_op, vector<Vector> &tmp_evol_coeff_arr, vector<Vector> &tmp_evol_correction)
{
	assert(false && "An abstract function was called.");
}
void Integration::col_ions_inner_edge_fMa_nl(const ElementX &ti, const EdgeV &te, const ElementV &tj0, const ElementV &tj1, const int &i_op, const int &sp_a_id, const vector<vector<int>> &sp_ab_ion_col_flag_arr, const Vector &fMa_coeff_part2, const Vector &del_fa_coeff_part2, const Vector &fhat_part, const vector<vector<Vector>> &gamma_aa_h_arr, const vector<vector<Vector>> &gamma_aa_g_arr, const vector<Vector> &col3_nUvthsq_fM, const vector<Vector> &col3_ccf_h0_fM, const int &hMgM_analytic_op, vector<Vector> &tmp_evol_coeff_arr, vector<Vector> &tmp_evol_correction)
{
	assert(false && "An abstract function was called.");
}

Vector Integration::col_edge_bd_mat_cal_Vec(const ElementX &, const EdgeV &, const ElementV &, const int &, const int &, const int &, const Vector &, const Vector &, const Vector &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
Vector Integration::col_edge_bd_mat_cal_Vec_test_ptl_col(const ElementX &, const EdgeV &, const ElementV &, const int &, const int &, const int &, const Vector &, const Vector &, const Vector &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
Vector Integration::col_edge_bd_mat_cal_Vec_nl(const ElementX &, const EdgeV &, const ElementV &, const int &, const int &, const int &, const int &, const Vector &, const Vector &, const Vector &, const Vector &, const Vector &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::col_ions_boundary_edge_nl(const ElementX &ti, const EdgeV &te, const ElementV &tj, const int &ele_index, const int &i_op, const int &sp_a_id, const vector<vector<int>> &sp_ab_ion_col_flag_arr, const Vector &coeff, const Vector &fhat_coeff, const vector<vector<Vector>> &gamma_aa_h_arr, const vector<vector<Vector>> &gamma_aa_g_arr, const vector<Vector> &col3_nUvthsq_fM, const vector<Vector> &col3_ccf_h0_fM, const int &hMgM_analytic_op, vector<Vector> &tmp_evol_coeff_arr, vector<Vector> &tmp_evol_correction, vector<vector<Vector>> &sum_consv)
{
	assert(false && "An abstract function was called.");
}
void Integration::col_electron_boundary_edge_nl(const ElementX &ti, const EdgeV &te, const ElementV &tj, const int &ele_index, const int &i_op, const int &sp_a_id, const vector<vector<int>> &sp_ab_ion_col_flag_arr, const Vector &coeff, const Vector &fhat_coeff, const vector<vector<Vector>> &gamma_aa_h_arr, const vector<vector<Vector>> &gamma_aa_g_arr, const vector<Vector> &col3_nUvthsq_fM, const vector<Vector> &col3_ccf_h0_fM, const int &hMgM_analytic_op, vector<Vector> &tmp_evol_coeff_arr, vector<Vector> &tmp_evol_correction, vector<vector<Vector>> &sum_consv)
{
	assert(false && "An abstract function was called.");
}
void Integration::col3_bc_mat_cal(const ElementX &, const ElementV &, const int &, const int &, const vector<value_type> &, const vector<value_type> &, const int &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::col3_hg_vol_mat_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &, vector<value_type> &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::col3_hg_b_to_hg_a_vol_mat_cal(const ElementX &, const ElementV &, const int &, const int &, vector<value_type> &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::col3_small_hg_b_to_hg_a_vol_mat_cal(const ElementX &, const ElementV &, const int &, const int &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
#ifdef SEO_COL_IMPLICIT 
void Integration::col3_slow_moment_to_fast_hg_cal(const ElementX &, const ElementV &, const int &, const int &, const value_type &, const int &, Vector &)
{
	assert(false && "An abstract function was called.");
}
void Integration::col3_slow_moment_to_fast_hg_adj_n_cal(const ElementX &, const ElementV &, const int &, const int &, const value_type &, const int &, Vector &, int &)
{
	assert(false && "An abstract function was called.");
}
#endif

void Integration::col3_small_fb_to_big_h_a_source_mat_cal(const ElementX &, const ElementV &, const int &, const int &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}
void Integration::col3_edge_bd_hg_stiffness(const ElementX &, const EdgeV &, const int &, const int &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
} 
void Integration::col_nUT_spatial_mat_cal(const ElementX &, const ElementV &, const int &, const int, const Vector &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col_diag_test_cal(const ElementX &, const ElementV &, const int &, const int &, const Vector &, const Vector &, const Vector &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}

Vector Integration::col_fM_coeff_vol_cal(const ElementX &, const ElementV &, const int &, const Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col_fM_hg_vol_source_cal(const ElementX &, const ElementV &, const int &, const Vector &, Vector &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_hg_to_Qv_vol_cal(const ElementX &, const ElementV &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_hMgM_to_Qv_vol_cal(const int &, const ElementX &, const ElementV &, const Vector &, const int &, Vector &)
{
	assert(false && "An abstract function was called.");
}
		
void Integration::col3_hMgM_to_Qv_vol_adj_n_cal(const int &, const ElementX &, const ElementV &, const Vector &, const int &, Vector &, int &, value_type)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_hg_to_Qx_cal(const ElementX &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_hg_S_LF_Vmat_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &, Vector &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_hMgM_S_Vmat_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_hMgM_S_Vmat_adj_n_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &, Vector &, int &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_g_S_Vmat_bd_cal(const ElementX &,  const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_gM_S_Vmat_bd_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_gM_S_Vmat_bd_adj_n_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &, int &)
{
	assert(false && "An abstract function was called.");
}
void Integration::col3_h_flux_qd_cal(const ElementX &,  const EdgeV &, const int &, const Vector &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_hM_flux_qd_cal(const int &, const ElementX &, const EdgeV &, const int &, const Vector &, const Vector &, const int &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_f_to_consv_vol_cal(const ElementX &, const ElementV &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_f_to_mom_consv_surf_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_f_to_en_consv_inner_surf_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &)
{
	assert(false && "An abstract function was called.");
}

void Integration::col3_f_to_en_consv_bd_surf_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &)
{
	assert(false && "An abstract function was called.");
}


#endif

void Integration::init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const EQ_reader &, const int)
{
	assert(false && "An abstract function was called.");
}

void Integration::source_cell_diag_mat_cal(const ElementX &, const ElementV &, vector<value_type> &)
{
	assert(false && "An abstract function was called.");
}

