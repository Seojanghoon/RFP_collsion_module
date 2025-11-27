#pragma once

#include <cmath>
#include <vector>
#include <iostream>
#include <functional>
#include <cstdarg>
#include <numeric>
#include <unordered_map>
#include <iterator>

//#include <Eigen/Dense>
//#include <Eigen/Sparse>
//#include <Eigen/PardisoSupport>

//#include <linalg/linalg.hpp>

#include "type.h"
#include "basis.h"
#include "functor.h"
#include "eq_reader.h"
#include "transform.h"
#include "index.h"
#include "flux.h"
#include "quadrature.h"
#include "indexer.h"
#include "debug.h"
#include "efld_data.h"
#include "elliptic_integral.hpp"
using namespace std;
using namespace gyrokinetic_type;



class Integration
{
	public:
		Integration() = default;
		virtual ~Integration() = default;
		virtual void ME_mat_cal_f_init(const ElementX &, const ElementV &, const int &, vector<value_type> &, vector<value_type> &, vector<value_type> &);


		virtual void init(const Mesh &, const Basis &, const Quadrature &, const Flux &, const Species_data &, const EQ_reader &, const int, const int &, const value_type &);
		virtual void col3_init(const int &, const int &, const vector<int> &, const vector<int> &, const vector<value_type> &, const vector<value_type> &, const vector<value_type> &, const vector<value_type> &, const int &, const int &, const value_type &);
			
		virtual void col_cell_mat_cal(const ElementX &, const ElementV &, vector<value_type> &);
		virtual void col_cell_consv_mat_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &, const int &);


		virtual void col_fhat_mat_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, vector<value_type> &, vector<int> &);
		virtual void col_fhat_mat_cal_bd(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, vector<value_type> &, vector<int> &);
		virtual void col_f_dg_to_valid_fhat_mat_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &);
		virtual void col3_bc_mat_cal(const ElementX &, const ElementV &, const int &, const int &, const vector<value_type> &, const vector<value_type> &, const int &, vector<value_type> &);
		virtual void col3_hg_vol_mat_cal(const ElementX &, const ElementV &, const int &, vector<value_type> &, vector<value_type> &, vector<value_type> &);
		virtual void col3_hg_b_to_hg_a_vol_mat_cal(const ElementX &, const ElementV &, const int &, const int &, vector<value_type> &, vector<value_type> &);
		virtual void col3_small_hg_b_to_hg_a_vol_mat_cal(const ElementX &, const ElementV &, const int &, const int &, vector<value_type> &);
		virtual void col3_slow_moment_to_fast_hg_cal(const ElementX &, const ElementV &, const int &, const int &, const value_type &, const int &, Vector &);
		virtual void col3_slow_moment_to_fast_hg_adj_n_cal(const ElementX &, const ElementV &, const int &, const int &, const value_type &, const int &, Vector &, int &);

		virtual void col3_small_fb_to_big_h_a_source_mat_cal(const ElementX &, const ElementV &, const int &, const int &, vector<value_type> &);
		virtual void col3_edge_bd_hg_stiffness(const ElementX &, const EdgeV &, const int &, const int &, vector<value_type> &);
		virtual Vector col_fM_coeff_vol_cal(const ElementX &, const ElementV &, const int &, const Vector &);
		

		virtual void col3_hg_to_Qx_cal(const ElementX &, Vector &);
		virtual void col3_hg_to_Qv_vol_cal(const ElementX &, const ElementV &, Vector &);
		virtual void col3_hMgM_to_Qv_vol_cal(const int &, const ElementX &, const ElementV &, const Vector &, const int &, Vector &);
		virtual void col3_hMgM_to_Qv_vol_adj_n_cal(const int &, const ElementX &, const ElementV &, const Vector &, const int &, Vector &, int &, value_type);


		virtual void col3_hg_S_LF_Vmat_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &, Vector &, Vector &);
		virtual void col3_hMgM_S_Vmat_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &, Vector &);
		virtual void col3_hMgM_S_Vmat_adj_n_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &, Vector &, int &);

		virtual void col3_g_S_Vmat_bd_cal(const ElementX &,  const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &);
		virtual void col3_gM_S_Vmat_bd_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &);
		virtual void col3_gM_S_Vmat_bd_adj_n_cal(const int &, const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, const Vector &, const int &, Vector &, int &);

		virtual void col3_h_flux_qd_cal(const ElementX &,  const EdgeV &, const int &, const Vector &, Vector &);
		virtual void col3_hM_flux_qd_cal(const int &, const ElementX &, const EdgeV &, const int &, const Vector &, const Vector &, const int &, Vector &);
		virtual void col3_f_to_consv_vol_cal(const ElementX &, const ElementV &, Vector &);
		virtual void col3_f_to_mom_consv_surf_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, Vector &);
		virtual void col3_f_to_en_consv_inner_surf_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &);
		virtual void col3_f_to_en_consv_bd_surf_cal(const ElementX &, const EdgeV &, const ElementV &, const ElementV &, const int &, Vector &);

	protected:
		function<value_type(const Point4 &, int, int)> _f0;
		const Basis *basis = nullptr;
		const Functor *fcm = nullptr, *bfield = nullptr, *bsfield = nullptr, *f = nullptr, *f0 = nullptr;
		const Flux *flux = nullptr;
		const EQ_reader *eq_reader = nullptr;
		const Mesh *mesh = nullptr;
		const Quadrature *quadrature = nullptr;
		const Indexer *indexer = nullptr;

		const Efld_data *efld_data = nullptr;
		function<value_type(const Point2 &)> rhs_func;
		const Species_data *species_data = nullptr;
		int dof, dof_S1, dof_S1_x, dof_x, dof_v;
		int flag_f0, flag_fcm, sp_id;

		unordered_map<int, Matrix> lookup_table_C[2];
		unordered_map<int, Matrix> lookup_table_W;
		unordered_map<int, vector<Matrix>> lookup_table_V;
		clock_t tot_dc1=0, tot_dc2=0, tot_dc3=0, tot_dc4=0;
		clock_t tot_dc5=0, tot_dc6=0, tot_dc7=0;
		clock_t tmp_t1, tmp_t2;
};
