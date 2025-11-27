#pragma once

#include <cmath>
#include <vector>
#include <tuple>

#include "type.h"
#include "quadrature.h"
#include "eq_reader.h"
#include "species.h"
#include "quadrature_dg.h"
#include "flux.h"

using namespace std;
using namespace gyrokinetic_type;

class Basis
{
	public:
		Basis() = default;
		virtual ~Basis() = default;

		virtual value_type operator()(const Point4 &, int, int = 0) const;
		virtual value_type operator()(const Point4 &, const Point4 &, const ElementX &, const ElementV &, int, int = 0) const;
		void init(const int, const Mesh &, const Species_data &, const EQ_reader &, const Flux &, const Quadrature &);
		value_type lookup(int s, int qid, int, int = 0) const;
		template<class T> Vector grad(const T&, int) const;
		template<class T> Vector grad(const T&, int, vector<value_type> &) const;
		template<class T> Vector grad(const T&, const Point4 &, const ElementX &, const ElementV &, int) const;

		virtual Vector ind_val(const Point4 &, int, int) const;
		virtual value_type ind_val_S1(const Point4 &, int) const;
		template<class T> Vector average_grad(const T &, int) const;
		template<typename T> value_type jump(const T &, int) const;
		auto operator()(const vector<Point4> &, const vector<int> &) const;
		int get_dofx() const;
		int get_dofv() const;
		int get_dof() const;
		int get_basis_index(int k, int j) const { return dof * k + j; }
		void col_fhat_init_setup(int);
		Vector fhat_vbasis_Vec(const Point4 &p, int bd_type, int ele_index) const;

		int fhat_org_vbasis_dim_out() const {return fhat_org_v_basis_dim;}
		int fhat_valid_vbasis_dim_out() const {return fhat_valid_v_basis_dim;}
		int fhat_valid_to_org_i(const int &, const int &) const;


		int get_dof_vp_1d() const;
		int get_dof_u_1d() const;
		void col_fhat_setup(int);
		value_type col_spatial_basis(const Point4 &, int) const;
		value_type col_velocity_org_basis(const Point4 &, int) const;
		value_type col_velocity_fhat_basis(const Point4 &, int, int, int) const;
		void col_spatial_basis_vec(const Point2 &, vector<value_type> &) const;
		void col_velocity_basis_vec(const Point2 &, const value_type, const value_type, vector<value_type> &) const;
		void col_velocity_basis_vec_edge(const int, const Point2 &, const value_type, const value_type, vector<value_type> &) const;
		void col_velocity_fhat_basis_info(int, int, int, value_type, value_type, int &, value_type &) const;
		void col_velocity_fhat_basis_bd_info(int, int, value_type, value_type, int &, value_type &) const;

		int get_col_edge_type_num() const;
		void get_col_fhat_data(const int, int &, vector<int> &, vector<int> &, vector<int> &, vector<int> &, vector<int> &) const;

		int get_dofv_fhat() const {return dofv_fhat;}
		int get_max_dof_v_1d_fhat() const {return max_dof_v_1d_fhat;}
		
	protected:
		int dofx = 0, dofv = 0, dof = 0;
		int dof_S1 = 0 , dof_S1_x = 0;
		int dim = 0;
		int dop = 0;

		//Matrix make_p(const Point2 &) const;
		virtual Vector make_p(const Point2 &, int) const;
		Vector make_pn(const Point2 &, const Vector &) const;
		Matrix make_Dp(const Point2 &) const;
		virtual Matrix make_coeff() const;

		Matrix coeff;
		vector<Vector> normal;
		vector<value_type> lookup_table;

		const Mesh *mesh = nullptr;
		const EQ_reader *eq_reader = nullptr;
		const Species_data *species_data = nullptr;
		const Flux *flux = nullptr;
		const Quadrature *quadrature = nullptr;

		int sp_id;
		value_type Zs, Ms, BPT_SIGN, rh0n, ccf2;

		int fhat_v_basis_group_id_default = -1;
		int fhat_v_basis_id, fhat_org_v_basis_dim, fhat_valid_v_basis_dim;
		vector<vector<int>> fhat_valid_to_org_vbasis_i;

		value_type mesh_dvp, mesh_du;

		int dof_vp_1d, dof_u_1d;
		int op_divide_v_basis = 0;
		int col_edge_type_num = 2;
		vector<int> col_fhat_mat_num;
		vector<vector<int>> col_fhat_mat_half_size, col_fhat_org_mat_group, col_fhat_org_mat_index, col_fhat_new_mat_group, col_fhat_new_mat_index;
		vector<vector<int>> col_fhat_fhat_basis_info1;
		vector<vector<value_type>> col_fhat_new_basis_info2;

		int x_basis_group_id = -1;
		int v_basis_group_id = -1;
		int fhat_v_basis_group_id = -1;

		int dofv_fhat, max_dof_v_1d_fhat;
		void col_fhat_dofv_setup();
};


#include "basis.hpp"
