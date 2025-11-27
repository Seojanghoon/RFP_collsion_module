#pragma once

#include <algorithm>
#include <utility>
#include <cmath>
#include <unordered_map>

#include "config.h"
#include "eq_reader.h"
#include "species.h"
#include "functor.h"
#include "type.h"
#include "debug.h"
#include "bspline.h"
#include "index.h"
#include "quadrature.h"

#include <Eigen/Core>

#include <mpi.h>

Vector cross(const Vector &, const Vector &);
Vector cross_2d(const Vector &, const Vector &);

class Sokuri;

class Flux
{
	public:
		Flux() = default;
		Flux(const Sokuri &, value_type eps, int);
		~Flux() = default;

		void init(const EQ_reader &, const Species_data &, value_type eps);

		Vector get_B(const Point2 &) const;
		Vector get_B(const Point4 &) const;
		Vector get_grad_B(const Point2 &, const Vector &B) const;
		Vector get_grad_B(const Point4 &, const Vector &B) const;
		Vector get_grad_B(const Point2 &) const;
		Vector get_grad_B(const Point4 &) const;
		Vector get_Bs(const Point4 &, const Vector &, int = 0) const;
		Vector get_Bs(const Point4 &, int = 0) const;
		Vector get_mod_curl_b(const Point2 &) const;
		value_type get_b_dot_Bs(const Point4 &, int = 0) const;
		Vector get_U0(const Point4 &, int = 0) const;
		Vector get_U1(const Point4 &, int = 0) const;
		Vector get_U(const Point4 &, int = 0) const;
		value_type get_divU(const Point4 &, int = 0) const;

		value_type get_fcm_SE_tot(const Point4 &, int ix, int species, int flag) const;
		value_type get_fcm_SE_X(const Point4 &, int ix, int species, int flag) const;
		value_type get_fcm_SE_V(const Point4 &, int ix, int species, int flag) const;

		value_type get_f0(const Point4 &, int species = 0, int flag = 0) const;
		value_type get_psi(const Point4 &) const;
		value_type get_psi(const Point2 &) const;
		value_type get_n(const Point4 &, int species) const;
		value_type get_n(const Point2 &, int species) const;
		Vector get_grad_n(const Point2 &, int species) const;
		Vector get_grad_n(const Point4 &, int species) const;
		value_type get_T(const Point2 &, int species) const;
		value_type get_T(const Point4 &, int species) const;
		Vector get_grad_T(const Point2 &, int species) const;
		Vector get_grad_T(const Point4 &, int species) const;

		value_type get_f0_nUT(const Point4 &, const int &species, const int &flag, const value_type, const value_type, const value_type) const;
		vector<value_type> get_b_dot_Bs_vec(const Point4 &, const Vector &, int = 0) const;

		vector<value_type> get_fcm_SE_vec(const Point4 &, int ix, int species = 0, int flag = 0) const;
		Vector get_b(const Point2 &) const;

	protected:
		const EQ_reader *eq_reader = nullptr;
		const Species_data *species_data = nullptr;
		const Quadrature *quadrature = nullptr;
		const Mesh *mesh = nullptr;
		value_type rh0n, B0, R0, vti0;
		value_type bp_sign, bt_sign;
		value_type eps = 0.0;

		Vector get_Ds(const Point2 &) const;
		Vector get_Ds(const Point4 &) const;
		Vector get_b(const Point4 &) const;
		Vector get_weighted_curl_b(const Point2 &) const;
		Vector get_weighted_curl_b(const Point4 &) const;
		Vector get_weighted_curl_b(const Point2 &, const Vector &) const;
		Vector get_weighted_curl_b(const Point4 &, const Vector &) const;
		vector<Vector> get_Vk(const Vector &) const;


		vector<value_type>  B_cen, n_cen, T_cen, U_cen;
};

