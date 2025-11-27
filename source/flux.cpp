#include "sokuri.h"
#include "flux.h"

Vector cross(const Vector &A, const Vector &B)
{
	Vector C(3);
	C[0] = A[1]*B[2] - A[2]*B[1];
	C[1] = A[2]*B[0] - A[0]*B[2];
	C[2] = A[0]*B[1] - A[1]*B[0];

	return C;
}


Vector cross_2d(const Vector &A, const Vector &B)
{
	Vector C(2);

	C[0] = A[1]*B[2] - A[2]*B[1];
	C[1] = A[2]*B[0] - A[0]*B[2];

	return C;
}



Flux::Flux(const Sokuri &sokuri, value_type eps, int sp_id)
{
	assert(sokuri.eq_reader != nullptr);
	assert(sokuri.species_data != nullptr);

	init(*sokuri.eq_reader, *sokuri.species_data, eps);

	assert(sokuri.quadrature_arr[sp_id] != nullptr);
	quadrature = sokuri.quadrature_arr[sp_id];

	assert(sokuri.mesh_arr[sp_id] != nullptr);
	mesh = sokuri.mesh_arr[sp_id];


	const int nx = mesh->size<ElementX>();
	B_cen.resize(nx, 0.0);
	n_cen.resize(nx, 0.0);
	T_cen.resize(nx, 0.0);
	U_cen.resize(nx, 0.0);
	for(int ix = 0; ix < nx; ix++)
	{
		auto tmp_vec = mesh->get_node_element(ElementX(ix));

		value_type R = 0.0, Z = 0.0;
		for (int i = 0; i < 3; i++) 
		{
			R += tmp_vec[i][0]/3.0;
			Z += tmp_vec[i][1]/3.0;
		}

		n_cen[ix] = species_data->prof2D_input(R, Z, sp_id, 0);
		T_cen[ix] = species_data->prof2D_input(R, Z, sp_id, 1);

		U_cen[ix] = bt_sign*species_data->prof2D_input(R, Z, sp_id, 2);
		B_cen[ix] = eq_reader->b_interpol(R, Z);
	}

}

void Flux::init(const EQ_reader &_eq_reader, const Species_data &_species_data, value_type _eps)
{
	eps = _eps;
	eq_reader = &_eq_reader;
	species_data = &_species_data;

	rh0n = eq_reader->ph_const.get_property(Ph_const::rh0n);
	B0 = eq_reader->ph_const.get_property(Ph_const::B0);
	R0 = eq_reader->ph_const.get_property(Ph_const::R0);
	vti0 = eq_reader->ph_const.get_property(Ph_const::vti0);
	bp_sign = eq_reader->get_property(EQ_reader::bp_sign);
	bt_sign = eq_reader->get_property(EQ_reader::bt_sign);
}


Vector Flux::get_B(const Point4 &p) const
{
	return get_B(Point2{p[0], p[1]});
}


Vector Flux::get_B(const Point2 &p) const
{
	assert(eq_reader != nullptr);

	value_type R = p[0], z = p[1];

	Vector B(3);
	eq_reader->bvec_interpol(R, z, B[0], B[1], B[2]);

	return B;
}


Vector Flux::get_Ds(const Point4 &p) const
{
	return get_Ds(Point2{p[0], p[1]});
}


Vector Flux::get_Ds(const Point2 &p) const
{
	//p = {R, z}
	value_type pv[12];

	eq_reader->intp_dpsi_RZ(p[0], p[1], pv);
	eq_reader->intp_dI_RZ(p[0], p[1], pv+6);

	//Swap Rz and z**2
	swap(pv[4], pv[5]);
	swap(pv[10], pv[11]);

	Vector Ds(12);

	copy(pv, pv+12, Ds.data());

	return Ds;
}


value_type Flux::get_psi(const Point4 &p) const
{
	Point2 q = p.segment(0,2);
	return get_psi(q);
}


value_type Flux::get_psi(const Point2 &p) const
{
	return eq_reader->psi_interpol(p[0], p[1]);
}

Vector Flux::get_b(const Point4 &p) const
{
	return get_b(Point2{p[0], p[1]});
}


Vector Flux::get_b(const Point2 &p) const
{ 
	auto B = get_B(p);

	B /= B.norm();

	return B;
}

Vector Flux::get_weighted_curl_b(const Point4 &p) const
{
	return get_weighted_curl_b(Point2{p[0], p[1]});
}


Vector Flux::get_weighted_curl_b(const Point2 &p) const
{
	auto B = get_B(p);

	return get_weighted_curl_b(p, B);
}


Vector Flux::get_weighted_curl_b(const Point4 &p, const Vector &B) const
{
	return get_weighted_curl_b(Point2{p[0], p[1]}, B);
}


Vector Flux::get_weighted_curl_b(const Point2 &p, const Vector &B) const
{
	auto R = p[0];
	auto lBl = B.norm();
	auto b = B / lBl;

	auto Ds = get_Ds(p);
	auto Vk = get_Vk(b);

	Vector curl_b(3);

	for(int i = 0; i < Vk.size(); i++)
		curl_b[i] = Vk[i].dot(Ds);

	curl_b[1] += B[2];

	return curl_b;
}


vector<Vector> Flux::get_Vk(const Vector &b) const
{
	vector<Vector> Vk(3, Vector::Zero(12));

	Vk[0][4] = -bp_sign*b[1]*b[2];
	Vk[0][5] = bp_sign*b[2]*b[0];
	Vk[0][8] = bt_sign*(b[2]*b[2]-1.0);

	Vk[1][3] = bp_sign*b[1]*b[2];
	Vk[1][4] = -bp_sign*b[2]*b[0];
	Vk[1][7] = bt_sign*(1.0-b[2]*b[2]);

	Vk[2][3] = bp_sign*(1.0-b[1]*b[1]);
	Vk[2][4] = bp_sign*2.0*b[1]*b[0];
	Vk[2][5] = bp_sign*(1.0-b[0]*b[0]);
	Vk[2][7] = bt_sign*b[1]*b[2];
	Vk[2][8] = -bt_sign*b[0]*b[2];

	return Vk;
}


Vector Flux::get_Bs(const Point4 &p, const Vector &B, int sp) const
{
	assert(eq_reader != nullptr && species_data != nullptr);

	/* sp 0 = electron
	 * sp 1 = ion1
	 * sp 2 = ion2 ...*/

	//p = {R, z, v, u}
	auto R = p[0], z = p[1], v = p[2];
	auto weighted_curl_b = get_weighted_curl_b(p, B);
	auto lBl = B.norm();

	auto Ms = species_data->normalized_mass(sp);
	//if (sp == 0) Ms = 0.0;
	auto Zs = species_data->normalized_charge(sp);
	auto ccf = (v*Ms) / (Zs*lBl*R);
	auto Bs = B + rh0n*ccf*weighted_curl_b;

	return Bs;
}

Vector Flux::get_mod_curl_b(const Point2 &p2_in) const
{
	assert(eq_reader != nullptr && species_data != nullptr);

	Point4 p;
	p[0] = p2_in[0];
	p[1] = p2_in[1];

	auto B = get_B(p);
	auto R = p[0], z = p[1], v = p[2];
	auto weighted_curl_b = get_weighted_curl_b(p, B);
	auto lBl = B.norm();

	auto ccf = 1.0 / (lBl*R);
	auto mod_curl_b = rh0n*ccf*weighted_curl_b;


	return mod_curl_b;
}

vector<value_type> Flux::get_b_dot_Bs_vec(const Point4 &p, const Vector &B, int sp) const
{
	assert(eq_reader != nullptr && species_data != nullptr);

	/* sp 0 = electron
	 * sp 1 = ion1
	 * sp 2 = ion2 ...*/

	vector<value_type> get_b_dot_Bs_arr(2);

	//p = {R, z, v, u}
	auto R = p[0], z = p[1], v = p[2];
	auto lBl = B.norm();
	auto weighted_curl_b_dot_b = get_weighted_curl_b(p, B).dot(B)/lBl;


	auto Ms = species_data->normalized_mass(sp);
	//if (sp == 0) Ms = 0.0;
	auto Zs = species_data->normalized_charge(sp);
	auto ccf0 = (Ms) / (Zs*lBl*R);
	auto ccf = v*ccf0;

	get_b_dot_Bs_arr[0] = lBl + rh0n*ccf*weighted_curl_b_dot_b;
	get_b_dot_Bs_arr[1] = rh0n*ccf0*weighted_curl_b_dot_b;

	return get_b_dot_Bs_arr;
}


Vector Flux::get_Bs(const Point4 &p, int sp) const
{
	return get_Bs(p, get_B(p), sp);
}



value_type Flux::get_b_dot_Bs(const Point4 &p, int sp) const
{
	auto b = get_b(p);
	auto val = get_Bs(p, get_B(p), sp).dot(b);
	//if(sp == 0) val = get_B(p).norm();

	return val;
}


Vector Flux::get_grad_B(const Point2 &p) const
{
	return get_grad_B(p, get_B(p));
}


Vector Flux::get_grad_B(const Point4 &p) const
{
	return get_grad_B(p, get_B(p));
}



Vector Flux::get_grad_B(const Point4 &p, const Vector &B) const
{
	return get_grad_B(Point2{p[0], p[1]}, B);
}


Vector Flux::get_grad_B(const Point2 &p, const Vector &B) const
{
	auto R = p[0];
	auto lBl = B.norm();
	auto b = B/lBl;
	auto Ds = get_Ds(p);
	Vector grad_B(3);
	Vector vk(12);

	//dBdR
	vk.setZero();
	vk[3] = -bp_sign*b[1];
	vk[4] = bp_sign*b[0];
	vk[7] = bt_sign*b[2];

	grad_B[0] = (-lBl + vk.dot(Ds))/R;

	//dBdz
	vk.setZero();
	vk[4] = -bp_sign*b[1];
	vk[5] = bp_sign*b[0];
	vk[8] = bt_sign*b[2];

	grad_B[1] = vk.dot(Ds)/R;

	//dBdvarphi
	grad_B[2] = 0.0;

	return grad_B;
}


Vector Flux::get_U(const Point4 &p, int sp) const
{
	/*
	   Vector U = Vector::Zero(4);
	//	U[1] = 1.0;
	//	U[0] = 4.0*p[1];

	return U;
	*/

	return get_U0(p, sp) + get_U1(p, sp);
}


value_type Flux::get_divU(const Point4 &p, int sp) const
{
	value_type eps = 1.0e-8;

	value_type div = 0.0;

	for(int d = 0; d < 4; d++)
	{
		Point4 a = p, b = p;
		b[d] += eps;
		a[d] -= eps;

		div += (get_U(b, sp)[d] - get_U(a, sp)[d]) / (2.0*eps);
	}

	return div;
}



Vector Flux::get_U1(const Point4 &p, int sp) const
{
	return Vector::Zero(4);
}


Vector Flux::get_U0(const Point4 &p, int sp) const
{
	auto Ms = species_data->normalized_mass(sp);
	auto Zs = species_data->normalized_charge(sp);

	auto R = p[0], z = p[1], v = p[2], u = p[3];
	auto mu = 0.5*Ms*u*u;
	auto B = get_B(p);
	auto Bs = get_Bs(p, B, sp);
	auto lBl = B.norm();
	auto b = B / lBl;
	auto Bs_ll = get_b_dot_Bs(p, sp);

	auto grad_B = get_grad_B(p, B);

	value_type ccf = mu/Zs;

	Vector bXgrad_B = -cross(b/Bs_ll, grad_B);
	bXgrad_B[2] /= R;	//Who are you?

	auto tmp = v*Bs/Bs_ll + rh0n*ccf*bXgrad_B;

	Vector U(4);
	U.segment(0, 2) = tmp.segment(0, 2);
	U[2] = -0.5*u*u*(Bs/Bs_ll).dot(grad_B);
	U[3] = 0.0;

	return U;
}

value_type Flux::get_n(const Point4 &pnt, int sp) const
{
	return get_n(Point2{pnt.segment(0, 2)}, sp);
}


value_type Flux::get_n(const Point2 &pnt, int sp) const
{
	return species_data->prof2D_input(pnt[0], pnt[1], sp, 0);
}


Vector Flux::get_grad_n(const Point4 &pnt, int sp) const
{
	return get_grad_n(Point2{pnt.segment(0, 2)}, sp);
}


Vector Flux::get_grad_n(const Point2 &pnt, int sp) const
{
	auto grad_n = species_data->prof2D_input_dfn(pnt[0], pnt[1], sp, 0);

	//grad = [1, dR, dz, dRdR, dzdz, dRdz]

	Vector tmp = operator&(Vector{grad_n.segment(1, 2)}, Vector::Zero(2));
	return tmp;
}


value_type Flux::get_T(const Point2 &pnt, int sp) const
{
	return species_data->prof2D_input(pnt[0], pnt[1], sp, 1);
}


Vector Flux::get_grad_T(const Point2 &pnt, int sp) const
{
	auto grad_T = species_data->prof2D_input_dfn(pnt[0], pnt[1], sp, 1);

	//grad = [1, dR, dz, dRdR, dzdz, dRdz]

	Vector tmp = operator&(Vector{grad_T.segment(1, 2)}, Vector::Zero(2));
	// return Vector{grad_T.segment(1, 2)} & Vector::Zero(2);
	return tmp;
}


Vector Flux::get_grad_T(const Point4 &pnt, int sp) const
{
	return get_grad_T(Point2{pnt.segment(0, 2)}, sp);
}


value_type Flux::get_f0(const Point4 &pnt, const int sp, const int flag) const
{
	//	 sp 0 = electron
	//	 sp 1 = ion1
	//	 sp 2 = ion2 ...

	//	 flag -1 = return 1.0;
	//	 flag 0= Maxwellian with no U_ll
	//	 flag 1= Maxwellian with input U_ll
	//	 flag 2= canoniacal Maxwellian with no U_ll modification
	//	 flag 3= canonical Maxwellian with U_ll modification

	if(flag == -1)
		return 1.0;

	auto R = pnt[0], Z = pnt[1], vll = pnt[2], u_perp = pnt[3];
	auto B = eq_reader->b_interpol(R, Z);

	value_type Ms = species_data->normalized_mass(sp);
	value_type ccf = 0.5*Ms*M_1_PI*sqrt(0.5*Ms*M_1_PI);

	value_type n_val, T_val, U_val = 0.0;

	switch(flag)
	{
		case -1:
			break;

		case 0:
			n_val = species_data->prof2D_input(R, Z, sp, 0);
			T_val = species_data->prof2D_input(R, Z, sp, 1);
			break;

		case 1:
			{
				n_val = species_data->prof2D_input(R, Z, sp, 0);
				T_val = species_data->prof2D_input(R, Z, sp, 1);

				value_type BT_SIGN = eq_reader->get_property(EQ_reader::bt_sign);
				U_val = BT_SIGN*species_data->prof2D_input(R, Z, sp, 2);
				break;
			}

		case 2:
			{
				if (sp == 0) 
				{
					cerr << "Canonical maxwellian for electron is not supported." << endl;
					exit(EXIT_FAILURE);
				}

				value_type I_val = eq_reader->I_interpol(R,Z);
				value_type Zs = species_data->normalized_charge(sp);
				value_type BPT_SIGN = eq_reader->get_property(EQ_reader::bpt_sign);
				value_type rh0n = eq_reader->ph_const.get_property(Ph_const::rh0n);
				value_type ccf2 = BPT_SIGN*rh0n*(Ms/Zs);

				value_type psi_h = eq_reader->psi_interpol(R,Z);
				psi_h -= ccf2*I_val/B*vll;

				value_type En = 0.5*Ms*(vll*vll + u_perp*u_perp*B);
				value_type E_ps = En - 0.5*Ms*u_perp*u_perp;

				n_val = species_data->prof1D_input(psi_h, sp, 0);
				T_val = species_data->prof1D_input(psi_h, sp, 1);
				value_type BT_SIGN = eq_reader->get_property(EQ_reader::bt_sign);
				U_val = BT_SIGN*species_data->prof1D_input(psi_h, sp, 2);
				break;
			}

		case 3:
			{
				if (sp == 0) 
				{
					cerr << "Canonical maxwellian for electron is not supported." << endl;
					exit(EXIT_FAILURE);
				}

				value_type I_val = eq_reader->I_interpol(R,Z);
				value_type Zs = species_data->normalized_charge(sp);
				value_type BPT_SIGN = eq_reader->get_property(EQ_reader::bpt_sign);
				value_type rh0n = eq_reader->ph_const.get_property(Ph_const::rh0n);
				value_type ccf2 = BPT_SIGN*rh0n*(Ms/Zs);

				auto psi_h = eq_reader->psi_interpol(R,Z);
				psi_h -= ccf2*I_val/B*vll;

				value_type En = 0.5*Ms*(vll*vll + u_perp*u_perp*B);
				value_type E_ps = En - 0.5*Ms*u_perp*u_perp;

				if (E_ps > 0.0)
					psi_h += vll > 0.0 ? ccf2*sqrt(2.0*E_ps/Ms) : -ccf2*sqrt(2.0*E_ps/Ms);

				n_val = species_data->prof1D_input(psi_h, sp, 0);
				T_val = species_data->prof1D_input(psi_h, sp, 1);
				value_type BT_SIGN = eq_reader->get_property(EQ_reader::bt_sign);
				U_val = BT_SIGN*species_data->prof1D_input(psi_h, sp, 2);
			}
			break;

		default:
			cerr << "Unknown flag in F0 : " << flag << endl;
			exit(EXIT_FAILURE);
			break;
	}

	value_type f_val = ccf*n_val/(T_val*sqrt(T_val))*exp(-0.5*Ms*((vll-U_val)*(vll-U_val)+u_perp*u_perp*B)/T_val);

	int tmp_anisotro_T_op = species_data->aniso_T_op(sp);
	if(tmp_anisotro_T_op == 1)
	{
		//anisotropic T test
		value_type mod_T_alpha = species_data->aniso_T_perp_T_para_ratio(sp);;
		//mod_T_alpha = 1.0;
		f_val = ccf*n_val/mod_T_alpha/(T_val*sqrt(T_val))*exp(-0.5*Ms*((vll-U_val)*(vll-U_val)+u_perp*u_perp*B/mod_T_alpha)/T_val);
	}
	return f_val;
}

value_type Flux::get_f0_nUT(const Point4 &pnt, const int &sp, const int &flag, const value_type n_in, const value_type T_in, const value_type U_in) const
{

	if(flag == -1)
		return 1.0;

	auto R = pnt[0], Z = pnt[1], vll = pnt[2], u_perp = pnt[3];
	auto B = eq_reader->b_interpol(R, Z);

	value_type Ms = species_data->normalized_mass(sp);
	value_type ccf = 0.5*Ms*M_1_PI*sqrt(0.5*Ms*M_1_PI);

	value_type n_val, T_val, U_val = 0.0;

	n_val = n_in;
	T_val = T_in;
	U_val = U_in;

	return ccf*n_val/(T_val*sqrt(T_val))*exp(-0.5*Ms*((vll-U_val)*(vll-U_val)+u_perp*u_perp*B)/T_val);

}


vector<value_type> Flux::get_fcm_SE_vec(const Point4 &pnt, const int ix, const int sp, const int flag) const
{
	vector<value_type> f_val_arr(3);
	value_type f_val;

	if (flag == -1 || ix < 0) 
	{
		f_val_arr[0] = 1.0;
		f_val_arr[1] = 0.0;
		f_val_arr[2] = 0.0;


		return f_val_arr;
	}

	auto vll = pnt[2], u_perp = pnt[3];
	auto B = B_cen[ix];

	value_type Ms = species_data->normalized_mass(sp);
	value_type ccf = 0.5*Ms*M_1_PI*sqrt(0.5*Ms*M_1_PI);

	value_type R = pnt[0], Z = pnt[1];
	value_type n_val, T_val, T_cen_val, U_cen_val;
	n_val = species_data->prof2D_input(R, Z, sp, 0);
	T_val = species_data->prof2D_input(R, Z, sp, 1);
	//T_val = T_cen[ix];
	T_cen_val = T_cen[ix];
	U_cen_val = U_cen[ix];

	if(flag == 0) U_cen_val = 0.0;

	if(flag == 2 || flag == 3) 
	{
		cout << "get_fcm_SE_vec is not compatilbe with canonical Maxwellian weight. " << endl;

		abort();
	}
	f_val = ccf*n_val/(T_val*sqrt(T_val))*exp(-0.5*Ms*((vll-U_cen_val)*(vll-U_cen_val)+u_perp*u_perp*B)/T_cen_val);


	f_val_arr[0] = f_val;
	f_val_arr[1] = -f_val*Ms*(vll-U_cen_val)/T_cen_val;
	f_val_arr[2] = -f_val*Ms*u_perp*B/T_cen_val;
	//f_val_arr[1] = -Ms*(vll-U_val)/T_val;
	//f_val_arr[2] = -Ms*u_perp*B/T_val;

	//cout << R << ", " << Z << ", " << vll << ", " << u_perp << ", " << B << endl;
	//cout << eq.psi_ov_psix_interpol(R,Z) << ", " << eq.psi_interpol(R,Z) << ", " <<psi_h << ", " << n_val << ", " << T_val << ", " << U_val << endl;

	return f_val_arr;
}


value_type Flux::get_fcm_SE_tot(const Point4 &pnt, const int ix, const int sp, const int flag) const
{
	//	 sp 0 = electron
	//	 sp 1 = ion1
	//	 sp 2 = ion2 ...

	//	 flag -1 = return 1.0;
	//	 flag 0= Maxwellian with no U_ll
	//	 flag 1= Maxwellian with input U_ll
	//	 flag 2= canoniacal Maxwellian with no U_ll modification
	//	 flag 3= canonical Maxwellian with U_ll modification

	if(flag == -1)
		return 1.0;

	if(ix < 0)
		return 1.0;


	auto vll = pnt[2], u_perp = pnt[3];
	auto B = B_cen[ix];

	value_type Ms = species_data->normalized_mass(sp);
	value_type ccf = 0.5*Ms*M_1_PI*sqrt(0.5*Ms*M_1_PI);

	value_type R = pnt[0], Z = pnt[1];
	value_type n_val, T_val, T_cen_val, U_cen_val;
	n_val = species_data->prof2D_input(R, Z, sp, 0);
	T_val = species_data->prof2D_input(R, Z, sp, 1);
	//T_val = T_cen[ix];
	T_cen_val = T_cen[ix];
	U_cen_val = U_cen[ix];

	//value_type n_val = n_cen[ix], T_val = T_cen[ix], U_val = U_cen[ix];

	if(flag == 0) U_cen_val = 0.0;

	if(flag == 2 || flag == 3) 
	{
		cout << "get_fcm_SE_tot is not compatilbe with canonical Maxwellian weight. " << endl;

		abort();
	}
	return ccf*n_val/(T_val*sqrt(T_val))*exp(-0.5*Ms*((vll-U_cen_val)*(vll-U_cen_val)+u_perp*u_perp*B)/T_cen_val);
}

value_type Flux::get_fcm_SE_X(const Point4 &pnt, const int ix, const int sp, const int flag) const
{
	//	 sp 0 = electron
	//	 sp 1 = ion1
	//	 sp 2 = ion2 ...

	//	 flag -1 = return 1.0;
	//	 flag 0= Maxwellian with no U_ll
	//	 flag 1= Maxwellian with input U_ll
	//	 flag 2= canoniacal Maxwellian with no U_ll modification
	//	 flag 3= canonical Maxwellian with U_ll modification

	if(flag == -1)
		return 1.0;

	if(ix < 0)
		return 1.0;

	value_type Ms = species_data->normalized_mass(sp);
	value_type ccf = 0.5*Ms*M_1_PI*sqrt(0.5*Ms*M_1_PI);


	value_type R = pnt[0], Z = pnt[1];
	value_type n_val, T_val, U_val;
	n_val = species_data->prof2D_input(R, Z, sp, 0);
	T_val = species_data->prof2D_input(R, Z, sp, 1);
	//T_val = T_cen[ix];
	U_val = bt_sign*species_data->prof2D_input(R, Z, sp, 2);


	//value_type n_val = n_cen[ix], T_val = T_cen[ix], U_val = U_cen[ix];

	if(flag == 0) U_val = 0.0;

	if(flag == 2 || flag == 3) 
	{
		cout << "get_fcm_SE_tot is not compatilbe with canonical Maxwellian weight. " << endl;

		abort();
	}
	return ccf*n_val/(T_val*sqrt(T_val));
}


value_type Flux::get_fcm_SE_V(const Point4 &pnt, const int ix, const int sp, const int flag) const
{
	//	 sp 0 = electron
	//	 sp 1 = ion1
	//	 sp 2 = ion2 ...

	//	 flag -1 = return 1.0;
	//	 flag 0= Maxwellian with no U_ll
	//	 flag 1= Maxwellian with input U_ll
	//	 flag 2= canoniacal Maxwellian with no U_ll modification
	//	 flag 3= canonical Maxwellian with U_ll modification

	if(flag == -1)
		return 1.0;

	if(ix < 0)
		return 1.0;


	auto vll = pnt[2], u_perp = pnt[3];
	auto B = B_cen[ix];

	value_type Ms = species_data->normalized_mass(sp);
	value_type ccf = 0.5*Ms*M_1_PI*sqrt(0.5*Ms*M_1_PI);

	value_type T_val = T_cen[ix], U_val = U_cen[ix];

	if(flag == 0) U_val = 0.0;

	if(flag == 2 || flag == 3) 
	{
		cout << "get_fcm_SE_tot is not compatilbe with canonical Maxwellian weight. " << endl;

		abort();
	}
	return exp(-0.5*Ms*((vll-U_val)*(vll-U_val)+u_perp*u_perp*B)/T_val);
}

