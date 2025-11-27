#include "basis.h"


value_type Basis::operator()(const Point4 &, int, int) const
{
	assert(false && "This is a dummy function.");
}

value_type Basis::operator()(const Point4 &, const Point4 &, const ElementX &, const ElementV &, int, int) const
{
	assert(false && "This is a dummy function.");
}
		
void Basis::init(const int _sp_id, const Mesh &_mesh, const Species_data &_species_data, const EQ_reader &_eq_reader, const Flux &_flux, const Quadrature &_quadrature)
{
	mesh = &_mesh;
	species_data = &_species_data;
	eq_reader = &_eq_reader;
	flux = &_flux;
	quadrature = &_quadrature;

	sp_id = _sp_id;
	Zs = species_data->normalized_charge(sp_id); 
	Ms = species_data->normalized_mass(sp_id);

	BPT_SIGN = eq_reader->get_property(EQ_reader::bpt_sign);
	rh0n = eq_reader->ph_const.get_property(Ph_const::rh0n);
	ccf2 = BPT_SIGN*rh0n*(Ms/Zs);
}

Vector Basis::ind_val(const Point4 &, int, int) const
{
	assert(false && "This is a dummy function.");
}


value_type Basis::ind_val_S1(const Point4 &, int) const
{
	assert(false && "This is a dummy function.");
}


int Basis::get_dofx() const
{
	return dofx;
}


int Basis::get_dofv() const
{
	return dofv;
}


int Basis::get_dof() const
{
	return dof;
}


Vector Basis::make_p(const Point2 &, int) const
{
	assert(false && "This is a dummy function.");
}


Vector Basis::make_pn(const Point2 &pt, const Vector &normal) const
{
	auto Dp = make_Dp(pt);
	Vector pn = normal.transpose() * Dp;

	return pn;
}


Matrix Basis::make_Dp(const Point2 &pt) const
{
	auto row0 = make_p(pt, 1);
	auto row1 = make_p(pt, 2);

	Matrix Dp(2, row0.size());

	Dp.row(0) = row0;
	Dp.row(1) = row1;

   	return Dp;
}

int Basis::get_dof_vp_1d() const
{
	return dof_vp_1d;
}

int Basis::get_dof_u_1d() const
{
	return dof_u_1d;
}

void Basis::col_fhat_init_setup(int fhat_v_basis_group_id_input)
{
	if(fhat_v_basis_group_id_input > -1)
	{
		fhat_v_basis_id = fhat_v_basis_group_id_input;
	}
	else
	{
		fhat_v_basis_id = fhat_v_basis_group_id_default;
	}

	fhat_valid_to_org_vbasis_i.resize(2);

	Point2 v_sp_node0, v_sp_node1;
	v_sp_node0 = mesh->get_node_element(ElementV(0))[0];
	v_sp_node1 = mesh->get_node_element(ElementV(0))[2];

	mesh_dvp = v_sp_node1[0] - v_sp_node0[0];
	mesh_du = v_sp_node1[1] - v_sp_node0[1];

	Point4 test_dummy_p = Point4({0.0, 0.0, 1.0, M_1_PI}); //: ad hoc sample point to find out fhat_valid_v_basis_dim
	Vector dummy_basis = fhat_vbasis_Vec(test_dummy_p, 0, 0);

	fhat_org_v_basis_dim = dummy_basis.size()/3;

	for(int k = 0; k < fhat_org_v_basis_dim; k++)
	{
		int valid_flag = 0;
		for(int j = 0; j < 3; j++)
		{
			int loc_index = j*fhat_org_v_basis_dim + k;
			value_type org_val = dummy_basis[loc_index];
			if(org_val != 0.0) valid_flag = 1;
		}

		if(valid_flag == 1) fhat_valid_to_org_vbasis_i[0].push_back(k);
	}

	fhat_valid_v_basis_dim = fhat_valid_to_org_vbasis_i[0].size();

	test_dummy_p = Point4({0.0, 0.0, M_1_PI, 1.0});
	dummy_basis = fhat_vbasis_Vec(test_dummy_p, 1, 0);

	for(int k = 0; k < fhat_org_v_basis_dim; k++)
	{
		int valid_flag = 0;
		for(int j = 0; j < 3; j++)
		{
			int loc_index = j*fhat_org_v_basis_dim + k;
			value_type org_val = dummy_basis[loc_index];
			if(org_val != 0.0) valid_flag = 1;
		}

		if(valid_flag == 1) fhat_valid_to_org_vbasis_i[1].push_back(k);
	}

	if(fhat_valid_v_basis_dim != fhat_valid_to_org_vbasis_i[1].size())
	{
		cout << "fhat_valid_to_org_vbasis_i[0].size() should be same to fhat_valid_to_org_vbasis_i[1].size() : " << fhat_valid_to_org_vbasis_i[0].size() << " " << fhat_valid_to_org_vbasis_i[1].size() << endl;
		abort();
	}

	if(fhat_valid_v_basis_dim < dofv)
	{
		cout << "fhat_valid_v_basis_dim should be equal or bigger than dofv : " << fhat_valid_v_basis_dim << " " << dofv << endl;

		abort();
	}

}

Vector Basis::fhat_vbasis_Vec(const Point4 &p, int bd_type, int ele_index) const
{

	int fhat_org_v_basis_dim_loc;
	vector<value_type> vf[3];
	vector<value_type> tmp_norm_fac = {1.0, 2.0/mesh_dvp, 2.0/mesh_du};
	value_type mod_vp, mod_u;

	switch(bd_type)
	{
		case 0:
				if (ele_index == 0) mod_vp = p[2] - 1.0;
				else mod_vp = p[2] + 1.0;
				break;

		case 1:
				if (ele_index == 0) mod_u = p[3] - 1.0;
				else mod_u = p[3] + 1.0;
				break;

		case 2:
				break;

		case 3:
				break;

		default:
			cerr << "Unknown bd_type : " << bd_type << endl;
			exit(EXIT_FAILURE);
	}

	switch(fhat_v_basis_id)
	{
		case 0:
			fhat_org_v_basis_dim_loc = 6;
			for(int j = 0; j < 3; j++) vf[j].resize(fhat_org_v_basis_dim_loc);

			if(bd_type == 0)
			{

				vf[0] = {1.0, mod_vp, pow(mod_vp,2.0), pow(mod_vp,3.0), p[3], p[3]*mod_vp};
				vf[1] = {0.0, 1.0, 2.0*pow(mod_vp,1.0), 3.0*pow(mod_vp,2.0), 0.0, p[3]};
				vf[2] = {0.0, 0.0, 0.0, 0.0, 1.0, mod_vp};
			}
			else if(bd_type == 1)
			{

				vf[0] = {1.0, mod_u, pow(mod_u,2.0), pow(mod_u,3.0), p[2], p[2]*mod_u};
				vf[1] = {0.0, 0.0, 0.0, 0.0, 1.0, mod_u};
				vf[2] = {0.0, 1.0, 2.0*pow(mod_u,1.0), 3.0*pow(mod_u,2.0), 0.0, p[2]};
			}
			break;

		case 1:
			fhat_org_v_basis_dim_loc = 12;
			for(int j = 0; j < 3; j++) vf[j].resize(fhat_org_v_basis_dim_loc);
			if(bd_type == 0)
			{
				vf[0] = {1.0, mod_vp, pow(mod_vp,2.0), pow(mod_vp,3.0), pow(mod_vp,4.0), pow(mod_vp,5.0), p[3], p[3]*mod_vp, p[3]*pow(mod_vp,2.0), p[3]*pow(mod_vp,3.0), p[3]*p[3],  p[3]*p[3]*mod_vp};
				vf[1] = {0.0, 1.0, 2.0*pow(mod_vp,1.0), 3.0*pow(mod_vp,2.0), 4.0*pow(mod_vp,3.0), 5.0*pow(mod_vp,4.0), 0.0, p[3], p[3]*2.0*pow(mod_vp,1.0), p[3]*3.0*pow(mod_vp,2.0), 0.0,  p[3]*p[3]};
				vf[2] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, mod_vp, pow(mod_vp,2.0), pow(mod_vp,3.0), 2.0*p[3],  2.0*p[3]*mod_vp};
			}
			else if(bd_type == 1)
			{
				vf[0] = {1.0, mod_u, pow(mod_u,2.0), pow(mod_u,3.0), pow(mod_u,4.0), pow(mod_u,5.0), p[2], p[2]*mod_u, p[2]*pow(mod_u,2.0), p[2]*pow(mod_u,3.0), p[2]*p[2],  p[2]*p[2]*mod_u};
				vf[1] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, mod_u, pow(mod_u,2.0), pow(mod_u,3.0), 2.0*p[2],  2.0*p[2]*mod_u};
				vf[2] = {0.0, 1.0, 2.0*pow(mod_u,1.0), 3.0*pow(mod_u,2.0), 4.0*pow(mod_u,3.0), 5.0*pow(mod_u,4.0), 0.0, p[2], p[2]*2.0*pow(mod_u,1.0), p[2]*3.0*pow(mod_u,2.0), 0.0,  p[2]*p[2]};
			}

			break;

		default:
			cerr << "Unknown fhat_v_basis_id : " << fhat_v_basis_id << endl;
			exit(EXIT_FAILURE);
	}


	Vector fhat_vbasis_Vec_out = Vector::Zero(fhat_org_v_basis_dim_loc*3);
	if(bd_type == 0 || bd_type == 1)
	{
		for(int j = 0; j < 3; j++)
		{
			for(int k = 0; k < fhat_org_v_basis_dim_loc; k++)
			{
				int loc_index = j*fhat_org_v_basis_dim_loc + k;
				fhat_vbasis_Vec_out[loc_index] = tmp_norm_fac[j]*vf[j][k];
			}
		}
	}
	else
	{
		//if bd_type ==2 or 3
		//dg vbasis val
		vector<value_type> basisV(dofv*6);
		Point2 qb = Point2({p[2], p[3]});
		col_velocity_basis_vec(qb, mesh_dvp, mesh_du, basisV);

		for(int j = 0; j < 3; j++)
		{
			for(int k = 0; k < dofv; k++)
			{
				int loc_index = j*fhat_org_v_basis_dim_loc + k;
				fhat_vbasis_Vec_out[loc_index] = basisV[k*6 + j];
			}
		}
	}
	return fhat_vbasis_Vec_out;
}

int Basis::fhat_valid_to_org_i(const int &fhat_valid_index, const int &i_op) const
{
	if(i_op != 0 && i_op != 1)
	{
		cout << "Invalid i_op in fhat_valid_to_org_i : " << i_op << endl;
		abort();
	}

	int fhat_org_index = fhat_valid_to_org_vbasis_i[i_op][fhat_valid_index];

	return fhat_org_index;
}




void Basis::col_fhat_setup(int fhat_v_basis_group_id_input)
{

	if(fhat_v_basis_group_id_input > -1)
	{
		fhat_v_basis_group_id = fhat_v_basis_group_id_input;
	}
	else
	{
		fhat_v_basis_group_id = fhat_v_basis_group_id_default;
	}

	col_fhat_dofv_setup();

	dofv_fhat = get_dofv_fhat();
	max_dof_v_1d_fhat = get_max_dof_v_1d_fhat();

	dof_vp_1d = dof_u_1d = max_dof_v_1d_fhat;

	col_fhat_mat_num.resize(col_edge_type_num);
	if (op_divide_v_basis == 0)
	{
		col_fhat_mat_num[0] = 1;
		col_fhat_mat_num[1] = 1;
	}
	else
	{
		col_fhat_mat_num[0] = dof_u_1d;
		col_fhat_mat_num[1] = dof_vp_1d;
	}

	col_fhat_mat_half_size.resize(col_edge_type_num);

	col_fhat_org_mat_group.resize(col_edge_type_num);
	col_fhat_org_mat_index.resize(col_edge_type_num);
	col_fhat_new_mat_group.resize(col_edge_type_num);
	col_fhat_new_mat_index.resize(col_edge_type_num);

	for(int i = 0; i < col_edge_type_num; i++)
	{
		col_fhat_mat_half_size[i].resize(col_fhat_mat_num[i]);
		col_fhat_org_mat_group[i].resize(dofv);
		col_fhat_org_mat_index[i].resize(dofv);
		col_fhat_new_mat_group[i].resize(2*dofv);
		col_fhat_new_mat_index[i].resize(2*dofv);

	}

	if (op_divide_v_basis == 0)
	{
		for(int i = 0; i < col_edge_type_num; i++)
		{
			col_fhat_mat_half_size[i][0] = dofv;
			col_fhat_org_mat_group[i].assign(col_fhat_org_mat_group[i].size(), 0.0);

		}
	}
	else
	{
		if (dofv == 3)
		{
			col_fhat_mat_half_size[0][0] = 2;
			col_fhat_mat_half_size[0][1] = 1;
			col_fhat_org_mat_group[0] = {0, 0, 1};

			col_fhat_mat_half_size[1][0] = 2;
			col_fhat_mat_half_size[1][1] = 1;
			col_fhat_org_mat_group[1] = {0, 1, 0};
		}
		else if (dofv == 6)
		{
			col_fhat_mat_half_size[0][0] = 3;
			col_fhat_mat_half_size[0][1] = 2;
			col_fhat_mat_half_size[0][2] = 1;
			col_fhat_org_mat_group[0] = {0, 0, 1, 0, 1, 2};

			col_fhat_mat_half_size[1][0] = 3;
			col_fhat_mat_half_size[1][1] = 2;
			col_fhat_mat_half_size[1][2] = 1;
			col_fhat_org_mat_group[1] = {0, 1, 0, 2, 1, 0};
		}
	}

	int tmp_index;
	for(int i = 0; i < col_edge_type_num; i++)
	{
		vector<int> tmp_mat_index(col_fhat_mat_num[i]);
		memset(&tmp_mat_index[0], 0, tmp_mat_index.size() * sizeof tmp_mat_index[0]);
		memset(&col_fhat_org_mat_index[i][0], 0, col_fhat_org_mat_index[i].size() * sizeof col_fhat_org_mat_index[i][0]);

		for(int j = 0; j < dofv; j++)
		{
			tmp_index = col_fhat_org_mat_group[i][j];

			col_fhat_org_mat_index[i][j] = tmp_mat_index[tmp_index];
			tmp_mat_index[tmp_index] += 1;
		}

		tmp_index = 0;
		memset(&col_fhat_new_mat_group[i][0], 0, col_fhat_new_mat_group[i].size() * sizeof col_fhat_new_mat_group[i][0]);
		for(int j = 0; j < col_fhat_mat_num[i]; j++)
		{
			for(int k = 0; k < col_fhat_mat_half_size[i][j]*2; k++)
			{
				col_fhat_new_mat_group[i][tmp_index] = j;
				tmp_index += 1;
			}
		}

		memset(&tmp_mat_index[0], 0, tmp_mat_index.size() * sizeof tmp_mat_index[0]);
		memset(&col_fhat_new_mat_index[i][0], 0, col_fhat_new_mat_index[i].size() * sizeof col_fhat_new_mat_index[i][0]);

		for(int j = 0; j < 2*dofv; j++)
		{
			tmp_index = col_fhat_new_mat_group[i][j];

			col_fhat_new_mat_index[i][j] = tmp_mat_index[tmp_index];
			tmp_mat_index[tmp_index] += 1;
		}
	}
}

int Basis::get_col_edge_type_num() const
{
	return col_edge_type_num;
}



void Basis::get_col_fhat_data(const int i, int &_col_fhat_mat_num, vector<int> &_col_fhat_mat_half_size, vector<int> &_col_fhat_org_mat_group, vector<int> &_col_fhat_org_mat_index, vector<int> &_col_fhat_new_mat_group, vector<int> &_col_fhat_new_mat_index) const
{

	_col_fhat_mat_num = col_fhat_mat_num[i];

	_col_fhat_mat_half_size.resize(_col_fhat_mat_num);
	for(int j = 0; j < _col_fhat_mat_num; j++)
	{
		_col_fhat_mat_half_size[j] = col_fhat_mat_half_size[i][j];
	}

	_col_fhat_org_mat_group.resize(dofv);
	_col_fhat_org_mat_index.resize(dofv);
	for(int j = 0; j < dofv; j++)
	{
		_col_fhat_org_mat_group[j] = col_fhat_org_mat_group[i][j];
		_col_fhat_org_mat_index[j] = col_fhat_org_mat_index[i][j];
	}

	_col_fhat_new_mat_group.resize(2*dofv);
	_col_fhat_new_mat_index.resize(2*dofv);
	for(int j = 0; j < 2*dofv; j++)
	{
		_col_fhat_new_mat_group[j] = col_fhat_new_mat_group[i][j];
		_col_fhat_new_mat_index[j] = col_fhat_new_mat_index[i][j];
	}
}

//original basis part
//col_spatial_basis_vec
//col_velocity_basis_vec
//col_velocity_fhat_basis_bd_info
void Basis::col_spatial_basis_vec(const Point2 &p, vector<value_type> &basisX) const
{
	vector<value_type> xf(dofx);
	switch(x_basis_group_id)
	{
		case 0:
			xf = {1.0 - p[0] - p[1], p[0], p[1]};
			break;

		case 1:
			xf = {1.0 - p[0] - p[1], p[0], p[1], (1 - p[0] - p[1])*(1 - p[0] - p[1]), p[0]*p[0], p[1]*p[1]};
			break;

		case 2:
			xf = {1.0};
			break;


		default:
			cerr << "Unknown x_basis_group_id in col_spatial_basis_vec : " << x_basis_group_id << endl;
			exit(EXIT_FAILURE);
	}

	for(int i = 0; i < dofx; i++) basisX[i] = xf[i];
}

void Basis::col_velocity_basis_vec(const Point2 &p, const value_type dvp, const value_type du, vector<value_type> &basisV) const
{
	vector<value_type> tmp_norm_fac = {1.0, 2.0/dvp, 2.0/du, 4.0/dvp/dvp, 4.0/dvp/du, 4.0/du/du};
	vector<vector<value_type>> vf;

	vf.resize(dofv);
	switch(v_basis_group_id)
	{
		case 0:
			vf[0] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			vf[1] = {p[0], 1.0, 0.0, 0.0, 0.0, 0.0};
			vf[2] = {p[1], 0.0, 1.0, 0.0, 0.0, 0.0};
			break;

		case 1:
			vf[0] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			vf[1] = {p[0], 1.0, 0.0, 0.0, 0.0, 0.0};
			vf[2] = {p[1], 0.0, 1.0, 0.0, 0.0, 0.0};
			vf[3] = {p[0]*p[0], 2.0*p[0], 0.0, 2.0, 0.0, 0.0};
			vf[4] = {p[0]*p[1], p[1], p[0], 0.0, 1.0, 0.0};
			vf[5] = {p[1]*p[1], 0.0, 2.0*p[1], 0.0, 0.0, 2.0};
			break;

		default:
			cerr << "Unknown v_basis_group_id in col_velocity_basis_vec : " << v_basis_group_id << endl;
			exit(EXIT_FAILURE);
	}

	for(int i = 0; i < dofv; i++)
		for(int j = 0; j < 6; j++)
			basisV[i*6 + j] = vf[i][j]*tmp_norm_fac[j];
}

//connect org v basis to the fhat values at the last v surface
void Basis::col_velocity_fhat_basis_bd_info(int q, int bd_type, value_type vp_val, value_type u_val, int &fhat_out, value_type &fhat_factor) const 
{
	assert(q >= 0 && q < dofv);

	vector<int> tmp_fhat_out;
	vector<value_type> tmp_fhat_factor;

	tmp_fhat_out.resize(dofv);
	tmp_fhat_factor.resize(dofv);

	//for bd_type == 2 -> tmp_fhat_out[j] : x of u^x dependence at vp = vp_val boundary when j = org v basis index
	//for bd_type == 3 -> tmp_fhat_out[j] : x of v_p^x dependence at u = u_val boundary when j = org v basis index
	//for bd_type == 2 -> tmp_fhat_factor[j] -> value of [org v basis with index j]/u^x at vp = vp_val boundary
	//for bd_type == 3 -> tmp_fhat_factor[j] -> value of [org v basis with index j]/v_p^x at u = u_val boundary
	switch(v_basis_group_id)
	{
		case 0:
			if(bd_type == 2)
			{
				tmp_fhat_out[0] = 0; //1
				tmp_fhat_out[1] = 0; //vp
				tmp_fhat_out[2] = 1; //u

				tmp_fhat_factor[0] = 1.0;
				tmp_fhat_factor[1] = vp_val;
				tmp_fhat_factor[2] = 1.0;
			}
			else if(bd_type == 3)
			{
				tmp_fhat_out[0] = 0; //1
				tmp_fhat_out[1] = 1; //vp
				tmp_fhat_out[2] = 0; //u

				tmp_fhat_factor[0] = 1.0;
				tmp_fhat_factor[1] = 1.0;
				tmp_fhat_factor[2] = u_val;
			}
			else
			{
				cerr << "Unknown bd_type : " << bd_type << endl;
				exit(EXIT_FAILURE);
			}
			break;

		case 1:
			if(bd_type == 2)
			{
				tmp_fhat_out[0] = 0; //1
				tmp_fhat_out[1] = 0; //vp
				tmp_fhat_out[2] = 1; //u
				tmp_fhat_out[3] = 0; //vp^2
				tmp_fhat_out[4] = 1; //u*vp
				tmp_fhat_out[5] = 2; //u^2

				tmp_fhat_factor[0] = 1.0;
				tmp_fhat_factor[1] = vp_val;
				tmp_fhat_factor[2] = 1.0;
				tmp_fhat_factor[3] = vp_val*vp_val;
				tmp_fhat_factor[4] = vp_val;
				tmp_fhat_factor[5] = 1.0;
			}
			else if(bd_type == 3)
			{
				tmp_fhat_out[0] = 0; //1
				tmp_fhat_out[1] = 1; //vp
				tmp_fhat_out[2] = 0; //u
				tmp_fhat_out[3] = 2; //vp^2
				tmp_fhat_out[4] = 1; //u*vp
				tmp_fhat_out[5] = 0; //u^2

				tmp_fhat_factor[0] = 1.0;
				tmp_fhat_factor[1] = 1.0;
				tmp_fhat_factor[2] = u_val;
				tmp_fhat_factor[3] = 1.0;
				tmp_fhat_factor[4] = u_val;
				tmp_fhat_factor[5] = u_val*u_val;
			}
			else
			{
				cerr << "Unknown bd_type : " << bd_type << endl;
				exit(EXIT_FAILURE);
			}
			break;

		default:
			cerr << "Unknown v_basis_group_id in col_velocity_fhat_basis_bd_info : " << v_basis_group_id << endl;
			exit(EXIT_FAILURE);
	}

	fhat_out = tmp_fhat_out[q];
	fhat_factor = tmp_fhat_factor[q];

	if(fhat_out >= max_dof_v_1d_fhat)
	{
		cerr << "fhat_out should be smaller than max_dof_v_1d_fhat in col_velocity_fhat_basis_bd_info : " << fhat_out << " " << max_dof_v_1d_fhat << endl;
		exit(EXIT_FAILURE);
	}

}

//fhat part
//col_fhat_dofv_setup
//col_velocity_fhat_basis
//col_velocity_fhat_basis_info
void Basis::col_fhat_dofv_setup()
{
	switch(fhat_v_basis_group_id)
	{
		case 0:
			dofv_fhat = 6;
			max_dof_v_1d_fhat = 2;
			break;

		case 1:
			dofv_fhat = 12;
			max_dof_v_1d_fhat = 3;
			break;

		case 2:
			dofv_fhat = 5;
			max_dof_v_1d_fhat = 2;
			break;

		case 3:
			dofv_fhat = 9;
			max_dof_v_1d_fhat = 3;
			break;

		case 4:
			dofv_fhat = 4;
			max_dof_v_1d_fhat = 2;
			break;

		case 5:
			dofv_fhat = 6;
			max_dof_v_1d_fhat = 3;
			break;

		case 6:
			dofv_fhat = 5;
			max_dof_v_1d_fhat = 2;
			break;

		case 7:
			dofv_fhat = 10;
			max_dof_v_1d_fhat = 3;
			break;

		default:
			cerr << "Unknown fhat_v_basis_group_id in col_fhat_dofv_setup : " << fhat_v_basis_group_id << endl;
			exit(EXIT_FAILURE);

	}
}

value_type Basis::col_velocity_fhat_basis(const Point4 &p, int q, int bd_type, int ele_index) const
{
	assert(q >= 0 && q < dofv_fhat);

	if(bd_type != 0 && bd_type != 1)
	{
		cerr << "Unknown bd_type : " << bd_type << endl;
		exit(EXIT_FAILURE);
	}

	vector<value_type> vf(dofv_fhat);
	value_type mod_vp, mod_u;
	if (ele_index == 0) 
	{
		mod_vp = p[2] - 1.0;
		mod_u = p[3] - 1.0;
	}
	else 
	{
		mod_vp = p[2] + 1.0;
		mod_u = p[3] + 1.0;
	}

	switch(fhat_v_basis_group_id)
	{
		case 0:
			if(bd_type == 0)
				vf = {1.0, mod_vp, pow(mod_vp,2.0), pow(mod_vp,3.0), p[3], p[3]*mod_vp};
			else if(bd_type == 1)
				vf = {1.0, mod_u, pow(mod_u,2.0), pow(mod_u,3.0), p[2], p[2]*mod_u};
			break;

		case 1:
			if(bd_type == 0)
				vf = {1.0, mod_vp, pow(mod_vp,2.0), pow(mod_vp,3.0), pow(mod_vp,4.0), pow(mod_vp,5.0), p[3], p[3]*mod_vp, p[3]*pow(mod_vp,2.0), p[3]*pow(mod_vp,3.0), p[3]*p[3],  p[3]*p[3]*mod_vp};
			else if(bd_type == 1)
				vf = {1.0, mod_u, pow(mod_u,2.0), pow(mod_u,3.0), pow(mod_u,4.0), pow(mod_u,5.0), p[2], p[2]*mod_u, p[2]*pow(mod_u,2.0), p[2]*pow(mod_u,3.0), p[2]*p[2],  p[2]*p[2]*mod_u};
			break;

		case 2:
			if(bd_type == 0)
				vf = {1.0, mod_vp, pow(mod_vp,2.0), p[3], p[3]*mod_vp};
			else if(bd_type == 1)
				vf = {1.0, mod_u, pow(mod_u,2.0), p[2], p[2]*mod_u};
			break;

		case 3:
			if(bd_type == 0)
				vf = {1.0, mod_vp, pow(mod_vp,2.0), pow(mod_vp,3.0), p[3], p[3]*mod_vp, p[3]*pow(mod_vp,2.0), p[3]*p[3],  p[3]*p[3]*mod_vp};
			else if(bd_type == 1)
				vf = {1.0, mod_u, pow(mod_u,2.0), pow(mod_u,3.0), p[2], p[2]*mod_u, p[2]*pow(mod_u,2.0), p[2]*p[2],  p[2]*p[2]*mod_u};
			break;

		case 4:
			if(bd_type == 0)
				vf = {1.0, mod_vp, p[3], p[3]*mod_vp};
			else if(bd_type == 1)
				vf = {1.0, mod_u, p[2], p[2]*mod_u};
			break;

		case 5:
			if(bd_type == 0)
				vf = {1.0, mod_vp, p[3], p[3]*mod_vp, p[3]*p[3], p[3]*p[3]*mod_vp};
			else if(bd_type == 1)
				vf = {1.0, mod_u, p[2], p[2]*mod_u, p[2]*p[2], p[2]*p[2]*mod_u};
			break;

		case 6:
			if(bd_type == 0)
				vf = {1.0, mod_vp, pow(mod_vp,2.0), p[3], p[3]*mod_vp};
			else if(bd_type == 1)
				vf = {1.0, mod_u, pow(mod_u,2.0), p[2], p[2]*mod_u};
			break;

		case 7:
			if(bd_type == 0)
				vf = {1.0, mod_vp, pow(mod_vp,2.0), pow(mod_vp,3.0), p[3], p[3]*mod_vp, p[3]*pow(mod_vp,2.0), p[3]*pow(mod_vp,3.0), p[3]*p[3],  p[3]*p[3]*mod_vp};
			else if(bd_type == 1)
				vf = {1.0, mod_u, pow(mod_u,2.0), pow(mod_u,3.0), p[2], p[2]*mod_u, p[2]*pow(mod_u,2.0), p[2]*pow(mod_u,3.0), p[2]*p[2],  p[2]*p[2]*mod_u};
			break;

		default:
			cerr << "Unknown fhat_v_basis_group_id in col_velocity_fhat_basis: " << fhat_v_basis_group_id << endl;
			exit(EXIT_FAILURE);
	}

	return vf[q];
}

void Basis::col_velocity_fhat_basis_info(int q, int bd_type, int i_derivs, value_type dvp, value_type du, int &fhat_out, value_type &fhat_factor) const 
{
	assert(q >= 0 && q < dofv_fhat);

	if(bd_type != 0 && bd_type != 1)
	{
		cerr << "Unknown bd_type : " << bd_type << endl;
		exit(EXIT_FAILURE);
	}

	vector<vector<int>> tmp_fhat_out;
	vector<vector<value_type>> tmp_fhat_factor;

	tmp_fhat_out.resize(3);
	tmp_fhat_factor.resize(3);
	for(int i = 0; i < 3; i++)
	{
		tmp_fhat_out[i].resize(dofv_fhat);
		tmp_fhat_factor[i].resize(dofv_fhat);
		tmp_fhat_out[i].assign(tmp_fhat_out[i].size(), -1);

		value_type tmp_norm_fac;

		switch(i)
		{
			case 0:
				tmp_norm_fac = 1.0;
				break;

			case 1:
				tmp_norm_fac = 2.0/dvp;
				break;

			case 2:
				tmp_norm_fac = 2.0/du;
				break;
		}
		tmp_fhat_factor[i].assign(tmp_fhat_factor[i].size(), tmp_norm_fac);
#ifdef SEO_COL_TEST2
		cout << "fhat_factor testa : " << i_derivs << " " << q << " " << fhat_factor << " " << tmp_fhat_factor[i][q] << " " << sizeof(tmp_fhat_factor[i][0]) << " " << tmp_fhat_factor[i].size() << endl;
#endif
	}

#ifdef SEO_COL_TEST2
	/*
	   cout << "fhat_factor testa : " << i_derivs << " " << q << " " << fhat_factor << " " << tmp_fhat_factor[i_derivs][q] << " " << sizeof(tmp_fhat_factor[0][0]) << endl;
	   cout << "fhat_factor testa : " << i_derivs << " " << q << " " << fhat_factor << " " << tmp_fhat_factor[i_derivs][q] << " " << sizeof(tmp_fhat_factor[1][0]) << endl;
	   cout << "fhat_factor testa : " << i_derivs << " " << q << " " << fhat_factor << " " << tmp_fhat_factor[i_derivs][q] << " " << sizeof(tmp_fhat_factor[2][0]) << endl;
	   */
	abort();
#endif

	//for bd_type == 0 -> tmp_fhat_out[i][j] : x of u^x dependence at mod_vp == 0 when i = fhat derivative index, j = fhat v basis index
	//for bd_type == 1 -> tmp_fhat_out[i][j] : x of v_p^x dependence at mod_u == 0 when i = fhat derivative index, j = fhat v basis index
	switch(fhat_v_basis_group_id)
	{
		case 0:
			if(bd_type == 0)
			{
				//REf : vf = {1.0, mod_vp, pow(mod_vp,2.0), pow(mod_vp,3.0), p[3], p[3]*mod_vp};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][4] = 1; //u

				tmp_fhat_out[1][1] = 0; //v_p
				tmp_fhat_out[1][5] = 1; //u*v_p

				tmp_fhat_out[2][4] = 0; //u
			}
			else if(bd_type == 1)
			{
				//REF : vf = {1.0, mod_u, pow(mod_u,2.0), pow(mod_u,3.0), p[2], p[2]*mod_u};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][4] = 1; //v_p

				tmp_fhat_out[1][4] = 0; //v_p
										//tmp_fhat_out[1][5] = 1; //u*v_p

				tmp_fhat_out[2][1] = 0; //u
				tmp_fhat_out[2][5] = 1; //u*v_p
			}
			break;

		case 1:
			if(bd_type == 0)
			{
				//REF : vf = {1.0, mod_vp, pow(mod_vp,2.0), pow(mod_vp,3.0), pow(mod_vp,4.0), pow(mod_vp,5.0), p[3], p[3]*mod_vp, p[3]*pow(mod_vp,2.0), p[3]*pow(mod_vp,3.0), p[3]*p[3],  p[3]*p[3]*mod_vp};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][6] = 1; //u
				tmp_fhat_out[0][10] = 2; //u^2

				tmp_fhat_out[1][1] = 0; //v_p
				tmp_fhat_out[1][7] = 1; //u*v_p
				tmp_fhat_out[1][11] = 2; //u^2*v_p

				tmp_fhat_out[2][6] = 0; //u
				tmp_fhat_out[2][10] = 1; //u^2
				tmp_fhat_factor[2][10] *= 2.0; // due to d u^2/du
			}
			else if(bd_type == 1)
			{
				//REF : vf = {1.0, mod_u, pow(mod_u,2.0), pow(mod_u,3.0), pow(mod_u,4.0), pow(mod_u,5.0), p[2], p[2]*mod_u, p[2]*pow(mod_u,2.0), p[2]*pow(mod_u,3.0), p[2]*p[2],  p[2]*p[2]*mod_u};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][6] = 1; //v_p
				tmp_fhat_out[0][10] = 2; //v_p^2

				tmp_fhat_out[1][6] = 0; //v_p
				tmp_fhat_out[1][10] = 1; //v_p^2
				tmp_fhat_factor[1][10] *= 2.0; // due to d v_p^2/dv_p

				tmp_fhat_out[2][1] = 0; //u
				tmp_fhat_out[2][7] = 1; //v_p*u
				tmp_fhat_out[2][11] = 2; //v_p^2*u
			}
			break;

		case 2:
			if(bd_type == 0)
			{
				//REf : vf = {1.0, mod_vp, pow(mod_vp,2.0), p[3], p[3]*mod_vp};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][3] = 1; //u

				tmp_fhat_out[1][1] = 0; //v_p
				tmp_fhat_out[1][4] = 1; //u*v_p

				tmp_fhat_out[2][3] = 0; //u
			}
			else if(bd_type == 1)
			{
				//REF : vf = {1.0, mod_u, pow(mod_u,2.0), p[2], p[2]*mod_u};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][3] = 1; //v_p

				tmp_fhat_out[1][3] = 0; //v_p
										//tmp_fhat_out[1][4] = 1; //u*v_p

				tmp_fhat_out[2][1] = 0; //u
				tmp_fhat_out[2][4] = 1; //u*v_p
			}
			break;

		case 3:
			if(bd_type == 0)
			{
				//REF : vf = {1.0, mod_vp, pow(mod_vp,2.0), pow(mod_vp,3.0), p[3], p[3]*mod_vp, p[3]*pow(mod_vp,2.0), p[3]*p[3],  p[3]*p[3]*mod_vp};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][4] = 1; //u
				tmp_fhat_out[0][7] = 2; //u^2

				tmp_fhat_out[1][1] = 0; //v_p
				tmp_fhat_out[1][5] = 1; //u*v_p
				tmp_fhat_out[1][8] = 2; //u^2*v_p

				tmp_fhat_out[2][4] = 0; //u
				tmp_fhat_out[2][7] = 1; //u^2
				tmp_fhat_factor[2][7] *= 2.0; // due to d u^2/du
			}
			else if(bd_type == 1)
			{
				//REF : vf = {1.0, mod_u, pow(mod_u,2.0), pow(mod_u,3.0), p[2], p[2]*mod_u, p[2]*pow(mod_u,2.0), p[2]*p[2],  p[2]*p[2]*mod_u};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][4] = 1; //v_p
				tmp_fhat_out[0][7] = 2; //v_p^2

				tmp_fhat_out[1][4] = 0; //v_p
				tmp_fhat_out[1][7] = 1; //v_p^2
				tmp_fhat_factor[1][7] *= 2.0; // due to d v_p^2/dv_p

				tmp_fhat_out[2][1] = 0; //u
				tmp_fhat_out[2][5] = 1; //v_p*u
				tmp_fhat_out[2][8] = 2; //v_p^2*u
			}
			break;

		case 4:
			if(bd_type == 0)
			{
				//REf : vf = {1.0, mod_vp, p[3], p[3]*mod_vp};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][2] = 1; //u

				tmp_fhat_out[1][1] = 0; //v_p
				tmp_fhat_out[1][3] = 1; //u*v_p

				tmp_fhat_out[2][2] = 0; //u
			}
			else if(bd_type == 1)
			{
				//REF : vf = {1.0, mod_u, p[2], p[2]*mod_u};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][2] = 1; //v_p

				tmp_fhat_out[1][2] = 0; //v_p
										//tmp_fhat_out[1][3] = 1; //u*v_p

				tmp_fhat_out[2][1] = 0; //u
				tmp_fhat_out[2][3] = 1; //u*v_p
			}
			break;


		case 5:
			if(bd_type == 0)
			{
				//REF : vf = {1.0, mod_vp, p[3], p[3]*mod_vp, p[3]*p[3], p[3]*p[3]*mod_vp};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][2] = 1; //u
				tmp_fhat_out[0][4] = 2; //u^2

				tmp_fhat_out[1][1] = 0; //v_p
				tmp_fhat_out[1][3] = 1; //u*v_p
				tmp_fhat_out[1][5] = 2; //u^2*v_p

				tmp_fhat_out[2][2] = 0; //u
				tmp_fhat_out[2][4] = 1; //u^2
				tmp_fhat_factor[2][4] *= 2.0; // due to d u^2/du
			}
			else if(bd_type == 1)
			{
				//REF : vf = {1.0, mod_u, p[2], p[2]*mod_u, p[2]*p[2], p[2]*p[2]*mod_u};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][2] = 1; //v_p
				tmp_fhat_out[0][4] = 2; //v_p^2

				tmp_fhat_out[1][2] = 0; //v_p
				tmp_fhat_out[1][4] = 1; //v_p^2
				tmp_fhat_factor[1][4] *= 2.0; // due to d v_p^2/dv_p

				tmp_fhat_out[2][1] = 0; //u
				tmp_fhat_out[2][3] = 1; //v_p*u
				tmp_fhat_out[2][5] = 2; //v_p^2*u
			}
			break;

		case 6:
			if(bd_type == 0)
			{
				//REf : vf = {1.0, mod_vp, pow(mod_vp,2.0), p[3], p[3]*mod_vp};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][4] = 1; //u

				tmp_fhat_out[1][1] = 0; //v_p
				tmp_fhat_out[1][5] = 1; //u*v_p

				tmp_fhat_out[2][4] = 0; //u
			}
			else if(bd_type == 1)
			{
				//REF : vf = {1.0, mod_u, pow(mod_u,2.0), p[2], p[2]*mod_u};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][4] = 1; //v_p

				tmp_fhat_out[1][4] = 0; //v_p
										//tmp_fhat_out[1][5] = 1; //u*v_p

				tmp_fhat_out[2][1] = 0; //u
				tmp_fhat_out[2][5] = 1; //u*v_p
			}
			break;

		case 7:
			if(bd_type == 0)
			{
				//REF : vf = {1.0, mod_vp, pow(mod_vp,2.0), pow(mod_vp,3.0), p[3], p[3]*mod_vp, p[3]*pow(mod_vp,2.0), p[3]*pow(mod_vp,3.0), p[3]*p[3],  p[3]*p[3]*mod_vp};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][4] = 1; //u
				tmp_fhat_out[0][8] = 2; //u^2

				tmp_fhat_out[1][1] = 0; //v_p
				tmp_fhat_out[1][5] = 1; //u*v_p
				tmp_fhat_out[1][9] = 2; //u^2*v_p

				tmp_fhat_out[2][4] = 0; //u
				tmp_fhat_out[2][8] = 1; //u^2
				tmp_fhat_factor[2][8] *= 2.0; // due to d u^2/du
			}
			else if(bd_type == 1)
			{
				//REF : vf = {1.0, mod_u, pow(mod_u,2.0), pow(mod_u,3.0), p[2], p[2]*mod_u, p[2]*pow(mod_u,2.0), p[2]*pow(mod_u,3.0), p[2]*p[2],  p[2]*p[2]*mod_u};
				tmp_fhat_out[0][0] = 0; //1
				tmp_fhat_out[0][4] = 1; //v_p
				tmp_fhat_out[0][8] = 2; //v_p^2

				tmp_fhat_out[1][4] = 0; //v_p
				tmp_fhat_out[1][8] = 1; //v_p^2
				tmp_fhat_factor[1][8] *= 2.0; // due to d v_p^2/dv_p

				tmp_fhat_out[2][1] = 0; //u
				tmp_fhat_out[2][5] = 1; //v_p*u
				tmp_fhat_out[2][9] = 2; //v_p^2*u
			}
			break;

		default:
			cerr << "Unknown fhat_v_basis_group_id : " << fhat_v_basis_group_id << endl;
			exit(EXIT_FAILURE);
	}

	fhat_out = tmp_fhat_out[i_derivs][q];
	fhat_factor = tmp_fhat_factor[i_derivs][q];


}

void Basis::col_velocity_basis_vec_edge(const int i_op, const Point2 &p2, const value_type dvp, const value_type du, vector<value_type> &basisVE) const
{
	value_type p = p2[1];
	vector<value_type> vp_arr(2), u_arr(2);


	if (i_op == 0)
	{
		vp_arr[0] = 1.0;
		vp_arr[1] = -1.0;
		u_arr[0] = u_arr[1] = p;
	}
	else
	{
		vp_arr[0] = vp_arr[1] = p;
		u_arr[0] = 1.0;
		u_arr[1] = -1.0;
	}

	for(int l = 0; l < 2; l++)
	{
		Point2 qa = {vp_arr[l], u_arr[l]};

		vector<value_type> basisV(dofv*6);
		col_velocity_basis_vec(qa, dvp, du, basisV);

		for(int i = 0; i < dofv; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				basisVE[l*dofv*3 + i*3 + j] = basisV[i*6 + j];
			}
		}
	}
}

value_type Basis::col_spatial_basis(const Point4 &p, int r) const
{
	assert(r >= 0 && r < dofx);
	vector<value_type> xf(dofx);

	Point2 qa ={p[0], p[1]};
	col_spatial_basis_vec(qa, xf);

	return xf[r];
}

value_type Basis::col_velocity_org_basis(const Point4 &p, int q) const
{
	assert(q >= 0 && q < dofv);

	Point2 qb = {p[2], p[3]};
	//use dummy dvp, du because derivatives are not needed here
	value_type dvp_dum = 1.0, du_dum = 1.0;
	vector<value_type> basisV(dofv*6);

	col_velocity_basis_vec(qb, dvp_dum, du_dum, basisV);

	return basisV[q*6];
}


