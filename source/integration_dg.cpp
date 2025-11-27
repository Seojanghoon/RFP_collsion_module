
void Integration_MSE_mat_cal::init(const Mesh &_mesh, const Basis &_basis, const Quadrature &_quadrature, const Flux &_flux, const Species_data &_species_data, const EQ_reader &_eq_reader, const Functor &_f0, const int sp_id)
{
	mesh = &_mesh;
	basis = &_basis;
	quadrature = &_quadrature;
	flux = &_flux;
	f0 = &_f0;
	species_data = &_species_data;
	eq_reader = &_eq_reader;
	dof = _basis.get_dof();
	dof_x = _basis.get_dofx();
	dof_v = _basis.get_dofv();
	if (dof_x ==3) dof_S1_x = 2;
	else if (dof_x ==6) dof_S1_x = 4;
	else if (dof_x ==4) dof_S1_x = 2;
	else if (dof_x ==1) dof_S1_x = 1;
	else
	{
		cout << "Only dof_x = 3 or 6 case is ready. dof_x = " << dof_x << endl;
		abort();
	}

	dof_S1 = dof_S1_x*dof_v;

}

void Integration_MSE_mat_cal::ME_mat_cal_f_init(const ElementX &ti, const ElementV &tj, const int &sp_id, vector<value_type> &M_arr, vector<value_type> &E_arr, vector<value_type> &finit_arr)
{
	value_type w, wa, wb;
	Point4 p, q;
	Point2 pa, pb, qa, qb;
	flag_f0 = species_data->f0_type(sp_id);
	flag_fcm = species_data->fcm_type(sp_id);
	auto HT = quadrature->get_HT(ti, tj);
	value_type F[dof], tmp_basis_val[dof];
	memset(&M_arr[0], 0.0, M_arr.size() * sizeof M_arr[0]);
	memset(&E_arr[0], 0.0, E_arr.size() * sizeof E_arr[0]);
	memset(&finit_arr[0], 0.0, finit_arr.size() * sizeof finit_arr[0]);
	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		for(int b = 0; b < quadrature->size(1); b++)
		{
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

			w = wa * wb;
			p = pa & pb;
			q = qa & qb;

			value_type R = p[0];
			value_type u = p[3];

			auto bs_val = flux->get_b_dot_Bs(p, sp_id);
			auto f0_val = flux->get_f0(p, sp_id, flag_f0);
			auto fcm_val = flux->get_fcm_SE_tot(p, int(ti), sp_id, flag_fcm);
			value_type tmp_w1 = w * bs_val * R * u;
			value_type tmp_wf0 = tmp_w1 * fcm_val;
			value_type tmp_wfinit = tmp_w1 * f0_val;

			auto U = flux->get_U(p, sp_id);

			for(int i = 0; i < dof; i++) F[i] = U.dot(HT*basis->grad(q, p, ti, tj, i));
			for(int j = 0; j < dof; j++) tmp_basis_val[j] = (*basis)(q, p, ti, tj, j);

			for(int j = 0; j < dof; j++)
			{
				value_type tmp_w2 = tmp_wf0 * tmp_basis_val[j];
				finit_arr[j] += tmp_wfinit * tmp_basis_val[j];

				for(int i = 0; i < dof; i++)
				{
					if (i < j+1) M_arr[j*dof + i] += tmp_w2 * tmp_basis_val[i];
					E_arr[j*dof + i] += tmp_w2 * F[i]; 
				}
			}
		}
	}
	for(int j = 0; j < dof; j++) for(int i = 0; i < j; i++) M_arr[i*dof + j] = M_arr[j*dof + i];
}

void Integration_Col::init(const Mesh &_mesh, const Basis &_basis, const Quadrature &_quadrature, const Flux &_flux, const Species_data &_species_data, const EQ_reader &_eq_reader, int _sp_id, const int &_diag_num, const value_type &_fhat_vol_wid)
{
	mesh = &_mesh;
	basis = &_basis;
	quadrature = &_quadrature;
	flux = &_flux;
	species_data = &_species_data;
	eq_reader = &_eq_reader;
	sp_id = _sp_id;

	diag_num = _diag_num;
	fhat_vol_wid = _fhat_vol_wid;

	flag_f0 = species_data->f0_type(sp_id);
	flag_fcm = species_data->fcm_type(sp_id);

	dof = basis->get_dof();

	dof_x = basis->get_dofx();
	dof_v = basis->get_dofv();

	dof_v_1d = basis->get_max_dof_v_1d_fhat();

	fhat_dim_per_edge = dof_x*3*dof_v_1d;

	Point2 v_sp_node0, v_sp_node1;
	v_sp_node0 = mesh->get_node_element(ElementV(0))[0];
	v_sp_node1 = mesh->get_node_element(ElementV(0))[2];

	dvp = v_sp_node1[0] - v_sp_node0[0];
	du = v_sp_node1[1] - v_sp_node0[1];


	int qd_sizeX = quadrature->size(0), qd_sizeV = quadrature->size(1);
	int qd_sizeVE = quadrature->size(3);

	basisX_qd.resize(qd_sizeX);
	basisV_qd.resize(qd_sizeV, Vector::Zero(dof_v));

	//basisXV_qd.resize(qd_sizeX);
	for(int a = 0; a < qd_sizeX; a++)
	{
		value_type wa;
		Point2 pa, qa;
		vector<value_type> basisX(dof_x);

		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ElementX(0), a); 
		basis->col_spatial_basis_vec(qa, basisX);

		basisX_qd[a].resize(dof_x);
		for(int i = 0; i < dof_x; i++) basisX_qd[a][i] = basisX[i];

		for(int b = 0; b < qd_sizeV; b++)
		{
			value_type wb;
			Point2 pb, qb;
			vector<value_type> basisV(dof_v*6);

			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ElementV(0), b);
			basis->col_velocity_basis_vec(qb, dvp, du, basisV);

			for(int k = 0; k < 6; k++)
			{
				for(int i = 0; i < dof_x; i++)
				{
					value_type tmp_basisX = basisX[i];
					for(int j = 0; j < dof_v; j++)
					{

						if(a== 0 && k == 0 && i == 0)
							basisV_qd[b][j] = basisV[j*6 + k];
					}
				}
			}
		}
	}

	int nx_loc = mesh->own_size<ElementX>();
	int nv_loc = mesh->size<ElementV>();
	int nve_loc = mesh->size<EdgeV>();
	col_edge_terms_qd.resize(nx_loc);

	for(int ix = 0; ix < nx_loc; ix++)
	{
		for(int iv = 0; iv < nv_loc; iv++)
		{
			for(int a = 0; a < qd_sizeX; a++)
			{
				value_type wa;
				Point2 pa, qa;
				Point4 p;

				tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ElementX(ix), a); 

				p = pa & Point2({0, 0});
				Vector B = flux->get_B(p);
				value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));
				value_type R = pa[0], Z = pa[1];

				for(int b = 0; b < qd_sizeV; b++)
				{

					value_type wb;
					Point2 pb, qb;
					tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ElementV(iv), b);
					value_type w = wa * wb;
					p = pa & pb;

					value_type v_para = p[2], u = p[3];

					vector<value_type> bs_val_arr(2);
					bs_val_arr = flux->get_b_dot_Bs_vec(p, B, sp_id);
					value_type bs_val = bs_val_arr[0];
					value_type dbs_dvp_ov_bs = bs_val_arr[1]/bs_val;

					auto fcm_val = flux->get_fcm_SE_tot(p, ix, sp_id, flag_fcm);
					value_type ccf = w*R*u*bs_val*fcm_val;
				}
			}
		}

		col_edge_terms_qd[ix].resize(nve_loc);

		for(int ie = 0; ie < nve_loc; ie++)
		{
			col_edge_terms_qd[ix][ie].resize(qd_sizeX);

			for(int a = 0; a < qd_sizeX; a++)
			{
				value_type wa;
				Point2 pa, qa;
				Point4 p;

				tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ElementX(ix), a); 

				p = pa & Point2({0, 0});
				Vector B = flux->get_B(p);
				value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));
				value_type R = pa[0], Z = pa[1];

				col_edge_terms_qd[ix][ie][a].resize(qd_sizeVE);

				for(int b = 0; b < qd_sizeVE; b++)
				{
					value_type wb;
					Point2 pb, qb;
					tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(EdgeV(ie), b);
					value_type w = wa * wb;
					p = pa & pb;

					value_type v_para = p[2], u = p[3];

					vector<value_type> fcm_val_arr(3), bs_val_arr(2), ccf_arr(3);

					bs_val_arr = flux->get_b_dot_Bs_vec(p, B, sp_id);
					fcm_val_arr = flux->get_fcm_SE_vec(p, ix, sp_id, flag_fcm);
					ccf_arr[0] = w*R*u*bs_val_arr[0]*fcm_val_arr[0];
					ccf_arr[1] = w*R*u*bs_val_arr[0]*fcm_val_arr[1];
					ccf_arr[2] = w*R*u*bs_val_arr[0]*fcm_val_arr[2];

					col_edge_terms_qd[ix][ie][a][b].resize(3);
					col_edge_terms_qd[ix][ie][a][b][0] = ccf_arr[0];
					col_edge_terms_qd[ix][ie][a][b][1] = ccf_arr[1];
					col_edge_terms_qd[ix][ie][a][b][2] = ccf_arr[2];
				}
			}
		}

	}


	basisXVE_qd.resize(2);
	for(int i_op = 0; i_op < 2; i_op++)
	{
		basisXVE_qd[i_op].resize(qd_sizeX);
		for(int a = 0; a < qd_sizeX; a++)
		{
			value_type wa;
			Point2 pa, qa;
			vector<value_type> basisX(dof_x);

			tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ElementX(0), a); 
			basis->col_spatial_basis_vec(qa, basisX);

			basisXVE_qd[i_op][a].resize(qd_sizeVE);

			for(int b = 0; b < qd_sizeVE; b++)
			{
				value_type wb;
				Point2 pb, qb;
				vector<value_type> basisVE(dof_v*2*3);

				tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(EdgeV(0), b);
				basis->col_velocity_basis_vec_edge(i_op, qb, dvp, du, basisVE);

				basisXVE_qd[i_op][a][b].resize(2);
				for(int l = 0; l < 2; l++)
				{
					basisXVE_qd[i_op][a][b][l].resize(3);
					for(int k = 0; k < 3; k++)
					{
						basisXVE_qd[i_op][a][b][l][k].resize(dof);
						for(int i = 0; i < dof_x; i++)
						{
							value_type tmp_basisX = basisX[i];
							for(int j = 0; j < dof_v; j++)
							{
								int loc_index = l*3*dof_v + j*3 +k;
								basisXVE_qd[i_op][a][b][l][k][j*dof_x + i] = tmp_basisX*basisVE[loc_index];
							}
						}
					}
				}
			}
		}
	}


	col_edge_type_num = basis->get_col_edge_type_num();

	col_fhat_mat_num.resize(col_edge_type_num);
	col_fhat_mat_half_size.resize(col_edge_type_num);
	col_fhat_org_mat_group.resize(col_edge_type_num);
	col_fhat_org_mat_index.resize(col_edge_type_num);
	col_fhat_new_mat_group.resize(col_edge_type_num);
	col_fhat_new_mat_index.resize(col_edge_type_num);

	col_fhat_org_mat_inv_info.resize(col_edge_type_num);
	col_fhat_new_mat_inv_info.resize(col_edge_type_num);

	for(int i = 0; i < col_edge_type_num; i++)
	{
		basis->get_col_fhat_data(i, col_fhat_mat_num[i], col_fhat_mat_half_size[i], col_fhat_org_mat_group[i], col_fhat_org_mat_index[i], col_fhat_new_mat_group[i], col_fhat_new_mat_index[i]);

		col_fhat_org_mat_inv_info[i].resize(col_fhat_mat_num[i]);
		col_fhat_new_mat_inv_info[i].resize(col_fhat_mat_num[i]);

		for(int j = 0; j < col_fhat_mat_num[i]; j++)
		{
			int mat_half_size = col_fhat_mat_half_size[i][j];
			col_fhat_org_mat_inv_info[i][j].resize(mat_half_size);

			for(int k = 0; k < mat_half_size; k++)
			{
				for(int l = 0; l < col_fhat_org_mat_group[i].size(); l++)
				{
					int tmp_g = col_fhat_org_mat_group[i][l];
					int tmp_i = col_fhat_org_mat_index[i][l];

					if (tmp_g == j && tmp_i == k)
						col_fhat_org_mat_inv_info[i][j][k] = l;
				}
			}


			col_fhat_new_mat_inv_info[i][j].resize(mat_half_size*2);
			for(int k = 0; k < mat_half_size*2; k++)
			{
				for(int l = 0; l < col_fhat_new_mat_group[i].size(); l++)
				{
					int tmp_g = col_fhat_new_mat_group[i][l];
					int tmp_i = col_fhat_new_mat_index[i][l];

					if (tmp_g == j && tmp_i == k)
						col_fhat_new_mat_inv_info[i][j][k] = l;
				}
			}
		}
	}
}

void Integration_Col::vol_cal_weight(const ElementX &ti, const int j, const int a, vector<value_type> &val_arr)
{

	value_type w, wa, loc_R, loc_Z, tmp_val=0.0;
	Point4 q;
	Point2 pa, qa;

	tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

	q = qa & Point2({0, 0});
	auto tmp_vec = (*basis).ind_val(q, j, 0);
	w = wa * tmp_vec[0];
	loc_R = pa[0];
	loc_Z = pa[1];

	val_arr[0] = w;
	val_arr[1] = loc_R;
	val_arr[2] = loc_Z;
}

void Integration_Col::col_cell_mat_cal(const ElementX &ti, const ElementV &tj, vector<value_type> &val_arr)
{
	value_type val = 0.0;
	value_type w, wa, wb, R, Z, v_para, u, bs_val, dbs_dvp_ov_bs, b_val, bphi_val;
	Point4 p, q;
	Point2 pa, pb, qa, qb;
	value_type tmp_val_arr[diag_num];
	value_type Ms = species_data->normalized_mass(sp_id);
	value_type Zs = species_data->normalized_charge(sp_id);
	value_type BPT_SIGN = eq_reader->get_property(EQ_reader::bpt_sign);
	value_type rh0n = eq_reader->ph_const.get_property(Ph_const::rh0n);
	value_type ccf2 = BPT_SIGN*rh0n*(Ms/Zs);

	Vector B(3);
	value_type B_val, Bphi_val, tmp_w1;

	val_arr.assign(val_arr.size(), 0.0);
	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		p = pa & Point2({0, 0});
		B = flux->get_B(p);
		B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));
		Bphi_val = B(2);

		R = p[0];
		Z = p[1];

		auto psi_val = eq_reader->psi_interpol(R,Z);
		value_type I_val = eq_reader->I_interpol(R,Z);

		for(int b = 0; b < quadrature->size(1); b++)
		{
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

			w = wa * wb;
			p = pa & pb;
			q = qa & qb;

			v_para = p[2];
			u = p[3];

			vector<value_type> bs_val_arr(2);
			bs_val_arr = flux->get_b_dot_Bs_vec(p, B, sp_id);
			bs_val = bs_val_arr[0];
			dbs_dvp_ov_bs = bs_val_arr[1]/bs_val;

			tmp_val_arr[0] = 1.0;
			tmp_val_arr[1] = v_para;
			tmp_val_arr[2] = 0.5*v_para*v_para + 0.5*u*u*B_val;
			tmp_val_arr[3] = dbs_dvp_ov_bs;
			tmp_val_arr[4] = 3.0 + v_para*dbs_dvp_ov_bs;
			tmp_val_arr[5] = 0.5*v_para*v_para;
			tmp_val_arr[6] = 0.5*u*u*B_val;

			auto fcm_val = flux->get_fcm_SE_tot(p, int(ti), sp_id, flag_fcm);

			tmp_w1 = w * bs_val * R * u * fcm_val;
			tmp_w1 *= 4.0*M_PI*M_PI;

			for(int j = 0; j < dof; j++)
			{
				val = tmp_w1*(*basis)(q, j);
				for(int val_id = 0; val_id < diag_num; val_id++) val_arr[j*diag_num + val_id] += val*tmp_val_arr[val_id];
			}
		}
	}
}

void Integration_Col::col_cell_consv_mat_cal(const ElementX &ti, const ElementV &tj, const int &col_method, vector<value_type> &val_arr, const int &col_consv_can_ang_mom_op)
{
	value_type w, wa, wb, R, Z, v_para, u, bs_val, dbs_dvp_ov_bs, b_val, bphi_val;
	Point4 p, q;
	Point2 pa, pb, qa, qb;
	value_type Ms = species_data->normalized_mass(sp_id);
	value_type Zs = species_data->normalized_charge(sp_id);
	value_type BPT_SIGN = eq_reader->get_property(EQ_reader::bpt_sign);
	value_type rh0n = eq_reader->ph_const.get_property(Ph_const::rh0n);
	value_type ccf2 = BPT_SIGN*rh0n*(Ms/Zs);

	int col_consv_num;
	if (col_method == 2) col_consv_num = 3;
	else col_consv_num = 2;
	value_type tmp_val_arr[col_consv_num];

	val_arr.assign(val_arr.size(), 0.0);
	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));
		value_type Bphi_val = B(2);

		R = p[0];
		Z = p[1];

		auto psi_val = eq_reader->psi_interpol(R,Z);
		value_type I_val = eq_reader->I_interpol(R,Z);

		for(int b = 0; b < quadrature->size(1); b++)
		{
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

			w = wa * wb;
			p = pa & pb;
			q = qa & qb;

			v_para = p[2];
			u = p[3];
			vector<value_type> bs_val_arr(2);
			bs_val_arr = flux->get_b_dot_Bs_vec(p, B, sp_id);
			bs_val = bs_val_arr[0];
			dbs_dvp_ov_bs = bs_val_arr[1]/bs_val;

			if (col_consv_can_ang_mom_op == 0) tmp_val_arr[0] = v_para;
			else
			{
				auto psi_h = psi_val - ccf2*I_val/B_val*v_para;
				tmp_val_arr[0] = psi_h;
			}


			tmp_val_arr[1] = 0.5*v_para*v_para + 0.5*u*u*B_val;
			if (col_method == 2) tmp_val_arr[2] = 1.0;
			auto fcm_val = flux->get_fcm_SE_tot(p, int(ti), sp_id, flag_fcm);

			value_type tmp_w1 = w * bs_val * R * u * fcm_val;
			tmp_w1 *= 4.0*M_PI*M_PI;

			for(int j = 0; j < dof; j++)
			{
				value_type val = tmp_w1*(*basis)(q, j);
				for(int val_id = 0; val_id < col_consv_num; val_id++) val_arr[j*col_consv_num + val_id] += val*tmp_val_arr[val_id]*Ms;
			}
		}
	}
}

void Integration_Col::col_f_dg_to_valid_fhat_mat_cal(const ElementX &ti, const EdgeV &te, const ElementV &tj0, const ElementV &tj1, const int &i_op, Vector &fhat_mat_arr)
{
	value_type w, R, u, bs_val, wa, wb;
	Point4 p, q;
	Point2 pa, pb, qa, qb;
	int flag_fcm = species_data->fcm_type(sp_id);
	value_type tmp_w1, tmp_wf0, tmp_w2, tmp_inner_product[dof], tmp_basis_val[dof];

	value_type tmp_x_basis[dof_x];
	value_type tmp_v_org_basis[dof_v];

	int dim_fhat_org_v_basis = basis->fhat_org_vbasis_dim_out();

	vector<Matrix> M_col, R_col;

	for (int i = 0; i < dof_x; i++)
	{
		M_col.push_back(Matrix(dim_fhat_org_v_basis, dim_fhat_org_v_basis));
		R_col.push_back(Matrix(dim_fhat_org_v_basis, dof_v*2));

		M_col[i].setZero();
		R_col[i].setZero();
	}


	//volume integration
	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 
		qb = {0.0, 0.0};
		q = qa & qb;
		for(int i = 0; i < dof_x; i++) 
		{
			tmp_x_basis[i] = basis->col_spatial_basis(q, i);
		}

		for(int iele = 0; iele < 2; iele++)
		{
			vector<Point2> ele_node_arr;
			if (iele == 0)
			{
				ele_node_arr = mesh->get_node_element(tj0);
			}
			else
			{
				ele_node_arr = mesh->get_node_element(tj1);
			}

			value_type min_vp = 1e10, max_vp = -1e10;
			value_type min_u = 1e10, max_u = -1e10;

			for(int t1 = 0; t1 < ele_node_arr.size(); t1++)
			{
				min_vp = min(min_vp, ele_node_arr[t1][0]);
				max_vp = max(max_vp, ele_node_arr[t1][0]);
				min_u = min(min_u, ele_node_arr[t1][1]);
				max_u = max(max_u, ele_node_arr[t1][1]);
			}
			value_type del_vp = max_vp - min_vp, del_u = max_u - min_u;

			for(int b = 0; b < quadrature->size(1); b++)
			{
				if (iele == 0) tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj0, b);
				else tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj1, b);

				w = wa * wb;

				if(i_op == 0)
				{
					if(iele == 0)
					{
						qb[0] = 1.0 + fhat_vol_wid*(qb[0] - 1.0);
					}
					else
					{
						qb[0] = -1.0 + fhat_vol_wid*(qb[0] + 1.0);
					}
					pb[0] = min_vp + 0.5*del_vp*(qb[0] + 1.0);
				}
				else
				{
					if(iele == 0)
					{
						qb[1] = 1.0 + fhat_vol_wid*(qb[1] - 1.0);
					}
					else
					{
						qb[1] = -1.0 + fhat_vol_wid*(qb[1] + 1.0);
					}
					pb[1] = min_u + 0.5*del_u*(qb[1] + 1.0);
				}

				p = pa & pb;
				q = qa & qb;

				R = p[0];
				u = p[3];
				bs_val = flux->get_b_dot_Bs(p, sp_id);
				tmp_w1 = w * bs_val * R * u;

				auto fcm_val = flux->get_fcm_SE_tot(p, int(ti), sp_id, flag_fcm);

				tmp_wf0 = tmp_w1 * fcm_val;

				value_type tmp_wf0_fhat = tmp_wf0;

				for(int j = 0; j < dof_v; j++) tmp_v_org_basis[j] = basis->col_velocity_org_basis(q, j);
				Vector fhat_basis_Vec = basis->fhat_vbasis_Vec(q, i_op, iele);

				int flag_vol_within = 0;

				flag_vol_within = 1;

				if (flag_vol_within == 1)
				{
					for(int i = 0; i < dof_x; i++) 
					{
						for(int k = 0; k < dim_fhat_org_v_basis; k++) 
						{
							for(int l = 0; l < dim_fhat_org_v_basis; l++)
							{
								M_col[i].coeffRef(k, l) += tmp_x_basis[i]*fhat_basis_Vec[k]*fhat_basis_Vec[l]*tmp_wf0_fhat;
							}

							for(int l = 0; l < dof_v; l++)
							{
								int tmp_l = iele*dof_v + l;
								R_col[i].coeffRef(k, tmp_l) += tmp_x_basis[i]*fhat_basis_Vec[k]*tmp_v_org_basis[l]*tmp_wf0;

							}
						}
					}
				}
			}
		}
	}

	int dim_fhat_valid_v_basis = basis->fhat_valid_vbasis_dim_out();
	//X basis index
	for (int i = 0; i < dof_x; i++)
	{
		Matrix M_inv = M_col[i].fullPivLu().inverse();
		Matrix M_inv_R = M_inv*R_col[i];

		for(int col = 0; col < 2*dof_v; col++) 
		{
			int f_dg_basis_index = col*dof_x + i;
			for(int row = 0; row < dim_fhat_valid_v_basis; row++) 
			{
				int tmp_index = (i*2*dof + f_dg_basis_index)*dim_fhat_valid_v_basis + row;

				int fhat_org_basis_i = basis->fhat_valid_to_org_i(row, i_op);
				fhat_mat_arr[tmp_index] = M_inv_R(fhat_org_basis_i, col);
			}
		}
	}
}

void Integration_Col::col_fhat_mat_cal(const ElementX &ti, const EdgeV &te, const ElementV &tj0, const ElementV &tj1, const int &i_op, vector<value_type> &fhat_mat_arr, vector<int> &fhat_mat_info)
{
	value_type w, R, u, bs_val, wa, wb;
	Point4 p, q;
	Point2 pa, pb, qa, qb;
	int flag_fcm = species_data->fcm_type(sp_id);
	value_type tmp_w1, tmp_wf0, tmp_w2, tmp_inner_product[dof], tmp_basis_val[dof];

	value_type tmp_x_basis[dof_x];
	value_type tmp_v_org_basis[dof_v];

	int dim_fhat_v_basis = basis->get_dofv_fhat();

	value_type tmp_v_new_basis[dim_fhat_v_basis];
	vector<SparseMatrix> M_col, R_col;


	for (int i = 0; i < dof_x; i++)
	{
		M_col.push_back(SparseMatrix(dim_fhat_v_basis, dim_fhat_v_basis));
		R_col.push_back(SparseMatrix(dim_fhat_v_basis, dof_v*2));

		M_col[i].setZero();
		R_col[i].setZero();
	}

	fhat_mat_arr.assign(fhat_mat_arr.size(), 0.0);
	memset(&fhat_mat_info[0], -1, fhat_mat_info.size() * sizeof fhat_mat_info[0]);

	//volume integration
	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 
		qb = {0.0, 0.0};
		q = qa & qb;
		for(int i = 0; i < dof_x; i++) 
		{
			tmp_x_basis[i] = basis->col_spatial_basis(q, i);
		}

		for(int iele = 0; iele < 2; iele++)
		{
			vector<Point2> ele_node_arr;
			if (iele == 0)
			{
				ele_node_arr = mesh->get_node_element(tj0);
			}
			else
			{
				ele_node_arr = mesh->get_node_element(tj1);
			}

			value_type min_vp = 1e10, max_vp = -1e10;
			value_type min_u = 1e10, max_u = -1e10;

			for(int t1 = 0; t1 < ele_node_arr.size(); t1++)
			{
				min_vp = min(min_vp, ele_node_arr[t1][0]);
				max_vp = max(max_vp, ele_node_arr[t1][0]);
				min_u = min(min_u, ele_node_arr[t1][1]);
				max_u = max(max_u, ele_node_arr[t1][1]);
			}
			value_type del_vp = max_vp - min_vp, del_u = max_u - min_u;

			for(int b = 0; b < quadrature->size(1); b++)
			{
				if (iele == 0) tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj0, b);
				else tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj1, b);

				w = wa * wb;

				if(i_op == 0)
				{
					if(iele == 0)
					{
						qb[0] = 1.0 + fhat_vol_wid*(qb[0] - 1.0);
					}
					else
					{
						qb[0] = -1.0 + fhat_vol_wid*(qb[0] + 1.0);
					}
					pb[0] = min_vp + 0.5*del_vp*(qb[0] + 1.0);
				}
				else
				{
					if(iele == 0)
					{
						qb[1] = 1.0 + fhat_vol_wid*(qb[1] - 1.0);
					}
					else
					{
						qb[1] = -1.0 + fhat_vol_wid*(qb[1] + 1.0);
					}
					pb[1] = min_u + 0.5*del_u*(qb[1] + 1.0);
				}
				p = pa & pb;
				q = qa & qb;

				R = p[0];
				u = p[3];
				bs_val = flux->get_b_dot_Bs(p, sp_id);
				tmp_w1 = w * bs_val * R * u;

				auto fcm_val = flux->get_fcm_SE_tot(p, int(ti), sp_id, flag_fcm);

				tmp_wf0 = tmp_w1 * fcm_val;

				value_type tmp_wf0_fhat = tmp_wf0;

				for(int j = 0; j < dof_v; j++) tmp_v_org_basis[j] = basis->col_velocity_org_basis(q, j);
				for(int j = 0; j < dim_fhat_v_basis; j++) tmp_v_new_basis[j] = basis->col_velocity_fhat_basis(q, j, i_op, iele);

				int flag_vol_within = 0;

				flag_vol_within = 1;

				if (flag_vol_within == 1)
				{
					for(int i = 0; i < dof_x; i++) 
					{
						for(int k = 0; k < dim_fhat_v_basis; k++) 
						{
							for(int l = 0; l < dim_fhat_v_basis; l++)
							{
								M_col[i].coeffRef(k, l) += tmp_x_basis[i]*tmp_v_new_basis[k]*tmp_v_new_basis[l]*tmp_wf0_fhat;
							}

							for(int l = 0; l < dof_v; l++)
							{
								int tmp_l = iele*dof_v + l;
								R_col[i].coeffRef(k, tmp_l) += tmp_x_basis[i]*tmp_v_new_basis[k]*tmp_v_org_basis[l]*tmp_wf0;
							}
						}
					}
				}
			}
		}
	}

	Eigen::PardisoLU<SparseMatrix> solver;

	//X basis index
	for (int i = 0; i < dof_x; i++)
	{
		solver.analyzePattern(M_col[i]);
		solver.factorize(M_col[i]);

		Vector X, F(dim_fhat_v_basis);

		//element index
		for(int k = 0; k < 2; k++) 
		{
			//org velocity basis index
			for(int l = 0; l < dof_v; l++)
			{
				for(int m = 0; m < dim_fhat_v_basis; m++)
				{
					F[m] = R_col[i].coeffRef(m, k*dof_v + l);
				}
				//org velocity basis l -> new fhat velocity basis
				X = solver.solve(F);


				//colume index of Matrix for org basis -> fhat value at the edge
				int col_pos = k*dof + dof_x*l + i;

				//new fhat velocity basis index
				for(int m = 0; m < dim_fhat_v_basis; m++)
				{
					//int tmp_index = col_fhat_new_mat_inv_info[i_op][j][m];

					//fhat derivative index -> n=0 : fhat, n=1 : dfhat/dvp, n=2 : dfhat/du
					for(int n = 0; n < 3; n++)
					{
						//when i_op == 0 -> fhat_out_pos : u dependence of new fhat velocity basis m
						//when i_op == 1 -> fhat_out_pos : vp dependence of new fhat velocity basis m
						int fhat_out_pos;
						value_type fhat_factor;
						basis->col_velocity_fhat_basis_info(m, i_op, n, dvp, du, fhat_out_pos, fhat_factor);

						if (fhat_out_pos > -1)
						{
							int row_pos, pos_1d;
							row_pos = i*3*dof_v_1d + n*dof_v_1d + fhat_out_pos;

							pos_1d = col_pos*fhat_dim_per_edge + row_pos;
							fhat_mat_info[pos_1d] = 1;
							fhat_mat_arr[pos_1d] = fhat_factor*X[m];
						}
					}
				}
			}
		}
	}
}

void Integration_Col::col_fhat_mat_cal_bd(const ElementX &ti, const EdgeV &te, const ElementV &tj0, const ElementV &tj1, const int &i_op, vector<value_type> &fhat_mat_arr, vector<int> &fhat_mat_info)
{

	int ele_pos;
	value_type vp_val, u_val;

	if (int(tj0) == -1) ele_pos = 1;
	else if (int(tj1) == -1) ele_pos = 0;
	else 
	{
		cout << "wrong tj0, tj1 in col_fhat_mat_cal_bd : " << int(tj0) << " " << int(tj1) << endl;
		abort();
	}

	fhat_mat_arr.assign(fhat_mat_arr.size(), 0.0);
	fhat_mat_info.assign(fhat_mat_info.size(), -1);

	if (i_op == 2)
	{
		if (ele_pos == 0) vp_val = 1.0;
		else vp_val = -1.0;
	}
	else if (i_op == 3)
	{
		if (ele_pos == 0) u_val = 1.0;
		else u_val = -1.0;
	}
	else 
	{
		cout << "wrong i_op in col_fhat_mat_cal_bd : " << i_op << endl;
		abort();
	}

	//X basis index
	for (int i = 0; i < dof_x; i++)
	{
		//org velocity basis index
		for(int l = 0; l < dof_v; l++)
		{
			int col_pos = ele_pos*dof + dof_x*l + i;

			int fhat_out_pos;
			value_type fhat_factor;

			basis->col_velocity_fhat_basis_bd_info(l, i_op, vp_val, u_val, fhat_out_pos, fhat_factor);

			//fhat derivative index -> n=0 : fhat, n=1 : dfhat/dvp, n=2 : dfhat/du
			//for last boundary, only n=0 is required.
			int n = 0; 
			int row_pos = i*3*dof_v_1d + n*dof_v_1d + fhat_out_pos;
			int pos_1d = col_pos*fhat_dim_per_edge + row_pos;

			fhat_mat_info[pos_1d] = 1;
			fhat_mat_arr[pos_1d] = fhat_factor;
		}
	}
}

void Integration_Col::col3_init(const int &_tot_species_num, const int &_sp_id, const vector<int> &_vp_n_arr, const vector<int> &_u_n_arr, const vector<value_type> &_vp_min_arr, const vector<value_type> &_vp_max_arr, const vector<value_type> &_u_min_arr, const vector<value_type> &_u_max_arr, const int &_ei_pitch_angle_col_op, const int &_ei_pitch_no_v_in_nu_ei_op, const value_type &_ei_pitch_no_v_in_nu_ei_v_e_min)
{
	tot_species_num = _tot_species_num;
	sp_own_id = _sp_id;
	ei_pitch_angle_col_op =_ei_pitch_angle_col_op;
	ei_pitch_no_v_in_nu_ei_op = _ei_pitch_no_v_in_nu_ei_op;
	ei_pitch_no_v_in_nu_ei_v_e_min = _ei_pitch_no_v_in_nu_ei_v_e_min;

	if (ei_pitch_no_v_in_nu_ei_op != 0 && ei_pitch_no_v_in_nu_ei_op != 1 && ei_pitch_no_v_in_nu_ei_op != 2)
	{
		cout << "Unknown ei_pitch_no_v_in_nu_ei_op : " << ei_pitch_no_v_in_nu_ei_op << endl;
		abort();
	}

	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
	{
		int vp_n = _vp_n_arr[sp_a_id];
		int u_n = _u_n_arr[sp_a_id];

		vp_n_arr.push_back(vp_n);
		u_n_arr.push_back(u_n);

		hg_inner_n.push_back((vp_n + 1)*(u_n + 1));
		hg_total_n.push_back((vp_n + 3)*(u_n + 3));

		value_type vp_min = _vp_min_arr[sp_a_id];
		value_type vp_max = _vp_max_arr[sp_a_id];
		value_type u_min = _u_min_arr[sp_a_id];
		value_type u_max = _u_max_arr[sp_a_id];

		vp_min_arr.push_back(vp_min);
		vp_max_arr.push_back(vp_max);
		u_min_arr.push_back(u_min);
		u_max_arr.push_back(u_max);

		dvp_arr.push_back((vp_max - vp_min)/vp_n);
		du_arr.push_back((u_max - u_min)/u_n);
		inv_dvp_arr.push_back(vp_n/(vp_max - vp_min));
		inv_du_arr.push_back(u_n/(u_max - u_min));

		int sp_b_id = sp_a_id;
		value_type Ma = species_data->normalized_mass(sp_own_id);
		value_type Mb = species_data->normalized_mass(sp_b_id);
		Mb_ov_Ma_arr.push_back(Mb/Ma);
	}

	col3_int_vp_n = vp_n_arr[sp_own_id];
	col3_int_u_n = u_n_arr[sp_own_id];
	col3_int_org_f_vp_n = species_data->n_vp(sp_own_id);
	col3_int_org_f_u_n = species_data->n_u(sp_own_id);

	col3_int_vp_min = vp_min_arr[sp_own_id];
	col3_int_vp_max = vp_max_arr[sp_own_id];
	col3_int_u_min = u_min_arr[sp_own_id];
	col3_int_u_max = u_max_arr[sp_own_id];

	col3_int_dvp = (col3_int_vp_max - col3_int_vp_min)/col3_int_vp_n;
	col3_int_du = (col3_int_u_max - col3_int_u_min)/col3_int_u_n;

	col3_int_vp_tot_min = col3_int_vp_min - col3_int_dvp;
	col3_int_u_tot_min = col3_int_u_min - col3_int_du;

	col3_int_inv_dvp = 1.0/col3_int_dvp;
	col3_int_inv_du = 1.0/col3_int_du;

	col3_int_hg_tot_n = (col3_int_vp_n + 3)*(col3_int_u_n + 3);
	col3_int_hg_inner_n = (col3_int_vp_n + 1)*(col3_int_u_n + 1);
	col3_int_spline_half_wid = 2.0;

	int nv_loc = mesh->size<ElementV>();
	int qd_sizeX = quadrature->size(0), qd_sizeV = quadrature->size(1);
	int qd_sizeVE = quadrature->size(3);

	col3_int_h_interpol_num = 2;
	col3_int_g_interpol_num = 4;

	col3_int_vol_hg_interpol_mat.resize(nv_loc);
	for(int iv = 0; iv < nv_loc; iv++)
	{
		vector<SparseMatrix_Triplet> tripletList_h, tripletList_g;
		tripletList_h.reserve(qd_sizeV*col3_int_h_interpol_num*100);
		tripletList_g.reserve(qd_sizeV*col3_int_g_interpol_num*100);

		for(int b = 0; b < qd_sizeV; b++)
		{
			value_type wb;
			Point2 pb, qb;
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ElementV(iv), b);

			value_type v_para = pb[0], u = pb[1];

			int loc_vp_tot_index = int((v_para - col3_int_vp_tot_min)*col3_int_inv_dvp);
			int loc_u_tot_index = int((u - col3_int_u_tot_min)*col3_int_inv_du);

			int vp_tot_min_index = max(0, loc_vp_tot_index - 1);
			int vp_tot_max_index = min(col3_int_vp_n + 2, loc_vp_tot_index + 2);
			int u_tot_min_index = max(0, loc_u_tot_index - 1);
			int u_tot_max_index = min(col3_int_u_n + 2, loc_u_tot_index + 2);

			for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
			{
				value_type t1_tot = (v_para - (col3_int_vp_tot_min + j3*col3_int_dvp))*col3_int_inv_dvp;
				Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_int_inv_dvp);

				for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
				{
					value_type t2_tot = (u - (col3_int_u_tot_min + j4*col3_int_du))*col3_int_inv_du;
					Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_int_inv_du);

					Matrix tmp_val = Matrix::Zero(3,3);

					for(int i1 = 0; i1 < 3; i1++)
					{
						for(int i2 = 0; i2 < 3; i2++)
						{
							tmp_val(i1, i2) = spline_vp_tot_vec(i1)*spline_u_tot_vec(i2);
						}
					}

					int l = j4*(col3_int_vp_n + 3) + j3;
					int kh = col3_int_h_interpol_num*b;
					int kg = col3_int_g_interpol_num*b;

					tripletList_h.push_back(SparseMatrix_Triplet(kh + 0, l, tmp_val(1,0)));
					tripletList_h.push_back(SparseMatrix_Triplet(kh + 1, l, tmp_val(0,1)));
					tripletList_g.push_back(SparseMatrix_Triplet(kg + 0, l, tmp_val(0,1)));
					tripletList_g.push_back(SparseMatrix_Triplet(kg + 1, l, tmp_val(2,0)));
					tripletList_g.push_back(SparseMatrix_Triplet(kg + 2, l, tmp_val(1,1)));
					tripletList_g.push_back(SparseMatrix_Triplet(kg + 3, l, tmp_val(0,2)));
				}
			}
		}

		col3_int_vol_hg_interpol_mat[iv].push_back(SparseMatrix(col3_int_h_interpol_num*qd_sizeV, col3_int_hg_tot_n));
		col3_int_vol_hg_interpol_mat[iv].push_back(SparseMatrix(col3_int_g_interpol_num*qd_sizeV, col3_int_hg_tot_n));

		col3_int_vol_hg_interpol_mat[iv][0].setFromTriplets(tripletList_h.begin(), tripletList_h.end());
		col3_int_vol_hg_interpol_mat[iv][1].setFromTriplets(tripletList_g.begin(), tripletList_g.end());

	}

	int nve_loc = mesh->size<EdgeV>();
	col3_int_edge_hg_interpol_mat.resize(nve_loc);
	for(int ie = 0; ie < nve_loc; ie++)
	{
		vector<SparseMatrix_Triplet> tripletList_h, tripletList_g;
		tripletList_h.reserve(qd_sizeVE*col3_int_h_interpol_num*100);
		tripletList_g.reserve(qd_sizeVE*col3_int_g_interpol_num*100);

		for(int b = 0; b < qd_sizeVE; b++)
		{
			value_type wb;
			Point2 pb, qb;
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(EdgeV(ie), b);

			value_type v_para = pb[0], u = pb[1];

			int loc_vp_tot_index = int((v_para - col3_int_vp_tot_min)*col3_int_inv_dvp);
			int loc_u_tot_index = int((u - col3_int_u_tot_min)*col3_int_inv_du);

			int vp_tot_min_index = max(0, loc_vp_tot_index - 1);
			int vp_tot_max_index = min(col3_int_vp_n + 2, loc_vp_tot_index + 2);
			int u_tot_min_index = max(0, loc_u_tot_index - 1);
			int u_tot_max_index = min(col3_int_u_n + 2, loc_u_tot_index + 2);

			for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
			{
				value_type t1_tot = (v_para - (col3_int_vp_tot_min + j3*col3_int_dvp))*col3_int_inv_dvp;
				Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_int_inv_dvp);

				for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
				{
					value_type t2_tot = (u - (col3_int_u_tot_min + j4*col3_int_du))*col3_int_inv_du;
					Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_int_inv_du);

					Matrix tmp_val = Matrix::Zero(3,3);

					for(int i1 = 0; i1 < 3; i1++)
					{
						for(int i2 = 0; i2 < 3; i2++)
						{
							tmp_val(i1, i2) = spline_vp_tot_vec(i1)*spline_u_tot_vec(i2);
						}
					}

					int l = j4*(col3_int_vp_n + 3) + j3;
					int kh = col3_int_h_interpol_num*b;
					int kg = col3_int_g_interpol_num*b;

					tripletList_h.push_back(SparseMatrix_Triplet(kh + 0, l, tmp_val(1,0)));
					tripletList_h.push_back(SparseMatrix_Triplet(kh + 1, l, tmp_val(0,1)));
					tripletList_g.push_back(SparseMatrix_Triplet(kg + 0, l, tmp_val(0,1)));
					tripletList_g.push_back(SparseMatrix_Triplet(kg + 1, l, tmp_val(2,0)));
					tripletList_g.push_back(SparseMatrix_Triplet(kg + 2, l, tmp_val(1,1)));
					tripletList_g.push_back(SparseMatrix_Triplet(kg + 3, l, tmp_val(0,2)));
				}
			}
		}

		col3_int_edge_hg_interpol_mat[ie].push_back(SparseMatrix(col3_int_h_interpol_num*qd_sizeVE, col3_int_hg_tot_n));
		col3_int_edge_hg_interpol_mat[ie].push_back(SparseMatrix(col3_int_g_interpol_num*qd_sizeVE, col3_int_hg_tot_n));

		col3_int_edge_hg_interpol_mat[ie][0].setFromTriplets(tripletList_h.begin(), tripletList_h.end());
		col3_int_edge_hg_interpol_mat[ie][1].setFromTriplets(tripletList_g.begin(), tripletList_g.end());

	}

}

void Integration_Col::col3_bc_mat_cal(const ElementX &ti, const ElementV &tj, const int &sp_id, const int &bc_n, const vector<value_type> &bc_vp_arr, const vector<value_type> &bc_u_arr, const int &hg_op, vector<value_type> &val_arr)
{
	value_type w, wa, wb;
	Point4 p, q;
	Point2 pa, pb, qa, qb;

	val_arr.assign(val_arr.size(), 0.0);
	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));

		value_type R = p[0], Z = p[1];

		for(int b = 0; b < quadrature->size(1); b++)
		{
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

			w = wa * wb;
			p = pa & pb;
			q = qa & qb;

			value_type v_para = p[2], u = p[3];

			auto fcm_val = flux->get_fcm_SE_tot(p, int(ti), sp_id, flag_fcm);

			value_type tmp_w1 = w * B_val * R * u * fcm_val;
			value_type tmp_w_arr[dof];
			for(int j = 0; j < dof; j++) tmp_w_arr[j] = tmp_w1*(*basis)(q, j); 

			for(int m = 0; m < bc_n; m++)
			{
				value_type vp_m = bc_vp_arr[m];
				value_type u_m = bc_u_arr[m];

				value_type tmp1 = B_val*(u + u_m)*(u + u_m) + (v_para - vp_m)*(v_para - vp_m);
				value_type k_val = sqrt(B_val*(4.0*u*u_m)/tmp1);

				value_type K_int_val = elliptic_fk(k_val);
				value_type E_int_val = elliptic_ek(k_val);

				for(int i = 0; i < dof; i++)
				{
					int mi_index = m*dof + i;
					if (hg_op == 0)
					{
						val_arr[mi_index] += 2.0*tmp_w_arr[i]*K_int_val/sqrt(tmp1);
					}
					else
					{
						val_arr[mi_index] += tmp_w_arr[i]*E_int_val*sqrt(tmp1);
					}
				}
			}
		}
	}
}

void Integration_Col::col3_hg_vol_mat_cal(const ElementX &ti, const ElementV &tj, const int &sp_id, vector<value_type> &val_arr, vector<value_type> &val_arr2, vector<value_type> &val_arr3)
{
	value_type w, wa, wb, v_para, u;
	Point4 p, q;
	Point2 pa, pb, qa, qb;

	vector<vector<int>> qd_to_hg_bin_inner_index, qd_to_hg_bin_tot_index;
	for(int b = 0; b < quadrature->size(1); b++)
	{
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

		v_para = pb[0];
		u = pb[1];

		int loc_vp_inner_index = min(col3_int_vp_n - 1, max(0, int((v_para - col3_int_vp_min)*col3_int_inv_dvp)));
		int loc_u_inner_index = min(col3_int_u_n - 1, max(0, int((u - col3_int_u_min)*col3_int_inv_du)));

		int loc_vp_tot_index = loc_vp_inner_index + 1;
		int loc_u_tot_index = loc_u_inner_index + 1;

		vector<int> tmp_qd_to_hg_inner_index(4), tmp_qd_to_hg_tot_index(4);
		tmp_qd_to_hg_inner_index[0] = max(0, loc_vp_inner_index - 1);
		tmp_qd_to_hg_inner_index[1] = min(col3_int_vp_n, loc_vp_inner_index + 2);
		tmp_qd_to_hg_inner_index[2] = max(0, loc_u_inner_index - 1);
		tmp_qd_to_hg_inner_index[3] = min(col3_int_u_n, loc_u_inner_index + 2);

		tmp_qd_to_hg_tot_index[0] = max(0, loc_vp_tot_index - 1);
		tmp_qd_to_hg_tot_index[1] = min(col3_int_vp_n + 2, loc_vp_tot_index + 2);
		tmp_qd_to_hg_tot_index[2] = max(0, loc_u_tot_index - 1);
		tmp_qd_to_hg_tot_index[3] = min(col3_int_u_n + 2, loc_u_tot_index + 2);

		qd_to_hg_bin_inner_index.push_back(tmp_qd_to_hg_inner_index);
		qd_to_hg_bin_tot_index.push_back(tmp_qd_to_hg_tot_index);
	}

	val_arr.assign(val_arr.size(), 0.0);
	val_arr2.assign(val_arr2.size(), 0.0);
	val_arr3.assign(val_arr3.size(), 0.0);
	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));

		value_type R = p[0], Z = p[1];

		for(int b = 0; b < quadrature->size(1); b++)
		{
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

			w = wa * wb;
			p = pa & pb;
			q = qa & qb;

			v_para = p[2];
			u = p[3];

			auto fcm_val = flux->get_fcm_SE_tot(p, int(ti), sp_id, flag_fcm);

			value_type tmp_w1 = 2.0 * M_PI * w * R * u;
			value_type tmp_w2 = tmp_w1 * B_val;
			value_type tmp_w = tmp_w2* fcm_val;

			value_type tmp_w_arr[dof];

			for(int i = 0; i < dof; i++) tmp_w_arr[i] = tmp_w*(*basis)(q, i); 

			for(int j1 = qd_to_hg_bin_inner_index[b][0]; j1 < qd_to_hg_bin_inner_index[b][1] + 1; j1++)
			{
				value_type t1_inner = (v_para - (col3_int_vp_min + j1*col3_int_dvp))*col3_int_inv_dvp;
				Vector spline_vp_inner_vec = Spline::cal_Cs_Vec(t1_inner, col3_int_inv_dvp);

				for(int j2 = qd_to_hg_bin_inner_index[b][2]; j2 < qd_to_hg_bin_inner_index[b][3] + 1; j2++)
				{
					value_type t2_inner = (u - (col3_int_u_min + j2*col3_int_du))*col3_int_inv_du;
					Vector spline_u_inner_vec = Spline::cal_Cs_Vec(t2_inner, col3_int_inv_du);

					int m = j2*(col3_int_vp_n + 1) + j1;

					//f to h source part
					for(int i = 0; i < dof; i++)
					{
						int mi_index = m*dof + i;

						val_arr[mi_index] += tmp_w_arr[i]*spline_vp_inner_vec(0)*spline_u_inner_vec(0);
					}

					for(int j3 = qd_to_hg_bin_tot_index[b][0]; j3 < qd_to_hg_bin_tot_index[b][1] + 1; j3++)
					{
						value_type t1_tot = (v_para - (col3_int_vp_tot_min + j3*col3_int_dvp))*col3_int_inv_dvp;
						Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_int_inv_dvp);

						for(int j4 = qd_to_hg_bin_tot_index[b][2]; j4 < qd_to_hg_bin_tot_index[b][3] + 1; j4++)
						{
							value_type t2_tot = (u - (col3_int_u_tot_min + j4*col3_int_du))*col3_int_inv_du;
							Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_int_inv_du);

							int l = j4*(col3_int_vp_n + 3) + j3;

							int ml_index1 = m*col3_int_hg_tot_n + l;
							int ml_index2 = ml_index1 + col3_int_hg_inner_n*col3_int_hg_tot_n;
							int ml_index3 = ml_index2 + col3_int_hg_inner_n*col3_int_hg_tot_n;

							//h to g source part
							val_arr2[ml_index1] += tmp_w2*spline_vp_inner_vec(0)*spline_u_inner_vec(0)*spline_vp_tot_vec(0)*spline_u_tot_vec(0);
							//h or g stiffness matrix part
							val_arr2[ml_index2] += tmp_w1*spline_vp_inner_vec(0)*spline_u_inner_vec(1)*spline_vp_tot_vec(0)*spline_u_tot_vec(1);
							val_arr2[ml_index3] += tmp_w2*spline_vp_inner_vec(1)*spline_u_inner_vec(0)*spline_vp_tot_vec(1)*spline_u_tot_vec(0);
						}
					}
				}
			}

			//conservation part
			vector<value_type> bs_val_arr(2);
			bs_val_arr = flux->get_b_dot_Bs_vec(p, B, sp_id);
			value_type bs_val = bs_val_arr[0];
			value_type dbs_dvp_ov_bs = bs_val_arr[1]/bs_val;

			tmp_w = tmp_w1*bs_val*fcm_val;
			for(int i = 0; i < dof; i++) tmp_w_arr[i] = tmp_w*(*basis)(q, i); 

			for(int j3 = qd_to_hg_bin_tot_index[b][0]; j3 < qd_to_hg_bin_tot_index[b][1] + 1; j3++)
			{
				value_type t1_tot = (v_para - (col3_int_vp_tot_min + j3*col3_int_dvp))*col3_int_inv_dvp;
				Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_int_inv_dvp);

				for(int j4 = qd_to_hg_bin_tot_index[b][2]; j4 < qd_to_hg_bin_tot_index[b][3] + 1; j4++)
				{
					value_type t2_tot = (u - (col3_int_u_tot_min + j4*col3_int_du))*col3_int_inv_du;
					Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_int_inv_du);

					int l = j4*(col3_int_vp_n + 3) + j3;

					for(int i = 0; i < dof; i++)
					{
						vector<value_type> tmp_val_out(4);
						tmp_val_out[0] = dbs_dvp_ov_bs*spline_vp_tot_vec(2)*spline_u_tot_vec(0);
						tmp_val_out[1] = -2.0*spline_vp_tot_vec(1)*spline_u_tot_vec(0);
						tmp_val_out[2] = (1.0 + v_para*dbs_dvp_ov_bs)*spline_vp_tot_vec(2)*spline_u_tot_vec(0);
						tmp_val_out[2] += u*dbs_dvp_ov_bs*spline_vp_tot_vec(1)*spline_u_tot_vec(1);
						tmp_val_out[2] += (spline_vp_tot_vec(0)*spline_u_tot_vec(1)/u + spline_vp_tot_vec(0)*spline_u_tot_vec(2))/B_val;
						tmp_val_out[3] = -2.0*v_para*spline_vp_tot_vec(1)*spline_u_tot_vec(0);
						tmp_val_out[3] += -2.0*u*spline_vp_tot_vec(0)*spline_u_tot_vec(1);

						for(int tmp_j = 0; tmp_j < 4; tmp_j++)
						{
							int li_index = l*dof + i + tmp_j*col3_int_hg_tot_n*dof;
							val_arr3[li_index] += tmp_w_arr[i]*tmp_val_out[tmp_j];
						}
					}
				}
			}
		}
	}
}

void Integration_Col::col3_hg_b_to_hg_a_vol_mat_cal(const ElementX &ti, const ElementV &tj, const int &sp_a_id, const int &sp_b_id, vector<value_type> &val_arr1, vector<value_type> &val_arr2)
{
	value_type w, wa, wb, v_para, u;
	Point4 p, q;
	Point2 pa, pb, qa, qb;

	int sp_a_int_vp_n = vp_n_arr[sp_a_id];
	int sp_a_int_u_n = u_n_arr[sp_a_id];
	int sp_a_hg_tot_n = (sp_a_int_vp_n + 3)*(sp_a_int_u_n + 3);
	int sp_b_int_vp_n = vp_n_arr[sp_b_id];
	int sp_b_int_u_n = u_n_arr[sp_b_id];
	int sp_b_hg_tot_n = (sp_b_int_vp_n + 3)*(sp_b_int_u_n + 3);

	value_type sp_a_int_vp_min = vp_min_arr[sp_a_id];
	value_type sp_a_int_u_min = u_min_arr[sp_a_id];
	value_type sp_a_int_dvp = dvp_arr[sp_a_id];
	value_type sp_a_int_du = du_arr[sp_a_id];
	value_type sp_a_int_inv_dvp = inv_dvp_arr[sp_a_id];
	value_type sp_a_int_inv_du = inv_du_arr[sp_a_id];
	value_type sp_a_tot_vp_min = sp_a_int_vp_min - sp_a_int_dvp;
	value_type sp_a_tot_u_min = sp_a_int_u_min - sp_a_int_du;


	value_type sp_b_int_vp_min = vp_min_arr[sp_b_id];
	value_type sp_b_int_u_min = u_min_arr[sp_b_id];
	value_type sp_b_int_dvp = dvp_arr[sp_b_id];
	value_type sp_b_int_du = du_arr[sp_b_id];
	value_type sp_b_int_inv_dvp = inv_dvp_arr[sp_b_id];
	value_type sp_b_int_inv_du = inv_du_arr[sp_b_id];
	value_type sp_b_tot_vp_min = sp_b_int_vp_min - sp_b_int_dvp;
	value_type sp_b_tot_u_min = sp_b_int_u_min - sp_b_int_du;

	vector<vector<int>> qd_to_hg_a_bin_inner_index, qd_to_hg_a_bin_tot_index;
	vector<vector<int>> qd_to_hg_b_bin_inner_index, qd_to_hg_b_bin_tot_index;
	for(int b = 0; b < quadrature->size(1); b++)
	{
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

		v_para = pb[0];
		u = pb[1];

		int loc_vp_inner_index = min(sp_a_int_vp_n - 1, max(0, int((v_para - sp_a_int_vp_min)*sp_a_int_inv_dvp)));
		int loc_u_inner_index = min(sp_a_int_u_n - 1, max(0, int((u - sp_a_int_u_min)*sp_a_int_inv_du)));
		int loc_vp_tot_index = loc_vp_inner_index + 1;
		int loc_u_tot_index = loc_u_inner_index + 1;

		vector<int> tmp_qd_to_hg_inner_index(4), tmp_qd_to_hg_tot_index(4);
		tmp_qd_to_hg_inner_index[0] = max(0, loc_vp_inner_index - 1);
		tmp_qd_to_hg_inner_index[1] = min(sp_a_int_vp_n, loc_vp_inner_index + 2);
		tmp_qd_to_hg_inner_index[2] = max(0, loc_u_inner_index - 1);
		tmp_qd_to_hg_inner_index[3] = min(sp_a_int_u_n, loc_u_inner_index + 2);

		tmp_qd_to_hg_tot_index[0] = max(0, loc_vp_tot_index - 1);
		tmp_qd_to_hg_tot_index[1] = min(sp_a_int_vp_n + 2, loc_vp_tot_index + 2);
		tmp_qd_to_hg_tot_index[2] = max(0, loc_u_tot_index - 1);
		tmp_qd_to_hg_tot_index[3] = min(sp_a_int_u_n + 2, loc_u_tot_index + 2);

		qd_to_hg_a_bin_inner_index.push_back(tmp_qd_to_hg_inner_index);
		qd_to_hg_a_bin_tot_index.push_back(tmp_qd_to_hg_tot_index);

		loc_vp_inner_index = min(sp_b_int_vp_n - 1, max(0, int((v_para - sp_b_int_vp_min)*sp_b_int_inv_dvp)));
		loc_u_inner_index = min(sp_b_int_u_n - 1, max(0, int((u - sp_b_int_u_min)*sp_b_int_inv_du)));

		loc_vp_tot_index = loc_vp_inner_index + 1;
		loc_u_tot_index = loc_u_inner_index + 1;

		tmp_qd_to_hg_inner_index[0] = max(0, loc_vp_inner_index - 1);
		tmp_qd_to_hg_inner_index[1] = min(sp_b_int_vp_n, loc_vp_inner_index + 2);
		tmp_qd_to_hg_inner_index[2] = max(0, loc_u_inner_index - 1);
		tmp_qd_to_hg_inner_index[3] = min(sp_b_int_u_n, loc_u_inner_index + 2);

		tmp_qd_to_hg_tot_index[0] = max(0, loc_vp_tot_index - 1);
		tmp_qd_to_hg_tot_index[1] = min(sp_b_int_vp_n + 2, loc_vp_tot_index + 2);
		tmp_qd_to_hg_tot_index[2] = max(0, loc_u_tot_index - 1);
		tmp_qd_to_hg_tot_index[3] = min(sp_b_int_u_n + 2, loc_u_tot_index + 2);

		qd_to_hg_b_bin_inner_index.push_back(tmp_qd_to_hg_inner_index);
		qd_to_hg_b_bin_tot_index.push_back(tmp_qd_to_hg_tot_index);
	}

	val_arr1.assign(val_arr1.size(), 0.0);
	val_arr2.assign(val_arr2.size(), 0.0);
	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));

		value_type R = p[0], Z = p[1];

		for(int b = 0; b < quadrature->size(1); b++)
		{
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

			w = wa * wb;
			p = pa & pb;
			q = qa & qb;

			v_para = p[2];
			u = p[3];

			value_type tmp_w1 = 2.0 * M_PI * w * R * u;
			value_type tmp_w2 = tmp_w1 * B_val;


			for(int j1 = qd_to_hg_a_bin_inner_index[b][0]; j1 < qd_to_hg_a_bin_inner_index[b][1] + 1; j1++)
			{
				value_type t1_inner = (v_para - (sp_a_int_vp_min + j1*sp_a_int_dvp))*sp_a_int_inv_dvp;
				Vector spline_vp_inner_vec = Spline::cal_Cs_Vec(t1_inner, sp_a_int_inv_dvp);

				for(int j2 = qd_to_hg_a_bin_inner_index[b][2]; j2 < qd_to_hg_a_bin_inner_index[b][3] + 1; j2++)
				{
					value_type t2_inner = (u - (sp_a_int_u_min + j2*sp_a_int_du))*sp_a_int_inv_du;
					Vector spline_u_inner_vec = Spline::cal_Cs_Vec(t2_inner, sp_a_int_inv_du);

					int m = j2*(sp_a_int_vp_n + 1) + j1;

					for(int j3 = qd_to_hg_b_bin_tot_index[b][0]; j3 < qd_to_hg_b_bin_tot_index[b][1] + 1; j3++)
					{
						value_type t1_tot = (v_para - (sp_b_tot_vp_min + j3*sp_b_int_dvp))*sp_b_int_inv_dvp;
						Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, sp_b_int_inv_dvp);

						for(int j4 = qd_to_hg_b_bin_tot_index[b][2]; j4 < qd_to_hg_b_bin_tot_index[b][3] + 1; j4++)
						{
							value_type t2_tot = (u - (sp_b_tot_u_min + j4*sp_b_int_du))*sp_b_int_inv_du;
							Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, sp_b_int_inv_du);

							int l = j4*(sp_b_int_vp_n + 3) + j3;

							int ml_index = m*sp_b_hg_tot_n + l;

							val_arr1[ml_index] += tmp_w2*spline_vp_inner_vec(0)*spline_u_inner_vec(0)*spline_vp_tot_vec(0)*spline_u_tot_vec(0);

						}
					}
				}
			}

			for(int j1 = qd_to_hg_b_bin_inner_index[b][0]; j1 < qd_to_hg_b_bin_inner_index[b][1] + 1; j1++)
			{
				value_type t1_inner = (v_para - (sp_b_int_vp_min + j1*sp_b_int_dvp))*sp_b_int_inv_dvp;
				Vector spline_vp_inner_vec = Spline::cal_Cs_Vec(t1_inner, sp_b_int_inv_dvp);

				for(int j2 = qd_to_hg_b_bin_inner_index[b][2]; j2 < qd_to_hg_b_bin_inner_index[b][3] + 1; j2++)
				{
					value_type t2_inner = (u - (sp_b_int_u_min + j2*sp_b_int_du))*sp_b_int_inv_du;
					Vector spline_u_inner_vec = Spline::cal_Cs_Vec(t2_inner, sp_b_int_inv_du);

					int m = j2*(sp_b_int_vp_n + 1) + j1;

					for(int j3 = qd_to_hg_a_bin_tot_index[b][0]; j3 < qd_to_hg_a_bin_tot_index[b][1] + 1; j3++)
					{
						value_type t1_tot = (v_para - (sp_a_tot_vp_min + j3*sp_a_int_dvp))*sp_a_int_inv_dvp;
						Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, sp_a_int_inv_dvp);

						for(int j4 = qd_to_hg_a_bin_tot_index[b][2]; j4 < qd_to_hg_a_bin_tot_index[b][3] + 1; j4++)
						{
							value_type t2_tot = (u - (sp_a_tot_u_min + j4*sp_a_int_du))*sp_a_int_inv_du;
							Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, sp_a_int_inv_du);

							int l = j4*(sp_a_int_vp_n + 3) + j3;

							int ml_index = m*sp_a_hg_tot_n + l;

							val_arr2[ml_index] += tmp_w2*spline_vp_inner_vec(0)*spline_u_inner_vec(0)*spline_vp_tot_vec(0)*spline_u_tot_vec(0);

						}
					}
				}
			}
		}
	}
}

void Integration_Col::col3_small_hg_b_to_hg_a_vol_mat_cal(const ElementX &ti, const ElementV &tj, const int &sp_a_id, const int &sp_b_id, vector<value_type> &val_arr1)
{
	value_type w, wa, wb, v_para, u;
	Point4 p, q;
	Point2 pa, pb, qa, qb;

	int sp_a_int_vp_n = vp_n_arr[sp_a_id];
	int sp_a_int_u_n = u_n_arr[sp_a_id];
	int sp_a_hg_tot_n = (sp_a_int_vp_n + 3)*(sp_a_int_u_n + 3);
	int sp_b_int_vp_n = vp_n_arr[sp_b_id];
	int sp_b_int_u_n = u_n_arr[sp_b_id];
	int sp_b_hg_tot_n = (sp_b_int_vp_n + 3)*(sp_b_int_u_n + 3);

	value_type sp_a_int_vp_min = vp_min_arr[sp_a_id];
	value_type sp_a_int_u_min = u_min_arr[sp_a_id];
	value_type sp_a_int_dvp = dvp_arr[sp_a_id];
	value_type sp_a_int_du = du_arr[sp_a_id];
	value_type sp_a_int_inv_dvp = inv_dvp_arr[sp_a_id];
	value_type sp_a_int_inv_du = inv_du_arr[sp_a_id];
	value_type sp_a_tot_vp_min = sp_a_int_vp_min - sp_a_int_dvp;
	value_type sp_a_tot_u_min = sp_a_int_u_min - sp_a_int_du;


	value_type sp_b_int_vp_min = vp_min_arr[sp_b_id];
	value_type sp_b_int_u_min = u_min_arr[sp_b_id];
	value_type sp_b_int_dvp = dvp_arr[sp_b_id];
	value_type sp_b_int_du = du_arr[sp_b_id];
	value_type sp_b_int_vp_max = sp_b_int_vp_min + sp_b_int_dvp*sp_b_int_vp_n;
	value_type sp_b_int_u_max = sp_b_int_u_min + sp_b_int_du*sp_b_int_u_n;


	value_type sp_b_int_inv_dvp = inv_dvp_arr[sp_b_id];
	value_type sp_b_int_inv_du = inv_du_arr[sp_b_id];
	value_type sp_b_tot_vp_min = sp_b_int_vp_min - sp_b_int_dvp;
	value_type sp_b_tot_u_min = sp_b_int_u_min - sp_b_int_du;

	vector<vector<int>> qd_to_hg_a_bin_inner_index, qd_to_hg_a_bin_tot_index;
	vector<vector<int>> qd_to_hg_b_bin_inner_index, qd_to_hg_b_bin_tot_index;
	for(int b = 0; b < quadrature->size(1); b++)
	{
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

		v_para = pb[0];
		u = pb[1];

		int loc_vp_inner_index = min(sp_a_int_vp_n - 1, max(0, int((v_para - sp_a_int_vp_min)*sp_a_int_inv_dvp)));
		int loc_u_inner_index = min(sp_a_int_u_n - 1, max(0, int((u - sp_a_int_u_min)*sp_a_int_inv_du)));
		int loc_vp_tot_index = loc_vp_inner_index + 1;
		int loc_u_tot_index = loc_u_inner_index + 1;

		vector<int> tmp_qd_to_hg_inner_index(4), tmp_qd_to_hg_tot_index(4);
		tmp_qd_to_hg_inner_index[0] = max(0, loc_vp_inner_index - 1);
		tmp_qd_to_hg_inner_index[1] = min(sp_a_int_vp_n, loc_vp_inner_index + 2);
		tmp_qd_to_hg_inner_index[2] = max(0, loc_u_inner_index - 1);
		tmp_qd_to_hg_inner_index[3] = min(sp_a_int_u_n, loc_u_inner_index + 2);

		tmp_qd_to_hg_tot_index[0] = max(0, loc_vp_tot_index - 1);
		tmp_qd_to_hg_tot_index[1] = min(sp_a_int_vp_n + 2, loc_vp_tot_index + 2);
		tmp_qd_to_hg_tot_index[2] = max(0, loc_u_tot_index - 1);
		tmp_qd_to_hg_tot_index[3] = min(sp_a_int_u_n + 2, loc_u_tot_index + 2);

		qd_to_hg_a_bin_inner_index.push_back(tmp_qd_to_hg_inner_index);
		qd_to_hg_a_bin_tot_index.push_back(tmp_qd_to_hg_tot_index);

		loc_vp_inner_index = min(sp_b_int_vp_n - 1, max(0, int((v_para - sp_b_int_vp_min)*sp_b_int_inv_dvp)));
		loc_u_inner_index = min(sp_b_int_u_n - 1, max(0, int((u - sp_b_int_u_min)*sp_b_int_inv_du)));

		loc_vp_tot_index = loc_vp_inner_index + 1;
		loc_u_tot_index = loc_u_inner_index + 1;

		tmp_qd_to_hg_inner_index[0] = max(0, loc_vp_inner_index - 1);
		tmp_qd_to_hg_inner_index[1] = min(sp_b_int_vp_n, loc_vp_inner_index + 2);
		tmp_qd_to_hg_inner_index[2] = max(0, loc_u_inner_index - 1);
		tmp_qd_to_hg_inner_index[3] = min(sp_b_int_u_n, loc_u_inner_index + 2);

		tmp_qd_to_hg_tot_index[0] = max(0, loc_vp_tot_index - 1);
		tmp_qd_to_hg_tot_index[1] = min(sp_b_int_vp_n + 2, loc_vp_tot_index + 2);
		tmp_qd_to_hg_tot_index[2] = max(0, loc_u_tot_index - 1);
		tmp_qd_to_hg_tot_index[3] = min(sp_b_int_u_n + 2, loc_u_tot_index + 2);

		qd_to_hg_b_bin_inner_index.push_back(tmp_qd_to_hg_inner_index);
		qd_to_hg_b_bin_tot_index.push_back(tmp_qd_to_hg_tot_index);
	}

	val_arr1.assign(val_arr1.size(), 0.0);
	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));

		value_type R = p[0], Z = p[1];

		for(int b = 0; b < quadrature->size(1); b++)
		{
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

			w = wa * wb;
			p = pa & pb;
			q = qa & qb;

			v_para = p[2];
			u = p[3];

			value_type tmp_w1 = 2.0 * M_PI * w * R * u;
			value_type tmp_w2 = tmp_w1 * B_val;

			for(int j1 = qd_to_hg_a_bin_inner_index[b][0]; j1 < qd_to_hg_a_bin_inner_index[b][1] + 1; j1++)
			{
				value_type t1_inner = (v_para - (sp_a_int_vp_min + j1*sp_a_int_dvp))*sp_a_int_inv_dvp;
				Vector spline_vp_inner_vec = Spline::cal_Cs_Vec(t1_inner, sp_a_int_inv_dvp);

				for(int j2 = qd_to_hg_a_bin_inner_index[b][2]; j2 < qd_to_hg_a_bin_inner_index[b][3] + 1; j2++)
				{
					value_type t2_inner = (u - (sp_a_int_u_min + j2*sp_a_int_du))*sp_a_int_inv_du;
					Vector spline_u_inner_vec = Spline::cal_Cs_Vec(t2_inner, sp_a_int_inv_du);

					int m = j2*(sp_a_int_vp_n + 1) + j1;

					if(v_para > sp_b_int_vp_min && v_para < sp_b_int_vp_max && u > sp_b_int_u_min && u < sp_b_int_u_max)
					{
						for(int j3 = qd_to_hg_b_bin_tot_index[b][0]; j3 < qd_to_hg_b_bin_tot_index[b][1] + 1; j3++)
						{
							value_type t1_tot = (v_para - (sp_b_tot_vp_min + j3*sp_b_int_dvp))*sp_b_int_inv_dvp;
							Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, sp_b_int_inv_dvp);

							for(int j4 = qd_to_hg_b_bin_tot_index[b][2]; j4 < qd_to_hg_b_bin_tot_index[b][3] + 1; j4++)
							{
								value_type t2_tot = (u - (sp_b_tot_u_min + j4*sp_b_int_du))*sp_b_int_inv_du;
								Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, sp_b_int_inv_du);

								int l = j4*(sp_b_int_vp_n + 3) + j3;

								int ml_index = m*sp_b_hg_tot_n + l;

								val_arr1[ml_index] += tmp_w2*spline_vp_inner_vec(0)*spline_u_inner_vec(0)*spline_vp_tot_vec(0)*spline_u_tot_vec(0);

							}
						}
					}
				}
			}
		}
	}
}


void Integration_Col::col3_slow_moment_to_fast_hg_cal(const ElementX &ti, const ElementV &tj, const int &sp_a_id, const int &sp_b_id, const value_type &U_b, const int &_num_of_bc_a, Vector &val_arr1)
{
	value_type w, wa, wb, v_para, u;
	Point4 p, q;
	Point2 pa, pb, qa, qb;

	int sp_a_int_vp_n = vp_n_arr[sp_a_id];
	int sp_a_int_u_n = u_n_arr[sp_a_id];
	int sp_a_hg_tot_n = (sp_a_int_vp_n + 3)*(sp_a_int_u_n + 3);
	int sp_b_int_vp_n = vp_n_arr[sp_b_id];
	int sp_b_int_u_n = u_n_arr[sp_b_id];
	int sp_b_hg_tot_n = (sp_b_int_vp_n + 3)*(sp_b_int_u_n + 3);

	value_type sp_a_int_vp_min = vp_min_arr[sp_a_id];
	value_type sp_a_int_u_min = u_min_arr[sp_a_id];
	value_type sp_a_int_dvp = dvp_arr[sp_a_id];
	value_type sp_a_int_du = du_arr[sp_a_id];
	value_type sp_a_int_inv_dvp = inv_dvp_arr[sp_a_id];
	value_type sp_a_int_inv_du = inv_du_arr[sp_a_id];
	value_type sp_a_tot_vp_min = sp_a_int_vp_min - sp_a_int_dvp;
	value_type sp_a_tot_u_min = sp_a_int_u_min - sp_a_int_du;


	value_type sp_b_int_vp_min = vp_min_arr[sp_b_id];
	value_type sp_b_int_u_min = u_min_arr[sp_b_id];
	value_type sp_b_int_dvp = dvp_arr[sp_b_id];
	value_type sp_b_int_du = du_arr[sp_b_id];
	value_type sp_b_int_vp_max = sp_b_int_vp_min + sp_b_int_dvp*sp_b_int_vp_n;
	value_type sp_b_int_u_max = sp_b_int_u_min + sp_b_int_du*sp_b_int_u_n;


	value_type sp_b_int_inv_dvp = inv_dvp_arr[sp_b_id];
	value_type sp_b_int_inv_du = inv_du_arr[sp_b_id];
	value_type sp_b_tot_vp_min = sp_b_int_vp_min - sp_b_int_dvp;
	value_type sp_b_tot_u_min = sp_b_int_u_min - sp_b_int_du;

	vector<vector<int>> qd_to_hg_a_bin_inner_index, qd_to_hg_a_bin_tot_index;
	vector<vector<int>> qd_to_hg_b_bin_inner_index, qd_to_hg_b_bin_tot_index;
	for(int b = 0; b < quadrature->size(1); b++)
	{
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

		v_para = pb[0];
		u = pb[1];

		int loc_vp_inner_index = min(sp_a_int_vp_n - 1, max(0, int((v_para - sp_a_int_vp_min)*sp_a_int_inv_dvp)));
		int loc_u_inner_index = min(sp_a_int_u_n - 1, max(0, int((u - sp_a_int_u_min)*sp_a_int_inv_du)));
		int loc_vp_tot_index = loc_vp_inner_index + 1;
		int loc_u_tot_index = loc_u_inner_index + 1;

		vector<int> tmp_qd_to_hg_inner_index(4), tmp_qd_to_hg_tot_index(4);
		tmp_qd_to_hg_inner_index[0] = max(0, loc_vp_inner_index - 1);
		tmp_qd_to_hg_inner_index[1] = min(sp_a_int_vp_n, loc_vp_inner_index + 2);
		tmp_qd_to_hg_inner_index[2] = max(0, loc_u_inner_index - 1);
		tmp_qd_to_hg_inner_index[3] = min(sp_a_int_u_n, loc_u_inner_index + 2);

		tmp_qd_to_hg_tot_index[0] = max(0, loc_vp_tot_index - 1);
		tmp_qd_to_hg_tot_index[1] = min(sp_a_int_vp_n + 2, loc_vp_tot_index + 2);
		tmp_qd_to_hg_tot_index[2] = max(0, loc_u_tot_index - 1);
		tmp_qd_to_hg_tot_index[3] = min(sp_a_int_u_n + 2, loc_u_tot_index + 2);

		qd_to_hg_a_bin_inner_index.push_back(tmp_qd_to_hg_inner_index);
		qd_to_hg_a_bin_tot_index.push_back(tmp_qd_to_hg_tot_index);

		loc_vp_inner_index = min(sp_b_int_vp_n - 1, max(0, int((v_para - sp_b_int_vp_min)*sp_b_int_inv_dvp)));
		loc_u_inner_index = min(sp_b_int_u_n - 1, max(0, int((u - sp_b_int_u_min)*sp_b_int_inv_du)));

		loc_vp_tot_index = loc_vp_inner_index + 1;
		loc_u_tot_index = loc_u_inner_index + 1;

		tmp_qd_to_hg_inner_index[0] = max(0, loc_vp_inner_index - 1);
		tmp_qd_to_hg_inner_index[1] = min(sp_b_int_vp_n, loc_vp_inner_index + 2);
		tmp_qd_to_hg_inner_index[2] = max(0, loc_u_inner_index - 1);
		tmp_qd_to_hg_inner_index[3] = min(sp_b_int_u_n, loc_u_inner_index + 2);

		tmp_qd_to_hg_tot_index[0] = max(0, loc_vp_tot_index - 1);
		tmp_qd_to_hg_tot_index[1] = min(sp_b_int_vp_n + 2, loc_vp_tot_index + 2);
		tmp_qd_to_hg_tot_index[2] = max(0, loc_u_tot_index - 1);
		tmp_qd_to_hg_tot_index[3] = min(sp_b_int_u_n + 2, loc_u_tot_index + 2);

		qd_to_hg_b_bin_inner_index.push_back(tmp_qd_to_hg_inner_index);
		qd_to_hg_b_bin_tot_index.push_back(tmp_qd_to_hg_tot_index);
	}

	//val_arr1.assign(val_arr1.size(), 0.0);
	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));

		value_type R = p[0], Z = p[1];

		for(int b = 0; b < quadrature->size(1); b++)
		{
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

			w = wa * wb;
			p = pa & pb;
			q = qa & qb;

			v_para = p[2];
			u = p[3];

			value_type tmp_w1 = 2.0 * M_PI * w * R * u;
			value_type tmp_w2 = tmp_w1 * B_val;

			for(int j1 = qd_to_hg_a_bin_inner_index[b][0]; j1 < qd_to_hg_a_bin_inner_index[b][1] + 1; j1++)
			{
				value_type t1_inner = (v_para - (sp_a_int_vp_min + j1*sp_a_int_dvp))*sp_a_int_inv_dvp;
				Vector spline_vp_inner_vec = Spline::cal_Cs_Vec(t1_inner, sp_a_int_inv_dvp);

				for(int j2 = qd_to_hg_a_bin_inner_index[b][2]; j2 < qd_to_hg_a_bin_inner_index[b][3] + 1; j2++)
				{
					value_type t2_inner = (u - (sp_a_int_u_min + j2*sp_a_int_du))*sp_a_int_inv_du;
					Vector spline_u_inner_vec = Spline::cal_Cs_Vec(t2_inner, sp_a_int_inv_du);

					int m = j2*(sp_a_int_vp_n + 1) + j1;

					if(v_para > sp_b_int_vp_min && v_para < sp_b_int_vp_max && u > sp_b_int_u_min && u < sp_b_int_u_max)
					{
					}
					else
					{
						value_type vp_mod = v_para - U_b;
						value_type vp_mod_sq = vp_mod*vp_mod;
						value_type u_sq = u*u*B_val;
						value_type v_sq = vp_mod_sq + u_sq;
						value_type v_5 = v_sq*v_sq*sqrt(v_sq);

						value_type loc_hb_from_moment = (vp_mod_sq - 0.5*u_sq)/(4.0*M_PI*v_5);
						val_arr1[2*_num_of_bc_a + m] += tmp_w2*spline_vp_inner_vec(0)*spline_u_inner_vec(0)*loc_hb_from_moment;

				
					}
				}
			}
		}
	}
}

void Integration_Col::col3_slow_moment_to_fast_hg_adj_n_cal(const ElementX &ti, const ElementV &tj, const int &sp_a_id, const int &sp_b_id, const value_type &U_b, const int &_num_of_bc_a, Vector &val_arr1, int &n_inout)
{
	value_type w, wa, v_para, u;
	Point4 p, q;
	Point2 pa, qa;

	int sp_a_int_vp_n = vp_n_arr[sp_a_id];
	int sp_a_int_u_n = u_n_arr[sp_a_id];
	int sp_a_hg_tot_n = (sp_a_int_vp_n + 3)*(sp_a_int_u_n + 3);
	int sp_b_int_vp_n = vp_n_arr[sp_b_id];
	int sp_b_int_u_n = u_n_arr[sp_b_id];
	int sp_b_hg_tot_n = (sp_b_int_vp_n + 3)*(sp_b_int_u_n + 3);

	value_type sp_a_int_vp_min = vp_min_arr[sp_a_id];
	value_type sp_a_int_u_min = u_min_arr[sp_a_id];
	value_type sp_a_int_dvp = dvp_arr[sp_a_id];
	value_type sp_a_int_du = du_arr[sp_a_id];
	value_type sp_a_int_inv_dvp = inv_dvp_arr[sp_a_id];
	value_type sp_a_int_inv_du = inv_du_arr[sp_a_id];
	value_type sp_a_tot_vp_min = sp_a_int_vp_min - sp_a_int_dvp;
	value_type sp_a_tot_u_min = sp_a_int_u_min - sp_a_int_du;


	value_type sp_b_int_vp_min = vp_min_arr[sp_b_id];
	value_type sp_b_int_u_min = u_min_arr[sp_b_id];
	value_type sp_b_int_dvp = dvp_arr[sp_b_id];
	value_type sp_b_int_du = du_arr[sp_b_id];
	value_type sp_b_int_vp_max = sp_b_int_vp_min + sp_b_int_dvp*sp_b_int_vp_n;
	value_type sp_b_int_u_max = sp_b_int_u_min + sp_b_int_du*sp_b_int_u_n;


	value_type sp_b_int_inv_dvp = inv_dvp_arr[sp_b_id];
	value_type sp_b_int_inv_du = inv_du_arr[sp_b_id];
	value_type sp_b_tot_vp_min = sp_b_int_vp_min - sp_b_int_dvp;
	value_type sp_b_tot_u_min = sp_b_int_u_min - sp_b_int_du;

	vector<vector<int>> qd_to_hg_a_bin_inner_index, qd_to_hg_a_bin_tot_index;
	vector<vector<int>> qd_to_hg_b_bin_inner_index, qd_to_hg_b_bin_tot_index;

	int qd_sizeV = n_inout;

	int u_index_loc = int(tj)/col3_int_org_f_vp_n;
	int vp_index_loc = int(tj) - u_index_loc*col3_int_org_f_vp_n;

	value_type wb = col3_int_du*col3_int_dvp/(qd_sizeV*qd_sizeV);
	value_type del_qb = 2.0/qd_sizeV;
	value_type qb1, qb2, pb1, pb2;

	for(int b1 = 0; b1 < qd_sizeV; b1++)
	{
		qb1 = -1.0 + (b1 + 0.5)*del_qb;
		pb1 = col3_int_vp_min + (vp_index_loc + (qb1 + 1.0)*0.5)*dvp;

		for(int b2 = 0; b2 < qd_sizeV; b2++)
		{
			int b = b1*qd_sizeV + b2;
			qb2 = -1.0 + (b2 + 0.5)*del_qb;
			pb2 = col3_int_u_min + (u_index_loc + (qb2 + 1.0)*0.5)*du;

			Point2 pb = Point2({pb1, pb2}), qb = Point2({qb1, qb2});

			v_para = pb[0];
			u = pb[1];

			int loc_vp_inner_index = min(sp_a_int_vp_n - 1, max(0, int((v_para - sp_a_int_vp_min)*sp_a_int_inv_dvp)));
			int loc_u_inner_index = min(sp_a_int_u_n - 1, max(0, int((u - sp_a_int_u_min)*sp_a_int_inv_du)));
			int loc_vp_tot_index = loc_vp_inner_index + 1;
			int loc_u_tot_index = loc_u_inner_index + 1;

			vector<int> tmp_qd_to_hg_inner_index(4), tmp_qd_to_hg_tot_index(4);
			tmp_qd_to_hg_inner_index[0] = max(0, loc_vp_inner_index - 1);
			tmp_qd_to_hg_inner_index[1] = min(sp_a_int_vp_n, loc_vp_inner_index + 2);
			tmp_qd_to_hg_inner_index[2] = max(0, loc_u_inner_index - 1);
			tmp_qd_to_hg_inner_index[3] = min(sp_a_int_u_n, loc_u_inner_index + 2);

			tmp_qd_to_hg_tot_index[0] = max(0, loc_vp_tot_index - 1);
			tmp_qd_to_hg_tot_index[1] = min(sp_a_int_vp_n + 2, loc_vp_tot_index + 2);
			tmp_qd_to_hg_tot_index[2] = max(0, loc_u_tot_index - 1);
			tmp_qd_to_hg_tot_index[3] = min(sp_a_int_u_n + 2, loc_u_tot_index + 2);

			qd_to_hg_a_bin_inner_index.push_back(tmp_qd_to_hg_inner_index);
			qd_to_hg_a_bin_tot_index.push_back(tmp_qd_to_hg_tot_index);

			loc_vp_inner_index = min(sp_b_int_vp_n - 1, max(0, int((v_para - sp_b_int_vp_min)*sp_b_int_inv_dvp)));
			loc_u_inner_index = min(sp_b_int_u_n - 1, max(0, int((u - sp_b_int_u_min)*sp_b_int_inv_du)));

			loc_vp_tot_index = loc_vp_inner_index + 1;
			loc_u_tot_index = loc_u_inner_index + 1;

			tmp_qd_to_hg_inner_index[0] = max(0, loc_vp_inner_index - 1);
			tmp_qd_to_hg_inner_index[1] = min(sp_b_int_vp_n, loc_vp_inner_index + 2);
			tmp_qd_to_hg_inner_index[2] = max(0, loc_u_inner_index - 1);
			tmp_qd_to_hg_inner_index[3] = min(sp_b_int_u_n, loc_u_inner_index + 2);

			tmp_qd_to_hg_tot_index[0] = max(0, loc_vp_tot_index - 1);
			tmp_qd_to_hg_tot_index[1] = min(sp_b_int_vp_n + 2, loc_vp_tot_index + 2);
			tmp_qd_to_hg_tot_index[2] = max(0, loc_u_tot_index - 1);
			tmp_qd_to_hg_tot_index[3] = min(sp_b_int_u_n + 2, loc_u_tot_index + 2);

			qd_to_hg_b_bin_inner_index.push_back(tmp_qd_to_hg_inner_index);
			qd_to_hg_b_bin_tot_index.push_back(tmp_qd_to_hg_tot_index);
		}
	}

	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));

		value_type R = p[0], Z = p[1];

		for(int b1 = 0; b1 < qd_sizeV; b1++)
		{
			qb1 = -1.0 + (b1 + 0.5)*del_qb;
			pb1 = col3_int_vp_min + (vp_index_loc + (qb1 + 1.0)*0.5)*dvp;

			for(int b2 = 0; b2 < qd_sizeV; b2++)
			{
				int b = b1*qd_sizeV + b2;

				qb2 = -1.0 + (b2 + 0.5)*del_qb;
				pb2 = col3_int_u_min + (u_index_loc + (qb2 + 1.0)*0.5)*du;

				Point2 pb = Point2({pb1, pb2}), qb = Point2({qb1, qb2});

				w = wa * wb;
				p = pa & pb;
				q = qa & qb;

				v_para = p[2];
				u = p[3];

				value_type tmp_w1 = 2.0 * M_PI * w * R * u;
				value_type tmp_w2 = tmp_w1 * B_val;

				for(int j1 = qd_to_hg_a_bin_inner_index[b][0]; j1 < qd_to_hg_a_bin_inner_index[b][1] + 1; j1++)
				{
					value_type t1_inner = (v_para - (sp_a_int_vp_min + j1*sp_a_int_dvp))*sp_a_int_inv_dvp;
					Vector spline_vp_inner_vec = Spline::cal_Cs_Vec(t1_inner, sp_a_int_inv_dvp);

					for(int j2 = qd_to_hg_a_bin_inner_index[b][2]; j2 < qd_to_hg_a_bin_inner_index[b][3] + 1; j2++)
					{
						value_type t2_inner = (u - (sp_a_int_u_min + j2*sp_a_int_du))*sp_a_int_inv_du;
						Vector spline_u_inner_vec = Spline::cal_Cs_Vec(t2_inner, sp_a_int_inv_du);

						int m = j2*(sp_a_int_vp_n + 1) + j1;

						if(v_para > sp_b_int_vp_min && v_para < sp_b_int_vp_max && u > sp_b_int_u_min && u < sp_b_int_u_max)
						{
						}
						else
						{
							value_type vp_mod = v_para - U_b;
							value_type vp_mod_sq = vp_mod*vp_mod;
							value_type u_sq = u*u*B_val;
							value_type v_sq = vp_mod_sq + u_sq;
							value_type v_5 = v_sq*v_sq*sqrt(v_sq);

							value_type loc_hb_from_moment = (vp_mod_sq - 0.5*u_sq)/(4.0*M_PI*v_5);
							val_arr1[2*_num_of_bc_a + m] += tmp_w2*spline_vp_inner_vec(0)*spline_u_inner_vec(0)*loc_hb_from_moment;
						}
					}
				}
			}
		}
	}
}



void Integration_Col::col3_small_fb_to_big_h_a_source_mat_cal(const ElementX &ti, const ElementV &tj, const int &sp_a_id, const int &sp_b_id, vector<value_type> &val_arr)
{
	value_type w, wa, wb, v_para, u;
	Point4 p, q;
	Point2 pa, pb, qa, qb;

	int sp_a_int_vp_n = vp_n_arr[sp_a_id];
	int sp_a_int_u_n = u_n_arr[sp_a_id];

	value_type sp_a_int_vp_min = vp_min_arr[sp_a_id];
	value_type sp_a_int_u_min = u_min_arr[sp_a_id];
	value_type sp_a_int_dvp = dvp_arr[sp_a_id];
	value_type sp_a_int_du = du_arr[sp_a_id];
	value_type sp_a_int_inv_dvp = inv_dvp_arr[sp_a_id];
	value_type sp_a_int_inv_du = inv_du_arr[sp_a_id];
	value_type sp_a_tot_vp_min = sp_a_int_vp_min - sp_a_int_dvp;
	value_type sp_a_tot_u_min = sp_a_int_u_min - sp_a_int_du;

	vector<vector<int>> qd_to_hg_a_bin_inner_index;

	for(int b = 0; b < quadrature->size(1); b++)
	{
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

		v_para = pb[0];
		u = pb[1];

		int loc_vp_inner_index = min(sp_a_int_vp_n - 1, max(0, int((v_para - sp_a_int_vp_min)*sp_a_int_inv_dvp)));
		int loc_u_inner_index = min(sp_a_int_u_n - 1, max(0, int((u - sp_a_int_u_min)*sp_a_int_inv_du)));

		vector<int> tmp_qd_to_hg_inner_index(4);
		tmp_qd_to_hg_inner_index[0] = max(0, loc_vp_inner_index - 1);
		tmp_qd_to_hg_inner_index[1] = min(sp_a_int_vp_n, loc_vp_inner_index + 2);
		tmp_qd_to_hg_inner_index[2] = max(0, loc_u_inner_index - 1);
		tmp_qd_to_hg_inner_index[3] = min(sp_a_int_u_n, loc_u_inner_index + 2);

		qd_to_hg_a_bin_inner_index.push_back(tmp_qd_to_hg_inner_index);
	}

	val_arr.assign(val_arr.size(), 0.0);

	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));

		value_type R = p[0], Z = p[1];

		for(int b = 0; b < quadrature->size(1); b++)
		{
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

			w = wa * wb;
			p = pa & pb;
			q = qa & qb;

			v_para = p[2];
			u = p[3];

			auto fcm_val = flux->get_fcm_SE_tot(p, int(ti), sp_id, flag_fcm);

			value_type tmp_w1 = 2.0 * M_PI * w * R * u;
			value_type tmp_w2 = tmp_w1 * B_val;
			value_type tmp_w = tmp_w2* fcm_val;

			value_type tmp_w_arr[dof];

			for(int i = 0; i < dof; i++) tmp_w_arr[i] = tmp_w*(*basis)(q, i); 

			for(int j1 = qd_to_hg_a_bin_inner_index[b][0]; j1 < qd_to_hg_a_bin_inner_index[b][1] + 1; j1++)
			{
				value_type t1_inner = (v_para - (sp_a_int_vp_min + j1*sp_a_int_dvp))*sp_a_int_inv_dvp;
				Vector spline_vp_inner_vec = Spline::cal_Cs_Vec(t1_inner, sp_a_int_inv_dvp);

				for(int j2 = qd_to_hg_a_bin_inner_index[b][2]; j2 < qd_to_hg_a_bin_inner_index[b][3] + 1; j2++)
				{
					value_type t2_inner = (u - (sp_a_int_u_min + j2*sp_a_int_du))*sp_a_int_inv_du;
					Vector spline_u_inner_vec = Spline::cal_Cs_Vec(t2_inner, sp_a_int_inv_du);

					int m = j2*(sp_a_int_vp_n + 1) + j1;

					//f to h source part
					for(int i = 0; i < dof; i++)
					{
						int mi_index = m*dof + i;

						val_arr[mi_index] += tmp_w_arr[i]*spline_vp_inner_vec(0)*spline_u_inner_vec(0);
					}

				}
			}


		}
	}
}



void Integration_Col::col3_edge_bd_hg_stiffness(const ElementX &ti, const EdgeV &te, const int &sp_id, const int &i_op, vector<value_type> &val_arr)
{
	value_type w, wa, wb, R, Z, v_para, u;
	Point4 p, q;
	Point2 pa, pb, qa, qb;


	vector<vector<int>> qd_to_hg_bin_inner_index, qd_to_hg_bin_tot_index;
	for(int b = 0; b < quadrature->size(3); b++)
	{
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(te, b);

		v_para = pb[0];
		u = pb[1];

		int loc_vp_inner_index = min(col3_int_vp_n - 1, max(0, int((v_para - col3_int_vp_min)*col3_int_inv_dvp)));
		int loc_u_inner_index = min(col3_int_u_n - 1, max(0, int((u - col3_int_u_min)*col3_int_inv_du)));

		int loc_vp_tot_index = loc_vp_inner_index + 1;
		int loc_u_tot_index = loc_u_inner_index + 1;

		vector<int> tmp_qd_to_hg_inner_index(4), tmp_qd_to_hg_tot_index(4);
		tmp_qd_to_hg_inner_index[0] = max(0, loc_vp_inner_index - 1);
		tmp_qd_to_hg_inner_index[1] = min(col3_int_vp_n, loc_vp_inner_index + 2);
		tmp_qd_to_hg_inner_index[2] = max(0, loc_u_inner_index - 1);
		tmp_qd_to_hg_inner_index[3] = min(col3_int_u_n, loc_u_inner_index + 2);

		tmp_qd_to_hg_tot_index[0] = max(0, loc_vp_tot_index - 1);
		tmp_qd_to_hg_tot_index[1] = min(col3_int_vp_n + 2, loc_vp_tot_index + 2);
		tmp_qd_to_hg_tot_index[2] = max(0, loc_u_tot_index - 1);
		tmp_qd_to_hg_tot_index[3] = min(col3_int_u_n + 2, loc_u_tot_index + 2);

		qd_to_hg_bin_inner_index.push_back(tmp_qd_to_hg_inner_index);
		qd_to_hg_bin_tot_index.push_back(tmp_qd_to_hg_tot_index);
	}

	val_arr.assign(val_arr.size(), 0.0);
	for(int a = 0; a < quadrature->size(0); a++)
	{
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));

		R = p[0];
		Z = p[1];

		for(int b = 0; b < quadrature->size(3); b++)
		{
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(te, b);

			w = wa * wb;
			p = pa & pb;
			q = qa & qb;

			v_para = p[2];
			u = p[3];

			value_type tmp_w = - 2.0 * M_PI * w * R * u;
			if(i_op == 2) tmp_w *= B_val;
			if(i_op == 2 && v_para < (col3_int_vp_min + col3_int_vp_max)*0.5) tmp_w *= -1.0;

			for(int j1 = qd_to_hg_bin_inner_index[b][0]; j1 < qd_to_hg_bin_inner_index[b][1] + 1; j1++)
			{
				value_type t1_inner = (v_para - (col3_int_vp_min + j1*col3_int_dvp))*col3_int_inv_dvp;
				Vector spline_vp_inner_vec = Spline::cal_Cs_Vec(t1_inner, col3_int_inv_dvp);

				for(int j2 = qd_to_hg_bin_inner_index[b][2]; j2 < qd_to_hg_bin_inner_index[b][3] + 1; j2++)
				{
					value_type t2_inner = (u - (col3_int_u_min + j2*col3_int_du))*col3_int_inv_du;
					Vector spline_u_inner_vec = Spline::cal_Cs_Vec(t2_inner, col3_int_inv_du);

					int m = j2*(col3_int_vp_n + 1) + j1;

					for(int j3 = qd_to_hg_bin_tot_index[b][0]; j3 < qd_to_hg_bin_tot_index[b][1] + 1; j3++)
					{
						value_type t1_tot = (v_para - (col3_int_vp_tot_min + j3*col3_int_dvp))*col3_int_inv_dvp;
						Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_int_inv_dvp);

						for(int j4 = qd_to_hg_bin_tot_index[b][2]; j4 < qd_to_hg_bin_tot_index[b][3] + 1; j4++)
						{
							value_type t2_tot = (u - (col3_int_u_tot_min + j4*col3_int_du))*col3_int_inv_du;
							Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_int_inv_du);

							int l = j4*(col3_int_vp_n + 3) + j3;

							int ml_index = m*col3_int_hg_tot_n + l;

							if(i_op == 2) 
							{
								val_arr[ml_index] += tmp_w*spline_vp_inner_vec(0)*spline_u_inner_vec(0)*spline_vp_tot_vec(1)*spline_u_tot_vec(0);
							}
							else
							{
								val_arr[ml_index] += tmp_w*spline_vp_inner_vec(0)*spline_u_inner_vec(0)*spline_vp_tot_vec(0)*spline_u_tot_vec(1);
							}
						}
					}
				}
			}
		}
	}
}





Vector Integration_Col::col_fM_coeff_vol_cal(const ElementX &ti, const ElementV &tj, const int &sp_id, const Vector &avged_Q)
{
	value_type den = avged_Q[0], U_b_para = avged_Q[1], vth_b_sq = avged_Q[2];

	Vector tmp_sum_w_evol = Vector::Zero(dof);
	value_type ccf_fM0 = 0.5*M_1_PI*sqrt(0.5*M_1_PI)*den/(vth_b_sq*sqrt(vth_b_sq));
	int ix = int(ti), iv = int(tj);
	int qd_sizeX = quadrature->size(0), qd_sizeV = quadrature->size(1);
	vector<value_type> wb_arr(qd_sizeV);
	vector<Point2> pb_arr(qd_sizeV), qb_arr(qd_sizeV);
	dof_x = basis->get_dofx();
	dof_v = basis->get_dofv();

	for(int b = 0; b < qd_sizeV; b++)
		tie(wb_arr[b], pb_arr[b], qb_arr[b]) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

	for(int a = 0; a < qd_sizeX; a++)
	{
		value_type wa;
		Point2 pa, qa;
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		vector<value_type> basisX(dof_x);
		for(int i = 0; i < dof_x; i++) basisX[i] = basisX_qd[a][i];

		Point4 p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));
		value_type R = p[0];

		for(int b = 0; b < qd_sizeV; b++)
		{
			p = pa & pb_arr[b];
			value_type w = wa * wb_arr[b];
			value_type v_para = p[2], u = p[3];

			auto bs_val = flux->get_b_dot_Bs(p, sp_id);

			Vector basisXV_loc = Vector::Zero(dof);


			for(int i = 0; i < dof_x; i++)
			{
				value_type tmp_basisX = basisX[i];
				for(int j = 0; j < dof_v; j++)
				{
					basisXV_loc[j*dof_x + i] = tmp_basisX*basisV_qd[b][j];
				}
			}
				value_type f0_val = ccf_fM0*exp(-0.5*((v_para-U_b_para)*(v_para-U_b_para) + u*u*B_val)/vth_b_sq);

			value_type tmp_w1 = w * bs_val * R * u;//*fcm_val;
			value_type tmp_wfinit = tmp_w1 * f0_val;

			tmp_sum_w_evol += tmp_wfinit*basisXV_loc;
		}
	}

	return tmp_sum_w_evol;
}





void Integration_Col::col3_hg_to_Qx_cal(const ElementX &ti, Vector &hg_to_Qx_Vec)
{
	int qd_sizeX = quadrature->size(0);

	for(int a = 0; a < qd_sizeX; a++)
	{
		value_type wa;
		Point2 pa, qa;
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 

		Point4 p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));
		value_type B_inv_val = 1.0/B_val;

		value_type R = pa[0], Z = pa[1];

		value_type fM_wX = flux->get_fcm_SE_X(p, int(ti), sp_id, flag_fcm);

		value_type tmp_wa1 = wa * R;
		value_type tmp_waf0 = tmp_wa1 * fM_wX;

		Vector Gj(3);

		Gj[0] = B_val;
		Gj[1] = 1.0;
		Gj[2] = B_inv_val;

		vector<value_type> basisX(dof_x);
		basis->col_spatial_basis_vec(qa, basisX);

		for(int lxp = 0; lxp < dof_x; lxp++)
		{
			value_type zeta_lxp_val = basisX[lxp];
			for(int lx = 0; lx < dof_x; lx++)
			{
				value_type zeta_lx_val = basisX[lx];

				for(int j = 0; j < 3; j++)
				{
					int tmp_index = (j*dof_x + lxp)*dof_x + lx;
					hg_to_Qx_Vec(tmp_index) += tmp_waf0*zeta_lxp_val*zeta_lx_val*Gj[j];
				}
			}
		}
	}
}



void Integration_Col::col3_hg_to_Qv_vol_cal(const ElementX &ti, const ElementV &tj, Vector &hg_to_Q_Vec)
{
	int qd_sizeV = quadrature->size(1);

	for(int b = 0; b < qd_sizeV; b++)
	{
		value_type wb;
		Point2 pb, qb;
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

		vector<value_type> basisV(dof_v*6);
		basis->col_velocity_basis_vec(qb, dvp, du, basisV);

		value_type v_para = pb[0], u = pb[1];
		Point4 p = Point2({0, 0}) & pb;
		value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_id, flag_fcm);

		value_type tmp_wb1 = wb * u;
		value_type tmp_wbf0 = tmp_wb1 * fM_wV;

		int loc_vp_tot_index = int((v_para - col3_int_vp_tot_min)*col3_int_inv_dvp);
		int loc_u_tot_index = int((u - col3_int_u_tot_min)*col3_int_inv_du);

		int vp_tot_min_index = max(0, loc_vp_tot_index - 1);
		int vp_tot_max_index = min(col3_int_vp_n + 2, loc_vp_tot_index + 2);
		int u_tot_min_index = max(0, loc_u_tot_index - 1);
		int u_tot_max_index = min(col3_int_u_n + 2, loc_u_tot_index + 2);

		for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
		{
			value_type t1_tot = (v_para - (col3_int_vp_tot_min + j3*col3_int_dvp))*col3_int_inv_dvp;
			Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_int_inv_dvp);

			for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
			{
				value_type t2_tot = (u - (col3_int_u_tot_min + j4*col3_int_du))*col3_int_inv_du;
				Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_int_inv_du);

				Matrix tmp_val = Matrix::Zero(3,3);

				for(int i1 = 0; i1 < 3; i1++)
				{
					for(int i2 = 0; i2 < 3; i2++)
					{
						tmp_val(i1, i2) = spline_vp_tot_vec(i1)*spline_u_tot_vec(i2);
					}
				}

				int l = j4*(col3_int_vp_n + 3) + j3;

				for(int lv = 0; lv < dof_v; lv++)
				{
					Vector Gj(5);

					Gj[0] = tmp_val(1,0)*basisV[lv*6 + 1];
					Gj[1] = tmp_val(0,1)*basisV[lv*6 + 2];
					Gj[2] = tmp_val(2,0)*basisV[lv*6 + 3];
					Gj[3] = 2.0*tmp_val(1,1)*basisV[lv*6 + 4];
					Gj[4] = tmp_val(0,1)*basisV[lv*6 + 2]/u/u + tmp_val(0,2)*basisV[lv*6 + 5];

					for(int j = 0; j < 5; j++)
					{
						for(int lvp = 0; lvp < dof_v; lvp++)
						{
							value_type zeta_lvp_val = basisV[lvp*6 + 0];

							int tmp_index = ((j*dof_v + lvp)*dof_v + lv)*col3_int_hg_tot_n + l;
							hg_to_Q_Vec(tmp_index) += tmp_wbf0*zeta_lvp_val*Gj[j];
						}
					}
				}
			}
		}
	}
}

void Integration_Col::col3_hMgM_to_Qv_vol_cal(const int &sp_b_id, const ElementX &ti, const ElementV &tj, const Vector &tmp_avged_Q, const int &col3_hMbgMb_with_UaTa_op, Vector &hMgM_to_Q_Vec)
{
	int sp_a_id = sp_id;
	int qd_sizeV = quadrature->size(1);

	vector<Vector> nUT_b_arr, ccf_fMb_arr;
	nUT_b_arr.resize(2, Vector::Zero(3));
	ccf_fMb_arr.resize(2, Vector::Zero(2));

	int k1_max = 1;
	if(col3_hMbgMb_with_UaTa_op == 1) k1_max = 2;

	int pitch_angle_ei_col_loc_op = 0;
	if(sp_a_id == 0 && sp_b_id != 0 && ei_pitch_angle_col_op == 1)
		pitch_angle_ei_col_loc_op = 1;

	value_type min_x_val_loc = 0.0, const_inv_x3_loc = 0.0;

	for(int k1 = 0; k1 < k1_max; k1++)
	{
		for(int k2 = 0; k2 < 3; k2++)
		{
			nUT_b_arr[k1][k2] = tmp_avged_Q[k1*3 + k2];

		}

		if(pitch_angle_ei_col_loc_op == 1)
			nUT_b_arr[k1][1] = 0.0;

		ccf_fMb_arr[k1][0] = M_1_PI/sqrt(32.0)*nUT_b_arr[k1][0]/sqrt(nUT_b_arr[k1][2]);
		ccf_fMb_arr[k1][1] = ccf_fMb_arr[k1][0]*nUT_b_arr[k1][2];

		if(pitch_angle_ei_col_loc_op == 1 && ei_pitch_no_v_in_nu_ei_op == 1)
		{
			min_x_val_loc = ei_pitch_no_v_in_nu_ei_v_e_min*sqrt(tmp_avged_Q[7 + 2]/tmp_avged_Q[k1*3 + 2]*0.5);
			value_type const_inv_x = 1.0/min_x_val_loc;
			const_inv_x3_loc = const_inv_x*const_inv_x*const_inv_x;
		}
	}

	value_type B_val = tmp_avged_Q[6];
	value_type pitch_angle_g_fac_mod;

	for(int b = 0; b < qd_sizeV; b++)
	{
		value_type wb;
		Point2 pb, qb;
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);

		vector<value_type> basisV(dof_v*6);
		basis->col_velocity_basis_vec(qb, dvp, du, basisV);

		value_type v_para = pb[0], u = pb[1];
		Point4 p = Point2({0, 0}) & pb;
		value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_a_id, flag_fcm);

		value_type tmp_wb1 = wb * u;
		value_type tmp_wbf0 = tmp_wb1 * fM_wV;

		Matrix h_tot_mat = Matrix::Zero(3,3);
		Matrix g_tot_mat = Matrix::Zero(3,3);

		for(int k1 = 0; k1 < k1_max; k1++)
		{
			value_type vp_mod = v_para - nUT_b_arr[k1][1];
			value_type v_mod_sq = vp_mod*vp_mod + u*u*B_val;
			value_type v_mod = sqrt(v_mod_sq);
			value_type inv_v_mod_sq = 1.0/v_mod_sq;

			value_type vth_b_sq = nUT_b_arr[k1][2];
			value_type x_val = v_mod/sqrt(2.0*vth_b_sq), x_val2, inv_x, inv_x2, inv_x3;
			x_val2 = x_val*x_val;
			inv_x = 1.0/x_val;
			inv_x2 = inv_x*inv_x;
			inv_x3 = inv_x2*inv_x;

			value_type inv_2_vth_b_sq_x_val = 1.0/(2.0*vth_b_sq*x_val);

			value_type erf_x = erf(x_val), derf_dx = 1.128379167095513*exp(-x_val2);
			value_type dx_du = u*B_val*inv_2_vth_b_sq_x_val;
			value_type dx_dvp = vp_mod*inv_2_vth_b_sq_x_val;

			value_type ccf_h0 = ccf_fMb_arr[k1][0];
			value_type ccf_g0 = ccf_fMb_arr[k1][1];

			value_type dh0_dx = ccf_h0*(-inv_x2*erf_x + inv_x*derf_dx);
			value_type dg0_dx = ccf_g0*((1.0-0.5*inv_x2)*erf_x + 0.5*inv_x*derf_dx);
			value_type d2g0_dx2 = ccf_g0*(inv_x3*erf_x - inv_x2*derf_dx);

			value_type d2x_dvp2 = (1.0 - (vp_mod*vp_mod)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;
			value_type d2x_dudvp = -B_val*vp_mod*u*inv_v_mod_sq*inv_2_vth_b_sq_x_val;
			value_type d2x_du2 = B_val*(1.0 - B_val*(u*u)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;

			h_tot_mat(1,0) = dh0_dx*dx_dvp;
			h_tot_mat(0,1)= dh0_dx*dx_du;

			g_tot_mat(0,1)= dg0_dx*dx_du;
			g_tot_mat(2,0) = dg0_dx*d2x_dvp2 + dx_dvp*dx_dvp*d2g0_dx2;
			g_tot_mat(1,1)= dg0_dx*d2x_dudvp + dx_dvp*dx_du*d2g0_dx2;
			g_tot_mat(0,2) = dg0_dx*d2x_du2 + dx_du*dx_du*d2g0_dx2;

			if(pitch_angle_ei_col_loc_op == 1)
			{ 
				h_tot_mat.setZero();
				g_tot_mat.setZero();

				pitch_angle_g_fac_mod = d2g0_dx2*x_val - dg0_dx;
				value_type mod_inv_x3 = inv_x3;

				if(ei_pitch_no_v_in_nu_ei_op == 1)
				{
					if(x_val < min_x_val_loc)
						mod_inv_x3 = const_inv_x3_loc;
				}

				pitch_angle_g_fac_mod *= mod_inv_x3/(4.0*vth_b_sq*vth_b_sq)*B_val;
			}

			for(int lv = 0; lv < dof_v; lv++)
			{
				Vector Gj(5);

				Gj[0] = h_tot_mat(1,0)*basisV[lv*6 + 1];
				Gj[1] = h_tot_mat(0,1)*basisV[lv*6 + 2];

				if(pitch_angle_ei_col_loc_op == 0)
				{
					Gj[2] = g_tot_mat(2,0)*basisV[lv*6 + 3];
					Gj[3] = 2.0*g_tot_mat(1,1)*basisV[lv*6 + 4];
					Gj[4] = g_tot_mat(0,1)*basisV[lv*6 + 2]/u/u + g_tot_mat(0,2)*basisV[lv*6 + 5];
				}
				else
				{
					Gj[2] = (-u*u)*basisV[lv*6 + 3];
					Gj[3] = (2.0*vp_mod)*basisV[lv*6 + 1] + (u)*basisV[lv*6 + 2] + (2*u*vp_mod)*basisV[lv*6 + 4];
					Gj[4] = (-vp_mod*vp_mod/u)*basisV[lv*6 + 2] + (-vp_mod*vp_mod)*basisV[lv*6 + 5];

					Gj *= pitch_angle_g_fac_mod;
				}

				for(int j = 0; j < 5; j++)
				{
					for(int lvp = 0; lvp < dof_v; lvp++)
					{
						value_type zeta_lvp_val = basisV[lvp*6 + 0];

						int tmp_index = ((k1*5 + j)*dof_v + lvp)*dof_v + lv;
						hMgM_to_Q_Vec(tmp_index) += tmp_wbf0*zeta_lvp_val*Gj[j];
					}
				}
			}
		}
	}
}

void Integration_Col::col3_hMgM_to_Qv_vol_adj_n_cal(const int &sp_b_id, const ElementX &ti, const ElementV &tj, const Vector &tmp_avged_Q, const int &col3_hMbgMb_with_UaTa_op, Vector &hMgM_to_Q_Vec, int &n_inout, value_type err_tor)
{

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);



	int sp_a_id = sp_id;
	int qd_sizeV = n_inout;

	vector<Vector> nUT_b_arr, ccf_fMb_arr;
	nUT_b_arr.resize(2, Vector::Zero(3));
	ccf_fMb_arr.resize(2, Vector::Zero(2));

	int k1_max = 1;
	if(col3_hMbgMb_with_UaTa_op == 1) k1_max = 2;

	int pitch_angle_ei_col_loc_op = 0;
	if(sp_a_id == 0 && sp_b_id != 0 && ei_pitch_angle_col_op == 1)
		pitch_angle_ei_col_loc_op = 1;

	value_type min_x_val_loc = 0.0, const_inv_x3_loc = 0.0;

	for(int k1 = 0; k1 < k1_max; k1++)
	{
		for(int k2 = 0; k2 < 3; k2++)
		{
			nUT_b_arr[k1][k2] = tmp_avged_Q[k1*3 + k2];

		}

		if(pitch_angle_ei_col_loc_op == 1)
			nUT_b_arr[k1][1] = 0.0;

		ccf_fMb_arr[k1][0] = M_1_PI/sqrt(32.0)*nUT_b_arr[k1][0]/sqrt(nUT_b_arr[k1][2]);
		ccf_fMb_arr[k1][1] = ccf_fMb_arr[k1][0]*nUT_b_arr[k1][2];

		if(pitch_angle_ei_col_loc_op == 1 && ei_pitch_no_v_in_nu_ei_op == 1)
		{
			min_x_val_loc = ei_pitch_no_v_in_nu_ei_v_e_min*sqrt(tmp_avged_Q[7 + 2]/tmp_avged_Q[k1*3 + 2]*0.5);
			value_type const_inv_x = 1.0/min_x_val_loc;
			const_inv_x3_loc = const_inv_x*const_inv_x*const_inv_x;
		}
	}

	int total_vec_size = 2*5*dof_v*dof_v;

	Vector tmp_out_0 = Vector::Zero(total_vec_size), tmp_out_1 = Vector::Zero(total_vec_size);

	value_type B_val = tmp_avged_Q[6];
	value_type pitch_angle_g_fac_mod;

	int u_index_loc = int(tj)/col3_int_org_f_vp_n;
	int vp_index_loc = int(tj) - u_index_loc*col3_int_org_f_vp_n;

	int iter = 0;
	value_type tor = 1e-3;
	tor = 5e-3;
	tor = err_tor;
	//tor = 1e-2;
	value_type err = 1.0;

	while(iter < 8 && err> tor)
	{

		tmp_out_1 = Vector::Zero(total_vec_size);

		value_type wb = col3_int_du*col3_int_dvp/(qd_sizeV*qd_sizeV);
		value_type del_qb = 2.0/qd_sizeV;
		value_type qb1, qb2, pb1, pb2;

		for(int b1 = 0; b1 < qd_sizeV; b1++)
		{
			qb1 = -1.0 + (b1 + 0.5)*del_qb;
			pb1 = col3_int_vp_min + (vp_index_loc + (qb1 + 1.0)*0.5)*dvp;

			for(int b2 = 0; b2 < qd_sizeV; b2++)
			{
				qb2 = -1.0 + (b2 + 0.5)*del_qb;
				pb2 = col3_int_u_min + (u_index_loc + (qb2 + 1.0)*0.5)*du;

				Point2 pb = Point2({pb1, pb2}), qb = Point2({qb1, qb2});

				vector<value_type> basisV(dof_v*6);
				basis->col_velocity_basis_vec(qb, dvp, du, basisV);

				value_type v_para = pb[0], u = pb[1];
				Point4 p = Point2({0, 0}) & pb;
				value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_a_id, flag_fcm);

				value_type tmp_wb1 = wb * u;
				value_type tmp_wbf0 = tmp_wb1 * fM_wV;

				Matrix h_tot_mat = Matrix::Zero(3,3);
				Matrix g_tot_mat = Matrix::Zero(3,3);

				for(int k1 = 0; k1 < k1_max; k1++)
				{
					value_type vp_mod = v_para - nUT_b_arr[k1][1];
					value_type v_mod_sq = vp_mod*vp_mod + u*u*B_val;
					value_type v_mod = sqrt(v_mod_sq);
					value_type inv_v_mod_sq = 1.0/v_mod_sq;

					value_type vth_b_sq = nUT_b_arr[k1][2];
					value_type x_val = v_mod/sqrt(2.0*vth_b_sq), x_val2, inv_x, inv_x2, inv_x3;
					x_val2 = x_val*x_val;
					inv_x = 1.0/x_val;
					inv_x2 = inv_x*inv_x;
					inv_x3 = inv_x2*inv_x;

					value_type inv_2_vth_b_sq_x_val = 1.0/(2.0*vth_b_sq*x_val);

					value_type erf_x = erf(x_val), derf_dx = 1.128379167095513*exp(-x_val2);
					value_type dx_du = u*B_val*inv_2_vth_b_sq_x_val;
					value_type dx_dvp = vp_mod*inv_2_vth_b_sq_x_val;

					value_type ccf_h0 = ccf_fMb_arr[k1][0];
					value_type ccf_g0 = ccf_fMb_arr[k1][1];

					value_type dh0_dx = ccf_h0*(-inv_x2*erf_x + inv_x*derf_dx);
					value_type dg0_dx = ccf_g0*((1.0-0.5*inv_x2)*erf_x + 0.5*inv_x*derf_dx);
					value_type d2g0_dx2 = ccf_g0*(inv_x3*erf_x - inv_x2*derf_dx);

					value_type d2x_dvp2 = (1.0 - (vp_mod*vp_mod)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;
					value_type d2x_dudvp = -B_val*vp_mod*u*inv_v_mod_sq*inv_2_vth_b_sq_x_val;
					value_type d2x_du2 = B_val*(1.0 - B_val*(u*u)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;

					h_tot_mat(1,0) = dh0_dx*dx_dvp;
					h_tot_mat(0,1)= dh0_dx*dx_du;

					g_tot_mat(0,1)= dg0_dx*dx_du;
					g_tot_mat(2,0) = dg0_dx*d2x_dvp2 + dx_dvp*dx_dvp*d2g0_dx2;
					g_tot_mat(1,1)= dg0_dx*d2x_dudvp + dx_dvp*dx_du*d2g0_dx2;
					g_tot_mat(0,2) = dg0_dx*d2x_du2 + dx_du*dx_du*d2g0_dx2;

					if(pitch_angle_ei_col_loc_op == 1)
					{ 
						h_tot_mat.setZero();
						g_tot_mat.setZero();

						pitch_angle_g_fac_mod = d2g0_dx2*x_val - dg0_dx;
						value_type mod_inv_x3 = inv_x3;

						if(ei_pitch_no_v_in_nu_ei_op == 1)
						{
							if(x_val < min_x_val_loc)
								mod_inv_x3 = const_inv_x3_loc;
						}

						pitch_angle_g_fac_mod *= mod_inv_x3/(4.0*vth_b_sq*vth_b_sq)*B_val;
					}

					for(int lv = 0; lv < dof_v; lv++)
					{
						Vector Gj(5);

						Gj[0] = h_tot_mat(1,0)*basisV[lv*6 + 1];
						Gj[1] = h_tot_mat(0,1)*basisV[lv*6 + 2];

						if(pitch_angle_ei_col_loc_op == 0)
						{
							Gj[2] = g_tot_mat(2,0)*basisV[lv*6 + 3];
							Gj[3] = 2.0*g_tot_mat(1,1)*basisV[lv*6 + 4];
							Gj[4] = g_tot_mat(0,1)*basisV[lv*6 + 2]/u/u + g_tot_mat(0,2)*basisV[lv*6 + 5];
						}
						else
						{
							Gj[2] = (-u*u)*basisV[lv*6 + 3];
							Gj[3] = (2.0*vp_mod)*basisV[lv*6 + 1] + (u)*basisV[lv*6 + 2] + (2*u*vp_mod)*basisV[lv*6 + 4];
							Gj[4] = (-vp_mod*vp_mod/u)*basisV[lv*6 + 2] + (-vp_mod*vp_mod)*basisV[lv*6 + 5];

							Gj *= pitch_angle_g_fac_mod;
						}

						for(int j = 0; j < 5; j++)
						{
							for(int lvp = 0; lvp < dof_v; lvp++)
							{
								value_type zeta_lvp_val = basisV[lvp*6 + 0];

								int tmp_index = ((k1*5 + j)*dof_v + lvp)*dof_v + lv;
								tmp_out_1(tmp_index) += tmp_wbf0*zeta_lvp_val*Gj[j];
							}
						}
					}
				}
			}
		}

		Vector del_tmp_out = tmp_out_1 - tmp_out_0;
		value_type err_sq = del_tmp_out.dot(del_tmp_out)/(tmp_out_1.dot(tmp_out_1));

		err = sqrt(err_sq);
		tmp_out_0 = tmp_out_1;
		qd_sizeV *= 2;

		iter += 1;
	}
	hMgM_to_Q_Vec = tmp_out_1;

	n_inout = qd_sizeV/4;
}


void Integration_Col::col3_hg_S_LF_Vmat_cal(const ElementX &ti,  const EdgeV &te, const ElementV &kl, const ElementV &kr, const int &i_op, Vector &tmp_h_LF_Vmat_Vec, Vector &tmp_h_S_Vmat_Vec, Vector &tmp_g_S_Vmat_Vec)
{
	const auto nb = quadrature->size(3);
	//const auto [kl, kr] = mesh->get_neighborhood(te);
	int dim_fhat_org_v_basis = basis->fhat_org_vbasis_dim_out();
	int dim_fhat_valid_v_basis = basis->fhat_valid_vbasis_dim_out();

	vector<Matrix> hg_basis_mat_arr[4];
	for(int i = 0; i < 4; i++)
		hg_basis_mat_arr[i].resize(4, Matrix::Zero(3,3));

	for(int b = 0; b < nb; b++)
	{
		value_type wb;
		Point2 pb, qb;
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(te, b); //qb[0] is a dummy variable

		Point4 p = Point2({1.0, 0}) & pb;
		value_type vp = pb[0];
		value_type u = pb[1];

		//fM_wV part
		//value_type fM_wV = 1.0;
		value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_id, flag_fcm);

		vector<value_type> dlnfMV_dv_arr(3);
		dlnfMV_dv_arr = flux->get_fcm_SE_vec(p, int(ti), sp_id, flag_fcm);

		dlnfMV_dv_arr[1] = dlnfMV_dv_arr[1]/dlnfMV_dv_arr[0];
		dlnfMV_dv_arr[2] = dlnfMV_dv_arr[2]/dlnfMV_dv_arr[0];

		value_type tmp_wb1 = wb * u;
		value_type tmp_wbf0 = tmp_wb1 * fM_wV;

		//hg_basis part
		int loc_vp_tot_index = int((vp - col3_int_vp_tot_min)*col3_int_inv_dvp);
		int loc_u_tot_index = int((u - col3_int_u_tot_min)*col3_int_inv_du);

		int vp_tot_min_index = max(0, loc_vp_tot_index - 1);
		int vp_tot_max_index = min(col3_int_vp_n + 2, loc_vp_tot_index + 2);
		int u_tot_min_index = max(0, loc_u_tot_index - 1);
		int u_tot_max_index = min(col3_int_u_n + 2, loc_u_tot_index + 2);

		for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
		{
			int vp_i_index = j3 - vp_tot_min_index;
			value_type t1_tot = (vp - (col3_int_vp_tot_min + j3*col3_int_dvp))*col3_int_inv_dvp;
			Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_int_inv_dvp);

			for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
			{
				int u_i_index = j4 - u_tot_min_index;
				value_type t2_tot = (u - (col3_int_u_tot_min + j4*col3_int_du))*col3_int_inv_du;
				Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_int_inv_du);

				for(int i1 = 0; i1 < 3; i1++)
				{
					for(int i2 = 0; i2 < 3; i2++)
					{
						hg_basis_mat_arr[vp_i_index][u_i_index](i1, i2) = spline_vp_tot_vec(i1)*spline_u_tot_vec(i2);
					}
				}
			}
		}

		//org basisV part of v elements adjacent to the edge te
		int a= 0;
		auto ql = quadrature->get_traced_point(ti, te, kl, a, b);
		auto qr = quadrature->get_traced_point(ti, te, kr, a, b);

		Point2 qb_l = Point2({ql[2], ql[3]});
		Point2 qb_r = Point2({qr[2], qr[3]});

		vector<value_type> basisV_lr[2];
		basisV_lr[0].resize(dof_v*6);
		basisV_lr[1].resize(dof_v*6);
		basis->col_velocity_basis_vec(qb_l, dvp, du, basisV_lr[0]);
		basis->col_velocity_basis_vec(qb_r, dvp, du, basisV_lr[1]);

		vector<value_type> zeta_lv_l_val_arr, zeta_lv_r_val_arr;

		for(int lv = 0; lv < dof_v; lv++)
		{
			zeta_lv_l_val_arr.push_back(basisV_lr[0][lv*6 + 0]);
			zeta_lv_r_val_arr.push_back(basisV_lr[1][lv*6 + 0]);
		}

		//h & LF part for mixed flux of upwind and LF
		for(int j = 0; j < 4; j++)
		{
			//int index_j_tmp = j*(4*dof_v*dof_v);
			for(int lvp = 0; lvp < dof_v; lvp++)
			{
				for(int lv = 0; lv < dof_v; lv++)
				{
					value_type basis_v_tot;
					switch(j)
					{
						case 0:
							basis_v_tot = zeta_lv_l_val_arr[lvp]*zeta_lv_l_val_arr[lv];
							break;

						case 1:
							basis_v_tot = -zeta_lv_l_val_arr[lvp]*zeta_lv_r_val_arr[lv];
							break;

						case 2:
							basis_v_tot = zeta_lv_r_val_arr[lvp]*zeta_lv_l_val_arr[lv];
							break;

						case 3:
							basis_v_tot = -zeta_lv_r_val_arr[lvp]*zeta_lv_r_val_arr[lv];
							break;

						default:
							cerr << "Unknown j type : " << j << endl;
							exit(EXIT_FAILURE);
					}

					int tmp_index = (j*dof_v + lvp)*dof_v + lv;

					tmp_h_LF_Vmat_Vec(tmp_index) += tmp_wbf0*basis_v_tot;


					for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
					{
						int vp_i_index = j3 - vp_tot_min_index;

						for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
						{
							int u_i_index = j4 - u_tot_min_index;

							Matrix tmp_hg_basis_mat = hg_basis_mat_arr[vp_i_index][u_i_index];

							int l = j4*(col3_int_vp_n + 3) + j3;

							Vector Gj(1);

							if(i_op == 0) Gj[0] = tmp_hg_basis_mat(1,0);
							else if(i_op == 1) Gj[0] = tmp_hg_basis_mat(0,1);

							int tmp_index2 = tmp_index*col3_int_hg_tot_n + l;
							tmp_h_S_Vmat_Vec[tmp_index2] += tmp_wbf0*basis_v_tot*Gj[0];

						}
					}
				}
			}
		}

		Vector fhat_basis_Vec = basis->fhat_vbasis_Vec(ql, i_op, 0); //fhat at edge with original fhat index
		for(int row = 0; row < dim_fhat_valid_v_basis; row++) 
		{
			int fhat_org_basis_i = basis->fhat_valid_to_org_i(row, i_op);
			//fhat_loc_arr[0] : fhat
			//fhat_loc_arr[1] : dfhat/dvp
			//fhat_loc_arr[2] : dfhat/du
			Vector fhat_loc_arr = Vector::Zero(3);
			fhat_loc_arr[0] = fhat_basis_Vec(fhat_org_basis_i);
			fhat_loc_arr[1] = fhat_basis_Vec(dim_fhat_org_v_basis + fhat_org_basis_i);
			fhat_loc_arr[2] = fhat_basis_Vec(2*dim_fhat_org_v_basis + fhat_org_basis_i);

			for(int i_ele = 0; i_ele < 2; i_ele++)
			{
				value_type n_dot_v = 1.0;
				if(i_ele == 1) n_dot_v = -1.0;

				for(int lv = 0; lv < dof_v; lv++)
				{
					Vector basis_loc_arr = Vector::Zero(3);
					basis_loc_arr[0] = basisV_lr[i_ele][lv*6 + 0];
					basis_loc_arr[1] = basisV_lr[i_ele][lv*6 + 1];
					basis_loc_arr[2] = basisV_lr[i_ele][lv*6 + 2];


					for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
					{
						int vp_i_index = j3 - vp_tot_min_index;

						for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
						{
							int u_i_index = j4 - u_tot_min_index;

							Matrix tmp_hg_basis_mat = hg_basis_mat_arr[vp_i_index][u_i_index];

							int l = j4*(col3_int_vp_n + 3) + j3;

							Vector Gj(2);
							if(i_op == 0)
							{
								Gj[0] = tmp_hg_basis_mat(2,0); 
								Gj[1] = tmp_hg_basis_mat(1,1); 
							}
							else
							{
								Gj[0] = tmp_hg_basis_mat(1,1); 
								Gj[1] = tmp_hg_basis_mat(0,2); 
							}

							for(int j = 0; j < 2; j++)
							{
								value_type tmp_val = 0.0;
								if(j == 0)
								{
									tmp_val += -basis_loc_arr[0]*fhat_loc_arr[1];
									tmp_val += -basis_loc_arr[0]*fhat_loc_arr[0]*dlnfMV_dv_arr[1];
									tmp_val += basis_loc_arr[1]*fhat_loc_arr[0];
								}
								else
								{
									tmp_val += -basis_loc_arr[0]*fhat_loc_arr[2];
									tmp_val += -basis_loc_arr[0]*fhat_loc_arr[0]*dlnfMV_dv_arr[2];
									tmp_val += basis_loc_arr[2]*fhat_loc_arr[0];
								}

								Gj[j] *= n_dot_v*tmp_val;

								int tmp_index = ((j*dim_fhat_valid_v_basis + row)*2 + i_ele)*dof_v + lv;
								int tmp_index2 = tmp_index*col3_int_hg_tot_n + l;
								tmp_g_S_Vmat_Vec[tmp_index2] += tmp_wbf0*Gj[j];

							}
						}
					}
				}
			}
		}
	}
}

void Integration_Col::col3_hMgM_S_Vmat_cal(const int &sp_b_id, const ElementX &ti,  const EdgeV &te, const ElementV &kl, const ElementV &kr, const int &i_op, const Vector &tmp_avged_Q, const int &col3_hMbgMb_with_UaTa_op, Vector &tmp_hM_S_Vmat_Vec, Vector &tmp_gM_S_Vmat_Vec)
{
	int sp_a_id = sp_id;
	const auto nb = quadrature->size(3);

	int dim_fhat_org_v_basis = basis->fhat_org_vbasis_dim_out();
	int dim_fhat_valid_v_basis = basis->fhat_valid_vbasis_dim_out();

	vector<Vector> nUT_b_arr, ccf_fMb_arr;
	nUT_b_arr.resize(2, Vector::Zero(3));
	ccf_fMb_arr.resize(2, Vector::Zero(2));

	int k1_max = 1;
	if(col3_hMbgMb_with_UaTa_op == 1) k1_max = 2;

	int pitch_angle_ei_col_loc_op = 0;
	if(sp_a_id == 0 && sp_b_id != 0 && ei_pitch_angle_col_op == 1)
		pitch_angle_ei_col_loc_op = 1;

	value_type min_x_val_loc = 0.0, const_inv_x3_loc = 0.0;

	for(int k1 = 0; k1 < k1_max; k1++)
	{
		for(int k2 = 0; k2 < 3; k2++)
		{
			nUT_b_arr[k1][k2] = tmp_avged_Q[k1*3 + k2];
		}

		if(pitch_angle_ei_col_loc_op == 1)
			nUT_b_arr[k1][1] = 0.0;

		ccf_fMb_arr[k1][0] = M_1_PI/sqrt(32.0)*nUT_b_arr[k1][0]/sqrt(nUT_b_arr[k1][2]);
		ccf_fMb_arr[k1][1] = ccf_fMb_arr[k1][0]*nUT_b_arr[k1][2];

		if(pitch_angle_ei_col_loc_op == 1 && ei_pitch_no_v_in_nu_ei_op == 1)
		{
			min_x_val_loc = ei_pitch_no_v_in_nu_ei_v_e_min*sqrt(tmp_avged_Q[7 + 2]/tmp_avged_Q[k1*3 + 2]*0.5);
			value_type const_inv_x = 1.0/min_x_val_loc;
			const_inv_x3_loc = const_inv_x*const_inv_x*const_inv_x;
		}
	}

	value_type B_val = tmp_avged_Q[6];
	value_type pitch_angle_g_fac_mod;

	for(int b = 0; b < nb; b++)
	{
		value_type wb;
		Point2 pb, qb;
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(te, b); //qb[0] is a dummy variable

		Point4 p = Point2({1.0, 0}) & pb;
		value_type vp = pb[0];
		value_type u = pb[1];

		//fM_wV part
		//value_type fM_wV = 1.0;
		value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_id, flag_fcm);

		vector<value_type> dlnfMV_dv_arr(3);
		dlnfMV_dv_arr = flux->get_fcm_SE_vec(p, int(ti), sp_id, flag_fcm);

		dlnfMV_dv_arr[1] = dlnfMV_dv_arr[1]/dlnfMV_dv_arr[0];
		dlnfMV_dv_arr[2] = dlnfMV_dv_arr[2]/dlnfMV_dv_arr[0];

		value_type tmp_wb1 = wb * u;
		value_type tmp_wbf0 = tmp_wb1 * fM_wV;

		//org basisV part of v elements adjacent to the edge te
		int a= 0;
		auto ql = quadrature->get_traced_point(ti, te, kl, a, b);
		auto qr = quadrature->get_traced_point(ti, te, kr, a, b);

		Point2 qb_l = Point2({ql[2], ql[3]});
		Point2 qb_r = Point2({qr[2], qr[3]});

		vector<value_type> basisV_lr[2];
		basisV_lr[0].resize(dof_v*6);
		basisV_lr[1].resize(dof_v*6);
		basis->col_velocity_basis_vec(qb_l, dvp, du, basisV_lr[0]);
		basis->col_velocity_basis_vec(qb_r, dvp, du, basisV_lr[1]);

		vector<value_type> zeta_lv_l_val_arr, zeta_lv_r_val_arr;

		for(int lv = 0; lv < dof_v; lv++)
		{
			zeta_lv_l_val_arr.push_back(basisV_lr[0][lv*6 + 0]);
			zeta_lv_r_val_arr.push_back(basisV_lr[1][lv*6 + 0]);
		}

		Matrix h_tot_mat = Matrix::Zero(3,3);
		Matrix g_tot_mat = Matrix::Zero(3,3);

		for(int k1 = 0; k1 < k1_max; k1++)
		{
			value_type vp_mod = vp - nUT_b_arr[k1][1];
			value_type v_mod_sq = vp_mod*vp_mod + u*u*B_val;
			value_type v_mod = sqrt(v_mod_sq);
			value_type inv_v_mod_sq = 1.0/v_mod_sq;

			value_type vth_b_sq = nUT_b_arr[k1][2];
			value_type x_val = v_mod/sqrt(2.0*vth_b_sq), x_val2, inv_x, inv_x2, inv_x3;
			x_val2 = x_val*x_val;
			inv_x = 1.0/x_val;
			inv_x2 = inv_x*inv_x;
			inv_x3 = inv_x2*inv_x;

			value_type inv_2_vth_b_sq_x_val = 1.0/(2.0*vth_b_sq*x_val);

			value_type erf_x = erf(x_val), derf_dx = 1.128379167095513*exp(-x_val2);
			value_type dx_du = u*B_val*inv_2_vth_b_sq_x_val;
			value_type dx_dvp = vp_mod*inv_2_vth_b_sq_x_val;

			value_type ccf_h0 = ccf_fMb_arr[k1][0];
			value_type ccf_g0 = ccf_fMb_arr[k1][1];

			value_type dh0_dx = ccf_h0*(-inv_x2*erf_x + inv_x*derf_dx);
			value_type dg0_dx = ccf_g0*((1.0-0.5*inv_x2)*erf_x + 0.5*inv_x*derf_dx);
			value_type d2g0_dx2 = ccf_g0*(inv_x3*erf_x - inv_x2*derf_dx);

			value_type d2x_dvp2 = (1.0 - (vp_mod*vp_mod)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;
			value_type d2x_dudvp = -B_val*vp_mod*u*inv_v_mod_sq*inv_2_vth_b_sq_x_val;
			value_type d2x_du2 = B_val*(1.0 - B_val*(u*u)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;

			h_tot_mat(1,0) = dh0_dx*dx_dvp;
			h_tot_mat(0,1)= dh0_dx*dx_du;

			g_tot_mat(2,0) = dg0_dx*d2x_dvp2 + dx_dvp*dx_dvp*d2g0_dx2;
			g_tot_mat(1,1)= dg0_dx*d2x_dudvp + dx_dvp*dx_du*d2g0_dx2;
			g_tot_mat(0,2) = dg0_dx*d2x_du2 + dx_du*dx_du*d2g0_dx2;

			if(pitch_angle_ei_col_loc_op == 1)
			{ 
				h_tot_mat.setZero();
				g_tot_mat.setZero();

				//pitch_angle_g_fac_mod = 0.0;
				pitch_angle_g_fac_mod = d2g0_dx2*x_val - dg0_dx;
				value_type mod_inv_x3 = inv_x3;

				if(ei_pitch_no_v_in_nu_ei_op == 1)
				{
					if(x_val < min_x_val_loc)
						mod_inv_x3 = const_inv_x3_loc;
				}

				pitch_angle_g_fac_mod *= mod_inv_x3/(4.0*vth_b_sq*vth_b_sq)*B_val;

				g_tot_mat(2,0) = -u*u*pitch_angle_g_fac_mod;
				g_tot_mat(1,1)= vp_mod*u*pitch_angle_g_fac_mod;
				g_tot_mat(0,2) = -vp_mod*vp_mod*pitch_angle_g_fac_mod;
			}

			//h part for mixed flux of upwind and LF
			for(int j = 0; j < 4; j++)
			{
				//int index_j_tmp = j*(4*dof_v*dof_v);
				for(int lvp = 0; lvp < dof_v; lvp++)
				{
					for(int lv = 0; lv < dof_v; lv++)
					{
						value_type basis_v_tot;
						switch(j)
						{
							case 0:
								basis_v_tot = zeta_lv_l_val_arr[lvp]*zeta_lv_l_val_arr[lv];
								break;

							case 1:
								basis_v_tot = -zeta_lv_l_val_arr[lvp]*zeta_lv_r_val_arr[lv];
								break;

							case 2:
								basis_v_tot = zeta_lv_r_val_arr[lvp]*zeta_lv_l_val_arr[lv];
								break;

							case 3:
								basis_v_tot = -zeta_lv_r_val_arr[lvp]*zeta_lv_r_val_arr[lv];
								break;

							default:
								cerr << "Unknown j type : " << j << endl;
								exit(EXIT_FAILURE);
						}

						Vector Gj(1);

						if(i_op == 0) Gj[0] = h_tot_mat(1,0);
						else if(i_op == 1) Gj[0] = h_tot_mat(0,1);

						int tmp_index = ((k1*4 + j)*dof_v + lvp)*dof_v + lv;
						tmp_hM_S_Vmat_Vec[tmp_index] += tmp_wbf0*basis_v_tot*Gj[0];
					}
				}
			}



			Vector fhat_basis_Vec = basis->fhat_vbasis_Vec(ql, i_op, 0); //fhat at edge with original fhat index
			for(int row = 0; row < dim_fhat_valid_v_basis; row++) 
			{
				int fhat_org_basis_i = basis->fhat_valid_to_org_i(row, i_op);
				//fhat_loc_arr[0] : fhat
				//fhat_loc_arr[1] : dfhat/dvp
				//fhat_loc_arr[2] : dfhat/du
				Vector fhat_loc_arr = Vector::Zero(3);
				fhat_loc_arr[0] = fhat_basis_Vec(fhat_org_basis_i);
				fhat_loc_arr[1] = fhat_basis_Vec(dim_fhat_org_v_basis + fhat_org_basis_i);
				fhat_loc_arr[2] = fhat_basis_Vec(2*dim_fhat_org_v_basis + fhat_org_basis_i);

				for(int i_ele = 0; i_ele < 2; i_ele++)
				{
					value_type n_dot_v = 1.0;
					if(i_ele == 1) n_dot_v = -1.0;

					for(int lv = 0; lv < dof_v; lv++)
					{
						Vector basis_loc_arr = Vector::Zero(3);
						basis_loc_arr[0] = basisV_lr[i_ele][lv*6 + 0];
						basis_loc_arr[1] = basisV_lr[i_ele][lv*6 + 1];
						basis_loc_arr[2] = basisV_lr[i_ele][lv*6 + 2];

						Vector Gj(2);
						if(i_op == 0)
						{
							Gj[0] = g_tot_mat(2,0); 
							Gj[1] = g_tot_mat(1,1); 
						}
						else
						{
							Gj[0] = g_tot_mat(1,1); 
							Gj[1] = g_tot_mat(0,2); 
						}

						for(int j = 0; j < 2; j++)
						{
							value_type tmp_val = 0.0;
							if(j == 0)
							{
								tmp_val += -basis_loc_arr[0]*fhat_loc_arr[1];
								tmp_val += -basis_loc_arr[0]*fhat_loc_arr[0]*dlnfMV_dv_arr[1];
								tmp_val += basis_loc_arr[1]*fhat_loc_arr[0];
							}
							else
							{
								tmp_val += -basis_loc_arr[0]*fhat_loc_arr[2];
								tmp_val += -basis_loc_arr[0]*fhat_loc_arr[0]*dlnfMV_dv_arr[2];
								tmp_val += basis_loc_arr[2]*fhat_loc_arr[0];
							}

							Gj[j] *= n_dot_v*tmp_val;

							int tmp_index = (((k1*2 + j)*dim_fhat_valid_v_basis + row)*2 + i_ele)*dof_v + lv;
							tmp_gM_S_Vmat_Vec[tmp_index] += tmp_wbf0*Gj[j];
						}
					}
				}
			}
		}
	}
}

void Integration_Col::col3_hMgM_S_Vmat_adj_n_cal(const int &sp_b_id, const ElementX &ti,  const EdgeV &te, const ElementV &kl, const ElementV &kr, const int &i_op, const Vector &tmp_avged_Q, const int &col3_hMbgMb_with_UaTa_op, Vector &tmp_hM_S_Vmat_Vec, Vector &tmp_gM_S_Vmat_Vec, int &n_inout)
{
	int sp_a_id = sp_id;
	int nb = n_inout;

	int dim_fhat_org_v_basis = basis->fhat_org_vbasis_dim_out();
	int dim_fhat_valid_v_basis = basis->fhat_valid_vbasis_dim_out();

	vector<Vector> nUT_b_arr, ccf_fMb_arr;
	nUT_b_arr.resize(2, Vector::Zero(3));
	ccf_fMb_arr.resize(2, Vector::Zero(2));

	int k1_max = 1;
	if(col3_hMbgMb_with_UaTa_op == 1) k1_max = 2;

	int pitch_angle_ei_col_loc_op = 0;
	if(sp_a_id == 0 && sp_b_id != 0 && ei_pitch_angle_col_op == 1)
		pitch_angle_ei_col_loc_op = 1;

	value_type min_x_val_loc = 0.0, const_inv_x3_loc = 0.0;

	for(int k1 = 0; k1 < k1_max; k1++)
	{
		for(int k2 = 0; k2 < 3; k2++)
		{
			nUT_b_arr[k1][k2] = tmp_avged_Q[k1*3 + k2];
		}

		if(pitch_angle_ei_col_loc_op == 1)
			nUT_b_arr[k1][1] = 0.0;

		ccf_fMb_arr[k1][0] = M_1_PI/sqrt(32.0)*nUT_b_arr[k1][0]/sqrt(nUT_b_arr[k1][2]);
		ccf_fMb_arr[k1][1] = ccf_fMb_arr[k1][0]*nUT_b_arr[k1][2];

		if(pitch_angle_ei_col_loc_op == 1 && ei_pitch_no_v_in_nu_ei_op == 1)
		{
			min_x_val_loc = ei_pitch_no_v_in_nu_ei_v_e_min*sqrt(tmp_avged_Q[7 + 2]/tmp_avged_Q[k1*3 + 2]*0.5);
			value_type const_inv_x = 1.0/min_x_val_loc;
			const_inv_x3_loc = const_inv_x*const_inv_x*const_inv_x;
		}
	}

	value_type B_val = tmp_avged_Q[6];
	value_type pitch_angle_g_fac_mod;

	int u_index_kl_loc = int(kl)/col3_int_org_f_vp_n;
	int vp_index_kl_loc = int(kl) - u_index_kl_loc*col3_int_org_f_vp_n;

	value_type wb;

	if(i_op == 0)
		wb = col3_int_du/nb;
	else
		wb = col3_int_dvp/nb;

	value_type del_qb = 2.0/nb;
	value_type qb1 = 0.0, qb2, pb1, pb2;


	for(int b = 0; b < nb; b++)
	{
		qb2 = -1.0 + (b + 0.5)*del_qb;
		if (i_op == 0)
		{
			pb1 = col3_int_vp_min + (vp_index_kl_loc + 1)*dvp;
			pb2 = col3_int_u_min + (u_index_kl_loc + (qb2 + 1.0)*0.5)*du;
		}
		else
		{
			pb1 = col3_int_vp_min + (vp_index_kl_loc + (qb2 + 1.0)*0.5)*dvp;
			pb2 = col3_int_u_min + (u_index_kl_loc + 1)*du;
		}

		Point2 pb = Point2({pb1, pb2}), qb = Point2({qb1, qb2});

		Point4 p = Point2({1.0, 0}) & pb;
		value_type vp = pb[0];
		value_type u = pb[1];

		//fM_wV part
		value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_id, flag_fcm);

		vector<value_type> dlnfMV_dv_arr(3);
		dlnfMV_dv_arr = flux->get_fcm_SE_vec(p, int(ti), sp_id, flag_fcm);

		dlnfMV_dv_arr[1] = dlnfMV_dv_arr[1]/dlnfMV_dv_arr[0];
		dlnfMV_dv_arr[2] = dlnfMV_dv_arr[2]/dlnfMV_dv_arr[0];

		value_type tmp_wb1 = wb * u;
		value_type tmp_wbf0 = tmp_wb1 * fM_wV;

		//org basisV part of v elements adjacent to the edge te
		int a= 0;
		auto ql = quadrature->get_traced_point(ti, te, kl, a, b);
		auto qr = quadrature->get_traced_point(ti, te, kr, a, b);

		if (i_op == 0)
		{
			ql[2] = 1.0;
			qr[2] = -1.0;
			ql[3] = qr[3] = qb2;
		}
		else
		{
			ql[2] = qr[2] = qb2;
			ql[3] = 1.0;
			qr[3] = -1.0;
		}

		Point2 qb_l = Point2({ql[2], ql[3]});
		Point2 qb_r = Point2({qr[2], qr[3]});

		vector<value_type> basisV_lr[2];
		basisV_lr[0].resize(dof_v*6);
		basisV_lr[1].resize(dof_v*6);
		basis->col_velocity_basis_vec(qb_l, dvp, du, basisV_lr[0]);
		basis->col_velocity_basis_vec(qb_r, dvp, du, basisV_lr[1]);

		vector<value_type> zeta_lv_l_val_arr, zeta_lv_r_val_arr;

		for(int lv = 0; lv < dof_v; lv++)
		{
			zeta_lv_l_val_arr.push_back(basisV_lr[0][lv*6 + 0]);
			zeta_lv_r_val_arr.push_back(basisV_lr[1][lv*6 + 0]);
		}

		Matrix h_tot_mat = Matrix::Zero(3,3);
		Matrix g_tot_mat = Matrix::Zero(3,3);

		for(int k1 = 0; k1 < k1_max; k1++)
		{
			value_type vp_mod = vp - nUT_b_arr[k1][1];
			value_type v_mod_sq = vp_mod*vp_mod + u*u*B_val;
			value_type v_mod = sqrt(v_mod_sq);
			value_type inv_v_mod_sq = 1.0/v_mod_sq;

			value_type vth_b_sq = nUT_b_arr[k1][2];
			value_type x_val = v_mod/sqrt(2.0*vth_b_sq), x_val2, inv_x, inv_x2, inv_x3;
			x_val2 = x_val*x_val;
			inv_x = 1.0/x_val;
			inv_x2 = inv_x*inv_x;
			inv_x3 = inv_x2*inv_x;

			value_type inv_2_vth_b_sq_x_val = 1.0/(2.0*vth_b_sq*x_val);

			value_type erf_x = erf(x_val), derf_dx = 1.128379167095513*exp(-x_val2);
			value_type dx_du = u*B_val*inv_2_vth_b_sq_x_val;
			value_type dx_dvp = vp_mod*inv_2_vth_b_sq_x_val;

			value_type ccf_h0 = ccf_fMb_arr[k1][0];
			value_type ccf_g0 = ccf_fMb_arr[k1][1];

			value_type dh0_dx = ccf_h0*(-inv_x2*erf_x + inv_x*derf_dx);
			value_type dg0_dx = ccf_g0*((1.0-0.5*inv_x2)*erf_x + 0.5*inv_x*derf_dx);
			value_type d2g0_dx2 = ccf_g0*(inv_x3*erf_x - inv_x2*derf_dx);

			value_type d2x_dvp2 = (1.0 - (vp_mod*vp_mod)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;
			value_type d2x_dudvp = -B_val*vp_mod*u*inv_v_mod_sq*inv_2_vth_b_sq_x_val;
			value_type d2x_du2 = B_val*(1.0 - B_val*(u*u)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;

			h_tot_mat(1,0) = dh0_dx*dx_dvp;
			h_tot_mat(0,1)= dh0_dx*dx_du;

			g_tot_mat(2,0) = dg0_dx*d2x_dvp2 + dx_dvp*dx_dvp*d2g0_dx2;
			g_tot_mat(1,1)= dg0_dx*d2x_dudvp + dx_dvp*dx_du*d2g0_dx2;
			g_tot_mat(0,2) = dg0_dx*d2x_du2 + dx_du*dx_du*d2g0_dx2;

			if(pitch_angle_ei_col_loc_op == 1)
			{ 
				h_tot_mat.setZero();
				g_tot_mat.setZero();

				//pitch_angle_g_fac_mod = 0.0;
				pitch_angle_g_fac_mod = d2g0_dx2*x_val - dg0_dx;
				value_type mod_inv_x3 = inv_x3;

				if(ei_pitch_no_v_in_nu_ei_op == 1)
				{
					if(x_val < min_x_val_loc)
						mod_inv_x3 = const_inv_x3_loc;
				}

				pitch_angle_g_fac_mod *= mod_inv_x3/(4.0*vth_b_sq*vth_b_sq)*B_val;

				g_tot_mat(2,0) = -u*u*pitch_angle_g_fac_mod;
				g_tot_mat(1,1)= vp_mod*u*pitch_angle_g_fac_mod;
				g_tot_mat(0,2) = -vp_mod*vp_mod*pitch_angle_g_fac_mod;
			}

			//h part for mixed flux of upwind and LF
			for(int j = 0; j < 4; j++)
			{
				//int index_j_tmp = j*(4*dof_v*dof_v);
				for(int lvp = 0; lvp < dof_v; lvp++)
				{
					for(int lv = 0; lv < dof_v; lv++)
					{
						value_type basis_v_tot;
						switch(j)
						{
							case 0:
								basis_v_tot = zeta_lv_l_val_arr[lvp]*zeta_lv_l_val_arr[lv];
								break;

							case 1:
								basis_v_tot = -zeta_lv_l_val_arr[lvp]*zeta_lv_r_val_arr[lv];
								break;

							case 2:
								basis_v_tot = zeta_lv_r_val_arr[lvp]*zeta_lv_l_val_arr[lv];
								break;

							case 3:
								basis_v_tot = -zeta_lv_r_val_arr[lvp]*zeta_lv_r_val_arr[lv];
								break;

							default:
								cerr << "Unknown j type : " << j << endl;
								exit(EXIT_FAILURE);
						}

						Vector Gj(1);

						if(i_op == 0) Gj[0] = h_tot_mat(1,0);
						else if(i_op == 1) Gj[0] = h_tot_mat(0,1);

						int tmp_index = ((k1*4 + j)*dof_v + lvp)*dof_v + lv;
						tmp_hM_S_Vmat_Vec[tmp_index] += tmp_wbf0*basis_v_tot*Gj[0];
					}
				}
			}



			Vector fhat_basis_Vec = basis->fhat_vbasis_Vec(ql, i_op, 0); //fhat at edge with original fhat index
			for(int row = 0; row < dim_fhat_valid_v_basis; row++) 
			{
				int fhat_org_basis_i = basis->fhat_valid_to_org_i(row, i_op);
				//fhat_loc_arr[0] : fhat
				//fhat_loc_arr[1] : dfhat/dvp
				//fhat_loc_arr[2] : dfhat/du
				Vector fhat_loc_arr = Vector::Zero(3);
				fhat_loc_arr[0] = fhat_basis_Vec(fhat_org_basis_i);
				fhat_loc_arr[1] = fhat_basis_Vec(dim_fhat_org_v_basis + fhat_org_basis_i);
				fhat_loc_arr[2] = fhat_basis_Vec(2*dim_fhat_org_v_basis + fhat_org_basis_i);

				for(int i_ele = 0; i_ele < 2; i_ele++)
				{
					value_type n_dot_v = 1.0;
					if(i_ele == 1) n_dot_v = -1.0;

					for(int lv = 0; lv < dof_v; lv++)
					{
						Vector basis_loc_arr = Vector::Zero(3);
						basis_loc_arr[0] = basisV_lr[i_ele][lv*6 + 0];
						basis_loc_arr[1] = basisV_lr[i_ele][lv*6 + 1];
						basis_loc_arr[2] = basisV_lr[i_ele][lv*6 + 2];

						Vector Gj(2);
						if(i_op == 0)
						{
							Gj[0] = g_tot_mat(2,0); 
							Gj[1] = g_tot_mat(1,1); 
						}
						else
						{
							Gj[0] = g_tot_mat(1,1); 
							Gj[1] = g_tot_mat(0,2); 
						}

						for(int j = 0; j < 2; j++)
						{
							value_type tmp_val = 0.0;
							if(j == 0)
							{
								tmp_val += -basis_loc_arr[0]*fhat_loc_arr[1];
								tmp_val += -basis_loc_arr[0]*fhat_loc_arr[0]*dlnfMV_dv_arr[1];
								tmp_val += basis_loc_arr[1]*fhat_loc_arr[0];
							}
							else
							{
								tmp_val += -basis_loc_arr[0]*fhat_loc_arr[2];
								tmp_val += -basis_loc_arr[0]*fhat_loc_arr[0]*dlnfMV_dv_arr[2];
								tmp_val += basis_loc_arr[2]*fhat_loc_arr[0];
							}

							Gj[j] *= n_dot_v*tmp_val;

							int tmp_index = (((k1*2 + j)*dim_fhat_valid_v_basis + row)*2 + i_ele)*dof_v + lv;
							tmp_gM_S_Vmat_Vec[tmp_index] += tmp_wbf0*Gj[j];
						}
					}
				}
			}
		}
	}
}


void Integration_Col::col3_g_S_Vmat_bd_cal(const ElementX &ti, const EdgeV &te, const ElementV &kl, const ElementV &kr, const int &i_op, Vector &tmp_g_S_Vmat_Vec)
{
	const auto nb = quadrature->size(3);
	int dim_fhat_valid_v_basis = basis->fhat_valid_vbasis_dim_out();

	int i_ele;
	if (int(kl) == -1) i_ele = 1;
	else if (int(kr) == -1) i_ele = 0;
	else 
	{
		cout << "wrong kl, kr in col3_g_S_Vmat_bd_cal : " << int(kl) << " " << int(kr) << endl;
		abort();
	}

	vector<Matrix> hg_basis_mat_arr[4];
	for(int i = 0; i < 4; i++)
		hg_basis_mat_arr[i].resize(4, Matrix::Zero(3,3));

	for(int b = 0; b < nb; b++)
	{
		value_type wb;
		Point2 pb, qb;
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(te, b); //qb[0] is a dummy variable


		Point4 p = Point2({1.0, 0}) & pb;
		value_type vp = pb[0];
		value_type u = pb[1];

		//fM_wV part
		value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_id, flag_fcm);

		value_type tmp_wb1 = wb * u;
		value_type tmp_wbf0 = tmp_wb1 * fM_wV;

		//hg_basis part
		int loc_vp_tot_index = int((vp - col3_int_vp_tot_min)*col3_int_inv_dvp);
		int loc_u_tot_index = int((u - col3_int_u_tot_min)*col3_int_inv_du);

		int vp_tot_min_index = max(0, loc_vp_tot_index - 1);
		int vp_tot_max_index = min(col3_int_vp_n + 2, loc_vp_tot_index + 2);
		int u_tot_min_index = max(0, loc_u_tot_index - 1);
		int u_tot_max_index = min(col3_int_u_n + 2, loc_u_tot_index + 2);

		for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
		{
			int vp_i_index = j3 - vp_tot_min_index;
			value_type t1_tot = (vp - (col3_int_vp_tot_min + j3*col3_int_dvp))*col3_int_inv_dvp;
			Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_int_inv_dvp);

			for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
			{
				int u_i_index = j4 - u_tot_min_index;
				value_type t2_tot = (u - (col3_int_u_tot_min + j4*col3_int_du))*col3_int_inv_du;
				Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_int_inv_du);

				for(int i1 = 0; i1 < 3; i1++)
				{
					for(int i2 = 0; i2 < 3; i2++)
					{
						hg_basis_mat_arr[vp_i_index][u_i_index](i1, i2) = spline_vp_tot_vec(i1)*spline_u_tot_vec(i2);
					}
				}
			}
		}


		int a= 0;
		Point4 q_t;
		if(i_ele == 0)
		{
			q_t = quadrature->get_traced_point(ti, te, kl, a, b);
		}
		else
		{
			q_t = quadrature->get_traced_point(ti, te, kr, a, b);
		}
		Point2 qb_t = Point2({q_t[2], q_t[3]});

		vector<value_type> basisV_t(dof_v*6);
		basis->col_velocity_basis_vec(qb_t, dvp, du, basisV_t);

		value_type n_dot_v = 1.0;
		if(i_ele == 1) n_dot_v = -1.0;

		for(int row = 0; row < dof_v; row++) 
		{
			Vector fhat_loc_arr = Vector::Zero(3);
			fhat_loc_arr[0] = basisV_t[row*6 + 0];

			for(int lv = 0; lv < dof_v; lv++)
			{

				Vector basis_loc_arr = Vector::Zero(3);
				basis_loc_arr[0] = basisV_t[lv*6 + 0];
				basis_loc_arr[1] = basisV_t[lv*6 + 1];
				basis_loc_arr[2] = basisV_t[lv*6 + 2];

				for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
				{
					int vp_i_index = j3 - vp_tot_min_index;

					for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
					{
						int u_i_index = j4 - u_tot_min_index;

						Matrix tmp_hg_basis_mat = hg_basis_mat_arr[vp_i_index][u_i_index];

						int l = j4*(col3_int_vp_n + 3) + j3;

						Vector Gj(2);
						if(i_op == 2)
						{
							Gj[0] = tmp_hg_basis_mat(2,0); 
							Gj[1] = tmp_hg_basis_mat(1,1); 
						}
						else
						{
							Gj[0] = tmp_hg_basis_mat(1,1); 
							Gj[1] = tmp_hg_basis_mat(0,2); 
						}

						for(int j = 0; j < 2; j++)
						{
							value_type tmp_val = 0.0;
							if(j == 0)
							{
								tmp_val += basis_loc_arr[1]*fhat_loc_arr[0];
							}
							else
							{
								tmp_val += basis_loc_arr[2]*fhat_loc_arr[0];
							}

							Gj[j] *= n_dot_v*tmp_val;

							int tmp_index = ((j*dim_fhat_valid_v_basis + row)*2 + i_ele)*dof_v + lv;
							int tmp_index2 = tmp_index*col3_int_hg_tot_n + l;
							tmp_g_S_Vmat_Vec[tmp_index2] += tmp_wbf0*Gj[j];
						}
					}
				}
			}
		}
	}
}

void Integration_Col::col3_gM_S_Vmat_bd_cal(const int &sp_b_id, const ElementX &ti, const EdgeV &te, const ElementV &kl, const ElementV &kr, const int &i_op, const Vector &tmp_avged_Q, const int &col3_hMbgMb_with_UaTa_op, Vector &tmp_gM_S_Vmat_Vec)
{
	int sp_a_id = sp_id; 
	const auto nb = quadrature->size(3);
	int dim_fhat_valid_v_basis = basis->fhat_valid_vbasis_dim_out();

	vector<Vector> nUT_b_arr, ccf_fMb_arr;
	nUT_b_arr.resize(2, Vector::Zero(3));
	ccf_fMb_arr.resize(2, Vector::Zero(2));

	int k1_max = 1;
	if(col3_hMbgMb_with_UaTa_op == 1) k1_max = 2;

	int pitch_angle_ei_col_loc_op = 0;
	if(sp_a_id == 0 && sp_b_id != 0 && ei_pitch_angle_col_op == 1)
		pitch_angle_ei_col_loc_op = 1;

	for(int k1 = 0; k1 < k1_max; k1++)
	{
		for(int k2 = 0; k2 < 3; k2++)
		{
			nUT_b_arr[k1][k2] = tmp_avged_Q[k1*3 + k2];
		}

		if(pitch_angle_ei_col_loc_op == 1)
			nUT_b_arr[k1][1] = 0.0;

		ccf_fMb_arr[k1][0] = M_1_PI/sqrt(32.0)*nUT_b_arr[k1][0]/sqrt(nUT_b_arr[k1][2]);
		ccf_fMb_arr[k1][1] = ccf_fMb_arr[k1][0]*nUT_b_arr[k1][2];
	}

	value_type B_val = tmp_avged_Q[6];
	value_type pitch_angle_g_fac_mod;

	int i_ele;
	if (int(kl) == -1) i_ele = 1;
	else if (int(kr) == -1) i_ele = 0;
	else 
	{
		cout << "wrong kl, kr in col3_g_S_Vmat_bd_cal : " << int(kl) << " " << int(kr) << endl;
		abort();
	}

	for(int b = 0; b < nb; b++)
	{
		value_type wb;
		Point2 pb, qb;
		tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(te, b); //qb[0] is a dummy variable


		Point4 p = Point2({1.0, 0}) & pb;
		value_type vp = pb[0];
		value_type u = pb[1];

		//fM_wV part
		value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_id, flag_fcm);

		value_type tmp_wb1 = wb * u;
		value_type tmp_wbf0 = tmp_wb1 * fM_wV;

		int a= 0;
		Point4 q_t;
		if(i_ele == 0)
		{
			q_t = quadrature->get_traced_point(ti, te, kl, a, b);
		}
		else
		{
			q_t = quadrature->get_traced_point(ti, te, kr, a, b);
		}
		Point2 qb_t = Point2({q_t[2], q_t[3]});

		vector<value_type> basisV_t(dof_v*6);
		basis->col_velocity_basis_vec(qb_t, dvp, du, basisV_t);

		value_type n_dot_v = 1.0;
		if(i_ele == 1) n_dot_v = -1.0;

		Matrix g_tot_mat = Matrix::Zero(3,3);

		for(int k1 = 0; k1 < k1_max; k1++)
		{
			value_type vp_mod = vp - nUT_b_arr[k1][1];
			value_type v_mod_sq = vp_mod*vp_mod + u*u*B_val;
			value_type v_mod = sqrt(v_mod_sq);
			value_type inv_v_mod_sq = 1.0/v_mod_sq;

			value_type vth_b_sq = nUT_b_arr[k1][2];
			value_type x_val = v_mod/sqrt(2.0*vth_b_sq), x_val2, inv_x, inv_x2, inv_x3;
			x_val2 = x_val*x_val;
			inv_x = 1.0/x_val;
			inv_x2 = inv_x*inv_x;
			inv_x3 = inv_x2*inv_x;

			value_type inv_2_vth_b_sq_x_val = 1.0/(2.0*vth_b_sq*x_val);

			value_type erf_x = erf(x_val), derf_dx = 1.128379167095513*exp(-x_val2);
			value_type dx_du = u*B_val*inv_2_vth_b_sq_x_val;
			value_type dx_dvp = vp_mod*inv_2_vth_b_sq_x_val;

			value_type ccf_h0 = ccf_fMb_arr[k1][0];
			value_type ccf_g0 = ccf_fMb_arr[k1][1];

			value_type dh0_dx = ccf_h0*(-inv_x2*erf_x + inv_x*derf_dx);
			value_type dg0_dx = ccf_g0*((1.0-0.5*inv_x2)*erf_x + 0.5*inv_x*derf_dx);
			value_type d2g0_dx2 = ccf_g0*(inv_x3*erf_x - inv_x2*derf_dx);

			value_type d2x_dvp2 = (1.0 - (vp_mod*vp_mod)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;
			value_type d2x_dudvp = -B_val*vp_mod*u*inv_v_mod_sq*inv_2_vth_b_sq_x_val;
			value_type d2x_du2 = B_val*(1.0 - B_val*(u*u)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;

			g_tot_mat(2,0) = dg0_dx*d2x_dvp2 + dx_dvp*dx_dvp*d2g0_dx2;
			g_tot_mat(1,1)= dg0_dx*d2x_dudvp + dx_dvp*dx_du*d2g0_dx2;
			g_tot_mat(0,2) = dg0_dx*d2x_du2 + dx_du*dx_du*d2g0_dx2;

			if(pitch_angle_ei_col_loc_op == 1)
			{ 
				g_tot_mat.setZero();

				pitch_angle_g_fac_mod = d2g0_dx2*x_val - dg0_dx;
				value_type mod_inv_x3 = inv_x3;

				pitch_angle_g_fac_mod *= mod_inv_x3/(4.0*vth_b_sq*vth_b_sq)*B_val;

				g_tot_mat(2,0) = -u*u*pitch_angle_g_fac_mod;
				g_tot_mat(1,1)= vp_mod*u*pitch_angle_g_fac_mod;
				g_tot_mat(0,2) = -vp_mod*vp_mod*pitch_angle_g_fac_mod;
			}

			for(int row = 0; row < dof_v; row++) 
			{
				Vector fhat_loc_arr = Vector::Zero(3);
				fhat_loc_arr[0] = basisV_t[row*6 + 0];

				for(int lv = 0; lv < dof_v; lv++)
				{

					Vector basis_loc_arr = Vector::Zero(3);
					basis_loc_arr[0] = basisV_t[lv*6 + 0];
					basis_loc_arr[1] = basisV_t[lv*6 + 1];
					basis_loc_arr[2] = basisV_t[lv*6 + 2];

					Vector Gj(2);
					if(i_op == 2)
					{
						Gj[0] = g_tot_mat(2,0); 
						Gj[1] = g_tot_mat(1,1); 
					}
					else
					{
						Gj[0] = g_tot_mat(1,1); 
						Gj[1] = g_tot_mat(0,2); 
					}

					for(int j = 0; j < 2; j++)
					{
						value_type tmp_val = 0.0;
						if(j == 0)
						{
							tmp_val += basis_loc_arr[1]*fhat_loc_arr[0];
						}
						else
						{
							tmp_val += basis_loc_arr[2]*fhat_loc_arr[0];
						}

						Gj[j] *= n_dot_v*tmp_val;

						int tmp_index = (((k1*2 + j)*dim_fhat_valid_v_basis + row)*2 + i_ele)*dof_v + lv;
						tmp_gM_S_Vmat_Vec[tmp_index] += tmp_wbf0*Gj[j];

					}
				}
			}
		}
	}
}

void Integration_Col::col3_gM_S_Vmat_bd_adj_n_cal(const int &sp_b_id, const ElementX &ti, const EdgeV &te, const ElementV &kl, const ElementV &kr, const int &i_op, const Vector &tmp_avged_Q, const int &col3_hMbgMb_with_UaTa_op, Vector &tmp_gM_S_Vmat_Vec, int &n_inout)
{
	int sp_a_id = sp_id; 
	int nb = n_inout;

	int dim_fhat_valid_v_basis = basis->fhat_valid_vbasis_dim_out();

	vector<Vector> nUT_b_arr, ccf_fMb_arr;
	nUT_b_arr.resize(2, Vector::Zero(3));
	ccf_fMb_arr.resize(2, Vector::Zero(2));

	int k1_max = 1;
	if(col3_hMbgMb_with_UaTa_op == 1) k1_max = 2;

	int pitch_angle_ei_col_loc_op = 0;
	if(sp_a_id == 0 && sp_b_id != 0 && ei_pitch_angle_col_op == 1)
		pitch_angle_ei_col_loc_op = 1;

	for(int k1 = 0; k1 < k1_max; k1++)
	{
		for(int k2 = 0; k2 < 3; k2++)
		{
			nUT_b_arr[k1][k2] = tmp_avged_Q[k1*3 + k2];
		}

		if(pitch_angle_ei_col_loc_op == 1)
			nUT_b_arr[k1][1] = 0.0;

		ccf_fMb_arr[k1][0] = M_1_PI/sqrt(32.0)*nUT_b_arr[k1][0]/sqrt(nUT_b_arr[k1][2]);
		ccf_fMb_arr[k1][1] = ccf_fMb_arr[k1][0]*nUT_b_arr[k1][2];
	}

	value_type B_val = tmp_avged_Q[6];
	value_type pitch_angle_g_fac_mod;

	int i_ele;
	if (int(kl) == -1) i_ele = 1;
	else if (int(kr) == -1) i_ele = 0;
	else 
	{
		cout << "wrong kl, kr in col3_g_S_Vmat_bd_cal : " << int(kl) << " " << int(kr) << endl;
		abort();
	}

	int u_index_kl_loc = int(kl)/col3_int_org_f_vp_n;
	int vp_index_kl_loc = int(kl) - u_index_kl_loc*col3_int_org_f_vp_n;
	int u_index_kr_loc = int(kr)/col3_int_org_f_vp_n;
	int vp_index_kr_loc = int(kr) - u_index_kr_loc*col3_int_org_f_vp_n;

	value_type wb;

	if(i_op == 2)
		wb = col3_int_du/nb;
	else
		wb = col3_int_dvp/nb;

	value_type del_qb = 2.0/nb;
	value_type qb1 = 0.0, qb2, pb1, pb2;

	for(int b = 0; b < nb; b++)
	{
		qb2 = -1.0 + (b + 0.5)*del_qb;
		if (i_op == 2)
		{
			if(i_ele == 0)
			{
				pb1 = col3_int_vp_min + (vp_index_kl_loc + 1)*dvp;
				pb2 = col3_int_u_min + (u_index_kl_loc + (qb2 + 1.0)*0.5)*du;
			}
			else
			{
				pb1 = col3_int_vp_min + (vp_index_kr_loc)*dvp;
				pb2 = col3_int_u_min + (u_index_kr_loc + (qb2 + 1.0)*0.5)*du;
			}
		}
		else
		{
			if(i_ele == 0)
			{
				pb1 = col3_int_vp_min + (vp_index_kl_loc + (qb2 + 1.0)*0.5)*dvp;
				pb2 = col3_int_u_min + (u_index_kl_loc + 1)*du;
			}
			else
			{
				pb1 = col3_int_vp_min + (vp_index_kr_loc + (qb2 + 1.0)*0.5)*dvp;
				pb2 = col3_int_u_min + (u_index_kr_loc)*du;
			}
		}

		Point2 pb = Point2({pb1, pb2}), qb = Point2({qb1, qb2});

		Point4 p = Point2({1.0, 0}) & pb;
		value_type vp = pb[0];
		value_type u = pb[1];

		//fM_wV part
		value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_id, flag_fcm);

		value_type tmp_wb1 = wb * u;
		value_type tmp_wbf0 = tmp_wb1 * fM_wV;

		int a= 0;
		Point4 q_t;
		if(i_ele == 0)
		{
			q_t = quadrature->get_traced_point(ti, te, kl, a, b);

			if (i_op == 2)
			{
				q_t[2] = 1.0;
				q_t[3] = qb2;
			}
			else
			{
				q_t[2] = qb2;
				q_t[3] = 1.0;
			}
		}
		else
		{
			q_t = quadrature->get_traced_point(ti, te, kr, a, b);

			if (i_op == 2)
			{
				q_t[2] = -1.0;
				q_t[3] = qb2;
			}
			else
			{
				q_t[2] = qb2;
				q_t[3] = -1.0;
			}
		}
		Point2 qb_t = Point2({q_t[2], q_t[3]});

		vector<value_type> basisV_t(dof_v*6);
		basis->col_velocity_basis_vec(qb_t, dvp, du, basisV_t);

		value_type n_dot_v = 1.0;
		if(i_ele == 1) n_dot_v = -1.0;

		Matrix g_tot_mat = Matrix::Zero(3,3);

		for(int k1 = 0; k1 < k1_max; k1++)
		{
			value_type vp_mod = vp - nUT_b_arr[k1][1];
			value_type v_mod_sq = vp_mod*vp_mod + u*u*B_val;
			value_type v_mod = sqrt(v_mod_sq);
			value_type inv_v_mod_sq = 1.0/v_mod_sq;

			value_type vth_b_sq = nUT_b_arr[k1][2];
			value_type x_val = v_mod/sqrt(2.0*vth_b_sq), x_val2, inv_x, inv_x2, inv_x3;
			x_val2 = x_val*x_val;
			inv_x = 1.0/x_val;
			inv_x2 = inv_x*inv_x;
			inv_x3 = inv_x2*inv_x;

			value_type inv_2_vth_b_sq_x_val = 1.0/(2.0*vth_b_sq*x_val);

			value_type erf_x = erf(x_val), derf_dx = 1.128379167095513*exp(-x_val2);
			value_type dx_du = u*B_val*inv_2_vth_b_sq_x_val;
			value_type dx_dvp = vp_mod*inv_2_vth_b_sq_x_val;

			value_type ccf_h0 = ccf_fMb_arr[k1][0];
			value_type ccf_g0 = ccf_fMb_arr[k1][1];

			value_type dh0_dx = ccf_h0*(-inv_x2*erf_x + inv_x*derf_dx);
			value_type dg0_dx = ccf_g0*((1.0-0.5*inv_x2)*erf_x + 0.5*inv_x*derf_dx);
			value_type d2g0_dx2 = ccf_g0*(inv_x3*erf_x - inv_x2*derf_dx);

			value_type d2x_dvp2 = (1.0 - (vp_mod*vp_mod)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;
			value_type d2x_dudvp = -B_val*vp_mod*u*inv_v_mod_sq*inv_2_vth_b_sq_x_val;
			value_type d2x_du2 = B_val*(1.0 - B_val*(u*u)*inv_v_mod_sq)*inv_2_vth_b_sq_x_val;

			g_tot_mat(2,0) = dg0_dx*d2x_dvp2 + dx_dvp*dx_dvp*d2g0_dx2;
			g_tot_mat(1,1)= dg0_dx*d2x_dudvp + dx_dvp*dx_du*d2g0_dx2;
			g_tot_mat(0,2) = dg0_dx*d2x_du2 + dx_du*dx_du*d2g0_dx2;

			if(pitch_angle_ei_col_loc_op == 1)
			{ 
				g_tot_mat.setZero();

				pitch_angle_g_fac_mod = d2g0_dx2*x_val - dg0_dx;
				value_type mod_inv_x3 = inv_x3;

				pitch_angle_g_fac_mod *= mod_inv_x3/(4.0*vth_b_sq*vth_b_sq)*B_val;

				g_tot_mat(2,0) = -u*u*pitch_angle_g_fac_mod;
				g_tot_mat(1,1)= vp_mod*u*pitch_angle_g_fac_mod;
				g_tot_mat(0,2) = -vp_mod*vp_mod*pitch_angle_g_fac_mod;
			}

			for(int row = 0; row < dof_v; row++) 
			{
				Vector fhat_loc_arr = Vector::Zero(3);
				fhat_loc_arr[0] = basisV_t[row*6 + 0];

				for(int lv = 0; lv < dof_v; lv++)
				{

					Vector basis_loc_arr = Vector::Zero(3);
					basis_loc_arr[0] = basisV_t[lv*6 + 0];
					basis_loc_arr[1] = basisV_t[lv*6 + 1];
					basis_loc_arr[2] = basisV_t[lv*6 + 2];

					Vector Gj(2);
					if(i_op == 2)
					{
						Gj[0] = g_tot_mat(2,0); 
						Gj[1] = g_tot_mat(1,1); 
					}
					else
					{
						Gj[0] = g_tot_mat(1,1); 
						Gj[1] = g_tot_mat(0,2); 
					}

					for(int j = 0; j < 2; j++)
					{
						value_type tmp_val = 0.0;
						if(j == 0)
						{
							tmp_val += basis_loc_arr[1]*fhat_loc_arr[0];
						}
						else
						{
							tmp_val += basis_loc_arr[2]*fhat_loc_arr[0];
						}

						Gj[j] *= n_dot_v*tmp_val;

						int tmp_index = (((k1*2 + j)*dim_fhat_valid_v_basis + row)*2 + i_ele)*dof_v + lv;
						tmp_gM_S_Vmat_Vec[tmp_index] += tmp_wbf0*Gj[j];

					}
				}
			}
		}
	}
}

void Integration_Col::col3_h_flux_qd_cal(const ElementX &ti, const EdgeV &te, const int &i_op, const Vector &qd_points, Vector &tmp_h_flux_qd)
{
	int nb_h_flux = qd_points.size();

	vector<Matrix> hg_basis_mat_arr[4];
	for(int i = 0; i < 4; i++)
		hg_basis_mat_arr[i].resize(4, Matrix::Zero(3,3));

	vector<Point2> pnt = mesh->get_node_edge(te);
	value_type del_vp = pnt[1][0] - pnt[0][0];
	value_type del_u = pnt[1][1] - pnt[0][1];

	for(int b = 0; b < nb_h_flux; b++)
	{
		value_type vp = pnt[0][0] + del_vp*(0.5*(qd_points[b] + 1.0));
		value_type u  = pnt[0][1] + del_u*(0.5*(qd_points[b] + 1.0));

		//hg_basis part
		int loc_vp_tot_index = int((vp - col3_int_vp_tot_min)*col3_int_inv_dvp);
		int loc_u_tot_index = int((u - col3_int_u_tot_min)*col3_int_inv_du);

		int vp_tot_min_index = max(0, loc_vp_tot_index - 1);
		int vp_tot_max_index = min(col3_int_vp_n + 2, loc_vp_tot_index + 2);
		int u_tot_min_index = max(0, loc_u_tot_index - 1);
		int u_tot_max_index = min(col3_int_u_n + 2, loc_u_tot_index + 2);

		for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
		{
			int vp_i_index = j3 - vp_tot_min_index;
			value_type t1_tot = (vp - (col3_int_vp_tot_min + j3*col3_int_dvp))*col3_int_inv_dvp;
			Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_int_inv_dvp);

			for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
			{
				int u_i_index = j4 - u_tot_min_index;
				value_type t2_tot = (u - (col3_int_u_tot_min + j4*col3_int_du))*col3_int_inv_du;
				Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_int_inv_du);

				for(int i1 = 0; i1 < 3; i1++)
				{
					for(int i2 = 0; i2 < 3; i2++)
					{
						hg_basis_mat_arr[vp_i_index][u_i_index](i1, i2) = spline_vp_tot_vec(i1)*spline_u_tot_vec(i2);
					}
				}
			}
		}

		for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
		{
			int vp_i_index = j3 - vp_tot_min_index;

			for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
			{
				int u_i_index = j4 - u_tot_min_index;

				Matrix tmp_hg_basis_mat = hg_basis_mat_arr[vp_i_index][u_i_index];

				int l = j4*(col3_int_vp_n + 3) + j3;
				Vector Gj(1);

				if(i_op == 0)
				{
					Gj[0] = tmp_hg_basis_mat(1,0);
				}
				else
				{
					Gj[0] = tmp_hg_basis_mat(0,1);
				}

				int tmp_index2 = b*col3_int_hg_tot_n + l;
				tmp_h_flux_qd[tmp_index2] = Gj[0];
			}
		}
	}
}

void Integration_Col::col3_hM_flux_qd_cal(const int &sp_b_id, const ElementX &ti, const EdgeV &te, const int &i_op, const Vector &qd_points, const Vector &tmp_avged_Q, const int &col3_hMbgMb_with_UaTa_op, Vector &tmp_hM_flux_qd)
{
	int sp_a_id = sp_id;
	int nb_h_flux = qd_points.size();

	vector<Matrix> hg_basis_mat_arr[4];
	for(int i = 0; i < 4; i++)
		hg_basis_mat_arr[i].resize(4, Matrix::Zero(3,3));

	vector<Point2> pnt = mesh->get_node_edge(te);
	value_type del_vp = pnt[1][0] - pnt[0][0];
	value_type del_u = pnt[1][1] - pnt[0][1];

	vector<Vector> nUT_b_arr, ccf_fMb_arr;
	nUT_b_arr.resize(2, Vector::Zero(3));
	ccf_fMb_arr.resize(2, Vector::Zero(2));

	int k1_max = 1;
	if(col3_hMbgMb_with_UaTa_op == 1) k1_max = 2;

	int pitch_angle_ei_col_loc_op = 0;
	if(sp_a_id == 0 && sp_b_id != 0 && ei_pitch_angle_col_op == 1)
		pitch_angle_ei_col_loc_op = 1;

	for(int k1 = 0; k1 < k1_max; k1++)
	{
		for(int k2 = 0; k2 < 3; k2++)
		{
			nUT_b_arr[k1][k2] = tmp_avged_Q[k1*3 + k2];
		}

		ccf_fMb_arr[k1][0] = M_1_PI/sqrt(32.0)*nUT_b_arr[k1][0]/sqrt(nUT_b_arr[k1][2]);
		ccf_fMb_arr[k1][1] = ccf_fMb_arr[k1][0]*nUT_b_arr[k1][2];
	}

	value_type B_val = tmp_avged_Q[6];

	for(int b = 0; b < nb_h_flux; b++)
	{
		value_type vp = pnt[0][0] + del_vp*(0.5*(qd_points[b] + 1.0));
		value_type u  = pnt[0][1] + del_u*(0.5*(qd_points[b] + 1.0));

		Matrix h_tot_mat = Matrix::Zero(3,3);
		Matrix g_tot_mat = Matrix::Zero(3,3);

		for(int k1 = 0; k1 < k1_max; k1++)
		{
			value_type vp_mod = vp - nUT_b_arr[k1][1];
			value_type v_mod = sqrt(vp_mod*vp_mod + u*u*B_val);

			value_type vth_b_sq = nUT_b_arr[k1][2];
			value_type x_val = v_mod/sqrt(2.0*vth_b_sq), x_val2, inv_x, inv_x2, inv_x3;
			x_val2 = x_val*x_val;
			inv_x = 1.0/x_val;
			inv_x2 = inv_x*inv_x;
			inv_x3 = inv_x2*inv_x;

			value_type inv_2_vth_b_sq_x_val = 1.0/(2.0*vth_b_sq*x_val);

			value_type erf_x = erf(x_val), derf_dx = 1.128379167095513*exp(-x_val2);
			value_type dx_du = u*B_val*inv_2_vth_b_sq_x_val;
			value_type dx_dvp = vp_mod*inv_2_vth_b_sq_x_val;

			value_type ccf_h0 = ccf_fMb_arr[k1][0];
			value_type dh0_dx = ccf_h0*(-inv_x2*erf_x + inv_x*derf_dx);

			h_tot_mat(1,0) = dh0_dx*dx_dvp;
			h_tot_mat(0,1)= dh0_dx*dx_du;

			if(pitch_angle_ei_col_loc_op == 1)
			{ 
				h_tot_mat.setZero();
			}

			Vector Gj(1);

			if(i_op == 0)
			{
				Gj[0] = h_tot_mat(1,0);
			}
			else
			{
				Gj[0] = h_tot_mat(0,1);
			}

			int tmp_index = k1*nb_h_flux + b;
			tmp_hM_flux_qd[tmp_index] = Gj[0];
		}
	}
}

void Integration_Col::col3_f_to_consv_vol_cal(const ElementX &ti, const ElementV &tj, Vector &f_to_consv_vol_Vec)
{
	int qd_sizeX = quadrature->size(0), qd_sizeV = quadrature->size(1);

	vector<value_type> wb_arr(qd_sizeV);
	vector<Point2> pb_arr(qd_sizeV), qb_arr(qd_sizeV);

	for(int b = 0; b < qd_sizeV; b++)
		tie(wb_arr[b], pb_arr[b], qb_arr[b]) = quadrature->get_quadrature_single_not_traced_1d_qd_point(tj, b);


	for(int a = 0; a < qd_sizeX; a++)
	{
		value_type wa;
		Point2 pa, qa;
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 
		Point4 p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));
		value_type B_inv_val = 1.0/B_val;

		value_type R = pa[0], Z = pa[1];

		value_type fM_wX = flux->get_fcm_SE_X(p, int(ti), sp_id, flag_fcm);

		value_type tmp_wa1 = wa * R;
		value_type tmp_waf0 = B_val * tmp_wa1 * fM_wX;

		vector<value_type> basisX(dof_x);
		basis->col_spatial_basis_vec(qa, basisX);

		for(int b = 0; b < qd_sizeV; b++)
		{
			p = pa & pb_arr[b];

			value_type v_para = p[2], u = p[3];
			value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_id, flag_fcm);

			value_type tmp_wb1 = wb_arr[b] * u;
			value_type tmp_wbf0 = tmp_wb1 * fM_wV;

			value_type tmp_wf0 = tmp_waf0*tmp_wbf0;

			vector<value_type> basisV(dof_v*6);
			basis->col_velocity_basis_vec(qb_arr[b], dvp, du, basisV);

			for(int lvp = 0; lvp < dof_v; lvp++)
			{
				value_type zeta_lvp_val = basisV[lvp*6 + 0];

				for(int lv = 0; lv < dof_v; lv++)
				{
					value_type zeta_dlv_dvp = basisV[lv*6 + 1];
					value_type zeta_dlv_du = basisV[lv*6 + 2];
					value_type zeta_d2lv_du2 = basisV[lv*6 + 5];

					for(int lxp = 0; lxp < dof_x; lxp++)
					{
						value_type zeta_lxp_val = basisX[lxp];
						value_type zeta_lp_val = zeta_lvp_val*zeta_lxp_val;

						for(int lx = 0; lx < dof_x; lx++)
						{
							value_type zeta_lx_val = basisX[lx];

							int l_index = lv*dof_x + lx;
							int lp_index = lvp*dof_x + lxp;

							int tmp_index = lp_index*dof + l_index;
							f_to_consv_vol_Vec[tmp_index] += tmp_wf0*zeta_lp_val*zeta_dlv_dvp*zeta_lx_val;


							value_type tmp_en_consv = zeta_dlv_du/u + zeta_d2lv_du2;
							tmp_en_consv *= B_inv_val*zeta_lx_val;
							tmp_index = dof*dof + lp_index*dof + l_index;
							f_to_consv_vol_Vec[tmp_index] += tmp_wf0*zeta_lp_val*tmp_en_consv;

						}
					}
				}
			}
		}
	}
}


void Integration_Col::col3_f_to_mom_consv_surf_cal(const ElementX &ti, const EdgeV &te, const ElementV &kl, const ElementV &kr, Vector &tmp_f_to_mom_consv_surf)
{
	int qd_sizeX = quadrature->size(0), qd_sizeV = quadrature->size(3);

	for(int a = 0; a < qd_sizeX; a++)
	{
		value_type wa;
		Point2 pa, qa;
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 
		Point4 p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));
		value_type B_inv_val = 1.0/B_val;

		value_type R = pa[0], Z = pa[1];

		value_type fM_wX = flux->get_fcm_SE_X(p, int(ti), sp_id, flag_fcm);

		value_type tmp_wa1 = wa * R;
		value_type tmp_waf0 = B_val * tmp_wa1 * fM_wX;

		vector<value_type> basisX(dof_x);
		basis->col_spatial_basis_vec(qa, basisX);

		for(int b = 0; b < qd_sizeV; b++)
		{
			value_type wb;
			Point2 pb, qb;
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(te, b); //qb[0] is a dummy variable

			p = pa & pb;
			value_type vp = pb[0];
			value_type u = pb[1];

			//fM_wV part
			value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_id, flag_fcm);

			value_type tmp_wb1 = wb * u;
			value_type tmp_wbf0 = tmp_wb1 * fM_wV;

			value_type tmp_wf0 = tmp_waf0*tmp_wbf0;

			auto ql = quadrature->get_traced_point(ti, te, kl, a, b);
			auto qr = quadrature->get_traced_point(ti, te, kr, a, b);

			Point2 qb_l = Point2({ql[2], ql[3]});
			Point2 qb_r = Point2({qr[2], qr[3]});

			vector<value_type> basisV_lr[2];
			basisV_lr[0].resize(dof_v*6);
			basisV_lr[1].resize(dof_v*6);
			basis->col_velocity_basis_vec(qb_l, dvp, du, basisV_lr[0]);
			basis->col_velocity_basis_vec(qb_r, dvp, du, basisV_lr[1]);

			for(int type_index = 0; type_index < 2; type_index++)
			{
				int ele_lp;
				if(type_index == 0) ele_lp = 0;
				else ele_lp = 1;

				for(int lp = 0; lp < dof; lp++)
				{
					int lp_tot = ele_lp*dof + lp;

					int lpX = lp % dof_x;
					int lpV = (lp - lpX) / dof_x;

					value_type basis_lp = basisX[lpX]*basisV_lr[ele_lp][lpV*6 + 0];

					for(int ele_l = 0; ele_l < 2; ele_l++)
					{
						value_type n_dot_vp_hat;

						if(ele_l == 0) n_dot_vp_hat = 1.0;
						else n_dot_vp_hat = -1.0;

						for(int l = 0; l < dof; l++)
						{
							int l_tot = ele_l*dof + l;

							int lX = l % dof_x;
							int lV = (l - lX) / dof_x;

							value_type basis_l = basisX[lX]*basisV_lr[ele_l][lV*6 + 0];

							int tmp_index = (type_index*2*dof + lp_tot)*2*dof + l_tot;

							tmp_f_to_mom_consv_surf[tmp_index] += -n_dot_vp_hat*basis_lp*basis_l*tmp_wf0;
						}
					}
				}
			}
		}
	}
}

void Integration_Col::col3_f_to_en_consv_inner_surf_cal(const ElementX &ti, const EdgeV &te, const ElementV &kl, const ElementV &kr, const int &i_op, Vector &tmp_f_to_en_consv_surf)
{
	int qd_sizeX = quadrature->size(0), qd_sizeV = quadrature->size(3);
	int dim_fhat_org_v_basis = basis->fhat_org_vbasis_dim_out();
	int dim_fhat_valid_v_basis = basis->fhat_valid_vbasis_dim_out();

	for(int a = 0; a < qd_sizeX; a++)
	{
		value_type wa;
		Point2 pa, qa;
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 
		Point4 p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));
		value_type B_inv_val = 1.0/B_val;

		value_type R = pa[0], Z = pa[1];

		value_type fM_wX = flux->get_fcm_SE_X(p, int(ti), sp_id, flag_fcm);

		value_type tmp_wa1 = wa * R;
		value_type tmp_waf0 = B_val * tmp_wa1 * fM_wX;

		vector<value_type> basisX(dof_x);
		basis->col_spatial_basis_vec(qa, basisX);

		for(int b = 0; b < qd_sizeV; b++)
		{
			value_type wb;
			Point2 pb, qb;
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(te, b); //qb[0] is a dummy variable

			p = pa & pb;
			value_type vp = pb[0];
			value_type u = pb[1];

			//fM_wV part
			value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_id, flag_fcm);

			vector<value_type> dlnfMV_dv_arr(3);
			dlnfMV_dv_arr = flux->get_fcm_SE_vec(p, int(ti), sp_id, flag_fcm);

			dlnfMV_dv_arr[1] = dlnfMV_dv_arr[1]/dlnfMV_dv_arr[0];
			dlnfMV_dv_arr[2] = dlnfMV_dv_arr[2]/dlnfMV_dv_arr[0];

			value_type tmp_wb1 = wb * u;
			value_type tmp_wbf0 = tmp_wb1 * fM_wV;

			value_type tmp_wf0 = tmp_waf0*tmp_wbf0;

			auto ql = quadrature->get_traced_point(ti, te, kl, a, b);
			auto qr = quadrature->get_traced_point(ti, te, kr, a, b);

			Point2 qb_l = Point2({ql[2], ql[3]});
			Point2 qb_r = Point2({qr[2], qr[3]});

			vector<value_type> basisV_lr[2];
			basisV_lr[0].resize(dof_v*6);
			basisV_lr[1].resize(dof_v*6);
			basis->col_velocity_basis_vec(qb_l, dvp, du, basisV_lr[0]);
			basis->col_velocity_basis_vec(qb_r, dvp, du, basisV_lr[1]);

			Vector fhat_basis_Vec = basis->fhat_vbasis_Vec(ql, i_op, 0); //fhat at edge with original fhat index

			vector<Vector> fhat_vbasis_loc_Vec;
			fhat_vbasis_loc_Vec.resize(dim_fhat_valid_v_basis, Vector::Zero(3));
			for(int lp_v = 0; lp_v < dim_fhat_valid_v_basis; lp_v++) 
			{
				int fhat_org_basis_i = basis->fhat_valid_to_org_i(lp_v, i_op);
				fhat_vbasis_loc_Vec[lp_v][0] = fhat_basis_Vec(fhat_org_basis_i);
				fhat_vbasis_loc_Vec[lp_v][1] = fhat_basis_Vec(dim_fhat_org_v_basis + fhat_org_basis_i);
				fhat_vbasis_loc_Vec[lp_v][2] = fhat_basis_Vec(2*dim_fhat_org_v_basis + fhat_org_basis_i);
			}

			for(int lp_x = 0; lp_x < dof_x; lp_x++)
			{
				for(int lp_v = 0; lp_v < dim_fhat_valid_v_basis; lp_v++) 
				{

					int lp_tot = lp_x*dim_fhat_valid_v_basis + lp_v;

					//fhat_loc_arr[0] : fhat
					//fhat_loc_arr[1] : dfhat/dvp
					//fhat_loc_arr[2] : dfhat/du
					Vector fhat_loc_arr = Vector::Zero(3);
					for(int tmp_i = 0; tmp_i < 3; tmp_i++)
						fhat_loc_arr[tmp_i] = basisX[lp_x]*fhat_vbasis_loc_Vec[lp_v][tmp_i];

					for(int ele_l = 0; ele_l < 2; ele_l++)
					{
						value_type n_dot_v_hat;

						if(ele_l == 0) n_dot_v_hat = 1.0;
						else n_dot_v_hat = -1.0;

						for(int l = 0; l < dof; l++)
						{
							int l_tot = ele_l*dof + l;

							int lX = l % dof_x;
							int lV = (l - lX) / dof_x;

							Vector basis_loc_arr = Vector::Zero(3);
							basis_loc_arr[0] = basisX[lX]*basisV_lr[ele_l][lV*6 + 0];
							basis_loc_arr[1] = basisX[lX]*basisV_lr[ele_l][lV*6 + 1];
							basis_loc_arr[2] = basisX[lX]*basisV_lr[ele_l][lV*6 + 2];

							value_type tmp_val = 0.0;
							if(i_op == 0)
							{
								tmp_val += basis_loc_arr[0]*fhat_loc_arr[1];
								tmp_val += basis_loc_arr[0]*fhat_loc_arr[0]*dlnfMV_dv_arr[1];
								tmp_val += -basis_loc_arr[1]*fhat_loc_arr[0];
							}
							else
							{
								tmp_val += basis_loc_arr[0]*fhat_loc_arr[2];
								tmp_val += basis_loc_arr[0]*fhat_loc_arr[0]*dlnfMV_dv_arr[2];
								tmp_val += -basis_loc_arr[2]*fhat_loc_arr[0];
								tmp_val *= B_inv_val;
							}

							int tmp_index = lp_tot*2*dof + l_tot;
							tmp_f_to_en_consv_surf[tmp_index] += n_dot_v_hat*tmp_val*tmp_wf0;

						}
					}
				}
			}
		}
	}
}

void Integration_Col::col3_f_to_en_consv_bd_surf_cal(const ElementX &ti, const EdgeV &te, const ElementV &kl_bd, const ElementV &kr_bd, const int &i_op, Vector &tmp_f_to_en_consv_surf_bd)
{
	int qd_sizeX = quadrature->size(0), qd_sizeV = quadrature->size(3);
	int dim_fhat_valid_v_basis = basis->fhat_valid_vbasis_dim_out();

	int valid_ele_index;
	ElementV valid_eleV;

	if(int(kl_bd) == -1)
	{
		valid_eleV = kr_bd;
		valid_ele_index = 1;
	}
	else if(int(kr_bd) == -1)
	{
		valid_eleV = kl_bd;
		valid_ele_index = 0;
	}
	else
	{
		cout << "Invalid ele for bd type 2 or 3 in col3_f_to_en_consv_bd_surf_cal : " << int(kl_bd) << " " << int(kr_bd) << " " << i_op << endl;
		abort();
	}

	for(int a = 0; a < qd_sizeX; a++)
	{
		value_type wa;
		Point2 pa, qa;
		tie(wa, pa, qa) = quadrature->get_quadrature_single_not_traced_1d_qd_point(ti, a); 
		Point4 p = pa & Point2({0, 0});
		Vector B = flux->get_B(p);
		value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));
		value_type B_inv_val = 1.0/B_val;

		value_type R = pa[0], Z = pa[1];

		value_type fM_wX = flux->get_fcm_SE_X(p, int(ti), sp_id, flag_fcm);

		value_type tmp_wa1 = wa * R;
		value_type tmp_waf0 = B_val * tmp_wa1 * fM_wX;

		vector<value_type> basisX(dof_x);
		basis->col_spatial_basis_vec(qa, basisX);

		for(int b = 0; b < qd_sizeV; b++)
		{
			value_type wb;
			Point2 pb, qb;
			tie(wb, pb, qb) = quadrature->get_quadrature_single_not_traced_1d_qd_point(te, b); //qb[0] is a dummy variable

			p = pa & pb;
			value_type vp = pb[0];
			value_type u = pb[1];

			//fM_wV part
			value_type fM_wV = flux->get_fcm_SE_V(p, int(ti), sp_id, flag_fcm);

			vector<value_type> dlnfMV_dv_arr(3);
			dlnfMV_dv_arr = flux->get_fcm_SE_vec(p, int(ti), sp_id, flag_fcm);

			dlnfMV_dv_arr[1] = dlnfMV_dv_arr[1]/dlnfMV_dv_arr[0];
			dlnfMV_dv_arr[2] = dlnfMV_dv_arr[2]/dlnfMV_dv_arr[0];

			value_type tmp_wb1 = wb * u;
			value_type tmp_wbf0 = tmp_wb1 * fM_wV;

			value_type tmp_wf0 = tmp_waf0*tmp_wbf0;

			auto q_valid = quadrature->get_traced_point(ti, te, valid_eleV, a, b);
			Point2 qb_valid = Point2({q_valid[2], q_valid[3]});


			vector<value_type> basisV_valid(dof_v*6);
			basis->col_velocity_basis_vec(qb_valid, dvp, du, basisV_valid);

			for(int lp_x = 0; lp_x < dof_x; lp_x++)
			{
				for(int lp_v = 0; lp_v < dof_v; lp_v++) 
				{
					int lp_tot = lp_x*dim_fhat_valid_v_basis + lp_v;

					//fhat_loc_arr[0] : fhat
					//fhat_loc_arr[1] : dfhat/dvp
					//fhat_loc_arr[2] : dfhat/du
					Vector fhat_loc_arr = Vector::Zero(3);
					for(int tmp_i = 0; tmp_i < 3; tmp_i++)
						fhat_loc_arr[tmp_i] = basisX[lp_x]*basisV_valid[lp_v*6 + tmp_i];

					value_type n_dot_v_hat;

					if(valid_ele_index == 0) n_dot_v_hat = 1.0;
					else n_dot_v_hat = -1.0;

					for(int l = 0; l < dof; l++)
					{
						int l_tot = valid_ele_index*dof + l;

						int lX = l % dof_x;
						int lV = (l - lX) / dof_x;

						Vector basis_loc_arr = Vector::Zero(3);
						for(int tmp_i = 0; tmp_i < 3; tmp_i++)
							basis_loc_arr[tmp_i] = basisX[lX]*basisV_valid[lV*6 + tmp_i];

						value_type tmp_val = 0.0;
						if(i_op == 2)
						{
							tmp_val += -basis_loc_arr[1]*fhat_loc_arr[0];
						}
						else
						{
							tmp_val += -basis_loc_arr[2]*fhat_loc_arr[0];
							tmp_val *= B_inv_val;
						}

						int tmp_index = lp_tot*2*dof + l_tot;
						tmp_f_to_en_consv_surf_bd[tmp_index] += n_dot_v_hat*tmp_val*tmp_wf0;

					}
				}
			}
		}
	}
}


