#include <collision_dg.h>

Collision_dg::Collision_dg(Mesh &_mesh, const EQ_reader &_eq_reader, const Config &_config, const Species_data &_species_data, const MPI_Comm &comm) : comm(comm)
{
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &nproc);

	mesh_sp1 = &_mesh;
	eq_reader = &_eq_reader;
	config = &_config;
	species_data = &_species_data;

	col_diag_quantity_num = 7;
	//normalization coefficient for temperature
	col3_norm_T0_eV = eq_reader->ph_const.get_property(Ph_const::Ti0);

	//numerical integration option
	col_adjustable_quad_op = config->get_option<int>("col_adjustable_quadrature_op");
	col_adjustable_quad_tor = config->get_option<double>("col_adjustable_quadrature_tor");

	col_method =  config->get_option<int>("col_method_op"); // dummy option
	col_ion_period =  config->get_option<int>("col_ion_period"); //collision period
	col_fM_update_period = config->get_option<int>("col_fm_update_period"); //lowest order part update period
	col3_T_eV_min_lim = config->get_option<double>("col_temp_ev_min"); //Minimum limit for valid collision
	col_diag_1d_period = config->get_option<int>("diag_1d_output_period"); //diagnostic period

	col_mult_fac = config->get_option<double>("col_mult_fac"); //multiplication factor for the collision frequency
	nstep = config->get_option<int>("sml_nstep"); //total time step for the simulation

	//implicit collision tolerance level
	col_implicit_tor = config->get_option<double>("col_implicit_tor");

	col_fast_to_slow_moment_op = config->get_option<int>("col_fast_to_slow_moment_op");
	col_fast_to_slow_moment_ratio_max = config->get_option<double>("col_fast_to_slow_moment_ratio_max");

	//delta t step size (in R/v_t)
	system_dt = config->get_option<double>("sml_dt");
	col_consv_can_ang_mom = 0; // 0 : v_para conservation, 1 : canonical toroidal angular momentum conservation (not completed yet)

	col_fhat_vol_wid = config->get_option<double>("col_fhat_interpolation_cell_width"); //width for fhat matrix calculation

	//spatial toroidal section # = 1
	tor_wedge_n = config->get_option<int>("eq_tor_wedge_n");

	col3_vdim_factor = 1;  //ratio between hg_v_space_size & f_dg_v space size. for higher value than 1, test are needed.. maybe many qd points required?
	col3_diff_sp_nl_col_op = 1; //different species nonlinear collision part option

	//compensate C[f_a(n_a, U_a, T_a), f_b(n_b, U_a, T_a)] option
	col3_hMbgMb_with_UaTa_op = config->get_option<int>("col_hmbgmb_with_uata_op");
	
	//conservation on-off option
	col_consv_onoff_op = config->get_option<int>("col_consv_onoff_op");

	//pitch angle collision only option : Not used
	col_ei_pitch_hMbgMb_with_UaTa_op = 0;
	col_ei_pitch_angle_op = config->get_option<int>("col_ei_pitch_angle_col_op");
	col_ei_pitch_no_v_in_nu_ei_op = config->get_option<int>("col_ei_pitch_angle_col_no_v_in_nu_ei_op");
	col_ei_pitch_no_v_in_nu_ei_v_e_min = config->get_option<double>("col_ei_pitch_no_v_in_nu_ei_v_e_min");

	//collision general setup
	collision_dg_general_setup();
}

Collision_dg::~Collision_dg()
{

}

//collision general setup
void Collision_dg::collision_dg_general_setup(void)
{
	tot_species_num = species_data->get_tot_sp_num();

	ni00_norm_mks = eq_reader->ph_const.get_property(Ph_const::ni00);
	ni00_norm_mks *= 1e6;

	nx = mesh_sp1->own_size<ElementX>();

	//spatial cell center (R, Z) calculation
	for(int ix = 0; ix < nx; ix++)
	{
		auto tmp_vec = mesh_sp1->get_node_element(ElementX(ix));

		value_type R_cen = 0.0, Z_cen = 0.0;
		for (int i = 0; i < 3; i++) 
		{
			R_cen += tmp_vec[i][0]/3.0;
			Z_cen += tmp_vec[i][1]/3.0;
		}
		R_cen_cell.push_back(R_cen);
		Z_cen_cell.push_back(Z_cen);
	}

	//Gather number_own_node from other computing nodes
	int number_own_ele_table[nproc];
	MPI_Allgather(&nx, 1, MPI_INT, number_own_ele_table, 1, MPI_INT, comm);

	//Define ilower and iupper
	/*
	   int ilower = mesh->get_global_index(NodeX(0));
	   int iupper = ilower + number_own_node - 1;
	   int number_own_node = iupper - ilower + 1;
	   */

	ilower_ele = 0;
	for(int i = 0; i < rank; i++)
		ilower_ele += number_own_ele_table[i]; 

	//col_implicit_linear_ion_op = config->get_option<int>("col_implicit_linear_ion_op");
	//col_implicit_linear_electron_op = config->get_option<int>("col_implicit_linear_electron_op");
	col_implicit_linear_ion_op = 0;
	col_implicit_linear_electron_op = 0;
	col_implicit_linear_sp_op.resize(tot_species_num, 0);

	col3_C_fMa_to_fMb_wUaTa_summed_arr.resize(tot_species_num);
	col3_C_fMa_to_fMb_diff_UbTb_UaTa_summed_arr.resize(tot_species_num);
	col_sp_Minv_arr.resize(tot_species_num);
	col_sp_M_arr.resize(tot_species_num);
	consv_fa_fb0_mat_sum_arr.resize(tot_species_num);
	col_implicit_linear_solver_arr.resize(tot_species_num);

	for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
	{
		col3_C_fMa_to_fMb_wUaTa_summed_arr[sp_id].resize(nx);
		col3_C_fMa_to_fMb_diff_UbTb_UaTa_summed_arr[sp_id].resize(nx);
	}

	nv.resize(tot_species_num);
	col_diag_cell_mat.resize(tot_species_num);
	col_diag_cell_consv_mat.resize(tot_species_num);
	col_fhat_cell_mat.resize(tot_species_num);
	col_f_dg_to_valid_fhat_mat_arr.resize(tot_species_num);
	vspace_edge_list_arr.resize(tot_species_num);
	vspace_edge_ele_arr.resize(tot_species_num);
	dof_arr.resize(tot_species_num);
	dof_x_arr.resize(tot_species_num);
	dof_v_arr.resize(tot_species_num);
	nve_arr.resize(tot_species_num);
	fhat_dim_per_edge_arr.resize(tot_species_num);

	tau_aa_init.resize(tot_species_num);

	col3_f_to_h_bc_mat.resize(tot_species_num);
	col3_f_to_g_bc_mat.resize(tot_species_num);
	col3_f_to_h_source_mat.resize(tot_species_num);
	col3_h_to_g_source_mat.resize(tot_species_num);

	col3_h_solver_stiffness.resize(tot_species_num);
	col3_h_stiffness_mat.resize(tot_species_num);
	col3_g_solver_stiffness.resize(tot_species_num);
	col3_g_stiffness_mat.resize(tot_species_num);

	col3_hg_formula_to_hg0_solver.resize(tot_species_num);

	col3_f_to_mom_consv_mat_arr.resize(tot_species_num);
	col3_f_to_en_consv_mat_arr.resize(tot_species_num);

	vp_tot_min_index_arr.resize(tot_species_num); 
	vp_tot_max_index_arr.resize(tot_species_num);
	u_tot_min_index_arr.resize(tot_species_num);
	u_tot_max_index_arr.resize(tot_species_num);


	col3_M_solver_arr.resize(tot_species_num);

	col3_hg_bc_points_vp_u_arr.resize(tot_species_num);
	den_arr.resize(tot_species_num);

	col3_adjustable_quad_n_arr.resize(tot_species_num);

	for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
	{
		col3_vp_n_arr.push_back(species_data->n_vp(sp_id)*col3_vdim_factor);
		col3_u_n_arr.push_back(species_data->n_u(sp_id)*col3_vdim_factor);

		sp_kinetic_op_arr.push_back(species_data->kinetic_op(sp_id));
		value_type vth_max = species_data->vmax(sp_id);

		col3_vp_min_arr.push_back(-vth_max);
		col3_vp_max_arr.push_back(vth_max);
		col3_u_min_arr.push_back(0.0);
		col3_u_max_arr.push_back(vth_max);

		vp_tot_min_index_arr[sp_id].resize(nx);
		vp_tot_max_index_arr[sp_id].resize(nx); 
		u_tot_min_index_arr[sp_id].resize(nx); 
		u_tot_max_index_arr[sp_id].resize(nx); 

		col3_dvp_arr.push_back((col3_vp_max_arr[sp_id] - col3_vp_min_arr[sp_id])/col3_vp_n_arr[sp_id]);
		col3_du_arr.push_back((col3_u_max_arr[sp_id] - col3_u_min_arr[sp_id])/col3_u_n_arr[sp_id]);
		col3_inv_dvp_arr.push_back(1.0/col3_dvp_arr[sp_id]);
		col3_inv_du_arr.push_back(1.0/col3_du_arr[sp_id]);

		col3_vp_tot_min_arr.push_back(col3_vp_min_arr[sp_id] - col3_dvp_arr[sp_id]);
		col3_u_tot_min_arr.push_back(col3_u_min_arr[sp_id] - col3_du_arr[sp_id]);

		Ms_arr.push_back(species_data->normalized_mass(sp_id));
		Zs_arr.push_back(species_data->normalized_charge(sp_id));
		Ms_ov_Mp_arr.push_back(species_data->mass_ov_mp(sp_id));
		Zs_ov_e_arr.push_back(species_data->charge_ov_e(sp_id));
		den_arr[sp_id].resize(nx, 0.0);

		int col3_vp_n = col3_vp_n_arr[sp_id];
		int col3_u_n = col3_u_n_arr[sp_id];
		col3_hg_inner_n_arr.push_back((col3_vp_n + 1)*(col3_u_n + 1));
		col3_hg_tot_n_arr.push_back((col3_vp_n + 3)*(col3_u_n + 3));

		value_type col3_vp_min = col3_vp_min_arr[sp_id];
		value_type col3_vp_max = col3_vp_max_arr[sp_id];
		value_type col3_u_min = col3_u_min_arr[sp_id];
		value_type col3_u_max = col3_u_max_arr[sp_id];

		value_type col3_dvp = col3_dvp_arr[sp_id];
		value_type col3_du = col3_du_arr[sp_id];
		value_type col3_inv_dvp = col3_inv_dvp_arr[sp_id];
		value_type col3_inv_du = col3_inv_du_arr[sp_id];

		value_type col3_vp_tot_min = col3_vp_tot_min_arr[sp_id];
		value_type col3_u_tot_min = col3_u_tot_min_arr[sp_id];

		if(sp_kinetic_op_arr[sp_id] == 1)
		{
			if(col_implicit_linear_ion_op == 1 && Zs_arr[sp_id] > 0.0) col_implicit_linear_sp_op[sp_id] = 1;
			if(col_implicit_linear_electron_op == 1 && sp_id == 0) col_implicit_linear_sp_op[sp_id] = 1;
		}

		//h boundary condition setting
		int num_of_additional_h_bc = 4;
		int col3_h_bc_n = col3_u_n*2 + col3_vp_n + 1 + num_of_additional_h_bc;
		col3_h_bc_n_arr.push_back(col3_h_bc_n);


		vector<value_type> col3_h_bc_vp_index(col3_h_bc_n), col3_h_bc_u_index(col3_h_bc_n), col3_h_bc_vp_val(col3_h_bc_n), col3_h_bc_u_val(col3_h_bc_n);
		col3_hg_bc_points_vp_u_arr[sp_id].resize(4);

		for(int i = 0; i < 2; i++)
		{
			for(int j = 0; j < col3_u_n; j++)
			{
				int k = i*col3_u_n + j;

				col3_h_bc_vp_index[k] = i*col3_vp_n;
				col3_h_bc_u_index[k] = j;
			}
		}
		for(int j = 0; j < col3_vp_n + 1; j++)
		{
			int k = 2*col3_u_n + j;

			col3_h_bc_vp_index[k] = j;
			col3_h_bc_u_index[k] = col3_u_n;
		}

		for(int k = 0; k < num_of_additional_h_bc; k++)
		{
			value_type tmp_vp_index, tmp_u_index;

			switch(k)
			{
				case 0:
					tmp_vp_index = 0.0;
					tmp_u_index = col3_u_n - 0.5;
					break;
				case 1:
					tmp_vp_index = col3_vp_n;
					tmp_u_index = col3_u_n - 0.5;
					break;
				case 2:
					tmp_vp_index = 0.5;
					tmp_u_index = col3_u_n;
					break;
				case 3:
					tmp_vp_index = col3_vp_n - 0.5;
					tmp_u_index = col3_u_n;
					break;
				default:
					cerr << "Unknown flag in h bc index : " << k << endl;
					exit(EXIT_FAILURE);
					break;
			}
			int tmp_k = col3_h_bc_n - num_of_additional_h_bc + k;
			col3_h_bc_vp_index[tmp_k] = tmp_vp_index;
			col3_h_bc_u_index[tmp_k] = tmp_u_index;
		}

		for(int k = 0; k < col3_h_bc_n; k++)
		{
			col3_h_bc_vp_val[k] = col3_vp_min + col3_h_bc_vp_index[k]*col3_dvp;
			col3_h_bc_u_val[k] = col3_u_min + col3_h_bc_u_index[k]*col3_du;
			col3_hg_bc_points_vp_u_arr[sp_id][0].push_back(col3_h_bc_vp_val[k]);
			col3_hg_bc_points_vp_u_arr[sp_id][1].push_back(col3_h_bc_u_val[k]);
		}

		//g boundary condition setting
		int num_of_additional_g_bc = num_of_additional_h_bc;
		int col3_g_bc_n = col3_h_bc_n;
		col3_g_bc_n_arr.push_back(col3_g_bc_n);

		vector<value_type> col3_g_bc_vp_index(col3_g_bc_n), col3_g_bc_u_index(col3_g_bc_n), col3_g_bc_vp_val(col3_g_bc_n), col3_g_bc_u_val(col3_g_bc_n);

		for(int i = 0; i < 2; i++)
		{
			for(int j = 0; j < col3_u_n; j++)
			{
				int k = i*col3_u_n + j;

				col3_g_bc_vp_index[k] = i*col3_vp_n;
				col3_g_bc_u_index[k] = j;
			}
		}
		for(int j = 0; j < col3_vp_n + 1; j++)
		{
			int k = 2*col3_u_n + j;

			col3_g_bc_vp_index[k] = j;
			col3_g_bc_u_index[k] = col3_u_n;
		}

		for(int k = 0; k < num_of_additional_g_bc; k++)
		{
			value_type tmp_vp_index, tmp_u_index;

			switch(k)
			{
				case 0:
					tmp_vp_index = 0.0;
					tmp_u_index = col3_u_n - 0.5;
					break;
				case 1:
					tmp_vp_index = col3_vp_n;
					tmp_u_index = col3_u_n - 0.5;
					break;
				case 2:
					tmp_vp_index = 0.5;
					tmp_u_index = col3_u_n;
					break;
				case 3:
					tmp_vp_index = col3_vp_n - 0.5;
					tmp_u_index = col3_u_n;
					break;
				default:
					cerr << "Unknown flag in g bc index : " << k << endl;
					exit(EXIT_FAILURE);
					break;
			}
			int tmp_k = col3_g_bc_n - num_of_additional_g_bc + k;
			col3_g_bc_vp_index[tmp_k] = tmp_vp_index;
			col3_g_bc_u_index[tmp_k] = tmp_u_index;
		}

		for(int k = 0; k < col3_g_bc_n; k++)
		{
			col3_g_bc_vp_val[k] = col3_vp_min + col3_g_bc_vp_index[k]*col3_dvp;
			col3_g_bc_u_val[k] = col3_u_min + col3_g_bc_u_index[k]*col3_du;
			col3_hg_bc_points_vp_u_arr[sp_id][2].push_back(col3_g_bc_vp_val[k]);
			col3_hg_bc_points_vp_u_arr[sp_id][3].push_back(col3_g_bc_u_val[k]);
		}
	}

	col3_fM_coeff_for_hg_arr.resize(tot_species_num);
	col3_h0_g0_avged_Q_arr.resize(tot_species_num); 
	col3_fM_coeff_for_hg_nx_arr.resize(tot_species_num);
	col3_nUT_arr.resize(tot_species_num);
	col3_nUT_before_arr.resize(tot_species_num);
	col3_fM_coeff_update_before_time_arr.resize(nx);
	col_fast_to_slow_moment_hg_Vec.resize(tot_species_num);
	col_f_dg_to_f_hat_mat_arr.resize(tot_species_num);
	col3_f0_h0g0_col_arr.resize(tot_species_num);
	col3_M_C_fMa_to_fMb_wUaTa_summed_arr.resize(tot_species_num);
	col3_ei_mom_en_transfer_arr.resize(tot_species_num);
	col3_ie_mom_en_transfer_arr.resize(tot_species_num);
	col3_f0_h0g0_col_flag_arr.resize(tot_species_num);

	col3_Qv_h0_ab_arr.resize(2);
	col3_Qv_g0_ab_arr.resize(2);
	for(int k1 = 0; k1 < 2; k1++)
	{
		col3_Qv_h0_ab_arr[k1].resize(nx);
		col3_Qv_g0_ab_arr[k1].resize(nx);
	}
	col3_fa_to_fMb_col_flag_arr.resize(nx);
	col3_stored_del_quantity_arr.resize(nx);
	col3_Qv_h0_ab_summed_arr.resize(nx);
	col3_Qv_g0_ab_summed_arr.resize(nx);

	col3_sp_ab_ix_col_flag_arr.resize(nx);
	consv_fa_fb0_mat_arr.resize(nx);
	consv_fa_del_hg_b_mat_arr.resize(nx);
	consv_fa_fb0_UaTa_mat_arr.resize(nx);
	col3_h_to_Qv_mat_whole.resize(nx);
	col3_g_to_Qv_mat.resize(nx);
	for(int ix = 0; ix < nx; ix++)
	{
		col3_h_to_Qv_mat_whole[ix].resize(tot_species_num);
		col3_g_to_Qv_mat[ix].resize(tot_species_num);
		for(int k1 = 0; k1 < 2; k1++)
		{
			col3_Qv_h0_ab_arr[k1][ix].resize(tot_species_num);
			col3_Qv_g0_ab_arr[k1][ix].resize(tot_species_num);
		}
		col3_fa_to_fMb_col_flag_arr[ix].resize(tot_species_num);
		col3_stored_del_quantity_arr[ix].resize(tot_species_num);
		col3_Qv_h0_ab_summed_arr[ix].resize(tot_species_num);
		col3_Qv_g0_ab_summed_arr[ix].resize(tot_species_num);
		col3_sp_ab_ix_col_flag_arr[ix].resize(tot_species_num);
		for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
			col3_sp_ab_ix_col_flag_arr[ix][sp_a_id].resize(tot_species_num, 0);
		consv_fa_fb0_mat_arr[ix].resize(tot_species_num);
		consv_fa_del_hg_b_mat_arr[ix].resize(tot_species_num);
		consv_fa_fb0_UaTa_mat_arr[ix].resize(tot_species_num);
		for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
		{
			consv_fa_fb0_mat_arr[ix][sp_a_id].resize(tot_species_num);
			consv_fa_fb0_UaTa_mat_arr[ix][sp_a_id].resize(tot_species_num, Vector::Zero(2));
		}
	}
	col3_valid_ix_flag_arr.resize(nx, 1);

	col3_gamma_ab_arr.resize(tot_species_num);
	col3_smaller_vspace_ab_arr.resize(tot_species_num);
	col3_vspace_domain_ratio_ab_arr.resize(tot_species_num);
	col3_hg_b_to_hg_a_source_arr.resize(tot_species_num);
	col3_hg_b_to_hg_a_bc_arr.resize(tot_species_num);

	col3_hg_Xmat_arr.resize(tot_species_num);
	col3_h_LF_Vmat_arr.resize(tot_species_num);
	col3_h_to_Qv_mat.resize(tot_species_num);


	col_edge_flux_qd_num.resize(tot_species_num);
	col_edge_flux_qd_points.resize(tot_species_num);

	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
	{
		col_fast_to_slow_moment_hg_Vec[sp_a_id].resize(tot_species_num);
		col_f_dg_to_f_hat_mat_arr[sp_a_id].resize(nx);
		col3_f0_h0g0_col_arr[sp_a_id].resize(tot_species_num);
		col3_M_C_fMa_to_fMb_wUaTa_summed_arr[sp_a_id].resize(nx, Vector::Zero(1));
		col3_f0_h0g0_col_flag_arr[sp_a_id].resize(nx, 1);

		col3_gamma_ab_arr[sp_a_id].resize(tot_species_num);
		col3_smaller_vspace_ab_arr[sp_a_id].resize(tot_species_num);
		col3_vspace_domain_ratio_ab_arr[sp_a_id].resize(tot_species_num);
		col3_hg_b_to_hg_a_source_arr[sp_a_id].resize(tot_species_num);
		col3_hg_b_to_hg_a_bc_arr[sp_a_id].resize(tot_species_num);

		for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
		{
			col_fast_to_slow_moment_hg_Vec[sp_a_id][sp_b_id].resize(nx, Vector::Zero(1));
			col3_f0_h0g0_col_arr[sp_a_id][sp_b_id].resize(nx, Vector::Zero(1));
			col3_hg_b_to_hg_a_source_arr[sp_a_id][sp_b_id].resize(nx);
			col3_hg_b_to_hg_a_bc_arr[sp_a_id][sp_b_id].resize(nx);
			col3_smaller_vspace_ab_arr[sp_a_id][sp_b_id].resize(nx, 0);
			col3_vspace_domain_ratio_ab_arr[sp_a_id][sp_b_id].resize(nx, 0.0);

			for(int ix = 0; ix < nx; ix++)
			{
				if(sp_b_id == 0) col3_fa_to_fMb_col_flag_arr[ix][sp_a_id].resize(tot_species_num, 0);
				if(sp_b_id == 0) col3_stored_del_quantity_arr[ix][sp_a_id].resize(tot_species_num, Vector::Zero(2));
				value_type R_cen = R_cen_cell[ix], Z_cen = Z_cen_cell[ix];

				value_type n_a = species_data->prof2D_input(R_cen, Z_cen, sp_a_id, 0);
				value_type n_b = species_data->prof2D_input(R_cen, Z_cen, sp_b_id, 0);
				value_type T_a = species_data->prof2D_input(R_cen, Z_cen, sp_a_id, 1);
				value_type T_b = species_data->prof2D_input(R_cen, Z_cen, sp_b_id, 1);

				col3_gamma_ab_arr[sp_a_id][sp_b_id].push_back(gamma_ab_gen_norm(sp_a_id, n_a, T_a, sp_b_id, n_b, T_b));

				if(sp_a_id <= sp_b_id)
				{
					int smaller_sp_id;
					if (col3_vp_max_arr[sp_a_id] > col3_vp_max_arr[sp_b_id])
						smaller_sp_id = sp_b_id;
					else smaller_sp_id = sp_a_id;
					col3_smaller_vspace_ab_arr[sp_a_id][sp_b_id][ix] = smaller_sp_id;
					col3_vspace_domain_ratio_ab_arr[sp_a_id][sp_b_id][ix] = (col3_vp_max_arr[sp_a_id]/col3_vp_max_arr[sp_b_id]);
				}
				else
				{
					col3_smaller_vspace_ab_arr[sp_a_id][sp_b_id][ix] = col3_smaller_vspace_ab_arr[sp_b_id][sp_a_id][ix];
					col3_vspace_domain_ratio_ab_arr[sp_a_id][sp_b_id][ix] = 1.0/col3_vspace_domain_ratio_ab_arr[sp_b_id][sp_a_id][ix];
				}
			}


			if (rank ==0 && sp_kinetic_op_arr[sp_a_id] == 1 && sp_kinetic_op_arr[sp_b_id] == 1)
			{
				cout << "sp_a, sp_b : " << sp_a_id << " " << sp_b_id << " gamma_ab : " << col3_gamma_ab_arr[sp_a_id][sp_b_id][0] << " smaller vth " << col3_smaller_vspace_ab_arr[sp_a_id][sp_b_id][0] << endl;

				if(sp_b_id == col3_smaller_vspace_ab_arr[sp_a_id][sp_b_id][0])
				{
					cout << "v_a_size/v_b_size, fast_to_slow_moment_ratio_max : " << col3_vspace_domain_ratio_ab_arr[sp_a_id][sp_b_id][0] << " " << col_fast_to_slow_moment_ratio_max << endl;
				}

			}
		}

		col3_ei_mom_en_transfer_arr[sp_a_id].resize(nx);
		col3_ie_mom_en_transfer_arr[sp_a_id].resize(nx);
		for(int ix = 0; ix < nx; ix++)
		{
			col3_ei_mom_en_transfer_arr[sp_a_id][ix].resize(3, Vector::Zero(2));
			col3_ie_mom_en_transfer_arr[sp_a_id][ix].resize(3, Vector::Zero(2));
		}
	}

	max_t_step_ratio.resize(tot_species_num, Vector::Zero(nx));
	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
	{
		int basis_order;
		value_type del_p, C_adv, C_dif;

		int basis_type_input =config->get_option<int>("sml_basis_type");
		switch(basis_type_input)
		{
			case 0:
				basis_order = 1;
				del_p = 2.512;
				C_adv = 1.0;
				C_dif = 0.94;

				break;

			case 2:
				basis_order = 2;
				del_p = 2.512;
				C_adv = 1.2;
				C_dif = 0.92;

				break;

			case 3:
				basis_order = 2;
				del_p = 2.512;
				C_adv = 1.2;
				C_dif = 0.92;

				break;

			case 5:
				basis_order = 2;
				del_p = 2.512;
				C_adv = 1.2;
				C_dif = 0.92;

				break;

			default:
				cout << "Unsupported type of basis for col3 : " << basis_type_input <<endl;
				abort();
				break;
		}

		value_type col3_vp_max = col3_vp_max_arr[sp_a_id];
		value_type col3_u_max = col3_u_max_arr[sp_a_id];
		value_type col3_dvp = col3_dvp_arr[sp_a_id];
		value_type col3_du = col3_du_arr[sp_a_id];

		for(int ix = 0; ix < nx; ix++)
		{
			value_type tmp_val = 0.0;
			int n_vp_loc = species_data->n_vp(sp_a_id);
			int u_loc = species_data->n_u(sp_a_id);

			value_type R_cen = R_cen_cell[ix], Z_cen = Z_cen_cell[ix];
			value_type T_a = species_data->prof2D_input(R_cen, Z_cen, sp_a_id, 1);

			value_type vt_sq = T_a/Ms_arr[sp_a_id];
			tmp_val += 2.0*C_adv*(2*basis_order + 1)*(col3_vp_max/col3_dvp + col3_u_max/col3_du);
			tmp_val += 4.0*C_dif*(basis_order + 1)*(basis_order + 1)*vt_sq*(1.0/(col3_dvp*col3_dvp) + 1.0/(col3_du*col3_du));

			max_t_step_ratio[sp_a_id][ix] = del_p/tmp_val;

			if(rank == 0 && ix == 0)
			{
				cout << "(sp_id, max_t_step_ratio (= del_t_max/tau_min) for sp_id) : " << sp_a_id << " " << max_t_step_ratio[sp_a_id][ix] << endl;
			}
		}

		tau_aa_init[sp_a_id].resize(nx);
	}

	mesh_arr.resize(tot_species_num);
	quadrature_arr.resize(tot_species_num);
	basis_arr.resize(tot_species_num);
	flux_arr.resize(tot_species_num);

	S_lp_f_index_arr.resize(4);
	for(int ej = 0; ej < 4; ej++)
	{
		if(ej == 0 || ej == 1) 
		{
			S_lp_f_index_arr[ej] = 0;
		}
		else
		{
			S_lp_f_index_arr[ej] = 1;
		}
	}

	Q_S_op_flags_arr.resize(4);
	for(int edge_flux_type = 0; edge_flux_type < 4; edge_flux_type++)
	{
		Q_S_op_flags_arr[edge_flux_type].resize(4);

		//lp : giver, l : receiver
		//case ej=0 : [lp = -, l = -] case : coeff_lp[0] used
		//case ej=1 : [lp = -, l = +] case : coeff_lp[1] used
		//case ej=2 : [lp = +, l = -] case : coeff_lp[2] used
		//case ej=3 : [lp = +, l = +] case : coeff_lp[3] used
		for(int ej = 0; ej < 4; ej++)
		{
			Q_S_op_flags_arr[edge_flux_type][ej].resize(2, 0);

			//lp : giver, l : receiver 
			if(edge_flux_type == 0) //pure Upwind from - to +
			{
				if(ej == 0 || ej == 1) 
				{
					Q_S_op_flags_arr[edge_flux_type][ej][0] = 2;
				}
			}
			else if(edge_flux_type == 1) //pure Upwind from + to -
			{
				if(ej == 2 || ej == 3) 
				{
					Q_S_op_flags_arr[edge_flux_type][ej][0] = 2;
				}
			}
			else //Mixed from - to +
			{
				Q_S_op_flags_arr[edge_flux_type][ej][0] = 1;
				Q_S_op_flags_arr[edge_flux_type][ej][1] = 1;
			}
		}
	}
}

//collision species setup
void Collision_dg::collision_dg_coeff_setup(Mesh &_mesh, Quadrature &_quadrature, Basis &_basis, Flux &_flux, int sp_id)
{

	mesh_arr[sp_id] = &_mesh;
	quadrature_arr[sp_id] = &_quadrature;
	basis_arr[sp_id] = &_basis;
	flux_arr[sp_id] = &_flux;

	Mesh *mesh_sp = mesh_arr[sp_id];
	Quadrature *quadrature_sp = quadrature_arr[sp_id];
	Basis *basis_sp = basis_arr[sp_id];
	Flux *flux_sp = flux_arr[sp_id];


	int fhat_v_basis_group_id_input = config->get_option<int>("col_fhat_v_basis_type");
	basis_arr[sp_id]->col_fhat_setup(fhat_v_basis_group_id_input);
	basis_arr[sp_id]->col_fhat_init_setup(fhat_v_basis_group_id_input);

	if(rank == 0)
	{
		cout << "sp_id, fhat_org_v_basis_dim, fhat_valid_v_basis_dim : " << sp_id << " " << basis_arr[sp_id]->fhat_org_vbasis_dim_out() << " " << basis_arr[sp_id]->fhat_valid_vbasis_dim_out() << endl;

	}

	vector<value_type> tmp_val_arr;

	int dof = dof_arr[sp_id] = basis_sp->get_dof();
	int dof_x = dof_x_arr[sp_id] = basis_sp->get_dofx();
	int dof_v = dof_v_arr[sp_id] = basis_sp->get_dofv();
	nve_arr[sp_id] = mesh_arr[sp_id]->size<EdgeV>();

	while (integration_col_arr.size() < sp_id)
	{
		integration_col_arr.push_back(new Integration_Col);
	}
	integration_col_arr.push_back(new Integration_Col);

	integration_col_arr[sp_id]->init(*mesh_sp, *basis_sp, *quadrature_sp, *flux_sp, *species_data, *eq_reader, sp_id, col_diag_quantity_num, col_fhat_vol_wid);

	if(vol_cell.size() < nx)
	{
		vol_cell.resize(nx);
		vol_B_cell.resize(nx);

		int dim_of_spatial_quadrature = quadrature_sp->size(0);

		std::fill(vol_cell.begin(), vol_cell.end(), 1e-99);
		std::fill(vol_B_cell.begin(), vol_B_cell.end(), 1e-99);

		// volume calculation
		for(int ix = 0; ix < nx; ix++)
		{
			for(int a = 0; a < dim_of_spatial_quadrature; a++)
			{
				value_type wa;
				Point2 pa, qa;

				tie(wa, pa, qa) = quadrature_sp->get_quadrature_single_not_traced_1d_qd_point(ElementX(ix), a);

				value_type loc_R = pa[0];
				Point4 p = pa & Point2({0, 0});

				Vector B(3);
				B = flux_sp->get_B(p);
				value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));

				vol_cell[ix] += wa*loc_R*2.0*M_PI/double(tor_wedge_n);
				vol_B_cell[ix] += wa*loc_R*B_val*2.0*M_PI/double(tor_wedge_n);
			}
		}
	}

	int nv_loc = nv[sp_id] = mesh_sp->size<ElementV>();
	int num_elm = nx*nv_loc;

	col3_adjustable_quad_n_arr[sp_id].resize(tot_species_num);

	int quad_n_init = config->get_option<int>("sml_qd_order_v_vol");
	for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
	{
		col3_adjustable_quad_n_arr[sp_id][sp_b_id].resize(nx);

		for(int ix = 0; ix < nx; ix++)
		{
			col3_adjustable_quad_n_arr[sp_id][sp_b_id][ix].resize(nv_loc, floor(quad_n_init*0.5)+1);
		}
	}



	col_sp_M_arr[sp_id].resize(nx, SparseMatrix(nv_loc*dof, nv_loc*dof));
	col_sp_Minv_arr[sp_id].resize(nx, SparseMatrix(nv_loc*dof, nv_loc*dof));
	consv_fa_fb0_mat_sum_arr[sp_id].resize(nx, SparseMatrix(nv_loc*dof, nv_loc*dof));

	if(col_implicit_linear_sp_op[sp_id] == 1)
	{
		col_implicit_linear_solver_arr[sp_id].resize(nx);
		for(int ix = 0; ix < nx; ix++)
		{
			col_implicit_linear_solver_arr[sp_id][ix].resize(3);

			for(int rk_sub_step = 0; rk_sub_step < 3; rk_sub_step++)
			{
				col_implicit_linear_solver_arr[sp_id][ix][rk_sub_step] = new Eigen::PardisoLU<SparseMatrix>;
			}
		}
	}


	int dim_of_velocity_quadrature = quadrature_sp->size(1);
	for(int ix = 0; ix < nx; ix++)
	{
		vp_tot_min_index_arr[sp_id][ix].resize(nv_loc, 100000);
		vp_tot_max_index_arr[sp_id][ix].resize(nv_loc, -100000);
		u_tot_min_index_arr[sp_id][ix].resize(nv_loc, 100000);
		u_tot_max_index_arr[sp_id][ix].resize(nv_loc, -100000);
		for(int iv = 0; iv < nv_loc; iv++)
		{
			for(int b = 0; b < dim_of_velocity_quadrature; b++)
			{
				value_type wb;
				Point2 pb, qb;

				tie(wb, pb, qb) = quadrature_sp->get_quadrature_single_not_traced_1d_qd_point(ElementV(iv), b);

				value_type v_para = pb[0], u = pb[1];

				int loc_vp_tot_index = int((v_para - col3_vp_tot_min_arr[sp_id])*col3_inv_dvp_arr[sp_id]);
				int loc_u_tot_index = int((u - col3_u_tot_min_arr[sp_id])*col3_inv_du_arr[sp_id]);

				int vp_tot_min_index = max(0, loc_vp_tot_index - 1);
				int vp_tot_max_index = min(col3_vp_n_arr[sp_id] + 2, loc_vp_tot_index + 2);
				int u_tot_min_index = max(0, loc_u_tot_index - 1);
				int u_tot_max_index = min(col3_u_n_arr[sp_id] + 2, loc_u_tot_index + 2);

				vp_tot_min_index_arr[sp_id][ix][iv] = min(vp_tot_min_index_arr[sp_id][ix][iv], vp_tot_min_index);
				vp_tot_max_index_arr[sp_id][ix][iv] = max(vp_tot_max_index_arr[sp_id][ix][iv], vp_tot_max_index);

				u_tot_min_index_arr[sp_id][ix][iv] = min(u_tot_min_index_arr[sp_id][ix][iv], u_tot_min_index);
				u_tot_max_index_arr[sp_id][ix][iv] = max(u_tot_max_index_arr[sp_id][ix][iv], u_tot_max_index);
			}
		}
	}

	tmp_val_arr.resize(col_diag_quantity_num*dof);


	for(int ix = 0; ix < nx; ix++)
	{
		vector<SparseMatrix_Triplet> tripletList_diag_cell;
		tripletList_diag_cell.reserve(col_diag_quantity_num*nv_loc*dof);

		col_diag_cell_mat[sp_id].push_back(SparseMatrix(col_diag_quantity_num,nv_loc*dof)); 
		for(int iv = 0; iv < nv_loc; iv++)
		{
			int k = ix*nv_loc + iv;

			integration_col_arr[sp_id]->col_cell_mat_cal(ElementX(ix), ElementV(iv), tmp_val_arr);

			for(int i = 0; i < dof; i++)
			{
				int loc_index = (k*dof + i)*col_diag_quantity_num;
				for(int j = 0; j < col_diag_quantity_num; j++)
				{
					tripletList_diag_cell.push_back(SparseMatrix_Triplet(j, iv*dof + i, tmp_val_arr[i*col_diag_quantity_num + j]/double(tor_wedge_n)));
				}
			}

		}

		col_diag_cell_mat[sp_id][ix].setFromTriplets(tripletList_diag_cell.begin(), tripletList_diag_cell.end());
	}

	int col_consv_num = 2;

	tmp_val_arr.resize(col_consv_num*dof);
	for(int ix = 0; ix < nx; ix++)
	{
		vector<SparseMatrix_Triplet> tripletList_diag_cell_consv;
		tripletList_diag_cell_consv.reserve(col_consv_num*nv_loc*dof);

		col_diag_cell_consv_mat[sp_id].push_back(SparseMatrix(col_consv_num,nv_loc*dof)); 

		for(int iv = 0; iv < nv_loc; iv++)
		{
			int k = ix*nv_loc + iv;

			integration_col_arr[sp_id]->col_cell_consv_mat_cal(ElementX(ix), ElementV(iv), col_method, tmp_val_arr, col_consv_can_ang_mom);
			for(int i = 0; i < dof; i++)
			{
				int loc_index = (k*dof + i)*col_consv_num;
				for(int j = 0; j < col_consv_num; j++)
				{

					tripletList_diag_cell_consv.push_back(SparseMatrix_Triplet(j, iv*dof + i, tmp_val_arr[i*col_consv_num + j]/double(tor_wedge_n)));
				}
			}

		}
		//Note :: col_diag_cell_consv_mat is not complete yet. M_inv should be multiplied at rhs in col_consv_mat_setup(M_A_col, sp_id, ix); 
		col_diag_cell_consv_mat[sp_id][ix].setFromTriplets(tripletList_diag_cell_consv.begin(), tripletList_diag_cell_consv.end());
	}

	auto nve = mesh_sp->size<EdgeV>();
	vspace_edge_list_arr[sp_id].resize(vspace_edge_num_of_type);
	vspace_edge_ele_arr[sp_id].resize(2*nve);
	for(int e = 0; e < nve; e++) //routine to make vertical v edge list [vspace_edge_list[sp_id]]
	{
		auto pnt_e_element = mesh_sp->get_neighborhood(EdgeV(e));
		int ele_0 = get<0>(pnt_e_element), ele_1 = get<1>(pnt_e_element);

		auto normal_vec = mesh_sp->get_normal(EdgeV(e));
		value_type hor_dist = abs(normal_vec[0]), ver_dist = abs(normal_vec[1]);

		if (ele_0 != -1 && ele_1 != -1)
		{
			auto ve_1st_pos = mesh_sp->get_midpoint(get<0>(pnt_e_element));
			auto ve_2nd_pos = mesh_sp->get_midpoint(get<1>(pnt_e_element));
			auto dist_vec = ve_2nd_pos - ve_1st_pos;
			if (hor_dist > ver_dist)
			{
				if (dist_vec[0] < 0.0)
				{
					int tmp_ele;
					tmp_ele = ele_0;
					ele_0 = ele_1;
					ele_1 = tmp_ele;
				}
				vspace_edge_list_arr[sp_id][0].push_back(e);
			}
			else
			{
				if (dist_vec[1] < 0.0)
				{
					int tmp_ele;
					tmp_ele = ele_0;
					ele_0 = ele_1;
					ele_1 = tmp_ele;
				}
				vspace_edge_list_arr[sp_id][1].push_back(e);
			}
		}
		else if (hor_dist > ver_dist)
		{
			int non_zero_ele;
			if (ele_0 != -1) non_zero_ele = ele_0;
			else non_zero_ele = ele_1;

			auto ve_1st_pos = mesh_sp->get_midpoint(ElementV(non_zero_ele));
			auto ve_2nd_pos = mesh_sp->get_midpoint(EdgeV(e));

			auto dist_vec = ve_2nd_pos - ve_1st_pos;
			if (dist_vec[0] > 0.0)
			{
				ele_0 = non_zero_ele;
				ele_1 = -1;
			}
			else
			{
				ele_0 = -1;
				ele_1 = non_zero_ele;
			}

			vspace_edge_list_arr[sp_id][2].push_back(e);
		}
		else
		{
			int non_zero_ele;
			if (ele_0 != -1) non_zero_ele = ele_0;
			else non_zero_ele = ele_1;

			auto ve_1st_pos = mesh_sp->get_midpoint(ElementV(non_zero_ele));
			auto ve_2nd_pos = mesh_sp->get_midpoint(EdgeV(e));

			auto dist_vec = ve_2nd_pos - ve_1st_pos;
			if (dist_vec[1] > 0.0)
			{
				ele_0 = non_zero_ele;
				ele_1 = -1;
			}
			else
			{
				ele_0 = -1;
				ele_1 = non_zero_ele;
			}
			vspace_edge_list_arr[sp_id][3].push_back(e);
		}

		vspace_edge_ele_arr[sp_id][2*e] = ele_0;
		vspace_edge_ele_arr[sp_id][2*e + 1] = ele_1;
	}

	int dof_v_1d = basis_sp->get_max_dof_v_1d_fhat();
	int fhat_dim_per_edge = fhat_dim_per_edge_arr[sp_id] = dof_x*3*dof_v_1d;

	vector<value_type> tmp_fhat_mat_arr(fhat_dim_per_edge*2*dof);
	vector<int> tmp_fhat_mat_info(fhat_dim_per_edge*2*dof);

	//org f_dg to fhat_matrix construction part
	for(int ix = 0; ix < nx; ix++)
	{
		vector<SparseMatrix_Triplet> tripletList_fhat_cell;
		tripletList_fhat_cell.reserve(nve*fhat_dim_per_edge*2*dof);

		col_fhat_cell_mat[sp_id].push_back(SparseMatrix(nve*fhat_dim_per_edge,nv_loc*dof)); 

		for(int i = 0; i < vspace_edge_num_of_type; i++)
		{
			int n_size = vspace_edge_list_arr[sp_id][i].size();
			for(int j = 0; j < n_size; j++)
			{
				int e = vspace_edge_list_arr[sp_id][i][j];
				int ele[2] = {vspace_edge_ele_arr[sp_id][2*e], vspace_edge_ele_arr[sp_id][2*e + 1]};


				tmp_fhat_mat_arr.assign(tmp_fhat_mat_arr.size(), 0.0);
				memset(&tmp_fhat_mat_info[0], -1, tmp_fhat_mat_info.size() * sizeof tmp_fhat_mat_info[0]);
				if (i == 0 || i == 1)
				{
					integration_col_arr[sp_id]->col_fhat_mat_cal(ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), i, tmp_fhat_mat_arr, tmp_fhat_mat_info);
				}
				else
				{
					integration_col_arr[sp_id]->col_fhat_mat_cal_bd(ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), i, tmp_fhat_mat_arr, tmp_fhat_mat_info);

				}

				//form tripletList_fhat_cell from tmp_fhat_mat_arr
				for(int ele_index = 0; ele_index < 2; ele_index++)
				{
					if(ele[ele_index] > -1)
					{
						for (int basisX_i = 0; basisX_i < dof_x; basisX_i++)
						{
							for(int org_basis_v_index = 0; org_basis_v_index < dof_v; org_basis_v_index++)
							{
								//int org_basis_index = basisX_i*dof_v + org_basis_v_index;
								int org_basis_index = org_basis_v_index*dof_x + basisX_i;
								int col_pos = ele_index*dof + org_basis_index;

								for(int n = 0; n < 3; n++)
								{
									for(int v_1d_i = 0; v_1d_i < dof_v_1d; v_1d_i++)
									{
										int row_pos = basisX_i*3*dof_v_1d + n*dof_v_1d + v_1d_i;
										int pos_1d = col_pos*fhat_dim_per_edge + row_pos;

										int loc_fhat_mat_info = tmp_fhat_mat_info[pos_1d];
										if (loc_fhat_mat_info > -1)
										{
											tripletList_fhat_cell.push_back(SparseMatrix_Triplet(e*fhat_dim_per_edge + row_pos, ele[ele_index]*dof + org_basis_index, tmp_fhat_mat_arr[pos_1d]));

										}
									}
								}
							}
						}
					}
				}
			}
		}
		col_fhat_cell_mat[sp_id][ix].setFromTriplets(tripletList_fhat_cell.begin(), tripletList_fhat_cell.end());
	}

	//org f_dg to fhat_matrix construction part
	int dim_fhat_valid_v_basis = basis_sp->fhat_valid_vbasis_dim_out();
	for(int ix = 0; ix < nx; ix++)
	{
		vector<SparseMatrix_Triplet> tripletList_fhat_cell;
		tripletList_fhat_cell.reserve(nve*dof_x*dim_fhat_valid_v_basis*2*dof);

		col_f_dg_to_valid_fhat_mat_arr[sp_id].push_back(SparseMatrix(nve*dof_x*dim_fhat_valid_v_basis, nv_loc*dof)); 

		col_f_dg_to_f_hat_mat_arr[sp_id][ix].resize(nve, Matrix::Zero(dof_x*dim_fhat_valid_v_basis, 2*dof));

		for(int i = 0; i < vspace_edge_num_of_type; i++)
		{
			int n_size = vspace_edge_list_arr[sp_id][i].size();
			for(int j = 0; j < n_size; j++)
			{
				//e : edge id
				int e = vspace_edge_list_arr[sp_id][i][j];
				//ele[2] : vspace elements adjacent to edge e
				int ele[2] = {vspace_edge_ele_arr[sp_id][2*e], vspace_edge_ele_arr[sp_id][2*e + 1]};

				//when e is one of inner edges
				if (i == 0 || i == 1)
				{
					Vector tmp_f_dg_to_fhat_mat_arr = Vector::Zero(dof_x*dim_fhat_valid_v_basis*2*dof);

					integration_col_arr[sp_id]->col_f_dg_to_valid_fhat_mat_cal(ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), i, tmp_f_dg_to_fhat_mat_arr);


					for (int ele_index = 0; ele_index < 2; ele_index++)
					{
						for(int dof_v_i = 0; dof_v_i < dof_v; dof_v_i++) 
						{
							for (int dof_x_j = 0; dof_x_j < dof_x; dof_x_j++)
							{
								int f_dg_basis_index2 = (ele[ele_index]*dof_v + dof_v_i)*dof_x + dof_x_j;
								int f_dg_basis_index = (ele_index*dof_v + dof_v_i)*dof_x + dof_x_j;

								for(int row = 0; row < dim_fhat_valid_v_basis; row++) 
								{

									int f_hat_row_index = (e*dof_x + dof_x_j)*dim_fhat_valid_v_basis + row;

									int tmp_index = (dof_x_j*2*dof + f_dg_basis_index)*dim_fhat_valid_v_basis + row;

									tripletList_fhat_cell.push_back(SparseMatrix_Triplet(f_hat_row_index, f_dg_basis_index2, tmp_f_dg_to_fhat_mat_arr[tmp_index]));

									col_f_dg_to_f_hat_mat_arr[sp_id][ix][e](dof_x_j*dim_fhat_valid_v_basis + row, f_dg_basis_index) += tmp_f_dg_to_fhat_mat_arr[tmp_index];

								}
							}
						}
					}
				}
				else //when e is a boundary edge
				{
					int valid_ele_i;
					if(ele[0] == -1)
						valid_ele_i = 1;
					else if(ele[1] == -1)
						valid_ele_i = 0;
					else
					{
						cout << "Invalid ele for bd type 2 or 3 : " << ele[0] << " " << ele[1] << " " << i << endl;
						abort();
					}


					//In this case, fhat is basically same to f_dg but the storage order is changed from [v_i*dof_x + x_j] to [x_j*dof_v + v_i]
					for(int dof_v_i = 0; dof_v_i < dof_v; dof_v_i++) 
					{
						for (int dof_x_j = 0; dof_x_j < dof_x; dof_x_j++)
						{
							int f_hat_row_index = (e*dof_x + dof_x_j)*dim_fhat_valid_v_basis + dof_v_i;
							int f_dg_basis_index2 = (ele[valid_ele_i]*dof_v + dof_v_i)*dof_x + dof_x_j;
							int f_dg_basis_index = (valid_ele_i*dof_v + dof_v_i)*dof_x + dof_x_j;
							tripletList_fhat_cell.push_back(SparseMatrix_Triplet(f_hat_row_index, f_dg_basis_index2, 1.0));

							col_f_dg_to_f_hat_mat_arr[sp_id][ix][e](dof_x_j*dim_fhat_valid_v_basis + dof_v_i, f_dg_basis_index) = 1.0;

						}
					}
				}
			}
		}

		col_f_dg_to_valid_fhat_mat_arr[sp_id][ix].setFromTriplets(tripletList_fhat_cell.begin(), tripletList_fhat_cell.end());
	}

	integration_col_arr[sp_id]->col3_init(tot_species_num, sp_id, col3_vp_n_arr, col3_u_n_arr, col3_vp_min_arr, col3_vp_max_arr, col3_u_min_arr, col3_u_max_arr, col_ei_pitch_angle_op, col_ei_pitch_no_v_in_nu_ei_op, col_ei_pitch_no_v_in_nu_ei_v_e_min);

	for(int ix = 0; ix < nx; ix++)
	{
		col3_fM_coeff_for_hg_nx_arr[sp_id].push_back(Vector::Zero(dof*nv_loc));
		col3_nUT_arr[sp_id].push_back(Vector::Zero(3));
		col3_nUT_before_arr[sp_id].push_back(Vector::Zero(3));
	}

	int col3_vp_n = col3_vp_n_arr[sp_id];
	int col3_u_n = col3_u_n_arr[sp_id];

	value_type col3_vp_min = col3_vp_min_arr[sp_id];
	value_type col3_vp_max = col3_vp_max_arr[sp_id];
	value_type col3_u_min = col3_u_min_arr[sp_id];
	value_type col3_u_max = col3_u_max_arr[sp_id];

	value_type col3_dvp = col3_dvp_arr[sp_id];
	value_type col3_du = col3_du_arr[sp_id];
	value_type col3_inv_dvp = col3_inv_dvp_arr[sp_id];
	value_type col3_inv_du = col3_inv_du_arr[sp_id];

	value_type col3_vp_tot_min = col3_vp_tot_min_arr[sp_id];
	value_type col3_u_tot_min = col3_u_tot_min_arr[sp_id];


	//hg boundary points data load
	int col3_h_bc_n = col3_h_bc_n_arr[sp_id];
	int col3_g_bc_n = col3_g_bc_n_arr[sp_id];
	vector<value_type> col3_h_bc_vp_val(col3_h_bc_n), col3_h_bc_u_val(col3_h_bc_n);
	for(int k = 0; k < col3_h_bc_n; k++)
	{
		col3_h_bc_vp_val[k] = col3_hg_bc_points_vp_u_arr[sp_id][0][k];
		col3_h_bc_u_val[k] = col3_hg_bc_points_vp_u_arr[sp_id][1][k];
	}

	vector<value_type> col3_g_bc_vp_val(col3_g_bc_n), col3_g_bc_u_val(col3_g_bc_n);
	for(int k = 0; k < col3_g_bc_n; k++)
	{
		col3_g_bc_vp_val[k] = col3_hg_bc_points_vp_u_arr[sp_id][2][k];
		col3_g_bc_u_val[k] = col3_hg_bc_points_vp_u_arr[sp_id][3][k];
	}

	//f to hg bc matrix setup
	for(int ix = 0; ix < nx; ix++)
	{
		vector<SparseMatrix_Triplet> tripletList_f_to_h;
		tripletList_f_to_h.reserve(col3_h_bc_n*nv_loc*dof);
		vector<SparseMatrix_Triplet> tripletList_f_to_g;
		tripletList_f_to_g.reserve(col3_g_bc_n*nv_loc*dof);

		col3_f_to_h_bc_mat[sp_id].push_back(SparseMatrix(col3_h_bc_n, nv_loc*dof)); 
		col3_f_to_g_bc_mat[sp_id].push_back(SparseMatrix(col3_g_bc_n, nv_loc*dof)); 

		for(int iv = 0; iv < nv_loc; iv++)
		{
			vector<value_type> tmp_h_bc_arr(col3_h_bc_n*dof); 
			integration_col_arr[sp_id]->col3_bc_mat_cal(ElementX(ix), ElementV(iv), sp_id, col3_h_bc_n, col3_h_bc_vp_val, col3_h_bc_u_val, 0, tmp_h_bc_arr);  

			for(int m = 0; m < col3_h_bc_n; m++)
			{
				for(int i = 0; i < dof; i++)
				{
					int mi_index = m*dof + i;
					tripletList_f_to_h.push_back(SparseMatrix_Triplet(m, iv*dof + i, tmp_h_bc_arr[mi_index]/vol_cell[ix]));
				}
			}

			vector<value_type> tmp_g_bc_arr(col3_g_bc_n*dof); 
			integration_col_arr[sp_id]->col3_bc_mat_cal(ElementX(ix), ElementV(iv), sp_id, col3_g_bc_n, col3_g_bc_vp_val, col3_g_bc_u_val, 1, tmp_g_bc_arr);  

			for(int m = 0; m < col3_g_bc_n; m++)
			{
				for(int i = 0; i < dof; i++)
				{
					int mi_index = m*dof + i;
					tripletList_f_to_g.push_back(SparseMatrix_Triplet(m, iv*dof + i, tmp_g_bc_arr[mi_index]/vol_cell[ix]));
				}
			}
		}
		col3_f_to_h_bc_mat[sp_id][ix].setFromTriplets(tripletList_f_to_h.begin(), tripletList_f_to_h.end());
		col3_f_to_g_bc_mat[sp_id][ix].setFromTriplets(tripletList_f_to_g.begin(), tripletList_f_to_g.end());
	}

	int col3_hg_inner_n = col3_hg_inner_n_arr[sp_id];
	int col3_hg_tot_n = col3_hg_tot_n_arr[sp_id];

	tmp_val_arr.resize(col3_hg_inner_n*dof);

	col3_h_solver_stiffness[sp_id].resize(nx);
	col3_g_solver_stiffness[sp_id].resize(nx);
	col3_hg_formula_to_hg0_solver[sp_id].resize(nx);

	for(int ix = 0; ix < nx; ix++)
	{
		vector<SparseMatrix_Triplet> tripletList_f_to_h_source, tripletList_h_to_g_source;
		tripletList_f_to_h_source.reserve(col3_vdim_factor*100*nv_loc*dof);
		tripletList_h_to_g_source.reserve(col3_hg_tot_n*100);
		col3_f_to_h_source_mat[sp_id].push_back(SparseMatrix(col3_hg_inner_n, nv_loc*dof)); 
		col3_h_to_g_source_mat[sp_id].push_back(SparseMatrix(col3_hg_inner_n, col3_hg_tot_n)); 

		//stiffness matrix -> total # of row = (n_vp + 3)*(n_u + 3)
		vector<SparseMatrix_Triplet> tripletList_h_stiffness, tripletList_g_stiffness, tripletList_hg_formula_to_hg0;
		tripletList_h_stiffness.reserve(col3_hg_tot_n*100);
		tripletList_g_stiffness.reserve(col3_hg_tot_n*100);
		tripletList_hg_formula_to_hg0.reserve(col3_hg_tot_n*100);


		int base_h_row_index = 0;
		int base_g_row_index = 0;
		//part 1 : Neumann b.c. at u = 0 -> # of row = n_vp + 3
		for(int k = 0; k < col3_vp_n + 3; k++)
		{
			int k2 = k + 2*(col3_vp_n + 3);
			tripletList_h_stiffness.push_back(SparseMatrix_Triplet(base_h_row_index + k, k, 1.0));
			tripletList_h_stiffness.push_back(SparseMatrix_Triplet(base_h_row_index + k, k2, -1.0));
			tripletList_g_stiffness.push_back(SparseMatrix_Triplet(base_g_row_index + k, k, 1.0));
			tripletList_g_stiffness.push_back(SparseMatrix_Triplet(base_g_row_index + k, k2, -1.0));
			tripletList_hg_formula_to_hg0.push_back(SparseMatrix_Triplet(base_h_row_index + k, k, 1.0));
			tripletList_hg_formula_to_hg0.push_back(SparseMatrix_Triplet(base_h_row_index + k, k2, -1.0));
		}
		base_h_row_index += col3_vp_n + 3;
		base_g_row_index += col3_vp_n + 3;

		//part 2 : FEM weak form integral part
		for(int iv = 0; iv < nv_loc; iv++)
		{
			vector<value_type> tmp_val_arr2(col3_hg_inner_n*col3_hg_tot_n*3);
			vector<value_type> tmp_val_arr3(dof*col3_hg_tot_n*4);
			integration_col_arr[sp_id]->col3_hg_vol_mat_cal(ElementX(ix), ElementV(iv), sp_id, tmp_val_arr, tmp_val_arr2, tmp_val_arr3);  

			for(int m = 0; m < col3_hg_inner_n; m++)
			{
				//part 2A : DG f to h volume source
				for(int i = 0; i < dof; i++)
				{
					int mi_index = m*dof + i;
					if (abs(tmp_val_arr[mi_index]) > 1e-99)
					{
						tripletList_f_to_h_source.push_back(SparseMatrix_Triplet(m, iv*dof + i, tmp_val_arr[mi_index]));
					}
				}

				for(int l = 0; l < col3_hg_tot_n; l++)
				{
					int ml_index1 = m*col3_hg_tot_n + l;
					int ml_index2 = ml_index1 + col3_hg_inner_n*col3_hg_tot_n;
					int ml_index3 = ml_index2 + col3_hg_inner_n*col3_hg_tot_n;

					//part 2B : h to g volume source
					if (abs(tmp_val_arr2[ml_index1]) > 1e-99)
					{
						//H to G source calculation
						tripletList_h_to_g_source.push_back(SparseMatrix_Triplet(m, l, -tmp_val_arr2[ml_index1]));
						//h_M, g_M formula to h0, g0 matrix solver
						tripletList_hg_formula_to_hg0.push_back(SparseMatrix_Triplet(base_h_row_index + m, l, tmp_val_arr2[ml_index1]));
					}

					//part 2C : FEM weak form integral part -> # of row = (n_vp + 1)*(n_u + 1) : d basis/dv_p dot d basis /dv_p
					if (abs(tmp_val_arr2[ml_index2]) > 1e-99)
					{
						tripletList_h_stiffness.push_back(SparseMatrix_Triplet(base_h_row_index + m, l, tmp_val_arr2[ml_index2]));
						tripletList_g_stiffness.push_back(SparseMatrix_Triplet(base_g_row_index + m, l, tmp_val_arr2[ml_index2]));
					}

					//part 2D : FEM weak form integral part -> # of row = (n_vp + 1)*(n_u + 1) : d basis/du dot d basis /du
					if (abs(tmp_val_arr2[ml_index3]) > 1e-99)
					{
						tripletList_h_stiffness.push_back(SparseMatrix_Triplet(base_h_row_index + m, l, tmp_val_arr2[ml_index3]));
						tripletList_g_stiffness.push_back(SparseMatrix_Triplet(base_g_row_index + m, l, tmp_val_arr2[ml_index3]));
					}
				}
			}
		}
		col3_f_to_h_source_mat[sp_id][ix].setFromTriplets(tripletList_f_to_h_source.begin(), tripletList_f_to_h_source.end());
		col3_h_to_g_source_mat[sp_id][ix].setFromTriplets(tripletList_h_to_g_source.begin(), tripletList_h_to_g_source.end());

		//part 2E : FEM weak form integral part -> surface terms part
		for(int i = 2; i < vspace_edge_num_of_type; i++)
		{
			int n_size = vspace_edge_list_arr[sp_id][i].size();
			for(int j = 0; j < n_size; j++)
			{
				int e = vspace_edge_list_arr[sp_id][i][j];
				vector<value_type> tmp_val_arr3(col3_hg_inner_n*col3_hg_tot_n);
				integration_col_arr[sp_id]->col3_edge_bd_hg_stiffness(ElementX(ix), EdgeV(e), sp_id, i, tmp_val_arr3);

				for(int m = 0; m < col3_hg_inner_n; m++)
				{
					for(int l = 0; l < col3_hg_tot_n; l++)
					{
						int ml_index = m*col3_hg_tot_n + l;

						//part 2 : FEM weak form integral part -> # of row = (n_vp + 1)*(n_u + 1)
						if (abs(tmp_val_arr3[ml_index]) > 1e-99)
						{
							tripletList_h_stiffness.push_back(SparseMatrix_Triplet(base_h_row_index + m, l, tmp_val_arr3[ml_index]));
							tripletList_g_stiffness.push_back(SparseMatrix_Triplet(base_g_row_index + m, l, tmp_val_arr3[ml_index]));

						}
					}
				}
			}
		}
		base_h_row_index += col3_hg_inner_n;
		base_g_row_index += col3_hg_inner_n;


		//part 3 : Dirichlet b.c. part -> # of row = 2*n_u + n_vp + 1 + alpha(= 4 if there is no other constraints)
		for(int k = 0; k < col3_h_bc_n; k++)
		{
			value_type v_para = col3_h_bc_vp_val[k];
			value_type u = col3_h_bc_u_val[k];

			int loc_vp_tot_index = int((v_para - col3_vp_tot_min)*col3_inv_dvp);
			int loc_u_tot_index = int((u - col3_u_tot_min)*col3_inv_du);
			int vp_tot_min_index = max(0, loc_vp_tot_index - 1);
			int vp_tot_max_index = min(col3_vp_n + 2, loc_vp_tot_index + 2);
			int u_tot_min_index = max(0, loc_u_tot_index - 1);
			int u_tot_max_index = min(col3_u_n + 2, loc_u_tot_index + 2);
			for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
			{
				value_type t1_tot = (v_para - (col3_vp_tot_min + j3*col3_dvp))*col3_inv_dvp;
				Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_inv_dvp);

				for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
				{
					value_type t2_tot = (u - (col3_u_tot_min + j4*col3_du))*col3_inv_du;
					Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_inv_du);

					value_type spline_tot_val = spline_vp_tot_vec(0)*spline_u_tot_vec(0);

					if (abs(spline_tot_val) > 1e-99)
					{
						int l = j4*(col3_vp_n + 3) + j3;
						tripletList_h_stiffness.push_back(SparseMatrix_Triplet(base_h_row_index + k, l, spline_tot_val));

						tripletList_hg_formula_to_hg0.push_back(SparseMatrix_Triplet(base_h_row_index + k, l, spline_tot_val));
					}
				}
			}			
		}
		base_h_row_index += col3_h_bc_n;

		for(int k = 0; k < col3_g_bc_n; k++)
		{
			value_type v_para = col3_g_bc_vp_val[k];
			value_type u = col3_g_bc_u_val[k];

			int loc_vp_tot_index = int((v_para - col3_vp_tot_min)*col3_inv_dvp);
			int loc_u_tot_index = int((u - col3_u_tot_min)*col3_inv_du);

			int vp_tot_min_index = max(0, loc_vp_tot_index - 1);
			int vp_tot_max_index = min(col3_vp_n + 2, loc_vp_tot_index + 2);
			int u_tot_min_index = max(0, loc_u_tot_index - 1);
			int u_tot_max_index = min(col3_u_n + 2, loc_u_tot_index + 2);

			for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
			{
				value_type t1_tot = (v_para - (col3_vp_tot_min + j3*col3_dvp))*col3_inv_dvp;
				Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_inv_dvp);

				for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
				{
					value_type t2_tot = (u - (col3_u_tot_min + j4*col3_du))*col3_inv_du;
					Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_inv_du);

					value_type spline_tot_val = spline_vp_tot_vec(0)*spline_u_tot_vec(0);

					if (abs(spline_tot_val) > 1e-99)
					{
						int l = j4*(col3_vp_n + 3) + j3;
						tripletList_g_stiffness.push_back(SparseMatrix_Triplet(base_g_row_index + k, l, spline_tot_val));
					}
				}
			}			
		}
		base_g_row_index += col3_g_bc_n;

		SparseMatrix H_stiffness_tmp(col3_hg_tot_n, col3_hg_tot_n);
		H_stiffness_tmp.setFromTriplets(tripletList_h_stiffness.begin(), tripletList_h_stiffness.end());

		col3_h_solver_stiffness[sp_id][ix] = new Eigen::PardisoLU<SparseMatrix>;
		col3_h_solver_stiffness[sp_id][ix]->analyzePattern(H_stiffness_tmp);
		col3_h_solver_stiffness[sp_id][ix]->factorize(H_stiffness_tmp);

		col3_g_stiffness_mat[sp_id].push_back(SparseMatrix(col3_hg_tot_n, col3_hg_tot_n));
		col3_g_stiffness_mat[sp_id][ix].setFromTriplets(tripletList_g_stiffness.begin(), tripletList_g_stiffness.end());
		col3_g_solver_stiffness[sp_id][ix] = new Eigen::PardisoLU<SparseMatrix>;
		col3_g_solver_stiffness[sp_id][ix]->analyzePattern(col3_g_stiffness_mat[sp_id][ix]);
		col3_g_solver_stiffness[sp_id][ix]->factorize(col3_g_stiffness_mat[sp_id][ix]);

		SparseMatrix hg_formula_to_hg0_tmp(col3_hg_tot_n, col3_hg_tot_n);
		hg_formula_to_hg0_tmp.setFromTriplets(tripletList_hg_formula_to_hg0.begin(), tripletList_hg_formula_to_hg0.end());

		col3_hg_formula_to_hg0_solver[sp_id][ix] = new Eigen::PardisoLU<SparseMatrix>;
		col3_hg_formula_to_hg0_solver[sp_id][ix]->analyzePattern(hg_formula_to_hg0_tmp);
		col3_hg_formula_to_hg0_solver[sp_id][ix]->factorize(hg_formula_to_hg0_tmp);

		//hg_b to hg_a source term
		for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
		{
			int sp_a_id = sp_id;
			int smaller_sp_id = col3_smaller_vspace_ab_arr[sp_a_id][sp_b_id][ix];
			// big h_b to small a only
			if (sp_b_id != sp_a_id && sp_kinetic_op_arr[sp_b_id] == 1 && sp_a_id == smaller_sp_id)
			{
				int hg_a_inner_n = col3_hg_inner_n_arr[sp_a_id];
				int hg_a_total_n = col3_hg_tot_n_arr[sp_a_id];
				int hg_b_inner_n = col3_hg_inner_n_arr[sp_b_id];
				int hg_b_total_n = col3_hg_tot_n_arr[sp_b_id];

				col3_hg_b_to_hg_a_source_arr[sp_a_id][sp_b_id][ix] = SparseMatrix(hg_a_inner_n, hg_b_total_n);
				//col3_hg_b_to_hg_a_source_arr[sp_b_id][sp_a_id][ix] = SparseMatrix(hg_b_inner_n, hg_a_total_n);

				vector<SparseMatrix_Triplet> tripletList_hg_b_to_hg_a_source1, tripletList_hg_b_to_hg_a_source2;
				vector<value_type> tmp_val_arr_bigger_to_smaller(hg_a_inner_n*hg_b_total_n);
				vector<value_type> tmp_val_arr_smaller_to_bigger(hg_b_inner_n*hg_a_total_n);
				for(int iv = 0; iv < nv_loc; iv++)
				{
					integration_col_arr[sp_id]->col3_hg_b_to_hg_a_vol_mat_cal(ElementX(ix), ElementV(iv), sp_a_id, sp_b_id, tmp_val_arr_bigger_to_smaller, tmp_val_arr_smaller_to_bigger);

					for(int m = 0; m < hg_a_inner_n; m++)
					{
						for(int l = 0; l < hg_b_total_n; l++)
						{
							int ml_index = m*hg_b_total_n + l;
							if (abs(tmp_val_arr_bigger_to_smaller[ml_index]) > 1e-99)
							{
								tripletList_hg_b_to_hg_a_source1.push_back(SparseMatrix_Triplet(m, l, tmp_val_arr_bigger_to_smaller[ml_index]));
							}
						}
					}
				}

				col3_hg_b_to_hg_a_source_arr[sp_a_id][sp_b_id][ix].setFromTriplets(tripletList_hg_b_to_hg_a_source1.begin(), tripletList_hg_b_to_hg_a_source1.end());

			}

			// small h_b to big a only
			if (sp_b_id != sp_a_id && sp_kinetic_op_arr[sp_b_id] == 1 && sp_b_id == smaller_sp_id)
			{
				int hg_a_inner_n = col3_hg_inner_n_arr[sp_a_id];
				int hg_b_total_n = col3_hg_tot_n_arr[sp_b_id];

				col3_hg_b_to_hg_a_source_arr[sp_a_id][sp_b_id][ix] = SparseMatrix(hg_a_inner_n, hg_b_total_n);

				vector<SparseMatrix_Triplet> tripletList_hg_b_to_hg_a_source;
				vector<value_type> tmp_val_arr_smaller_to_bigger(hg_a_inner_n*hg_b_total_n);
				for(int iv = 0; iv < nv_loc; iv++)
				{
					integration_col_arr[sp_id]->col3_small_hg_b_to_hg_a_vol_mat_cal(ElementX(ix), ElementV(iv), sp_a_id, sp_b_id, tmp_val_arr_smaller_to_bigger);

					for(int m = 0; m < hg_a_inner_n; m++)
					{
						for(int l = 0; l < hg_b_total_n; l++)
						{
							int ml_index = m*hg_b_total_n + l;
							if (abs(tmp_val_arr_smaller_to_bigger[ml_index]) > 1e-99)
							{
								tripletList_hg_b_to_hg_a_source.push_back(SparseMatrix_Triplet(m, l, tmp_val_arr_smaller_to_bigger[ml_index]));
							}
						}
					}
				}

				col3_hg_b_to_hg_a_source_arr[sp_a_id][sp_b_id][ix].setFromTriplets(tripletList_hg_b_to_hg_a_source.begin(), tripletList_hg_b_to_hg_a_source.end());
			}
		}

		//hg_b to hg_a b.c. term

		for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
		{
			int sp_b_id = sp_id;
			if (sp_a_id != sp_b_id && sp_kinetic_op_arr[sp_a_id] == 1 && sp_kinetic_op_arr[sp_b_id] == 1)
			{
				int hg_b_total_n = col3_hg_tot_n_arr[sp_b_id];
				int num_of_bc_a = col3_hg_bc_points_vp_u_arr[sp_a_id][0].size(); 
				int smaller_sp_id = col3_smaller_vspace_ab_arr[sp_a_id][sp_b_id][ix];

				if (sp_a_id == smaller_sp_id)
				{
					if(sp_a_id == 0)
					{
						cout << "electron can not be slower species : " << sp_a_id << " " << sp_b_id << endl;
						abort();
					}

					col3_hg_b_to_hg_a_bc_arr[sp_a_id][sp_b_id][ix] = SparseMatrix(num_of_bc_a, hg_b_total_n);
					vector<SparseMatrix_Triplet> tripletList_hg_b_to_hg_a_bc_source;

					for(int k = 0; k < num_of_bc_a; k++)
					{
						value_type v_para = col3_hg_bc_points_vp_u_arr[sp_a_id][0][k];
						value_type u = col3_hg_bc_points_vp_u_arr[sp_a_id][1][k];

						int loc_vp_tot_index = int((v_para - col3_vp_tot_min)*col3_inv_dvp);
						int loc_u_tot_index = int((u - col3_u_tot_min)*col3_inv_du);
						int vp_tot_min_index = max(0, loc_vp_tot_index - 1);
						int vp_tot_max_index = min(col3_vp_n + 2, loc_vp_tot_index + 2);
						int u_tot_min_index = max(0, loc_u_tot_index - 1);
						int u_tot_max_index = min(col3_u_n + 2, loc_u_tot_index + 2);
						for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
						{
							value_type t1_tot = (v_para - (col3_vp_tot_min + j3*col3_dvp))*col3_inv_dvp;
							Vector spline_vp_tot_vec = Spline::cal_Cs_Vec(t1_tot, col3_inv_dvp);

							for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
							{
								value_type t2_tot = (u - (col3_u_tot_min + j4*col3_du))*col3_inv_du;
								Vector spline_u_tot_vec = Spline::cal_Cs_Vec(t2_tot, col3_inv_du);

								value_type spline_tot_val = spline_vp_tot_vec(0)*spline_u_tot_vec(0);

								if (abs(spline_tot_val) > 1e-99)
								{
									int l = j4*(col3_vp_n + 3) + j3;
									tripletList_hg_b_to_hg_a_bc_source.push_back(SparseMatrix_Triplet(k, l, spline_tot_val));
								}
							}
						}			
					}

					col3_hg_b_to_hg_a_bc_arr[sp_a_id][sp_b_id][ix].setFromTriplets(tripletList_hg_b_to_hg_a_bc_source.begin(), tripletList_hg_b_to_hg_a_bc_source.end());


				}
				else
				{
					int num_of_nv_b_dof_b = dof*nv_loc;
					int hg_a_inner_n = col3_hg_inner_n_arr[sp_a_id];

					col3_hg_b_to_hg_a_bc_arr[sp_a_id][sp_b_id][ix] = SparseMatrix(2*num_of_bc_a + hg_a_inner_n, num_of_nv_b_dof_b);

					vector<SparseMatrix_Triplet> tripletList_f_to_hg;
					tripletList_f_to_hg.reserve(2*num_of_bc_a*num_of_nv_b_dof_b);

					for(int iv = 0; iv < nv_loc; iv++)
					{
						vector<value_type> tmp_h_bc_arr(num_of_bc_a*dof); 
						integration_col_arr[sp_b_id]->col3_bc_mat_cal(ElementX(ix), ElementV(iv), sp_b_id, num_of_bc_a, col3_hg_bc_points_vp_u_arr[sp_a_id][0], col3_hg_bc_points_vp_u_arr[sp_a_id][1], 0, tmp_h_bc_arr);  

						for(int m = 0; m < num_of_bc_a; m++)
						{
							for(int i = 0; i < dof; i++)
							{
								int mi_index = m*dof + i;
								tripletList_f_to_hg.push_back(SparseMatrix_Triplet(m, iv*dof + i, tmp_h_bc_arr[mi_index]/vol_cell[ix]));
							}
						}

						vector<value_type> tmp_g_bc_arr(num_of_bc_a*dof); 
						integration_col_arr[sp_b_id]->col3_bc_mat_cal(ElementX(ix), ElementV(iv), sp_b_id, num_of_bc_a, col3_hg_bc_points_vp_u_arr[sp_a_id][0], col3_hg_bc_points_vp_u_arr[sp_a_id][1], 1, tmp_g_bc_arr);  

						for(int m = 0; m < num_of_bc_a; m++)
						{
							for(int i = 0; i < dof; i++)
							{
								int mi_index = m*dof + i;
								tripletList_f_to_hg.push_back(SparseMatrix_Triplet(m + num_of_bc_a, iv*dof + i, tmp_g_bc_arr[mi_index]/vol_cell[ix]));
							}
						}

						vector<value_type> tmp_fb_to_ha_source_arr(hg_a_inner_n*dof); 
						integration_col_arr[sp_b_id]->col3_small_fb_to_big_h_a_source_mat_cal(ElementX(ix), ElementV(iv), sp_a_id, sp_b_id, tmp_fb_to_ha_source_arr);

						for(int m = 0; m < hg_a_inner_n; m++)
						{
							for(int i = 0; i < dof; i++)
							{
								int mi_index = m*dof + i;
								if (abs(tmp_fb_to_ha_source_arr[mi_index]) > 1e-99)
								{
									tripletList_f_to_hg.push_back(SparseMatrix_Triplet(m + 2*num_of_bc_a, iv*dof + i, tmp_fb_to_ha_source_arr[mi_index]));
								}
							}
						}
					}
					col3_hg_b_to_hg_a_bc_arr[sp_a_id][sp_b_id][ix].setFromTriplets(tripletList_f_to_hg.begin(), tripletList_f_to_hg.end());


				}
			}

		}

	}

	//volume term part
	int tot_Qv_h_size_per_ix = nv_loc*2*dof_v*dof_v;
	int tot_Qv_h_size_per_ix_vol = nv_loc*2*dof_v*dof_v;
	int tot_Qv_h_size_per_ix_tmp = nv_loc*2*dof_v*dof_v;
	//surface term part
	tot_Qv_h_size_per_ix += nve*4*dof_v*dof_v;
	tot_Qv_h_size_per_ix_tmp += nve*4*dof_v*dof_v;

	//flux from h at edge qd points
	//col_edge_flux_qd_num[sp_id] = 5;
	col_edge_flux_qd_num[sp_id] = quadrature_sp->size(3);
	col_edge_flux_qd_points[sp_id] = Vector::Zero(col_edge_flux_qd_num[sp_id]);
	for(int qd_i = 0; qd_i < col_edge_flux_qd_num[sp_id]; qd_i++)
	{
		value_type del_qd = 2.0/(col_edge_flux_qd_num[sp_id] + 1);
		col_edge_flux_qd_points[sp_id][qd_i] = -1.0 + del_qd*(qd_i+1);
	}
	tot_Qv_h_size_per_ix += nve*col_edge_flux_qd_num[sp_id];


	//volume term part
	int tot_Qv_g_size_per_ix = nv_loc*3*dof_v*dof_v;
	//surface term part
	tot_Qv_g_size_per_ix += nve*2*(dim_fhat_valid_v_basis)*(2*dof_v);


	for(int k1 = 0; k1 < 2; k1++)
	{
		for(int ix = 0; ix < nx; ix++)
		{
			col3_Qv_h0_ab_arr[k1][ix][sp_id].resize(tot_species_num, Vector::Zero(tot_Qv_h_size_per_ix));
			col3_Qv_g0_ab_arr[k1][ix][sp_id].resize(tot_species_num, Vector::Zero(tot_Qv_g_size_per_ix));
		}
	}

	col3_hg_Xmat_arr[sp_id].resize(nx);
	col3_h_LF_Vmat_arr[sp_id].resize(nx*nve);

	col3_h_to_Qv_mat[sp_id].resize(nx, SparseMatrix(tot_Qv_h_size_per_ix, col3_hg_tot_n));


	for(int ix = 0; ix < nx; ix++)
	{
		col3_h_to_Qv_mat_whole[ix][sp_id] = SparseMatrix(tot_Qv_h_size_per_ix, 2*col3_hg_tot_n);
		col3_g_to_Qv_mat[ix][sp_id]= SparseMatrix(tot_Qv_g_size_per_ix, col3_hg_tot_n);
		//common spatial part matrix setup
		col3_hg_Xmat_arr[sp_id][ix].resize(3, Vector::Zero(dof_x*dof_x));

		Vector tmp_hg_to_Qx_Vec = Vector::Zero(3*dof_x*dof_x);

		integration_col_arr[sp_id]->col3_hg_to_Qx_cal(ElementX(ix), tmp_hg_to_Qx_Vec);

		for(int j = 0; j < 3; j++)
		{
			for(int lxp = 0; lxp < dof_x; lxp++)
			{
				for(int lx = 0; lx < dof_x; lx++)
				{

					int tmp_index = (j*dof_x + lxp)*dof_x + lx;

					col3_hg_Xmat_arr[sp_id][ix][j][lxp*dof_x + lx] = tmp_hg_to_Qx_Vec[tmp_index];

				}
			}
		}

		int h_to_Qv_loc_index = 0, g_to_Qv_loc_index = 0;

		//volume term part
		vector<SparseMatrix_Triplet> tripletList_h_to_Qv, tripletList_g_to_Qv;
		vector<SparseMatrix_Triplet> tripletList_h_to_Qv_part1, tripletList_h_to_Qv_part2;
		vector<SparseMatrix_Triplet> tripletList_h_to_Qv_whole;
		for(int iv = 0; iv < nv_loc; iv++)
		{
			Vector tmp_hg_to_Q_Vec = Vector::Zero(5*dof_v*dof_v*col3_hg_tot_n);
			integration_col_arr[sp_id]->col3_hg_to_Qv_vol_cal(ElementX(ix), ElementV(iv), tmp_hg_to_Q_Vec);

			for(int j = 0; j < 5; j++)
			{
				for(int lvp = 0; lvp < dof_v; lvp++)
				{
					for(int lv = 0; lv < dof_v; lv++) 
					{
						for(int l = 0; l < col3_hg_tot_n; l++)
						{
							int tmp_index = ((j*dof_v + lvp)*dof_v + lv)*col3_hg_tot_n + l;
							if(tmp_hg_to_Q_Vec[tmp_index] != 0.0)
							{
								if(j < 2)
								{
									int row_index = ((j*nv_loc + iv)*dof_v + lvp)*dof_v + lv;
									tripletList_h_to_Qv.push_back(SparseMatrix_Triplet(row_index, l, tmp_hg_to_Q_Vec[tmp_index]));
									tripletList_h_to_Qv_part1.push_back(SparseMatrix_Triplet(row_index, l, tmp_hg_to_Q_Vec[tmp_index]));
									tripletList_h_to_Qv_whole.push_back(SparseMatrix_Triplet(row_index, l, tmp_hg_to_Q_Vec[tmp_index]));

								}
								else if(j < 5)
								{
									int row_index = (((j-2)*nv_loc + iv)*dof_v + lvp)*dof_v + lv;
									tripletList_g_to_Qv.push_back(SparseMatrix_Triplet(row_index, l, tmp_hg_to_Q_Vec[tmp_index]));
								}
							}
						}
					}
				}
			}
		}

		h_to_Qv_loc_index += nv_loc*2*dof_v*dof_v;
		g_to_Qv_loc_index += nv_loc*3*dof_v*dof_v;

		//surface term part 
		for(int i = 0; i < vspace_edge_num_of_type; i++)
		{
			int n_size = vspace_edge_list_arr[sp_id][i].size();
			for(int j = 0; j < n_size; j++)
			{
				int e = vspace_edge_list_arr[sp_id][i][j];
				int ele[2] = {vspace_edge_ele_arr[sp_id][2*e], vspace_edge_ele_arr[sp_id][2*e + 1]};

				Vector tmp_g_S_Vmat_Vec = Vector::Zero((2*(dim_fhat_valid_v_basis*2*dof_v))*col3_hg_tot_n);
				if (i == 0 || i == 1)
				{
					Vector tmp_h_LF_Vmat_Vec = Vector::Zero(4*dof_v*dof_v);
					Vector tmp_h_S_Vmat_Vec = Vector::Zero(4*dof_v*dof_v*col3_hg_tot_n);
					integration_col_arr[sp_id]->col3_hg_S_LF_Vmat_cal(ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), i, tmp_h_LF_Vmat_Vec, tmp_h_S_Vmat_Vec, tmp_g_S_Vmat_Vec);

					//h & LF part for mixed flux of upwind and LF
					col3_h_LF_Vmat_arr[sp_id][ix*nve + e].resize(4, Vector::Zero(dof_v*dof_v));
					for(int tmp_j = 0; tmp_j < 4; tmp_j++)
					{
						col3_h_LF_Vmat_arr[sp_id][ix*nve + e][tmp_j] = tmp_h_LF_Vmat_Vec.segment(tmp_j*dof_v*dof_v, dof_v*dof_v);

						for(int lvp = 0; lvp < dof_v; lvp++)
						{
							for(int lv = 0; lv < dof_v; lv++) 
							{
								for(int l = 0; l < col3_hg_tot_n; l++)
								{
									int tmp_index2 = ((tmp_j*dof_v + lvp)*dof_v + lv)*col3_hg_tot_n + l;
									if(tmp_h_S_Vmat_Vec[tmp_index2] != 0.0)
									{

										int row_index = h_to_Qv_loc_index + ((e*4 + tmp_j)*dof_v + lvp)*dof_v + lv;
										tripletList_h_to_Qv.push_back(SparseMatrix_Triplet(row_index, l, tmp_h_S_Vmat_Vec[tmp_index2]));
										tripletList_h_to_Qv_part2.push_back(SparseMatrix_Triplet(row_index - tot_Qv_h_size_per_ix_vol, l, tmp_h_S_Vmat_Vec[tmp_index2]));
										tripletList_h_to_Qv_whole.push_back(SparseMatrix_Triplet(row_index, l + col3_hg_tot_n, tmp_h_S_Vmat_Vec[tmp_index2]));
									}
								}
							}
						}
					}
				}
				else
				{
					integration_col_arr[sp_id]->col3_g_S_Vmat_bd_cal(ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), i, tmp_g_S_Vmat_Vec);
				}

				//g part
				for(int tmp_j = 0; tmp_j < 2; tmp_j++)
				{
					for(int row = 0; row < dim_fhat_valid_v_basis; row++) 
					{
						for(int i_ele = 0; i_ele < 2; i_ele++)
						{
							for(int lv = 0; lv < dof_v; lv++)
							{
								for(int l = 0; l < col3_hg_tot_n; l++)
								{
									int tmp_index = ((tmp_j*dim_fhat_valid_v_basis + row)*2 + i_ele)*dof_v + lv;
									int tmp_index2 = tmp_index*col3_hg_tot_n + l;

									if(tmp_g_S_Vmat_Vec[tmp_index2] != 0.0)
									{
										int row_index = g_to_Qv_loc_index + (((e*2 +tmp_j)*dim_fhat_valid_v_basis + row)*2 + i_ele)*dof_v + lv;

										tripletList_g_to_Qv.push_back(SparseMatrix_Triplet(row_index, l, tmp_g_S_Vmat_Vec[tmp_index2]));

									}
								}
							}
						}
					}
				}
			}
		}

		h_to_Qv_loc_index += nve*4*dof_v*dof_v;
		g_to_Qv_loc_index += nve*2*(dim_fhat_valid_v_basis)*(2*dof_v);

		//flux from h at edge qd points
		for(int i = 0; i < vspace_edge_num_of_type; i++)
		{
			int n_size = vspace_edge_list_arr[sp_id][i].size();
			for(int j = 0; j < n_size; j++)
			{
				int e = vspace_edge_list_arr[sp_id][i][j];
				int ele[2] = {vspace_edge_ele_arr[sp_id][2*e], vspace_edge_ele_arr[sp_id][2*e + 1]};

				if (i == 0 || i == 1)
				{
					Vector tmp_h_flux_qd = Vector::Zero(col_edge_flux_qd_num[sp_id]*col3_hg_tot_n);
					integration_col_arr[sp_id]->col3_h_flux_qd_cal(ElementX(ix), EdgeV(e), i, col_edge_flux_qd_points[sp_id], tmp_h_flux_qd);


					for(int tmp_j = 0; tmp_j < col_edge_flux_qd_num[sp_id]; tmp_j++)
					{
						for(int l = 0; l < col3_hg_tot_n; l++)
						{
							int tmp_index2 = tmp_j*col3_hg_tot_n + l;
							if(tmp_h_flux_qd[tmp_index2] != 0.0)
							{
								int row_index = h_to_Qv_loc_index + e*col_edge_flux_qd_num[sp_id] + tmp_j;
								tripletList_h_to_Qv.push_back(SparseMatrix_Triplet(row_index, l, tmp_h_flux_qd[tmp_index2]));
								tripletList_h_to_Qv_part2.push_back(SparseMatrix_Triplet(row_index - tot_Qv_h_size_per_ix_vol, l, tmp_h_flux_qd[tmp_index2]));
								tripletList_h_to_Qv_whole.push_back(SparseMatrix_Triplet(row_index, l + col3_hg_tot_n, tmp_h_flux_qd[tmp_index2]));
							}
						}
					}
				}
			}
		}

		col3_h_to_Qv_mat[sp_id][ix].setFromTriplets(tripletList_h_to_Qv.begin(), tripletList_h_to_Qv.end()); 
		col3_h_to_Qv_mat_whole[ix][sp_id].setFromTriplets(tripletList_h_to_Qv_whole.begin(), tripletList_h_to_Qv_whole.end()); 
		col3_g_to_Qv_mat[ix][sp_id].setFromTriplets(tripletList_g_to_Qv.begin(), tripletList_g_to_Qv.end());
	}

	//consv operator setup
	//col3_f_to_mom_consv_mat_arr is a parallel shift in V_parallel direction : change momentum mainly
	//col3_f_to_en_consv_mat_arr is a isotropic diffusion in v space : change energy mainly
	col3_f_to_mom_consv_mat_arr[sp_id].resize(nx);
	col3_f_to_en_consv_mat_arr[sp_id].resize(nx);
	for(int ix = 0; ix < nx; ix++)
	{
		col3_f_to_mom_consv_mat_arr[sp_id][ix].resize(2, SparseMatrix(nv_loc*dof, nv_loc*dof));
		col3_f_to_en_consv_mat_arr[sp_id][ix] = SparseMatrix(nv_loc*dof, nv_loc*dof);

		vector<SparseMatrix_Triplet> tripletList_f_to_mom_consv1, tripletList_f_to_mom_consv2; 
		vector<SparseMatrix_Triplet> tripletList_f_to_en_consv_vol, tripletList_f_to_en_consv_surf;

		//volume term part
		for(int iv = 0; iv < nv_loc; iv++)
		{
			Vector tmp_f_to_consv_vol = Vector::Zero(2*dof*dof);
			integration_col_arr[sp_id]->col3_f_to_consv_vol_cal(ElementX(ix), ElementV(iv), tmp_f_to_consv_vol);


			for(int type_index = 0; type_index < 2; type_index++)
			{

				for(int lp = 0; lp < dof; lp++)
				{
					for(int l = 0; l < dof; l++)
					{
						int tmp_index = (type_index*dof + lp)*dof + l;
						if(tmp_f_to_consv_vol[tmp_index] != 0.0)
						{
							if(type_index == 0)
							{
								tripletList_f_to_mom_consv1.push_back(SparseMatrix_Triplet(iv*dof + l, iv*dof + lp, tmp_f_to_consv_vol[tmp_index]));
								tripletList_f_to_mom_consv2.push_back(SparseMatrix_Triplet(iv*dof + l, iv*dof + lp, tmp_f_to_consv_vol[tmp_index]));

							}
							else
							{
								tripletList_f_to_en_consv_vol.push_back(SparseMatrix_Triplet(iv*dof + l, iv*dof + lp, tmp_f_to_consv_vol[tmp_index]));
							}
						}
					}
				}
			}
		}


		//surface term part 
		for(int i = 0; i < vspace_edge_num_of_type; i++)
		{
			int n_size = vspace_edge_list_arr[sp_id][i].size();
			for(int j = 0; j < n_size; j++)
			{
				int e = vspace_edge_list_arr[sp_id][i][j];
				int ele[2] = {vspace_edge_ele_arr[sp_id][2*e], vspace_edge_ele_arr[sp_id][2*e + 1]};

				if (i == 0)
				{
					Vector tmp_f_to_mom_consv_surf = Vector::Zero(2*(2*dof)*(2*dof));

					integration_col_arr[sp_id]->col3_f_to_mom_consv_surf_cal(ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), tmp_f_to_mom_consv_surf);

					for(int type_index = 0; type_index < 2; type_index++)
					{
						for(int ele_lp = 0; ele_lp < 2; ele_lp++)
						{
							for(int lp = 0; lp < dof; lp++)
							{
								int lp_tot = ele_lp*dof + lp;

								for(int ele_l = 0; ele_l < 2; ele_l++)
								{
									for(int l = 0; l < dof; l++)
									{
										int l_tot = ele_l*dof + l;
										int tmp_index = (type_index*2*dof + lp_tot)*2*dof + l_tot;
										if(tmp_f_to_mom_consv_surf[tmp_index] != 0.0)
										{
											if(type_index == 0)
											{
												tripletList_f_to_mom_consv1.push_back(SparseMatrix_Triplet(ele[ele_l]*dof + l, ele[ele_lp]*dof + lp, tmp_f_to_mom_consv_surf[tmp_index]));
											}
											else
											{
												tripletList_f_to_mom_consv2.push_back(SparseMatrix_Triplet(ele[ele_l]*dof + l, ele[ele_lp]*dof + lp, tmp_f_to_mom_consv_surf[tmp_index]));
											}
										}
									}
								}
							}
						}
					}
				}


				Vector tmp_f_to_en_consv_surf = Vector::Zero((dof_x*dim_fhat_valid_v_basis)*(2*dof));

				if (i == 0 || i == 1)
				{
					integration_col_arr[sp_id]->col3_f_to_en_consv_inner_surf_cal(ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), i, tmp_f_to_en_consv_surf);
				}
				else
				{
					integration_col_arr[sp_id]->col3_f_to_en_consv_bd_surf_cal(ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), i, tmp_f_to_en_consv_surf);
				}

				for(int lp_x = 0; lp_x < dof_x; lp_x++)
				{
					for(int lp_v = 0; lp_v < dim_fhat_valid_v_basis; lp_v++) 
					{
						int lp_tot = lp_x*dim_fhat_valid_v_basis + lp_v;

						for (int ele_index = 0; ele_index < 2; ele_index++)
						{
							for(int l = 0; l < dof; l++)
							{
								int l_tot = ele_index*dof + l;
								int tmp_index = lp_tot*2*dof + l_tot;
								if(tmp_f_to_en_consv_surf[tmp_index] != 0.0)
								{
									tripletList_f_to_en_consv_surf.push_back(SparseMatrix_Triplet(ele[ele_index]*dof + l, e*dof_x*dim_fhat_valid_v_basis + lp_tot, tmp_f_to_en_consv_surf[tmp_index]));


									if(rank == 0)
									{
										int t1 = ele[ele_index]*dof + l;
										int t2 = e*dof_x*dim_fhat_valid_v_basis + lp_tot;
										if(t1 < 0 || t1 >= nv_loc*dof)
										{
											cout << "error in t1 : " << t1 << " " << nv_loc*dof << " " << i << " " << ele_index << " " << ele[ele_index] << " " << tmp_f_to_en_consv_surf[tmp_index] << endl;
											cout << e << " " << ele[0] << " " << ele[1] << " " << lp_tot << " " << l << " " << endl;
										}

										if(t2 < 0 || t2 >= nve*dof_x*dim_fhat_valid_v_basis)
										{
											cout << "error in t2 : " << t2 << " " << nve*dof_x*dim_fhat_valid_v_basis << endl;
										}

									}
								}
							}
						}
					}
				}
			}
		}

		col3_f_to_mom_consv_mat_arr[sp_id][ix][0].setFromTriplets(tripletList_f_to_mom_consv1.begin(), tripletList_f_to_mom_consv1.end());
		col3_f_to_mom_consv_mat_arr[sp_id][ix][1].setFromTriplets(tripletList_f_to_mom_consv2.begin(), tripletList_f_to_mom_consv2.end());

		col3_f_to_en_consv_mat_arr[sp_id][ix].setFromTriplets(tripletList_f_to_en_consv_vol.begin(), tripletList_f_to_en_consv_vol.end());

		SparseMatrix tmp_mat = SparseMatrix(nv_loc*dof, nve*dof_x*dim_fhat_valid_v_basis);

		tmp_mat.setFromTriplets(tripletList_f_to_en_consv_surf.begin(), tripletList_f_to_en_consv_surf.end());
		col3_f_to_en_consv_mat_arr[sp_id][ix] += tmp_mat*col_f_dg_to_valid_fhat_mat_arr[sp_id][ix];

	}

	col3_M_solver_arr[sp_id].resize(nx);
}

//setup for col_diag_cell_consv_mat[sp_id][ix] : needed for conservation 
//[mom, en] for [sp_id, ix] = col_diag_cell_consv_mat[sp_id][ix] * M * f 
void Collision_dg::col_consv_mat_setup(SparseMatrix &org_mat, int sp_id, int ix)
{

	int dof = dof_arr[sp_id];

	col3_M_solver_arr[sp_id][ix] = new Eigen::PardisoLLT<SparseMatrix>;
	col3_M_solver_arr[sp_id][ix]->analyzePattern(org_mat);
	col3_M_solver_arr[sp_id][ix]->factorize(org_mat);

	int nv_loc = nv[sp_id];
	int n = nv_loc*dof;

	vector<SparseMatrix_Triplet> tripletList_new_mat;
	tripletList_new_mat.reserve(nv_loc*dof*dof);

	SparseMatrix new_mat(n,n);

	vector<SparseMatrix_Triplet> tripletList_Minv_mat;
	for(int iv = 0; iv < nv_loc; iv++)
	{
		Matrix M(dof, dof);
		int k = iv*dof;
		for(int i = 0; i < dof; i++)
		{
			for(int j = 0; j < dof; j++) M(i,j) = org_mat.coeffRef(k+i, k+j);
		}

		Matrix M_inv = M.fullPivLu().inverse();
	
		for(int i = 0; i < dof; i++)
		{
			for(int j = 0; j < dof; j++) 
			{
				tripletList_Minv_mat.push_back(SparseMatrix_Triplet(k+i, k+j, M_inv(i,j)));
			}
		}

		for(int i = 0; i < dof; i++)
		{
			for(int j = 0; j < dof; j++)
			{
				tripletList_new_mat.push_back(SparseMatrix_Triplet(k+i, k+j, M_inv(i,j)));
			}
		}
	}

	col_sp_M_arr[sp_id][ix] = org_mat;
	col_sp_Minv_arr[sp_id][ix].setFromTriplets(tripletList_Minv_mat.begin(), tripletList_Minv_mat.end());



	new_mat.setFromTriplets(tripletList_new_mat.begin(), tripletList_new_mat.end());
	col_diag_cell_consv_mat[sp_id][ix] = col_diag_cell_consv_mat[sp_id][ix]*new_mat;

	vector<SparseMatrix_Triplet> tripletList_MC_fa_del_hg_b;

	int dof_x = dof_x_arr[sp_id]; 
	int dof_v = dof_v_arr[sp_id]; 
	int col3_hg_tot_n = col3_hg_tot_n_arr[sp_id];

	for(int iv = 0; iv < nv_loc; iv++)
	{

		int vp_tot_min_index = vp_tot_min_index_arr[sp_id][ix][iv];
		int vp_tot_max_index = vp_tot_max_index_arr[sp_id][ix][iv];
		int u_tot_min_index = u_tot_min_index_arr[sp_id][ix][iv];
		int u_tot_max_index = u_tot_max_index_arr[sp_id][ix][iv];

		for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
		{
			for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
			{
				int l_hg = j4*(col3_vp_n_arr[sp_id] + 3) + j3;

				//for h part
				for(int j = 0; j < 2; j++)
				{
					for(int lvp = 0; lvp < dof_v; lvp++)
					{
						int index_lvp = ((j*nv_loc + iv)*dof_v + lvp)*dof_v;

						for(int lv = 0; lv < dof_v; lv++)
						{
							int loc_index1 = index_lvp + lv;
							value_type h_coeff = col3_h_to_Qv_mat[sp_id][ix].coeffRef(loc_index1, l_hg);

							if(h_coeff != 0.0)
							{
								for(int lxp = 0; lxp < dof_x; lxp++)
								{	
									int lp_index = lvp*dof_x + lxp;
									int col_index = iv*dof*2*col3_hg_tot_n;
									col_index += lp_index*2*col3_hg_tot_n + (0 + l_hg);

									for(int lx = 0; lx < dof_x; lx++)
									{
										int l_index = lv*dof_x + lx;
										value_type tmp_X_part = col3_hg_Xmat_arr[sp_id][ix][j][lxp*dof_x + lx];

										int row_index = iv*dof + l_index; 
										tripletList_MC_fa_del_hg_b.push_back(SparseMatrix_Triplet(row_index, col_index, tmp_X_part*h_coeff));
									}
								}
							}
						}
					}
				}

				//for g part
				for(int j = 0; j < 3; j++)
				{
					for(int lvp = 0; lvp < dof_v; lvp++)
					{
						int index_lvp = ((j*nv_loc + iv)*dof_v + lvp)*dof_v;

						for(int lv = 0; lv < dof_v; lv++)
						{
							int loc_index1 = index_lvp + lv;
							value_type g_coeff = col3_g_to_Qv_mat[ix][sp_id].coeffRef(loc_index1, l_hg);

							if(g_coeff != 0.0)
							{
								for(int lxp = 0; lxp < dof_x; lxp++)
								{	
									int lp_index = lvp*dof_x + lxp;
									int col_index = iv*dof*2*col3_hg_tot_n;
									col_index += lp_index*2*col3_hg_tot_n + (col3_hg_tot_n + l_hg);

									for(int lx = 0; lx < dof_x; lx++)
									{
										int l_index = lv*dof_x + lx;
										value_type tmp_X_part = col3_hg_Xmat_arr[sp_id][ix][j][lxp*dof_x + lx];

										int row_index = iv*dof + l_index; 
										tripletList_MC_fa_del_hg_b.push_back(SparseMatrix_Triplet(row_index, col_index ,tmp_X_part*g_coeff));
									}
								}
							}
						}
					}
				}
			}
		}
	}

	SparseMatrix MC_fa_del_hg_b_mat_part1 = SparseMatrix(nv_loc*dof, nv_loc*dof*2*col3_hg_tot_n);
	MC_fa_del_hg_b_mat_part1.setFromTriplets(tripletList_MC_fa_del_hg_b.begin(), tripletList_MC_fa_del_hg_b.end());

	SparseMatrix consv_MC_fa_del_hg_b_mat = col_diag_cell_consv_mat[sp_id][ix]*MC_fa_del_hg_b_mat_part1;

	vector<SparseMatrix_Triplet> tripletList_MC_fa_del_hg_b_part2;
	for(int iv = 0; iv < nv_loc; iv++)
	{
		int vp_tot_min_index = vp_tot_min_index_arr[sp_id][ix][iv];
		int vp_tot_max_index = vp_tot_max_index_arr[sp_id][ix][iv];
		int u_tot_min_index = u_tot_min_index_arr[sp_id][ix][iv];
		int u_tot_max_index = u_tot_max_index_arr[sp_id][ix][iv];

		for(int j3 = vp_tot_min_index; j3 < vp_tot_max_index + 1; j3++)
		{
			for(int j4 = u_tot_min_index; j4 < u_tot_max_index + 1; j4++)
			{
				for(int hg_index = 0; hg_index < 2; hg_index++)
				{
					int l_hg2 = hg_index*col3_hg_tot_n + j4*(col3_vp_n_arr[sp_id] + 3) + j3;

					for(int consv_i = 0; consv_i < 2; consv_i++)
					{
						for(int lp = 0; lp < dof; lp++)
						{
							int col_index = iv*dof*2*col3_hg_tot_n;
							col_index += lp*2*col3_hg_tot_n + l_hg2;

							value_type tmp_coeff = consv_MC_fa_del_hg_b_mat.coeffRef(consv_i, col_index);
							if(tmp_coeff != 0.0) 
							{
								int row_index2 = (iv*dof + lp)*2 + consv_i;
								tripletList_MC_fa_del_hg_b_part2.push_back(SparseMatrix_Triplet(row_index2, l_hg2, tmp_coeff));
							}
						}
					}


				}
			}
		}


	}

	consv_fa_del_hg_b_mat_arr[ix][sp_id] = SparseMatrix(2*nv_loc*dof, 2*col3_hg_tot_n);
	consv_fa_del_hg_b_mat_arr[ix][sp_id].setFromTriplets(tripletList_MC_fa_del_hg_b_part2.begin(), tripletList_MC_fa_del_hg_b_part2.end());

}

//lowest order part of f update module
//called when Gyrokinetic class is generated
void Collision_dg::initial_h0g0_setup(vector<vector<value_type>> &_coeff_tot)
{
	f_data_size_per_ix.resize(tot_species_num, 0);
	for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
		f_data_size_per_ix[sp_id] = nv[sp_id]*dof_arr[sp_id];

	for(int ix = 0; ix < nx; ix++)
	{
		value_type R_cen = R_cen_cell[ix], Z_cen = Z_cen_cell[ix];
		int loc_rgn = eq_reader->rgn_fn(R_cen, Z_cen);

		for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
		{
			for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
			{
				if(sp_kinetic_op_arr[sp_a_id] == 1 && sp_kinetic_op_arr[sp_b_id] == 1)
				{
					col3_sp_ab_ix_col_flag_arr[ix][sp_a_id][sp_b_id] = 1;
				}
				else
					col3_sp_ab_ix_col_flag_arr[ix][sp_a_id][sp_b_id] = 0;

			}
		}

		for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
		{
			if(sp_kinetic_op_arr[sp_id] == 1)
			{
				value_type n_val, U_val, T_val;
				Vector f_loc = Vector::Zero(f_data_size_per_ix[sp_id]);

				//copy f data from _coeff_tot to f for ix element
				for(int i = 0; i < f_data_size_per_ix[sp_id]; i++)
				{
					int loc_i = ix*f_data_size_per_ix[sp_id] + i;
					f_loc[i] = _coeff_tot[sp_id][loc_i];
				}

				Vector tmp_avg = col_diag_cell_mat[sp_id][ix]*f_loc;

				col3_nUT_arr[sp_id][ix][0] = n_val = tmp_avg[0]/vol_cell[ix];
				col3_nUT_arr[sp_id][ix][1] = U_val = tmp_avg[1]/tmp_avg[0];
				col3_nUT_arr[sp_id][ix][2] = T_val = Ms_arr[sp_id]*(tmp_avg[2] - 0.5*U_val*U_val*tmp_avg[0])/(1.5*tmp_avg[0]);
				col3_nUT_before_arr[sp_id][ix] = col3_nUT_arr[sp_id][ix];
				col3_fM_coeff_update_before_time_arr[ix] = 0;

				int h_size = col3_h_to_Qv_mat[sp_id][ix].rows();
				int g_size = col3_g_to_Qv_mat[ix][sp_id].rows();

				col3_Qv_h0_ab_summed_arr[ix][sp_id] = Vector::Zero(h_size);
				col3_Qv_g0_ab_summed_arr[ix][sp_id] = Vector::Zero(g_size);

				value_type T_para_eV = Ms_arr[sp_id]*(tmp_avg[5] - 0.5*U_val*U_val*tmp_avg[0])/(1.5*tmp_avg[0])*col3_norm_T0_eV/3.0;
				value_type T_perp_eV = Ms_arr[sp_id]*tmp_avg[6]/(1.5*tmp_avg[0])*col3_norm_T0_eV/3.0;

				if(n_val < 0.0 || T_para_eV < col3_T_eV_min_lim/3.0 || T_perp_eV < col3_T_eV_min_lim*2.0/3.0)
				{
					for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
					{
						col3_sp_ab_ix_col_flag_arr[ix][sp_id][sp_b_id] = 0;
						col3_sp_ab_ix_col_flag_arr[ix][sp_b_id][sp_id] = 0;
					}
				}
			}
		}

		value_type del_t_max_theo_tot = 1e5;
		vector<value_type> nu_sp_arr(tot_species_num, 0.0);

		for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
		{
			for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
			{
				if(sp_kinetic_op_arr[sp_a_id] == 1 && sp_kinetic_op_arr[sp_b_id] == 1)
				{
					value_type tau_ab = taui_ab_coll_time_norm(sp_a_id, col3_nUT_arr[sp_a_id][ix][0], col3_nUT_arr[sp_a_id][ix][2], sp_b_id, col3_nUT_arr[sp_b_id][ix][0], col3_nUT_arr[sp_b_id][ix][2]);
					nu_sp_arr[sp_a_id] = max(nu_sp_arr[sp_a_id], 1.0/tau_ab);

					if(sp_a_id == sp_b_id) tau_aa_init[sp_a_id][ix] = tau_ab; 
				}
			}
		}


		for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
		{	
			if(sp_kinetic_op_arr[sp_id] == 1)
			{
				sp_h0g0_SE_update(sp_id, ix);
			}
		}

		if(col3_hMbgMb_with_UaTa_op == 1)
		{
			for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
			{
				if(sp_kinetic_op_arr[sp_id] == 1)
				{
					sp_hMbgMb_with_UaTa_update(sp_id, ix);
				}
			}
		}


	}
}

//lowest order part of f update module
//called periodically
void Collision_dg::h0g0fM_update(vector<vector<value_type>> &_coeff_tot, const int &t_step, const int &restart_op) 
{

	if (t_step%col_fM_update_period == 0 || restart_op == 1)
	{
		clock_t testc_1, testc_2;
		if(rank == 0) testc_1 = clock();

		for(int ix = 0; ix < nx; ix++)
		{
			value_type R_cen = R_cen_cell[ix], Z_cen = Z_cen_cell[ix];
			int loc_rgn = eq_reader->rgn_fn(R_cen, Z_cen);

			col3_valid_ix_flag_arr[ix] = 1;
			for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
			{
				for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
				{
					if(sp_kinetic_op_arr[sp_a_id] == 1 && sp_kinetic_op_arr[sp_b_id] == 1)
					{
						col3_sp_ab_ix_col_flag_arr[ix][sp_a_id][sp_b_id] = 1;
					}
					else
						col3_sp_ab_ix_col_flag_arr[ix][sp_a_id][sp_b_id] = 0;


					
				}
			}

		}

		for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
		{
			if(sp_kinetic_op_arr[sp_id] == 1)
			{
				Vector f_loc = Vector::Zero(f_data_size_per_ix[sp_id]);

				int neg_nT_count_local = 0, neg_nT_count_total;
				int tot_nx_count_local = 0, tot_nx_count_total;
				int valid_nx_count_local = 0, valid_nx_count_total;

				for(int ix = 0; ix < nx; ix++)
				{
					tot_nx_count_local += 1;
					int valid_flag = 0;
					for(int sp_dum_id = 0; sp_dum_id < tot_species_num; sp_dum_id++)
					{
						if(col3_sp_ab_ix_col_flag_arr[ix][sp_id][sp_dum_id] == 1 || col3_sp_ab_ix_col_flag_arr[ix][sp_dum_id][sp_id] == 1)
							valid_flag = 1;
					}


					if(valid_flag == 1)
					{
						//copy f data from _coeff_tot to f_loc for ix element
						copy(_coeff_tot[sp_id].begin() + ix*f_data_size_per_ix[sp_id], _coeff_tot[sp_id].begin() + (ix + 1)*f_data_size_per_ix[sp_id], f_loc.data());


						Vector tmp_avg = col_diag_cell_mat[sp_id][ix]*f_loc;

						value_type n_val, U_val, T_val;
						col3_nUT_arr[sp_id][ix][0] = n_val = tmp_avg[0]/vol_cell[ix];
						col3_nUT_arr[sp_id][ix][1] = U_val = tmp_avg[1]/tmp_avg[0];
						col3_nUT_arr[sp_id][ix][2] = T_val = Ms_arr[sp_id]*(tmp_avg[2] - 0.5*U_val*U_val*tmp_avg[0])/(1.5*tmp_avg[0]);


						value_type T_para, T_perp;
						T_para = 2.0*Ms_arr[sp_id]*(tmp_avg[5] - 0.5*U_val*U_val*tmp_avg[0])/(tmp_avg[0]);
						T_perp = Ms_arr[sp_id]*tmp_avg[6]/(tmp_avg[0]);


						value_type T_para_eV = T_para*col3_norm_T0_eV;
						value_type T_perp_eV = T_perp*col3_norm_T0_eV;

						if(n_val < 0.0 || (T_para_eV/3.0 + T_perp_eV*2.0/3.0) < col3_T_eV_min_lim)
						{
							for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
							{
								col3_sp_ab_ix_col_flag_arr[ix][sp_id][sp_b_id] = 0;
								col3_sp_ab_ix_col_flag_arr[ix][sp_b_id][sp_id] = 0;
							}


							value_type R_cen = R_cen_cell[ix], Z_cen = Z_cen_cell[ix];
							int loc_rgn = eq_reader->rgn_fn(R_cen, Z_cen);
							value_type loc_xx = sqrt(eq_reader->psi_ov_psix_interpol(R_cen, Z_cen));

							cout << "negative n or small T found [rank, ix, global ix, sp_id, n (MKS), T_para (eV), T_perp (eV), R, Z, xx, rgn] : " << rank << " " << ix << " " << ilower_ele + ix << " " << sp_id << " " << n_val*ni00_norm_mks << " " << T_para_eV << " " << T_perp_eV << " " << R_cen << " " << Z_cen << " " << loc_xx << " " << loc_rgn << endl;
							neg_nT_count_local += 1;

						}
						else
							valid_nx_count_local += 1;
					}

				}

				MPI_Allreduce(&neg_nT_count_local, &neg_nT_count_total, 1, MPI_INT, MPI_SUM, comm);
				MPI_Allreduce(&tot_nx_count_local, &tot_nx_count_total, 1, MPI_INT, MPI_SUM, comm);
				MPI_Allreduce(&valid_nx_count_local, &valid_nx_count_total, 1, MPI_INT, MPI_SUM, comm);

				if(rank == 0 && neg_nT_count_total > 0) cout << "negative n or small T found [sp_id, total ix] : " << sp_id << " " << neg_nT_count_total << endl;
				value_type diag_ratio = double(nstep)/double(col_diag_1d_period*500);
				int diag_ratio_int = max(1, int(diag_ratio));
				if(rank == 0 && t_step%(col_diag_1d_period*diag_ratio_int) == 0)
					cout << "total valid nx [sp_id, valid nx, total nx] : " << sp_id << " " << valid_nx_count_total << " " << tot_nx_count_total << endl;

			}
		}

		for(int ix = 0; ix < nx; ix++)
		{
			value_type fM_update_1st_err = 5e-2;
			value_type fM_update_2nd_err = 1e-2;
			value_type fM_update_3rd_err = 1e-3;

			int t_2nd_err = 100, t_3rd_err = 500;


			value_type max_err = 0.0, err_tor = 1.0;

			value_type n0, U0, T0, n1, U1, T1;

			int max_err_sp, max_err_mom;
			value_type max_err_v0, max_err_v1;
			for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
			{
				if(sp_kinetic_op_arr[sp_id] == 1)
				{
					n0 = col3_nUT_before_arr[sp_id][ix][0];
					U0 = col3_nUT_before_arr[sp_id][ix][1];
					T0 = col3_nUT_before_arr[sp_id][ix][2];

					n1 = col3_nUT_arr[sp_id][ix][0];
					U1 = col3_nUT_arr[sp_id][ix][1];
					T1 = col3_nUT_arr[sp_id][ix][2];

					value_type n_diff, U_diff, T_diff, v_th;

					v_th = sqrt(T0/Ms_arr[sp_id]);

					n_diff = abs((n1 - n0)/n0);
					T_diff = abs((T1 - T0)/T0);
					U_diff = abs((U1 - U0)/(10.0*v_th));

					max_err = max(max_err, max(n_diff, max(T_diff, U_diff)));

					if(n_diff == max_err)
					{
						max_err_mom = 0;
						max_err_sp = sp_id;
						max_err_v0 = n0;
						max_err_v1 = n1;
					}
					if(T_diff == max_err)
					{
						max_err_mom = 1;
						max_err_sp = sp_id;
						max_err_v0 = T0;
						max_err_v1 = T1;
					}
					if(U_diff == max_err)
					{
						max_err_mom = 2;
						max_err_sp = sp_id;
						max_err_v0 = U0;
						max_err_v1 = U1;
					}
				}
			}

			int f0_update_flag = 0;
			int del_time_step = t_step - col3_fM_coeff_update_before_time_arr[ix]; 

			if(del_time_step >= (t_3rd_err*col_fM_update_period))
				err_tor = fM_update_3rd_err;
			else if(del_time_step >= (t_2nd_err*col_fM_update_period))
				err_tor = abs(fM_update_2nd_err - (fM_update_2nd_err - fM_update_3rd_err)/(t_3rd_err - t_2nd_err)*(del_time_step/double(col_fM_update_period) - t_2nd_err));
			else if(del_time_step >= (1*col_fM_update_period))
				err_tor = abs(fM_update_1st_err - (fM_update_1st_err - fM_update_2nd_err)/(t_2nd_err - 1)*(del_time_step/double(col_fM_update_period) - 1));


			if(restart_op == 1 || t_step == 0)
				f0_update_flag = 1;
			else if(max_err > err_tor)
				f0_update_flag = 1;

			//f0_update_flag = 1;



			if(f0_update_flag == 1)
			{
				if(rank == 0 && ix == 0) cout << "fM is updated at time step = " << t_step << ", max_err, sp, mom_type, v0, v1 = " << max_err << " " << max_err_sp << " " << max_err_mom << " " << max_err_v0 << " " << max_err_v1 << endl;
				col3_fM_coeff_update_before_time_arr[ix] = t_step;

				for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
				{
					if(sp_kinetic_op_arr[sp_a_id] == 1)
					{
						col3_nUT_before_arr[sp_a_id][ix] = col3_nUT_arr[sp_a_id][ix];

						//flag for the collision of fM_a to fM_b with n_b, U_a, T_a 
						col3_f0_h0g0_col_flag_arr[sp_a_id][ix] = 0;

						col3_Qv_h0_ab_summed_arr[ix][sp_a_id].setZero();
						col3_Qv_g0_ab_summed_arr[ix][sp_a_id].setZero();
					}
				}


				for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
				{

					int valid_flag = 0;
					for(int sp_dum_id = 0; sp_dum_id < tot_species_num; sp_dum_id++)
					{
						if(col3_sp_ab_ix_col_flag_arr[ix][sp_id][sp_dum_id] == 1 || col3_sp_ab_ix_col_flag_arr[ix][sp_dum_id][sp_id] == 1)
							valid_flag = 1;
					}

					if(sp_kinetic_op_arr[sp_id] == 1 && valid_flag == 1)
					{
						sp_h0g0_SE_update(sp_id, ix);

						
					}
				}

				for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
				{
					int valid_flag = 0;
					for(int sp_dum_id = 0; sp_dum_id < tot_species_num; sp_dum_id++)
					{
						if(col3_sp_ab_ix_col_flag_arr[ix][sp_id][sp_dum_id] == 1 || col3_sp_ab_ix_col_flag_arr[ix][sp_dum_id][sp_id] == 1)
							valid_flag = 1;
					}

					if(sp_kinetic_op_arr[sp_id] == 1 && col3_hMbgMb_with_UaTa_op == 1 && valid_flag == 1)
					{
						sp_hMbgMb_with_UaTa_update(sp_id, ix);
					}

					int sp_a_id = sp_id;
					int nv_loc = nv[sp_a_id];
					int dof = dof_arr[sp_a_id];

					//SparseMatrix consv_fa_fb0_mat_sum = SparseMatrix(nv_loc*dof, nv_loc*dof);
					for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
					{
						int smaller_sp_id = col3_smaller_vspace_ab_arr[sp_a_id][sp_b_id][ix];

						//if(sp_a_id == smaller_sp_id)
						if(sp_a_id != sp_b_id && sp_kinetic_op_arr[sp_a_id] == 1 && sp_kinetic_op_arr[sp_b_id] == 1)
						{
							//consv_fa_fb0_mat_arr update
							int dof_v = dof_v_arr[sp_a_id];
							Vector h_ab = col3_Qv_h0_ab_arr[0][ix][sp_a_id][sp_b_id];
							h_ab.segment(0, 2*nv_loc*dof_v*dof_v) *= 1.0 + Ms_arr[sp_b_id]/Ms_arr[sp_a_id];
							Vector g_ab = col3_Qv_g0_ab_arr[0][ix][sp_a_id][sp_b_id];

							SparseMatrix M_C_fa_to_fb0;
							M_C_mat_fa_to_hg_ab(sp_a_id, ix, h_ab, g_ab, M_C_fa_to_fb0, -1);

							consv_fa_fb0_mat_arr[ix][sp_a_id][sp_b_id] = col_diag_cell_consv_mat[sp_a_id][ix]*M_C_fa_to_fb0;

							//consv_fa_fb0_UaTa_mat_arr update
							if(col3_hMbgMb_with_UaTa_op == 1 && valid_flag == 1)
							{
								consv_fa_fb0_UaTa_mat_arr[ix][sp_a_id][sp_b_id] = -col_diag_cell_consv_mat[sp_a_id][ix]*col3_f0_h0g0_col_arr[sp_a_id][sp_b_id][ix];
							}
						}
					}

					if(col3_hMbgMb_with_UaTa_op == 1 && sp_kinetic_op_arr[sp_a_id] == 1) 
					{
						col3_C_fMa_to_fMb_wUaTa_summed_arr[sp_a_id][ix] = col_sp_Minv_arr[sp_a_id][ix]*col3_M_C_fMa_to_fMb_wUaTa_summed_arr[sp_a_id][ix];
					}
				}
			}
		}

		if(rank == 0) col_dc1 += clock() - testc_1; 
	}

	if(rank == 0) 
	{
		int t_sum_period = min(1000, max(1, int(nstep/1000.0)));

		if(t_step == 0 || restart_op == 1)
		{
			col_dc1 = 0.0;
			col_dc2 = 0.0;
			col_dc2_1 = 0.0;
			col_dc2_2 = 0.0;
			col_dc2_3_1_1 = 0.0;
			col_dc2_3_1_2 = 0.0;
			col_dc2_3_2 = 0.0;
			col_dc2_4 = 0.0;
			col_dc2_5 = 0.0;  
			col_dc3 = 0.0;
			col_dc4 = 0.0;
			col_dc5 = 0.0;
			col_dc6 = 0.0;
			col_dc7 = 0.0;
			col_dc7_1 = 0.0;
		}

		value_type diag_ratio = double(nstep)/double(t_sum_period*500);
		int diag_ratio_int = max(1, int(diag_ratio));
		if(t_step % (t_sum_period*diag_ratio_int) == 0 && t_step > 0)
		{

			ofstream fout;
			string fname = string("summary_col_simulation_time.txt"); 
			fout.open(fname.c_str(), ios_base::out | ios_base::app);
			fout.precision(3);

			clock_t tot_col_time = col_dc1 + col_dc2;

			fout << "Time step : " << t_step << endl;
			fout << "CPU time used : time collision total for rank0 : " << double(tot_col_time)/CLOCKS_PER_SEC << " s" << endl;

			fout << "CPU time used : hMgM update : " << double(col_dc1)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc1)/double(tot_col_time)*100 << " % " << endl;
			fout << "CPU time used : main part : " << double(col_dc2)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc2)/double(tot_col_time)*100 << " % " << endl;
			fout << "CPU time used : main loop part1 : " << double(col_dc3)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc3)/double(tot_col_time)*100 << " % " << endl;
			fout << "CPU time used : col total in main part : " << double(col_dc5)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc5)/double(tot_col_time)*100 << " % " << endl;
			fout << "CPU time used : self h,g cal in main part : " << double(col_dc2_1)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc2_1)/double(tot_col_time)*100 << " % " << endl;
			fout << "CPU time used : cross h,g cal in main part : " << double(col_dc2_2)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc2_2)/double(tot_col_time)*100 << " % " << endl;
			fout << "CPU time used : Qv construction in main part : " << double(col_dc7)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc7)/double(tot_col_time)*100 << " % " << endl;
			fout << "CPU time used : Qv construction part 1 in main part : " << double(col_dc7_1)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc7_1)/double(tot_col_time)*100 << " % " << endl;
			fout << "CPU time used : Qv construction part 2 in main part : " << double(col_dc7_2)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc7_2)/double(tot_col_time)*100 << " % " << endl;
				fout << "CPU time used : MC cal for implicit iteration in main part : " << double(col_dc2_3_1_1)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc2_3_1_1)/double(tot_col_time)*100 << " % " << endl;
				fout << "CPU time used : iter matrix solve for implicit iteration in main part : " << double(col_dc2_3_1_2)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc2_3_1_2)/double(tot_col_time)*100 << " % " << endl;

			fout << "CPU time used : org part in main part : " << double(col_dc2_3_2)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc2_3_2)/double(tot_col_time)*100 << " % " << endl;
			fout << "CPU time used : consv part1 in main part : " << double(col_dc2_6)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc2_6)/double(tot_col_time)*100 << " % " << endl;
			fout << "CPU time used : consv part2 in main part : " << double(col_dc2_4)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc2_4)/double(tot_col_time)*100 << " % " << endl;
			fout << "CPU time used : Lf solve in main part : " << double(col_dc2_5)/CLOCKS_PER_SEC << " s, percentage : " << double(col_dc2_5)/double(tot_col_time)*100 << " % " << endl;
			fout << endl;
			fout.close();
		}

	}

}

//core function for the lowest order part of f update module
void Collision_dg::sp_h0g0_SE_update(const int sp_b_id, const int ix)
{

	//reset col3_fa_to_fMb_col_flag_arr 
	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
		col3_fa_to_fMb_col_flag_arr[ix][sp_a_id][sp_b_id] = 0;

	value_type R_cen = R_cen_cell[ix], Z_cen = Z_cen_cell[ix];
	Point4 p = Point2({R_cen, Z_cen}) & Point2({0, 0});

	Flux *flux_sp = flux_arr[sp_b_id];
	Vector B(3);
	B = flux_sp->get_B(p);
	value_type B_val = sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2));

	value_type n_b_val = col3_nUT_arr[sp_b_id][ix][0];
	value_type org_U_b_val = col3_nUT_arr[sp_b_id][ix][1];
	value_type org_vth_sq_b_val = col3_nUT_arr[sp_b_id][ix][2]/Ms_arr[sp_b_id];

	Vector tmp_avged_Q = Vector::Zero(10);
	tmp_avged_Q[0] = n_b_val;
	tmp_avged_Q[1] = org_U_b_val;
	tmp_avged_Q[2] = org_vth_sq_b_val;
	tmp_avged_Q[6] = B_val;

	if(n_b_val > 0.0) den_arr[sp_b_id][ix] = n_b_val;

	if(tmp_avged_Q[2] < 0.0)
	{
		cout << "negative T found at rank, sp_id, ix, den (MKS), T, R, Z, B : " << rank << " " << sp_b_id << " " << ix << " " << n_b_val*ni00_norm_mks << " " << tmp_avged_Q[2] << " " << R_cen << " " << Z_cen << " " << B_val << endl;
		return;
	}
	else if(tmp_avged_Q[0] < 0.0)
	{
		int loc_rgn = eq_reader->rgn_fn(R_cen, Z_cen);
		value_type loc_xx = sqrt(eq_reader->psi_ov_psix_interpol(R_cen, Z_cen));

		cout << "negative density found at rank, sp_id, ix, den (MKS), R, Z, B, xx, rgn : " << rank << " " << sp_b_id << " " << ix << " " << n_b_val*ni00_norm_mks << " " << R_cen << " " << Z_cen << " " << B_val << " " << loc_xx << " " << loc_rgn << endl;
		return;
	}
	else
	{
		for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
		{
			if(sp_kinetic_op_arr[sp_a_id] == 1 && sp_kinetic_op_arr[sp_b_id] == 1)
				col3_fa_to_fMb_col_flag_arr[ix][sp_a_id][sp_b_id] = 1;
		}
	}

	//f_Mb update
	int nv_loc = nv[sp_b_id];
	int dof = dof_arr[sp_b_id];
	Vector fM_ab_coeff_vol_int_for_hg = Vector::Zero(nv_loc*dof);
	for(int iv = 0; iv < nv_loc; iv++)
	{
		Vector tmp_sum_w_evol = integration_col_arr[sp_b_id]->col_fM_coeff_vol_cal(ElementX(ix), ElementV(iv), sp_b_id, tmp_avged_Q);
		fM_ab_coeff_vol_int_for_hg.segment(dof*iv, dof) += tmp_sum_w_evol;
	}
	col3_fM_coeff_for_hg_nx_arr[sp_b_id][ix] = col3_M_solver_arr[sp_b_id][ix]->solve(fM_ab_coeff_vol_int_for_hg);

	col3_f0_h0g0_col_flag_arr[sp_b_id][ix] = 0;


	int k1_max = 1;
	if(col3_hMbgMb_with_UaTa_op == 1) k1_max = 2;

	//h0, g0 update
	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++) 
	{
		if(sp_kinetic_op_arr[sp_a_id] == 1)
		{

			for(int k1 = 0; k1 < k1_max; k1++)
			{
				col3_Qv_h0_ab_arr[k1][ix][sp_a_id][sp_b_id].setZero();
				col3_Qv_g0_ab_arr[k1][ix][sp_a_id][sp_b_id].setZero();
			}

			int dof_v = dof_v_arr[sp_a_id];
			int nv_loc = nv[sp_a_id];
			int nve = nve_arr[sp_a_id];
			int dim_fhat_valid_v_basis = basis_arr[sp_a_id]->fhat_valid_vbasis_dim_out();

			value_type gamma_ab = col3_gamma_ab_arr[sp_a_id][sp_b_id][ix];
			value_type mass_ratio = Ms_arr[sp_a_id]/Ms_arr[sp_b_id];

			//redefine tmp_avged_Q for fMb & fMb_with Ua, Ta
			tmp_avged_Q[3 + 0] = col3_nUT_arr[sp_b_id][ix][0];
			tmp_avged_Q[3 + 1] = col3_nUT_arr[sp_a_id][ix][1];
			tmp_avged_Q[3 + 2] = col3_nUT_arr[sp_a_id][ix][2]/Ms_arr[sp_b_id];

			//redefine tmp_avged_Q for pitch-angle scattering with ei_pitch_no_v_in_nu_ei_op = 1
			tmp_avged_Q[7 + 0] = col3_nUT_arr[sp_a_id][ix][0];
			tmp_avged_Q[7 + 1] = col3_nUT_arr[sp_a_id][ix][1];
			tmp_avged_Q[7 + 2] = col3_nUT_arr[sp_a_id][ix][2]/Ms_arr[sp_a_id];

			int tot_Qv_h_size_per_ix_tmp = 0;
			int tot_Qv_g_size_per_ix_tmp = 0;

			//volume term part
			for(int iv = 0; iv < nv_loc; iv++)
			{
				Vector tmp_hMgM_to_Q_Vec = Vector::Zero(2*5*dof_v*dof_v);
				if(col_adjustable_quad_op == 1)
					integration_col_arr[sp_a_id]->col3_hMgM_to_Qv_vol_adj_n_cal(sp_b_id, ElementX(ix), ElementV(iv), tmp_avged_Q, col3_hMbgMb_with_UaTa_op, tmp_hMgM_to_Q_Vec, col3_adjustable_quad_n_arr[sp_a_id][sp_b_id][ix][iv], col_adjustable_quad_tor);
				else
					integration_col_arr[sp_a_id]->col3_hMgM_to_Qv_vol_cal(sp_b_id, ElementX(ix), ElementV(iv), tmp_avged_Q, col3_hMbgMb_with_UaTa_op, tmp_hMgM_to_Q_Vec);


				for(int k1 = 0; k1 < k1_max; k1++)
				{
					for(int j = 0; j < 5; j++)
					{
						for(int lvp = 0; lvp < dof_v; lvp++)
						{
							for(int lv = 0; lv < dof_v; lv++) 
							{
								int tmp_index = ((k1*5 + j)*dof_v + lvp)*dof_v + lv;
								if(tmp_hMgM_to_Q_Vec[tmp_index] != 0.0)
								{
									if(j < 2)
									{
										int row_index = tot_Qv_h_size_per_ix_tmp + ((j*nv_loc + iv)*dof_v + lvp)*dof_v + lv;
										col3_Qv_h0_ab_arr[k1][ix][sp_a_id][sp_b_id][row_index] = tmp_hMgM_to_Q_Vec[tmp_index]*gamma_ab*mass_ratio;
									}
									else if(j < 5)
									{
										int row_index = tot_Qv_g_size_per_ix_tmp + (((j-2)*nv_loc + iv)*dof_v + lvp)*dof_v + lv;
										col3_Qv_g0_ab_arr[k1][ix][sp_a_id][sp_b_id][row_index] = tmp_hMgM_to_Q_Vec[tmp_index]*gamma_ab;
									}
								}
							}
						}
					}
				}
			}

			tot_Qv_h_size_per_ix_tmp += 2*nv_loc*dof_v*dof_v;
			tot_Qv_g_size_per_ix_tmp += 3*nv_loc*dof_v*dof_v;


			//surface term part
			for(int i = 0; i < vspace_edge_num_of_type; i++)
			{
				int n_size = vspace_edge_list_arr[sp_a_id][i].size();
				for(int j = 0; j < n_size; j++)
				{
					int e = vspace_edge_list_arr[sp_a_id][i][j];
					int ele[2] = {vspace_edge_ele_arr[sp_a_id][2*e], vspace_edge_ele_arr[sp_a_id][2*e + 1]};

					Vector tmp_gM_S_Vmat_Vec = Vector::Zero(2*2*(dim_fhat_valid_v_basis*2*dof_v));
					if (i == 0 || i == 1)
					{
						Vector tmp_hM_S_Vmat_Vec = Vector::Zero(2*4*dof_v*dof_v);
						if(col_adjustable_quad_op == 1)
						{
							int loc_n_in = max(col3_adjustable_quad_n_arr[sp_a_id][sp_b_id][ix][ele[0]], col3_adjustable_quad_n_arr[sp_a_id][sp_b_id][ix][ele[1]]);
							integration_col_arr[sp_a_id]->col3_hMgM_S_Vmat_adj_n_cal(sp_b_id, ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), i, tmp_avged_Q, col3_hMbgMb_with_UaTa_op, tmp_hM_S_Vmat_Vec, tmp_gM_S_Vmat_Vec, loc_n_in);
						}
						else
							integration_col_arr[sp_a_id]->col3_hMgM_S_Vmat_cal(sp_b_id, ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), i, tmp_avged_Q, col3_hMbgMb_with_UaTa_op, tmp_hM_S_Vmat_Vec, tmp_gM_S_Vmat_Vec);

						for(int k1 = 0; k1 < k1_max; k1++)
						{
							//h part for mixed flux of upwind and LF
							for(int tmp_j = 0; tmp_j < 4; tmp_j++)
							{
								for(int lvp = 0; lvp < dof_v; lvp++)
								{
									for(int lv = 0; lv < dof_v; lv++) 
									{
										int tmp_index = ((k1*4 + tmp_j)*dof_v + lvp)*dof_v + lv;
										if(tmp_hM_S_Vmat_Vec[tmp_index] != 0.0)
										{
											int row_index = tot_Qv_h_size_per_ix_tmp + ((e*4 + tmp_j)*dof_v + lvp)*dof_v + lv;
											col3_Qv_h0_ab_arr[k1][ix][sp_a_id][sp_b_id][row_index] = tmp_hM_S_Vmat_Vec[tmp_index]*gamma_ab*mass_ratio;
										}
									}
								}
							}
						}
					}
					else
					{
						if(col_adjustable_quad_op == 1)
						{
							int loc_n_in;
							if(ele[0] > -1)
								loc_n_in = col3_adjustable_quad_n_arr[sp_a_id][sp_b_id][ix][ele[0]];
							else
								loc_n_in = col3_adjustable_quad_n_arr[sp_a_id][sp_b_id][ix][ele[1]];

							integration_col_arr[sp_a_id]->col3_gM_S_Vmat_bd_adj_n_cal(sp_b_id, ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), i, tmp_avged_Q, col3_hMbgMb_with_UaTa_op, tmp_gM_S_Vmat_Vec, loc_n_in);
						}
						else
							integration_col_arr[sp_a_id]->col3_gM_S_Vmat_bd_cal(sp_b_id, ElementX(ix), EdgeV(e), ElementV(ele[0]), ElementV(ele[1]), i, tmp_avged_Q, col3_hMbgMb_with_UaTa_op, tmp_gM_S_Vmat_Vec);
					}

					//g part
					for(int k1 = 0; k1 < k1_max; k1++)
					{
						for(int tmp_j = 0; tmp_j < 2; tmp_j++)
						{
							for(int row = 0; row < dim_fhat_valid_v_basis; row++) 
							{
								for(int i_ele = 0; i_ele < 2; i_ele++)
								{
									for(int lv = 0; lv < dof_v; lv++)
									{
										int tmp_index = (((k1*2 + tmp_j)*dim_fhat_valid_v_basis + row)*2 + i_ele)*dof_v + lv;

										if(tmp_gM_S_Vmat_Vec[tmp_index] != 0.0)
										{
											int row_index = tot_Qv_g_size_per_ix_tmp + (((e*2 +tmp_j)*dim_fhat_valid_v_basis + row)*2 + i_ele)*dof_v + lv;

											col3_Qv_g0_ab_arr[k1][ix][sp_a_id][sp_b_id][row_index] = tmp_gM_S_Vmat_Vec[tmp_index]*gamma_ab;

										}
									}
								}
							}
						}
					}
				}
			}

			tot_Qv_h_size_per_ix_tmp += nve*4*dof_v*dof_v;
			tot_Qv_g_size_per_ix_tmp += nve*2*(dim_fhat_valid_v_basis)*(2*dof_v);


			//flux from h at edge qd points
			for(int i = 0; i < vspace_edge_num_of_type; i++)
			{
				int n_size = vspace_edge_list_arr[sp_a_id][i].size();
				for(int j = 0; j < n_size; j++)
				{
					int e = vspace_edge_list_arr[sp_a_id][i][j];
					int ele[2] = {vspace_edge_ele_arr[sp_a_id][2*e], vspace_edge_ele_arr[sp_a_id][2*e + 1]};

					if (i == 0 || i == 1)
					{
						Vector tmp_hM_flux_qd = Vector::Zero(2*col_edge_flux_qd_num[sp_a_id]);
						integration_col_arr[sp_a_id]->col3_hM_flux_qd_cal(sp_b_id, ElementX(ix), EdgeV(e), i, col_edge_flux_qd_points[sp_a_id], tmp_avged_Q, col3_hMbgMb_with_UaTa_op, tmp_hM_flux_qd);

						for(int k1 = 0; k1 < k1_max; k1++) 
						{
							for(int tmp_j = 0; tmp_j < col_edge_flux_qd_num[sp_a_id]; tmp_j++)
							{
								int tmp_index = k1*col_edge_flux_qd_num[sp_a_id] + tmp_j;
								if(tmp_hM_flux_qd[tmp_index] != 0.0)
								{
									int row_index = tot_Qv_h_size_per_ix_tmp + e*col_edge_flux_qd_num[sp_a_id] + tmp_j;

									col3_Qv_h0_ab_arr[k1][ix][sp_a_id][sp_b_id][row_index] = tmp_hM_flux_qd[tmp_index]*gamma_ab*mass_ratio;
								}
							}
						}
					}
				}
			}



			Vector tmp_Qv_h_ab = col3_Qv_h0_ab_arr[0][ix][sp_a_id][sp_b_id];
			tmp_Qv_h_ab.segment(0, 2*nv_loc*dof_v*dof_v) *= 1.0 + Ms_arr[sp_b_id]/Ms_arr[sp_a_id];

			col3_Qv_h0_ab_summed_arr[ix][sp_a_id] += tmp_Qv_h_ab;
			col3_Qv_g0_ab_summed_arr[ix][sp_a_id] += col3_Qv_g0_ab_arr[0][ix][sp_a_id][sp_b_id];

		}
	}


	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++) 
	{
		if(sp_kinetic_op_arr[sp_a_id] == 1 && sp_a_id != sp_b_id && col_fast_to_slow_moment_op == 1)
		{
			int smaller_sp_id = col3_smaller_vspace_ab_arr[sp_a_id][sp_b_id][ix];

			if(sp_b_id == smaller_sp_id)
			{
				int num_of_bc_a = col3_hg_bc_points_vp_u_arr[sp_a_id][0].size(); 
				int hg_a_inner_tmp_n = col3_hg_inner_n_arr[sp_a_id];

				col_fast_to_slow_moment_hg_Vec[sp_a_id][sp_b_id][ix] = Vector::Zero(2*num_of_bc_a + hg_a_inner_tmp_n);

				int nv_a_loc = nv[sp_a_id]; 

				value_type U_para_b = org_U_b_val;
				for(int iv = 0; iv < nv_a_loc; iv++)
				{
					if(col_adjustable_quad_op == 1)
					{
						integration_col_arr[sp_a_id]->col3_slow_moment_to_fast_hg_adj_n_cal(ElementX(ix), ElementV(iv), sp_a_id, sp_b_id, U_para_b, num_of_bc_a, col_fast_to_slow_moment_hg_Vec[sp_a_id][sp_b_id][ix], col3_adjustable_quad_n_arr[sp_a_id][sp_b_id][ix][iv]);
					}
					else
						integration_col_arr[sp_a_id]->col3_slow_moment_to_fast_hg_cal(ElementX(ix), ElementV(iv), sp_a_id, sp_b_id, U_para_b, num_of_bc_a, col_fast_to_slow_moment_hg_Vec[sp_a_id][sp_b_id][ix]);

				}

				for(int k = 0; k < num_of_bc_a; k++)
				{
					value_type v_para = col3_hg_bc_points_vp_u_arr[sp_a_id][0][k];
					value_type u = col3_hg_bc_points_vp_u_arr[sp_a_id][1][k];

					value_type vp_mod = v_para - U_para_b;
					value_type vp_mod_sq = vp_mod*vp_mod;
					value_type u_sq = u*u*B_val;
					value_type v_sq = vp_mod_sq + u_sq;
					value_type v_3 = v_sq*sqrt(v_sq);
					value_type v_5 = v_sq*v_sq*sqrt(v_sq);

					value_type loc_hb_from_moment = (vp_mod_sq - 0.5*u_sq)/(4.0*M_PI*v_5);
					value_type loc_gb_from_moment = -(vp_mod_sq - 0.5*u_sq)/(24.0*M_PI*v_3);

					col_fast_to_slow_moment_hg_Vec[sp_a_id][sp_b_id][ix][k] = loc_hb_from_moment;
					col_fast_to_slow_moment_hg_Vec[sp_a_id][sp_b_id][ix][num_of_bc_a + k] = loc_gb_from_moment;
				}
			}
		}
	}




}

//Correction term for non-zero equilibrium collision effect
void Collision_dg::sp_hMbgMb_with_UaTa_update(const int sp_a_id, const int ix)
{

	//collsion
	int dof = dof_arr[sp_a_id];
	int dof_sq = dof*dof;
	int dof_x = dof_x_arr[sp_a_id];
	int dof_v = dof_v_arr[sp_a_id];
	int nv_loc = nv[sp_a_id];
	int nve = nve_arr[sp_a_id];
	int dim_fhat_valid_v_basis = basis_arr[sp_a_id]->fhat_valid_vbasis_dim_out();

	col3_M_C_fMa_to_fMb_wUaTa_summed_arr[sp_a_id][ix] = Vector::Zero(nv_loc*dof);

	for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
	{
		if(sp_kinetic_op_arr[sp_b_id] == 1)
		{
			col3_f0_h0g0_col_arr[sp_a_id][sp_b_id][ix] = Vector::Zero(nv_loc*dof);
		}
	}

	int max_sp_b_id = tot_species_num;
	//if(sp_a_id == 0 && col_ei_pitch_angle_op == 1) max_sp_b_id = 1;

	int tot_Qv_h_size_per_ix = nv_loc*2*dof_v*dof_v;
	int tot_Qv_h_size_vol_part_per_ix = nv_loc*2*dof_v*dof_v;
	tot_Qv_h_size_per_ix += nve*4*dof_v*dof_v;
	tot_Qv_h_size_per_ix += nve*col_edge_flux_qd_num[sp_a_id];

	//volume term part
	int tot_Qv_g_size_per_ix = nv_loc*3*dof_v*dof_v;
	//surface term part
	tot_Qv_g_size_per_ix += nve*2*(dim_fhat_valid_v_basis)*(2*dof_v);

	Vector Qv_h_ab_UaTa_sum = Vector::Zero(tot_Qv_h_size_per_ix); 
	Vector Qv_g_ab_UaTa_sum = Vector::Zero(tot_Qv_g_size_per_ix); 

	for(int sp_b_id = 0; sp_b_id < max_sp_b_id; sp_b_id++)
	{
		if(sp_kinetic_op_arr[sp_b_id] == 1)
		{

			int k1_ref = 1; //0 : f_a-> f_{M,b}(n_b, U_b, T_b), 1 : f_a-> f_{M,b}(n_b, U_a, T_a)
			Vector Qv_h_ab = col3_Qv_h0_ab_arr[1][ix][sp_a_id][sp_b_id];
			Vector Qv_g_ab = col3_Qv_g0_ab_arr[1][ix][sp_a_id][sp_b_id];

			Qv_h_ab_UaTa_sum += Qv_h_ab;
			Qv_h_ab_UaTa_sum.segment(0, tot_Qv_h_size_vol_part_per_ix) += (Ms_arr[sp_b_id]/Ms_arr[sp_a_id])*Qv_h_ab.segment(0, tot_Qv_h_size_vol_part_per_ix);

			Qv_g_ab_UaTa_sum += Qv_g_ab;

			Vector E_mat_col_1d = Vector::Zero(nv_loc*dof_sq);

			//Edge terms
			int Qv_h_index = nv_loc*2*dof_v*dof_v;
			int Qv_h_edge_flux_qd_index = nv_loc*2*dof_v*dof_v + nve*4*dof_v*dof_v;
			int h_edge_qd_num = col_edge_flux_qd_num[sp_a_id];

			int Qv_g_edge_index = nv_loc*3*dof_v*dof_v;

			Vector fhat_arr = col_f_dg_to_valid_fhat_mat_arr[sp_a_id][ix]*col3_fM_coeff_for_hg_nx_arr[sp_a_id][ix];

			vector<int> lp_elev_index(4), l_elev_index(4);
			vector<vector<value_type>> coeff_lp;
			coeff_lp.resize(4);

			for(int ej = 0; ej < 4; ej++)
			{
				coeff_lp[ej].resize(2);
			}

			//for self term calculation
			vector<int> S2_MEL_self_coeff_flag_arr;
			vector<vector<Vector>> S2_MEL_self_coeff_vpart_arr;

			S2_MEL_self_coeff_flag_arr.resize(nv_loc, 0);
			S2_MEL_self_coeff_vpart_arr.resize(nv_loc);
			for(int iv = 0; iv < nv_loc; iv++)
				S2_MEL_self_coeff_vpart_arr[iv].resize(2, Vector::Zero(dof_v*dof_v));

			for(int i = 0; i < vspace_edge_num_of_type; i++)
			{
				int n_size = vspace_edge_list_arr[sp_a_id][i].size();
				for(int j = 0; j < n_size; j++)
				{
					int e = vspace_edge_list_arr[sp_a_id][i][j];
					int ele[2] = {vspace_edge_ele_arr[sp_a_id][2*e], vspace_edge_ele_arr[sp_a_id][2*e + 1]};

					//h part
					if (i == 0 || i == 1)
					{
						Vector Qv_h_edge_flux_qd_arr = Qv_h_ab.segment(Qv_h_edge_flux_qd_index + e*h_edge_qd_num, h_edge_qd_num);

						value_type u_tot_max = Qv_h_edge_flux_qd_arr.maxCoeff();
						value_type u_tot_min = Qv_h_edge_flux_qd_arr.minCoeff();

						value_type x_c = 0.5*(u_tot_max + u_tot_min);
						value_type del_x = 0.5*(u_tot_max - u_tot_min);
						if(abs(del_x) < 1e-50) del_x = 1e-50;
						value_type w_plus, w_minus, w_LF_half;
						value_type lambda_LF = max(abs(u_tot_max), abs(u_tot_min));

						int edge_flux_type;
						if(x_c > del_x) //pure Upwind from - to +
						{
							edge_flux_type = 0;
							w_minus = 1.0;
							w_plus = 0.0;
						}
						else if(x_c < - del_x) //pure Upwind from + to -
						{
							edge_flux_type = 1;
							w_minus = 0.0;
							w_plus = 1.0;

						}
						else if(x_c > 0.0)
						{
							edge_flux_type = 2;
							w_minus = get_w_loc_mod(x_c/del_x);
							w_plus = 0.0;
						}
						else
						{
							edge_flux_type = 3;
							w_minus = 0.0;
							w_plus = get_w_loc_mod(x_c/del_x);
						}
						w_LF_half = (1.0 - w_minus - w_plus)*0.5;


						int kl = ele[0];
						int kr = ele[1];

						vector<Vector> f_sp_loc_lr(2);
						f_sp_loc_lr[0] = col3_fM_coeff_for_hg_nx_arr[sp_a_id][ix].segment(kl*dof, dof);
						f_sp_loc_lr[1] = col3_fM_coeff_for_hg_nx_arr[sp_a_id][ix].segment(kr*dof, dof);



						//lp : giver, l : receiver
						//case ej=0 : [lp = -, l = -] case : coeff_lp[0] used
						//case ej=1 : [lp = -, l = +] case : coeff_lp[1] used
						//case ej=2 : [lp = +, l = -] case : coeff_lp[2] used
						//case ej=3 : [lp = +, l = +] case : coeff_lp[3] used
						//

						for(int ej = 0; ej < 4; ej++)
						{
							if(ej == 0 || ej == 1)
							{
								coeff_lp[ej][0] = w_minus + w_LF_half;
								coeff_lp[ej][1] = w_LF_half*lambda_LF;

								lp_elev_index[ej] = kl;
							}
							else
							{
								coeff_lp[ej][0] = w_plus + w_LF_half;
								coeff_lp[ej][1] = -w_LF_half*lambda_LF;
								lp_elev_index[ej] = kr;
							}

							if(ej == 0 || ej == 2)
							{
								l_elev_index[ej] = kl;
							}
							else
							{
								l_elev_index[ej] = kr;
							}

							if(Q_S_op_flags_arr[edge_flux_type][ej][0] > 0)
							{
								int l_iv_index = l_elev_index[ej];
								if(ej == 0 || ej == 3) //self part
								{
									S2_MEL_self_coeff_flag_arr[l_iv_index] = 1;

									int Qv_h_edge_flux_mat_index = nv_loc*2*dof_v*dof_v + (e*4 + ej)*dof_v*dof_v;

									S2_MEL_self_coeff_vpart_arr[l_iv_index][i] += coeff_lp[ej][0]*Qv_h_ab.segment(Qv_h_edge_flux_mat_index, dof_v*dof_v);

									//U_LF part flag
									if(Q_S_op_flags_arr[edge_flux_type][ej][1] == 1)
									{
										S2_MEL_self_coeff_vpart_arr[l_iv_index][i] += coeff_lp[ej][1]*col3_h_LF_Vmat_arr[sp_a_id][ix*nve + e][ej];
									}
								}
								else //cross part
								{
									Vector U_tot_cross_mat_tmp = Vector::Zero(dof_sq);

									int Qv_h_edge_flux_mat_index = nv_loc*2*dof_v*dof_v + (e*4 + ej)*dof_v*dof_v;
									Vector Vmat_h = coeff_lp[ej][0]*Qv_h_ab.segment(Qv_h_edge_flux_mat_index, dof_v*dof_v);

									//U_LF part flag
									if(Q_S_op_flags_arr[edge_flux_type][ej][1] == 1)
									{
										Vmat_h += coeff_lp[ej][1]*col3_h_LF_Vmat_arr[sp_a_id][ix*nve + e][ej];
									}

									for(int lvp = 0; lvp < dof_v; lvp++)
									{
										for(int lv = 0; lv < dof_v; lv++)
										{
											value_type tmp_V_part = Vmat_h[lvp*dof_v + lv];

											if(tmp_V_part != 0.0)
											{
												for(int lxp = 0; lxp < dof_x; lxp++)
												{	
													int lp_index_mod = (lvp*dof_x + lxp)*dof;
													for(int lx = 0; lx < dof_x; lx++)
													{
														int l_index = lv*dof_x + lx;

														value_type tmp_X_part = col3_hg_Xmat_arr[sp_a_id][ix][i][lxp*dof_x + lx];
														U_tot_cross_mat_tmp(lp_index_mod + l_index) += tmp_X_part*tmp_V_part;
													}
												}
											}
										}
									}

									col3_f0_h0g0_col_arr[sp_a_id][sp_b_id][ix].segment(l_iv_index*dof, dof) -= U_tot_cross_mat_tmp.reshaped(dof, dof)*f_sp_loc_lr[S_lp_f_index_arr[ej]];

								}
							}
						}
					}


					//g part
					Matrix Q_g_2d_edge_loc_mat = Matrix::Zero(2*dof_v*dof_x, dof_x*dim_fhat_valid_v_basis);
					for(int tmp_j = 0; tmp_j < 2; tmp_j++)
					{
						int Qv_g_edge_index_loc = (e*2 + tmp_j)*(dim_fhat_valid_v_basis)*(2*dof_v);
						//Vector Qv_g_1d_edge_mat_loc = Qv_g_1d_edge_mat.segment(Qv_g_edge_index_loc, (dim_fhat_valid_v_basis)*(2*dof_v));

						Vector Qv_g_1d_edge_mat_loc = Qv_g_ab.segment(Qv_g_edge_index + Qv_g_edge_index_loc, (dim_fhat_valid_v_basis)*(2*dof_v));
						for(int lxp = 0; lxp < dof_x; lxp++)
						{	
							for(int lvp = 0; lvp < dim_fhat_valid_v_basis; lvp++)
							{
								int lp_index = lxp*dim_fhat_valid_v_basis + lvp;

								int S_mat_index_lvp = lvp*(2*dof_v);
								int S_mat_index1 = lp_index*(2*dof_v*dof_x);
								for(int lv = 0; lv < 2*dof_v; lv++)
								{
									value_type tmp_V_part = Qv_g_1d_edge_mat_loc[S_mat_index_lvp + lv];
									if(tmp_V_part != 0.0)
									{
										for(int lx = 0; lx < dof_x; lx++)
										{
											int l_index = lv*dof_x + lx;
											value_type tmp_X_part;

											if(i == 0 || i == 2) 
											{
												tmp_X_part = col3_hg_Xmat_arr[sp_a_id][ix][tmp_j][lxp*dof_x + lx];
											}
											else if(i == 1 || i == 3)
											{
												tmp_X_part = col3_hg_Xmat_arr[sp_a_id][ix][tmp_j + 1][lxp*dof_x + lx];
											}

											Q_g_2d_edge_loc_mat(l_index, lp_index) += tmp_X_part*tmp_V_part;
										}
									}
								}
							}
						}
					}

					Vector S2_g_1d_loc = Q_g_2d_edge_loc_mat*fhat_arr.segment(e*dof_x*dim_fhat_valid_v_basis, dof_x*dim_fhat_valid_v_basis);

					if (i == 0 || i == 1)
					{
						col3_f0_h0g0_col_arr[sp_a_id][sp_b_id][ix].segment(ele[0]*dof, dof) -= S2_g_1d_loc.segment(0, dof);
						col3_f0_h0g0_col_arr[sp_a_id][sp_b_id][ix].segment(ele[1]*dof, dof) -= S2_g_1d_loc.segment(dof, dof);
					}
					else
					{
						int valid_ele_index;
						if(ele[0] == -1) valid_ele_index = 1;
						else valid_ele_index = 0;

						col3_f0_h0g0_col_arr[sp_a_id][sp_b_id][ix].segment(ele[valid_ele_index]*dof, dof) -= S2_g_1d_loc.segment(valid_ele_index*dof, dof);
					}
				}
			}

			for(int iv = 0; iv < nv_loc; iv++)
			{
				int l_iv_index = iv;

				if(S2_MEL_self_coeff_flag_arr[l_iv_index] == 1)
				{
					Vector U_tot_self_mat_tmp = Vector::Zero(dof_sq);

					//U_ExB & LF mat construct
					for(int j_tmp = 0; j_tmp < 2; j_tmp++)
					{
						for(int lvp = 0; lvp < dof_v; lvp++)
						{
							for(int lv = 0; lv < dof_v; lv++)
							{
								value_type tmp_V_part = S2_MEL_self_coeff_vpart_arr[l_iv_index][j_tmp][lvp*dof_v + lv];
								if(tmp_V_part != 0.0)
								{
									for(int lxp = 0; lxp < dof_x; lxp++)
									{	
										int lp_index_mod = (lvp*dof_x + lxp)*dof;
										for(int lx = 0; lx < dof_x; lx++)
										{
											int l_index = lv*dof_x + lx;

											value_type tmp_X_part = col3_hg_Xmat_arr[sp_a_id][ix][j_tmp][lxp*dof_x + lx];

											U_tot_self_mat_tmp(lp_index_mod + l_index) += tmp_X_part*tmp_V_part;
										}
									}
								}
							}
						}

					}
					E_mat_col_1d.segment(iv*dof_sq, dof_sq) -= U_tot_self_mat_tmp;
				}
			}


			//Volume terms
			Vector Qv_1d_mat = Qv_g_ab.segment(0, 3*nv_loc*dof_v*dof_v);
			Qv_1d_mat.segment(0, 2*nv_loc*dof_v*dof_v) += (1.0 + Ms_arr[sp_b_id]/Ms_arr[sp_a_id])*Qv_h_ab.segment(0, 2*nv_loc*dof_v*dof_v);
			for(int j = 0; j < 3; j++)
			{
				for(int iv = 0; iv < nv_loc; iv++)
				{
					int E_mat_index0 = iv*dof_sq;
					for(int lvp = 0; lvp < dof_v; lvp++)
					{
						int E_mat_index_lvp = ((j*nv_loc + iv)*dof_v + lvp)*dof_v;

						for(int lxp = 0; lxp < dof_x; lxp++)
						{	
							int lp_index = lvp*dof_x + lxp;
							int E_mat_index1 = E_mat_index0 + lp_index*dof;
							for(int lv = 0; lv < dof_v; lv++)
							{
								value_type tmp_V_part = Qv_1d_mat[E_mat_index_lvp + lv];
								for(int lx = 0; lx < dof_x; lx++)
								{
									int l_index = lv*dof_x + lx;
									value_type tmp_X_part = col3_hg_Xmat_arr[sp_a_id][ix][j][lxp*dof_x + lx];

									E_mat_col_1d[E_mat_index1 + l_index] += tmp_X_part*tmp_V_part;
								}
							}
						}
					}
				}
			}

			for(int iv = 0; iv < nv_loc; iv++)
			{
				col3_f0_h0g0_col_arr[sp_a_id][sp_b_id][ix].segment(iv*dof, dof) += E_mat_col_1d.segment(iv*dof_sq, dof_sq).reshaped(dof, dof)*col3_fM_coeff_for_hg_nx_arr[sp_a_id][ix].segment(iv*dof, dof);
			}

		}
	}

	SparseMatrix M_C_fa_to_fb_UaTa;
	M_C_mat_fa_to_hg_ab(sp_a_id, ix, Qv_h_ab_UaTa_sum, Qv_g_ab_UaTa_sum, M_C_fa_to_fb_UaTa, -1);
	col3_M_C_fMa_to_fMb_wUaTa_summed_arr[sp_a_id][ix] = M_C_fa_to_fb_UaTa*col3_fM_coeff_for_hg_nx_arr[sp_a_id][ix];
}

//Big matrix construction for the collision operation
void Collision_dg::M_C_mat_fa_to_hg_ab(const int &sp_a_id, const int &ix, const Vector &Qv_h_ab, const Vector &Qv_g_ab, SparseMatrix &M_C_mat, int flag)
{
	//collsion
	int dof = dof_arr[sp_a_id];
	int dof_sq = dof*dof;
	int dof_x = dof_x_arr[sp_a_id];
	int dof_v = dof_v_arr[sp_a_id];
	int nv_loc = nv[sp_a_id];
	int nve = nve_arr[sp_a_id];
	int dim_fhat_valid_v_basis = basis_arr[sp_a_id]->fhat_valid_vbasis_dim_out();

	//1. Qv_h_ab, Qv_g_ab sum
	//2. iv loop
	//2.1. surf loop -> define 5 neighbor iv and make [dof*dof] size matrix for each iv
	//2.2. volume term + surf self term
	//2.3. mult M^(-1) for each neighbor : use M_inv from col_consv_mat_setup
	//2.4. update tripletList for iv cell

	Vector E_mat_col_1d = Vector::Zero(nv_loc*dof_sq);

	vector<SparseMatrix_Triplet> tripletList_M_C_cell;
	int tmp_cell_size = nv_loc*dof*(8*2*dof + dof);
	tripletList_M_C_cell.reserve(tmp_cell_size);

	//Edge terms
	int Qv_h_edge_flux_qd_index = nv_loc*2*dof_v*dof_v + nve*4*dof_v*dof_v;
	int h_edge_qd_num = col_edge_flux_qd_num[sp_a_id];

	int Qv_g_edge_index = nv_loc*3*dof_v*dof_v;

	vector<int> lp_elev_index(4), l_elev_index(4);
	vector<vector<value_type>> coeff_lp(4);

	for(int ej = 0; ej < 4; ej++)
		coeff_lp[ej].resize(2);

	//for self term calculation
	vector<int> S2_MEL_self_coeff_flag_arr;
	vector<vector<Vector>> S2_MEL_self_coeff_vpart_arr;

	S2_MEL_self_coeff_flag_arr.resize(nv_loc, 0);
	S2_MEL_self_coeff_vpart_arr.resize(nv_loc);
	for(int iv = 0; iv < nv_loc; iv++)
		S2_MEL_self_coeff_vpart_arr[iv].resize(2, Vector::Zero(dof_v*dof_v));

	//Edge terms
	for(int i = 0; i < vspace_edge_num_of_type; i++)
	{
		int n_size = vspace_edge_list_arr[sp_a_id][i].size();
		for(int j = 0; j < n_size; j++)
		{
			int e = vspace_edge_list_arr[sp_a_id][i][j];
			int ele[2] = {vspace_edge_ele_arr[sp_a_id][2*e], vspace_edge_ele_arr[sp_a_id][2*e + 1]};

			//h part : HLL flux
			if (i == 0 || i == 1)
			{
				Vector Qv_h_edge_flux_qd_arr = Qv_h_ab.segment(Qv_h_edge_flux_qd_index + e*h_edge_qd_num, h_edge_qd_num);

				value_type u_tot_max = Qv_h_edge_flux_qd_arr.maxCoeff();
				value_type u_tot_min = Qv_h_edge_flux_qd_arr.minCoeff();

				value_type x_c = 0.5*(u_tot_max + u_tot_min);
				value_type del_x = 0.5*(u_tot_max - u_tot_min);
				if(abs(del_x) < 1e-50) del_x = 1e-50;
				value_type w_plus, w_minus, w_LF_half;
				value_type lambda_LF = max(abs(u_tot_max), abs(u_tot_min));

				int edge_flux_type;
				if(x_c > del_x) //pure Upwind from - to +
				{
					edge_flux_type = 0;
					w_minus = 1.0;
					w_plus = 0.0;
				}
				else if(x_c < - del_x) //pure Upwind from + to -
				{
					edge_flux_type = 1;
					w_minus = 0.0;
					w_plus = 1.0;

				}
				else if(x_c > 0.0)
				{
					edge_flux_type = 2;
					w_minus = get_w_loc_mod(x_c/del_x);
					w_plus = 0.0;
				}
				else
				{
					edge_flux_type = 3;
					w_minus = 0.0;
					w_plus = get_w_loc_mod(x_c/del_x);
				}
				w_LF_half = (1.0 - w_minus - w_plus)*0.5;

				int kl = ele[0];
				int kr = ele[1];

				//lp : giver, l : receiver
				//case ej=0 : [lp = -, l = -] case : coeff_lp[0] used
				//case ej=1 : [lp = -, l = +] case : coeff_lp[1] used
				//case ej=2 : [lp = +, l = -] case : coeff_lp[2] used
				//case ej=3 : [lp = +, l = +] case : coeff_lp[3] used
				//

				for(int ej = 0; ej < 4; ej++)
				{
					if(ej == 0 || ej == 1)
					{
						coeff_lp[ej][0] = w_minus + w_LF_half;
						coeff_lp[ej][1] = w_LF_half*lambda_LF;

						lp_elev_index[ej] = kl;
					}
					else
					{
						coeff_lp[ej][0] = w_plus + w_LF_half;
						coeff_lp[ej][1] = -w_LF_half*lambda_LF;
						lp_elev_index[ej] = kr;
					}

					if(ej == 0 || ej == 2)
					{
						l_elev_index[ej] = kl;
					}
					else
					{
						l_elev_index[ej] = kr;
					}

					if(Q_S_op_flags_arr[edge_flux_type][ej][0] > 0)
					{
						int l_iv_index = l_elev_index[ej];
						if(ej == 0 || ej == 3) //self part
						{
							S2_MEL_self_coeff_flag_arr[l_iv_index] = 1;

							int Qv_h_edge_flux_mat_index = nv_loc*2*dof_v*dof_v + (e*4 + ej)*dof_v*dof_v;

							S2_MEL_self_coeff_vpart_arr[l_iv_index][i] += coeff_lp[ej][0]*Qv_h_ab.segment(Qv_h_edge_flux_mat_index, dof_v*dof_v);

							//U_LF part flag
							if(Q_S_op_flags_arr[edge_flux_type][ej][1] == 1)
							{
								S2_MEL_self_coeff_vpart_arr[l_iv_index][i] += coeff_lp[ej][1]*col3_h_LF_Vmat_arr[sp_a_id][ix*nve + e][ej];
							}
						}
						else //cross part
						{

							int Qv_h_edge_flux_mat_index = nv_loc*2*dof_v*dof_v + (e*4 + ej)*dof_v*dof_v;
							Vector Vmat_h = coeff_lp[ej][0]*Qv_h_ab.segment(Qv_h_edge_flux_mat_index, dof_v*dof_v);

							//U_LF part flag
							if(Q_S_op_flags_arr[edge_flux_type][ej][1] == 1)
							{
								Vmat_h += coeff_lp[ej][1]*col3_h_LF_Vmat_arr[sp_a_id][ix*nve + e][ej];
							}

							for(int lvp = 0; lvp < dof_v; lvp++)
							{
								for(int lv = 0; lv < dof_v; lv++)
								{
									value_type tmp_V_part = Vmat_h[lvp*dof_v + lv];

									if(tmp_V_part != 0.0)
									{
										for(int lxp = 0; lxp < dof_x; lxp++)
										{	
											int lp_index = lvp*dof_x + lxp;
											int lp_tot_index = lp_elev_index[ej]*dof + lp_index;
											for(int lx = 0; lx < dof_x; lx++)
											{
												int l_index = lv*dof_x + lx;
												int l_tot_index = l_elev_index[ej]*dof + l_index;

												value_type tmp_X_part = col3_hg_Xmat_arr[sp_a_id][ix][i][lxp*dof_x + lx];

												tripletList_M_C_cell.push_back(SparseMatrix_Triplet(l_tot_index, lp_tot_index, -tmp_X_part*tmp_V_part));

											}
										}
									}
								}
							}
						}
					}
				}
			}


			//g part
			Matrix Q_g_2d_edge_loc_mat = Matrix::Zero(2*dof_v*dof_x, dof_x*dim_fhat_valid_v_basis);
			for(int tmp_j = 0; tmp_j < 2; tmp_j++)
			{
				int Qv_g_edge_index_loc = (e*2 + tmp_j)*(dim_fhat_valid_v_basis)*(2*dof_v);

				Vector Qv_g_1d_edge_mat_loc = Qv_g_ab.segment(Qv_g_edge_index + Qv_g_edge_index_loc, (dim_fhat_valid_v_basis)*(2*dof_v));
				for(int lxp = 0; lxp < dof_x; lxp++)
				{	
					for(int lvp = 0; lvp < dim_fhat_valid_v_basis; lvp++)
					{
						int lp_index = lxp*dim_fhat_valid_v_basis + lvp;

						int S_mat_index_lvp = lvp*(2*dof_v);
						int S_mat_index1 = lp_index*(2*dof_v*dof_x);
						for(int lv = 0; lv < 2*dof_v; lv++)
						{
							value_type tmp_V_part = Qv_g_1d_edge_mat_loc[S_mat_index_lvp + lv];
							if(tmp_V_part != 0.0)
							{
								for(int lx = 0; lx < dof_x; lx++)
								{
									int l_index = lv*dof_x + lx;
									value_type tmp_X_part;

									if(i == 0 || i == 2) 
									{
										tmp_X_part = col3_hg_Xmat_arr[sp_a_id][ix][tmp_j][lxp*dof_x + lx];
									}
									else if(i == 1 || i == 3)
									{
										tmp_X_part = col3_hg_Xmat_arr[sp_a_id][ix][tmp_j + 1][lxp*dof_x + lx];
									}

									Q_g_2d_edge_loc_mat(l_index, lp_index) += tmp_X_part*tmp_V_part;
								}
							}
						}
					}
				}
			}

			Matrix S2_g_loc_mat = Q_g_2d_edge_loc_mat*col_f_dg_to_f_hat_mat_arr[sp_a_id][ix][e];

			if (i == 0 || i == 1)
			{
				for(int local_lp_ele_i = 0; local_lp_ele_i < 2; local_lp_ele_i++)
				{
					for(int lp_index = 0; lp_index < dof; lp_index++)
					{	
						int lp_tot_index = ele[local_lp_ele_i]*dof + lp_index;
						for(int local_l_ele_i = 0; local_l_ele_i < 2; local_l_ele_i++)
						{
							for(int l_index = 0; l_index < dof; l_index++)
							{	
								int l_tot_index = ele[local_l_ele_i]*dof + l_index;

								tripletList_M_C_cell.push_back(SparseMatrix_Triplet(l_tot_index, lp_tot_index, -S2_g_loc_mat(local_l_ele_i*dof + l_index ,local_lp_ele_i*dof + lp_index)));
							}
						}
					}
				}
			}
			else
			{
				int valid_ele_index;
				if(ele[0] == -1) valid_ele_index = 1;
				else valid_ele_index = 0;

				int local_lp_ele_i = valid_ele_index;
				int local_l_ele_i = valid_ele_index;

				for(int lp_index = 0; lp_index < dof; lp_index++)
				{	
					int lp_tot_index = ele[local_lp_ele_i]*dof + lp_index;
					for(int l_index = 0; l_index < dof; l_index++)
					{	
						int l_tot_index = ele[local_l_ele_i]*dof + l_index;

						tripletList_M_C_cell.push_back(SparseMatrix_Triplet(l_tot_index, lp_tot_index, -S2_g_loc_mat(local_l_ele_i*dof + l_index ,local_lp_ele_i*dof + lp_index)));
					}
				}
			}
		}
	}

	for(int iv = 0; iv < nv_loc; iv++)
	{
		int l_iv_index = iv;
		int E_mat_index0 = iv*dof_sq;

		if(S2_MEL_self_coeff_flag_arr[l_iv_index] == 1)
		{

			//U_ExB & LF mat construct
			for(int j_tmp = 0; j_tmp < 2; j_tmp++)
			{
				for(int lvp = 0; lvp < dof_v; lvp++)
				{
					for(int lv = 0; lv < dof_v; lv++)
					{
						value_type tmp_V_part = S2_MEL_self_coeff_vpart_arr[l_iv_index][j_tmp][lvp*dof_v + lv];
						if(tmp_V_part != 0.0)
						{
							for(int lxp = 0; lxp < dof_x; lxp++)
							{	
								int lp_index = lvp*dof_x + lxp;
								int E_mat_index1 = E_mat_index0 + lp_index*dof;

								for(int lx = 0; lx < dof_x; lx++)
								{
									int l_index = lv*dof_x + lx;

									value_type tmp_X_part = col3_hg_Xmat_arr[sp_a_id][ix][j_tmp][lxp*dof_x + lx];

									int lp_tot_index = iv*dof + lp_index;
									int l_tot_index = iv*dof + l_index;

									//Sum up edge terms
									E_mat_col_1d[E_mat_index1 + l_index] -= tmp_X_part*tmp_V_part;
								}
							}
						}
					}
				}

			}
		}
	}

	//Volume terms
	Vector Qv_1d_mat = Qv_g_ab.segment(0, 3*nv_loc*dof_v*dof_v);
	Qv_1d_mat.segment(0, 2*nv_loc*dof_v*dof_v) += Qv_h_ab.segment(0, 2*nv_loc*dof_v*dof_v);
	for(int j = 0; j < 3; j++)
	{
		for(int iv = 0; iv < nv_loc; iv++)
		{
			int E_mat_index0 = iv*dof_sq;
			for(int lvp = 0; lvp < dof_v; lvp++)
			{
				int E_mat_index_lvp = ((j*nv_loc + iv)*dof_v + lvp)*dof_v;

				for(int lxp = 0; lxp < dof_x; lxp++)
				{	
					int lp_index = lvp*dof_x + lxp;
					int lp_tot_index = iv*dof + lp_index;
					int E_mat_index1 = E_mat_index0 + lp_index*dof;
					for(int lv = 0; lv < dof_v; lv++)
					{
						value_type tmp_V_part = Qv_1d_mat[E_mat_index_lvp + lv];
						for(int lx = 0; lx < dof_x; lx++)
						{
							int l_index = lv*dof_x + lx;
							int l_tot_index = iv*dof + l_index;

							value_type tmp_X_part = col3_hg_Xmat_arr[sp_a_id][ix][j][lxp*dof_x + lx];

							//Add volume terms
							E_mat_col_1d[E_mat_index1 + l_index] += tmp_X_part*tmp_V_part;
						}
					}
				}
			}
		}
	}

	for(int iv = 0; iv < nv_loc; iv++)
	{
		int E_mat_index0 = iv*dof_sq;
		for(int lp_index = 0; lp_index < dof; lp_index++)
		{	
			int lp_tot_index = iv*dof + lp_index;
			int E_mat_index1 = E_mat_index0 + lp_index*dof;
			for(int l_index = 0; l_index < dof; l_index++)
			{
				int l_tot_index = iv*dof + l_index;
				tripletList_M_C_cell.push_back(SparseMatrix_Triplet(l_tot_index, lp_tot_index, E_mat_col_1d[E_mat_index1 + l_index]));

			}
		}
	}

	M_C_mat = SparseMatrix(nv_loc*dof, nv_loc*dof);
	M_C_mat.setFromTriplets(tripletList_M_C_cell.begin(), tripletList_M_C_cell.end());

}

//implicit time scheme collision module
void Collision_dg::RK_implicit_col(vector<vector<value_type>> &_coeff_tot, const int &t_step) 
{
	if ((t_step)%col_ion_period == 0)
	{
		clock_t testc_1, testc_2;
		if(rank == 0) testc_1 = clock();

		vector<Vector> f[4], Lf_arr(tot_species_num);

		for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
		{
			if(sp_kinetic_op_arr[sp_id] == 1)
				Lf_arr[sp_id] = Vector::Zero(f_data_size_per_ix[sp_id]);
			else
				Lf_arr[sp_id] = Vector::Zero(1);

			for(int rk_sub_step = 0; rk_sub_step < 4; rk_sub_step++)
			{
				if(sp_kinetic_op_arr[sp_id] == 1)
					f[rk_sub_step].push_back(Vector::Zero(f_data_size_per_ix[sp_id]));
				else
					f[rk_sub_step].push_back(Vector::Zero(1));
			}
		}

		int ix_max = nx;


		for(int ix = 0; ix < ix_max; ix++)
		{
			for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
			{
				if(sp_kinetic_op_arr[sp_id] == 1)
				{
					//copy f data from _coeff_tot to f for ix element
					copy(_coeff_tot[sp_id].begin() + ix*f_data_size_per_ix[sp_id], _coeff_tot[sp_id].begin() + (ix + 1)*f_data_size_per_ix[sp_id], f[3][sp_id].data());

					value_type n_val, U_val, T_val;
					Vector tmp_avg = col_diag_cell_mat[sp_id][ix]*f[3][sp_id];

					col3_nUT_arr[sp_id][ix][0] = n_val = tmp_avg[0]/vol_cell[ix];
					col3_nUT_arr[sp_id][ix][1] = U_val = tmp_avg[1]/tmp_avg[0];
					col3_nUT_arr[sp_id][ix][2] = T_val = Ms_arr[sp_id]*(tmp_avg[2] - 0.5*U_val*U_val*tmp_avg[0])/(1.5*tmp_avg[0]);

					value_type T_para_eV = Ms_arr[sp_id]*(tmp_avg[5] - 0.5*U_val*U_val*tmp_avg[0])/(1.5*tmp_avg[0])*col3_norm_T0_eV/3.0;
					value_type T_perp_eV = Ms_arr[sp_id]*tmp_avg[6]/(1.5*tmp_avg[0])*col3_norm_T0_eV/3.0;

					if(n_val < 0.0 || T_para_eV < col3_T_eV_min_lim/3.0 || T_perp_eV < col3_T_eV_min_lim*2.0/3.0)
					{
						for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
						{
							col3_sp_ab_ix_col_flag_arr[ix][sp_id][sp_b_id] = 0;
							col3_sp_ab_ix_col_flag_arr[ix][sp_b_id][sp_id] = 0;
						}
					}
				}
			}

			int n_iter = 1;

			value_type del_t_in = double(col_ion_period)*system_dt;
			value_type sub_dt = del_t_in/double(n_iter);

			if(col3_valid_ix_flag_arr[ix] == 1 && n_iter > 0)
			{
				if(rank == 0) testc_2 = clock();

				RK_implicit_col_single(f[3], ix, t_step, 0, 0, sub_dt, Lf_arr);

				for(int sp_id = 0; sp_id < tot_species_num; sp_id++)
				{
					int valid_flag = 0;
					for(int sp_dum_id = 0; sp_dum_id < tot_species_num; sp_dum_id++)
					{
						if(col3_sp_ab_ix_col_flag_arr[ix][sp_id][sp_dum_id] == 1 || col3_sp_ab_ix_col_flag_arr[ix][sp_dum_id][sp_id] == 1)
							valid_flag = 1;
					}

					if(sp_kinetic_op_arr[sp_id] == 1 && valid_flag == 1)
					{
						f[3][sp_id] += sub_dt*Lf_arr[sp_id];
						//copy f data to _coeff_tot from f for ix element
						copy(f[3][sp_id].data(), f[3][sp_id].data() + f_data_size_per_ix[sp_id],_coeff_tot[sp_id].begin() + ix*f_data_size_per_ix[sp_id]);
					}
				}

				if(rank == 0) col_dc3 += clock() - testc_2;

			}
		}

		if(rank == 0) col_dc2 += clock() - testc_1;
	}
}

//core function for the collision operation
void Collision_dg::RK_implicit_col_single(const vector<Vector> &coeff_arr, const int &ix, const int &t_step, const int &iter_step, const int &sub_step, const value_type &sub_dt, vector<Vector> &Lf_out_arr) 
{
	clock_t testc_3, testc_4, testc_5, testc_6, testc_7;
	if(rank == 0) testc_5 = clock();

	vector<Vector> local_coeff_arr, local_coeff_prev_arr, rhs_implicit_arr(tot_species_num);
	vector<vector<Vector>> evol_coeff_arr(tot_species_num);
	vector<vector<Vector>> gamma_aa_h_arr, gamma_aa_g_arr, gamma_aa_h_init_arr, gamma_aa_g_init_arr;
	local_coeff_arr.resize(tot_species_num);
	local_coeff_prev_arr.resize(tot_species_num);
	gamma_aa_h_arr.resize(tot_species_num);
	gamma_aa_g_arr.resize(tot_species_num);
	gamma_aa_h_init_arr.resize(tot_species_num);
	gamma_aa_g_init_arr.resize(tot_species_num);

	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
	{
		evol_coeff_arr[sp_a_id].resize(tot_species_num);
		gamma_aa_h_arr[sp_a_id].resize(tot_species_num);
		gamma_aa_g_arr[sp_a_id].resize(tot_species_num);
		gamma_aa_h_init_arr[sp_a_id].resize(tot_species_num);
		gamma_aa_g_init_arr[sp_a_id].resize(tot_species_num);

		//local_coeff_arr[sp_a_id] : DG coefficients of f_a for species a
		local_coeff_arr[sp_a_id] = coeff_arr[sp_a_id];
	}

	vector<Vector> gamma_aa_g_sum_arr;
	gamma_aa_g_sum_arr.resize(tot_species_num);
	vector<Vector> gamma_aa_h_whole_sum_arr;
	gamma_aa_h_whole_sum_arr.resize(tot_species_num);

	vector<Vector> coeff_hg_source_arr(tot_species_num);

	int flag_under_tor = 0;
	int iter_num = 0;
	int iter_min_num = 2;
	int iter_max_num = 40;
	//tolerance lvel for implicit cycle
	value_type tor_level = col_implicit_tor;

	//implicit collision iteration
	while (((iter_num < iter_min_num) || (flag_under_tor == 0) ) && iter_num < iter_max_num)
	{
		if(rank == 0) testc_3 = clock();

		//self-collision h, g calculation
		for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
		{
			int sp_b_id = sp_a_id;

			if(sp_kinetic_op_arr[sp_a_id] == 1)
			{
				//coeff_hg_source_arr[sp_a_id] : DG coefficients for [delta f_a = f_a - f_M]
				coeff_hg_source_arr[sp_a_id] = local_coeff_arr[sp_a_id] - col3_fM_coeff_for_hg_nx_arr[sp_a_id][ix];

				int col3_vp_n = col3_vp_n_arr[sp_a_id];
				int col3_u_n = col3_u_n_arr[sp_a_id];

				int col3_hg_tot_n = col3_hg_tot_n_arr[sp_a_id];
				int col3_hg_inner_n = col3_hg_inner_n_arr[sp_a_id];

				Vector tot_h_source = Vector::Zero(col3_hg_tot_n);
				int base_h_row_index = col3_vp_n + 3;
				tot_h_source.segment(base_h_row_index, col3_hg_inner_n) = col3_f_to_h_source_mat[sp_a_id][ix]*coeff_hg_source_arr[sp_a_id];
				base_h_row_index += col3_hg_inner_n;
				int col3_h_bc_n = col3_h_bc_n_arr[sp_a_id];
				tot_h_source.segment(base_h_row_index, col3_h_bc_n) = col3_f_to_h_bc_mat[sp_a_id][ix]*coeff_hg_source_arr[sp_a_id];
				base_h_row_index += col3_h_bc_n;


				//Rosenbluth potential h for self-collision
				gamma_aa_h_arr[sp_a_id][sp_b_id] = col3_h_solver_stiffness[sp_a_id][ix]->solve(tot_h_source);

				Vector tot_g_source = Vector::Zero(col3_hg_tot_n);
				int base_g_row_index = col3_vp_n + 3;
				tot_g_source.segment(base_g_row_index, col3_hg_inner_n) = col3_h_to_g_source_mat[sp_a_id][ix]*gamma_aa_h_arr[sp_a_id][sp_b_id];
				base_g_row_index += col3_hg_inner_n;
				int col3_g_bc_n = col3_g_bc_n_arr[sp_a_id];
				tot_g_source.segment(base_g_row_index, col3_g_bc_n) = col3_f_to_g_bc_mat[sp_a_id][ix]*coeff_hg_source_arr[sp_a_id];
				base_g_row_index += col3_g_bc_n;

				//Rosenbluth potential g for self-collision
				gamma_aa_g_arr[sp_a_id][sp_b_id] = col3_g_solver_stiffness[sp_a_id][ix]->solve(tot_g_source);

				gamma_aa_h_whole_sum_arr[sp_a_id] = Vector::Zero(2*col3_hg_tot_n);
				gamma_aa_g_sum_arr[sp_a_id] = Vector::Zero(col3_hg_tot_n);
			}
		}

		if(rank == 0) 
		{
			testc_4 = clock();
			col_dc2_1 += testc_4 - testc_3;
			testc_3 = testc_4;
		}

		//up to this point, 
		//gamma_aa_h_arr[sp_a_id][sp_a_id] is delta_h_a
		//gamma_aa_g_arr[sp_a_id][sp_a_id] is delta_g_a

		//inter ion species del_h, del_g calculation
		for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
		{
			int col_fast_sp_a_to_slow_sp_b_moment_flag = 0;
			value_type col_fast_sp_a_to_slow_sp_b_moment_val = 0.0;
			for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
			{
				if(sp_a_id != sp_b_id)
				{
					if(sp_kinetic_op_arr[sp_a_id] == 1 && sp_kinetic_op_arr[sp_b_id] == 1 && col3_diff_sp_nl_col_op == 1)
					{
						value_type gamma_aa = col3_gamma_ab_arr[sp_a_id][sp_a_id][ix];
						value_type gamma_ab = col3_gamma_ab_arr[sp_a_id][sp_b_id][ix];
						value_type gamma_ratio = gamma_ab/gamma_aa;
						value_type mass_ratio = Ms_arr[sp_a_id]/Ms_arr[sp_b_id];

						//smaller_sp_id : which has smaller v space size between a & b
						int smaller_sp_id = col3_smaller_vspace_ab_arr[sp_a_id][sp_b_id][ix];

						int hg_a_inner_tmp_n = col3_hg_inner_n_arr[sp_a_id];
						int hg_a_total_tmp_n = col3_hg_tot_n_arr[sp_a_id];

						Vector tot_h_ba_source_tot = Vector::Zero(hg_a_total_tmp_n);
						Vector tot_g_ba_source_tot = Vector::Zero(hg_a_total_tmp_n);


						if(sp_a_id == smaller_sp_id)
						{
							//when a has smaller v space size
							//source & bc for hg_ab is interpolated from hg_bb

							//h,g volume source
							int base_h_row_index = col3_vp_n_arr[sp_a_id] + 3; //# of Neumann b.c. [dh/du = 0] points at u=0

							tot_h_ba_source_tot.segment(base_h_row_index, hg_a_inner_tmp_n) = col3_hg_b_to_hg_a_source_arr[sp_a_id][sp_b_id][ix]*gamma_aa_h_arr[sp_b_id][sp_b_id];
							tot_g_ba_source_tot.segment(base_h_row_index, hg_a_inner_tmp_n) = col3_hg_b_to_hg_a_source_arr[sp_a_id][sp_b_id][ix]*gamma_aa_g_arr[sp_b_id][sp_b_id];
							base_h_row_index += hg_a_inner_tmp_n; //# of volume source data

							//h,g b.c.
							int num_of_bc_a = col3_hg_bc_points_vp_u_arr[sp_a_id][0].size(); 
							tot_h_ba_source_tot.segment(base_h_row_index, num_of_bc_a) = col3_hg_b_to_hg_a_bc_arr[sp_a_id][sp_b_id][ix]*gamma_aa_h_arr[sp_b_id][sp_b_id];
							tot_g_ba_source_tot.segment(base_h_row_index, num_of_bc_a) = col3_hg_b_to_hg_a_bc_arr[sp_a_id][sp_b_id][ix]*gamma_aa_g_arr[sp_b_id][sp_b_id];

							gamma_aa_h_arr[sp_a_id][sp_b_id] = gamma_ab*mass_ratio*col3_hg_formula_to_hg0_solver[sp_a_id][ix]->solve(tot_h_ba_source_tot);
							gamma_aa_g_arr[sp_a_id][sp_b_id] = gamma_ab*col3_hg_formula_to_hg0_solver[sp_a_id][ix]->solve(tot_g_ba_source_tot);

						}
						else
						{
							//when a has bigger v space size
							//source & bc for hg_ab is interpolated from f_b

							int base_h_row_index = col3_vp_n_arr[sp_a_id] + 3; //# of Neumann b.c. [dh/du = 0] points at u=0
							int num_of_bc_a = col3_hg_bc_points_vp_u_arr[sp_a_id][0].size(); 


							Vector tmp_f_to_hg_bc_and_source = Vector::Zero(2*num_of_bc_a + hg_a_inner_tmp_n);

							Vector tmp_del_h;

							if(col_fast_to_slow_moment_op == 1 && col3_vspace_domain_ratio_ab_arr[sp_a_id][sp_b_id][ix] > col_fast_to_slow_moment_ratio_max) //Moments of slow species are used
							{
								Vector tmp_avg = col_diag_cell_mat[sp_b_id][ix]*coeff_arr[sp_b_id];
								value_type Ub_loc_sq_n = 0.5*tmp_avg[1]*tmp_avg[1]/tmp_avg[0];
								col_fast_sp_a_to_slow_sp_b_moment_val = (2.0*(tmp_avg[5] - Ub_loc_sq_n) - tmp_avg[6])/vol_cell[ix];

								tmp_f_to_hg_bc_and_source.segment(2*num_of_bc_a, hg_a_inner_tmp_n) = col3_hg_b_to_hg_a_source_arr[sp_a_id][sp_b_id][ix]*gamma_aa_h_arr[sp_b_id][sp_b_id];
								tmp_f_to_hg_bc_and_source += col_fast_to_slow_moment_hg_Vec[sp_a_id][sp_b_id][ix]*col_fast_sp_a_to_slow_sp_b_moment_val; 

								//modify here if needed for the correction of h_b boundary condition.. maybe not necessary
								// h_(ab)->0 anyway if v_a/v_(th,b) is big enough
								//1. calculate 1 point h_b value from coeff_hg_source_arr[sp_b_id] at (v_para, u) = (0, max u_a) 
								//2. calculate 1 point h_b value from col_fast_sp_a_to_slow_sp_b_moment_val at (v_para, u) = (0, max u_a)
								//3. find del(h_b)
								//4. add tmp_f_to_hg_bc_and_source += col_fast_to_slow_moment_del_h_Vec[sp_a_id][sp_b_id][ix]*[del(h_b)]

								tot_h_ba_source_tot.segment(base_h_row_index, hg_a_inner_tmp_n) = tmp_f_to_hg_bc_and_source.segment(2*num_of_bc_a, hg_a_inner_tmp_n); //volume source for h_ab from f_b
								base_h_row_index += hg_a_inner_tmp_n; //# of volume source data

								tot_h_ba_source_tot.segment(base_h_row_index, num_of_bc_a) = tmp_f_to_hg_bc_and_source.segment(0, num_of_bc_a); //b.c. for h_ab from f_b

								tmp_del_h = col3_hg_formula_to_hg0_solver[sp_a_id][ix]->solve(tot_h_ba_source_tot);
							}
							else //Brute force approach without moment information of slow species
							{
								tmp_f_to_hg_bc_and_source = col3_hg_b_to_hg_a_bc_arr[sp_a_id][sp_b_id][ix]*coeff_hg_source_arr[sp_b_id];

								tot_h_ba_source_tot.segment(base_h_row_index, hg_a_inner_tmp_n) = tmp_f_to_hg_bc_and_source.segment(2*num_of_bc_a, hg_a_inner_tmp_n); //volume source for h_ab from f_b
								base_h_row_index += hg_a_inner_tmp_n; //# of volume source data

								tot_h_ba_source_tot.segment(base_h_row_index, num_of_bc_a) = tmp_f_to_hg_bc_and_source.segment(0, num_of_bc_a); //b.c. for h_ab from f_b

								tmp_del_h = col3_h_solver_stiffness[sp_a_id][ix]->solve(tot_h_ba_source_tot);
							}

							base_h_row_index = col3_vp_n_arr[sp_a_id] + 3; //# of Neumann b.c. [dh/du = 0] points at u=0

							tot_g_ba_source_tot.segment(base_h_row_index, hg_a_inner_tmp_n) = col3_h_to_g_source_mat[sp_a_id][ix]*tmp_del_h; //volume source for g_ab from h_b
							base_h_row_index += hg_a_inner_tmp_n; //# of volume source data
							tot_g_ba_source_tot.segment(base_h_row_index, num_of_bc_a) = tmp_f_to_hg_bc_and_source.segment(num_of_bc_a, num_of_bc_a); //b.c. for g_ab from f_b

							gamma_aa_h_arr[sp_a_id][sp_b_id] = gamma_ab*mass_ratio*tmp_del_h;
							gamma_aa_g_arr[sp_a_id][sp_b_id] = gamma_ab*col3_g_solver_stiffness[sp_a_id][ix]->solve(tot_g_ba_source_tot);
						}
					}
				}
			}
		}

		if(rank == 0) 
		{
			testc_4 = clock();
			col_dc2_2 += testc_4 - testc_3;
			testc_3 = testc_4;
			testc_6 = testc_4;
		}

		//up to this point, 
		//a == b
		//gamma_aa_h_arr[sp_a_id][sp_b_id] is delta_h_b
		//gamma_aa_g_arr[sp_a_id][sp_b_id] is delta_g_b
		//a !=b
		//gamma_aa_h_arr[sp_a_id][sp_b_id] is gamma_ab*m_a/m_b*delta_h_b
		//gamma_aa_g_arr[sp_a_id][sp_b_id] is gamma_ab*delta_g_b

		vector<Vector> summed_Qv_h_ab(tot_species_num), summed_Qv_g_ab(tot_species_num);
		for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
		{
			int dof_v = dof_v_arr[sp_a_id];
			int nv_loc = nv[sp_a_id];
			int nve = nve_arr[sp_a_id];

			int hg_a_total_tmp_n = col3_hg_tot_n_arr[sp_a_id];
			value_type gamma_aa = col3_gamma_ab_arr[sp_a_id][sp_a_id][ix];

			if(rank == 0) 
			{
				testc_6 = clock();
			}

			//sum all of (delta h_ab, delta g_ab) for fixed species a
			for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
			{
				if(sp_kinetic_op_arr[sp_a_id] == 1 && sp_kinetic_op_arr[sp_b_id] == 1)
				{
					if(sp_a_id == sp_b_id)
					{
						gamma_aa_h_arr[sp_a_id][sp_b_id] *= gamma_aa;
						gamma_aa_g_arr[sp_a_id][sp_b_id] *= gamma_aa;
					}
					else if(iter_num == 0)
					{
						gamma_aa_h_init_arr[sp_a_id][sp_b_id] = gamma_aa_h_arr[sp_a_id][sp_b_id];
						gamma_aa_g_init_arr[sp_a_id][sp_b_id] = gamma_aa_g_arr[sp_a_id][sp_b_id];
					}

					gamma_aa_h_whole_sum_arr[sp_a_id].segment(0, hg_a_total_tmp_n)  += (1.0 + Ms_arr[sp_b_id]/Ms_arr[sp_a_id])*gamma_aa_h_arr[sp_a_id][sp_b_id];
					gamma_aa_h_whole_sum_arr[sp_a_id].segment(hg_a_total_tmp_n, hg_a_total_tmp_n) += gamma_aa_h_arr[sp_a_id][sp_b_id];
					gamma_aa_g_sum_arr[sp_a_id] += gamma_aa_g_arr[sp_a_id][sp_b_id];
				}
			}

			if(rank == 0) 
			{
				testc_7 = clock();
				col_dc7_1 += testc_7 - testc_6;
				testc_6 = testc_7;
			}

			//add lowest part of (h, g) to (delta h, delta g)
			if(sp_kinetic_op_arr[sp_a_id] == 1)
			{
				summed_Qv_h_ab[sp_a_id] = col3_Qv_h0_ab_summed_arr[ix][sp_a_id] + col3_h_to_Qv_mat_whole[ix][sp_a_id]*gamma_aa_h_whole_sum_arr[sp_a_id];
				summed_Qv_g_ab[sp_a_id] = col3_Qv_g0_ab_summed_arr[ix][sp_a_id] + col3_g_to_Qv_mat[ix][sp_a_id]*gamma_aa_g_sum_arr[sp_a_id];
			}

			if(rank == 0) 
			{
				testc_7 = clock();
				col_dc7_2 += testc_7 - testc_6;
				testc_6 = testc_7;
			}

		}

		if(rank == 0) 
		{
			testc_4 = clock();
			col_dc7 += testc_4 - testc_3;
			testc_3 = testc_4;
		}

		int max_err_sp;
		value_type max_err_val = 0.0, max_err1;

		flag_under_tor = 1;
			
		//the trapezoidal method with fixed point iterations for each species
		for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
		{
			if(sp_kinetic_op_arr[sp_a_id] == 1)
			{
				int sp_b_id = sp_a_id;

				if(rank == 0) testc_3 = clock();

				SparseMatrix M_C_fa_to_fb;
				//Generate M_C_fa_to_fb which is [M C(f_a, f_b) / f_a] matrix
				M_C_mat_fa_to_hg_ab(sp_a_id, ix, summed_Qv_h_ab[sp_a_id], summed_Qv_g_ab[sp_a_id], M_C_fa_to_fb, 1);

				SparseMatrix C_fa_to_fb = col_sp_Minv_arr[sp_a_id][ix]*M_C_fa_to_fb;

				if(rank == 0) 
				{
					testc_4 = clock();
					col_dc2_3_1_1 += testc_4 - testc_3;
					testc_3 = testc_4;
				}

				//SparseMatrix 
				value_type loc_tot_dt = sub_dt;
				value_type loc_half_dt = 0.5*sub_dt;

				int dof = dof_arr[sp_a_id];
				int nv_loc = nv[sp_a_id];

				vector<SparseMatrix_Triplet> tripletList_identity;

				for(int i = 0; i < nv_loc*dof; i++)
					tripletList_identity.push_back(SparseMatrix_Triplet(i, i, 1.0));

				//Identity matrix construction
				SparseMatrix Mat_identity = SparseMatrix(nv_loc*dof, nv_loc*dof);
				Mat_identity.setFromTriplets(tripletList_identity.begin(), tripletList_identity.end());

				//LFS matrix construction
				SparseMatrix im_loc_mat = Mat_identity - loc_half_dt*C_fa_to_fb;

				//Solver construction
				Eigen::PardisoLU<SparseMatrix> col_im_loc_solver;
				col_im_loc_solver.compute(im_loc_mat);

				//RHS Vector construction 
				if(iter_num == 0)
				{
					rhs_implicit_arr[sp_a_id] = coeff_arr[sp_a_id] + loc_half_dt*(C_fa_to_fb*coeff_arr[sp_a_id]);
					//apply compensation with -C[f_a(n_a, U_a, T_a), f_b(n_b, U_a, T_a)]
					if(col3_hMbgMb_with_UaTa_op == 1)
					{
						rhs_implicit_arr[sp_a_id] -= col3_C_fMa_to_fMb_wUaTa_summed_arr[sp_a_id][ix]*loc_tot_dt;
					}
				}

				local_coeff_prev_arr[sp_a_id] = local_coeff_arr[sp_a_id];
				//Update DG coefficient 
				local_coeff_arr[sp_a_id] = col_im_loc_solver.solve(rhs_implicit_arr[sp_a_id]);

				//Check the difference between the previous solution and the current one
				Vector del_coeff_diag = local_coeff_arr[sp_a_id] - local_coeff_prev_arr[sp_a_id];
				Vector tmp_avg = col_diag_cell_mat[sp_a_id][ix]*del_coeff_diag.cwiseAbs();

				value_type err_val = tmp_avg[0]/vol_cell[ix]/col3_nUT_arr[sp_a_id][ix][0];

				//If err_val > tolerance level, another iteration is needed. 
				if (err_val > tor_level || err_val < 0.0) flag_under_tor = 0;

				if(err_val > max_err_val)
				{
					max_err_val = err_val;
					max_err_sp = sp_a_id;
					max_err1 = tmp_avg[0]/vol_cell[ix];
				}

				if(rank == 0) col_dc2_3_1_2 += clock() - testc_3;

			}
		}

		value_type diag_ratio = double(nstep)/double(col_diag_1d_period*500);
		int diag_ratio_int = max(1, int(diag_ratio));
		if(rank == 0 && t_step%(col_diag_1d_period*diag_ratio_int) == 0 && ix == 0)
		{
			cout << "implicit iter col (iter_num, sp, del n, n0, del n/n0) : " << iter_num << " " << max_err_sp << " " << max_err1 << " " << col3_nUT_arr[max_err_sp][ix][0] << " " << max_err_val << endl;
		}

		iter_num += 1;
	}

	//generate evol_coeff_arr which is [M (df/dt)]
	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
	{
		if(sp_kinetic_op_arr[sp_a_id] == 1)
		{
			evol_coeff_arr[sp_a_id][sp_a_id] = col_sp_M_arr[sp_a_id][ix]*((local_coeff_arr[sp_a_id] - coeff_arr[sp_a_id])/sub_dt);
		}
	}

	//conservation part
	//col to fMb(Ua, Ta) + conserved quantity exchange rate calculation

	//initialize del_quantity_arr
	vector<vector<Vector>> del_quantity_arr(tot_species_num);
	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
	{
		del_quantity_arr[sp_a_id].resize(tot_species_num, Vector::Zero(2));
	}

	if(rank == 0) testc_3 = clock();
	//calculate del [Mom, En] due to C_RFP(f_a, f_b)
	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
	{
		if(sp_kinetic_op_arr[sp_a_id] == 1)
		{
			Vector sum_consv_w_evol = col_diag_cell_consv_mat[sp_a_id][ix]*evol_coeff_arr[sp_a_id][sp_a_id];
			del_quantity_arr[sp_a_id][sp_a_id] += sum_consv_w_evol;
		}
	}

	//calculate del_quantity_arr which is the (mom, en) change from collision
	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
	{
		if(sp_kinetic_op_arr[sp_a_id] == 1)
		{
			for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
			{
				if(sp_kinetic_op_arr[sp_b_id] == 1 && col3_sp_ab_ix_col_flag_arr[ix][sp_a_id][sp_b_id] == 1 && sp_b_id != sp_a_id)
				{
					//calculate col3_stored_del_quantity_arr[ix][sp_a_id][sp_b_id] here..
					int col3_hg_tot_n = col3_hg_tot_n_arr[sp_a_id];
					int dof = dof_arr[sp_a_id];
					int nv_loc = nv[sp_a_id];


					Vector hg_b_Vec = Vector::Zero(2*col3_hg_tot_n);

					for(int j = 0; j < 2; j++)
					{
						if(j== 0) 
						{
							hg_b_Vec.segment(0, col3_hg_tot_n) = (1.0 + Ms_arr[sp_b_id]/Ms_arr[sp_a_id])*gamma_aa_h_init_arr[sp_a_id][sp_b_id];
							hg_b_Vec.segment(col3_hg_tot_n, col3_hg_tot_n) = gamma_aa_g_init_arr[sp_a_id][sp_b_id];
						}
						else
						{
							hg_b_Vec.segment(0, col3_hg_tot_n) = (1.0 + Ms_arr[sp_b_id]/Ms_arr[sp_a_id])*gamma_aa_h_arr[sp_a_id][sp_b_id];
							hg_b_Vec.segment(col3_hg_tot_n, col3_hg_tot_n) = gamma_aa_g_arr[sp_a_id][sp_b_id];
						}

						Matrix f_to_consv = (consv_fa_del_hg_b_mat_arr[ix][sp_a_id]*hg_b_Vec).reshaped(2, nv_loc*dof);
						f_to_consv += consv_fa_fb0_mat_arr[ix][sp_a_id][sp_b_id];

						if(j== 0) 
						{
							col3_stored_del_quantity_arr[ix][sp_a_id][sp_b_id] = f_to_consv*coeff_arr[sp_a_id];
						}
						else
						{
							col3_stored_del_quantity_arr[ix][sp_a_id][sp_b_id] += f_to_consv*local_coeff_arr[sp_a_id];

							col3_stored_del_quantity_arr[ix][sp_a_id][sp_b_id] *= 0.5;
						}
					}

					if(col3_hMbgMb_with_UaTa_op == 1) col3_stored_del_quantity_arr[ix][sp_a_id][sp_b_id] += consv_fa_fb0_UaTa_mat_arr[ix][sp_a_id][sp_b_id];

					del_quantity_arr[sp_a_id][sp_b_id] = col3_stored_del_quantity_arr[ix][sp_a_id][sp_b_id];
					del_quantity_arr[sp_a_id][sp_a_id] -= del_quantity_arr[sp_a_id][sp_b_id];
				}
			}
		}
	}

	if(rank == 0) 
	{
		testc_4 = clock();
		col_dc2_6 += testc_4 - testc_3;
	}

	//Perform the conservation correction
	for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
	{
		if(sp_kinetic_op_arr[sp_a_id] == 1)
		{
			if(rank == 0) testc_3 = clock();
			int dof = dof_arr[sp_a_id];
			int nv_loc = nv[sp_a_id];
			Vector tot_rhs_arr = Vector::Zero(nv_loc*dof);
			Vector F2 = Vector::Zero(2);

			//set up F2 which is target (momentum, energy) via conservation from del_quantity_arr in the previous section
			for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
			{
				if(col3_sp_ab_ix_col_flag_arr[ix][sp_a_id][sp_b_id] == 1)
				{

					Vector del_tot_quantity = del_quantity_arr[sp_a_id][sp_b_id] + del_quantity_arr[sp_b_id][sp_a_id];

					if(sp_a_id == sp_b_id) F2 -= del_tot_quantity*0.5;
					else
					{
						if(sp_b_id == col3_smaller_vspace_ab_arr[sp_a_id][sp_b_id][ix])
						{
							//In this case, sp_a has a bigger vspace and Ma is simply assumed to be lighter than Mb.
							F2 -= del_tot_quantity;
						}
					}

					if(sp_b_id == sp_a_id) tot_rhs_arr += evol_coeff_arr[sp_a_id][sp_b_id];
				}
			}

			//apply conservation routine
			if (col_consv_onoff_op == 1) 
			{
					vector<Vector> evol_correction(2);
					int mom_consv_direction;

					if(F2[0] > 0.0) mom_consv_direction = 0;
					else mom_consv_direction = 1;

	
					//col3_f_to_mom_consv_mat_arr is a parallel shift in V_parallel direction : change momentum mainly
					//col3_f_to_en_consv_mat_arr is a isotropic diffusion in v space : change energy mainly
					evol_correction[0] = col3_f_to_mom_consv_mat_arr[sp_a_id][ix][mom_consv_direction]*local_coeff_arr[sp_a_id];
					evol_correction[1] = col3_f_to_en_consv_mat_arr[sp_a_id][ix]*local_coeff_arr[sp_a_id];

					Matrix col_consv_mat2 = Matrix::Zero(2,2);
					for(int k = 0; k < 2; k++)
					{
						Vector tmp_avg_w_evol = col_diag_cell_consv_mat[sp_a_id][ix]*evol_correction[k];
						for(int k2 = 0; k2 < 2; k2++) col_consv_mat2(k2, k) = tmp_avg_w_evol[k2];
					}

					//solve 2x2 matrix to determine the correction coefficients
					Vector X2 = col_consv_mat2.fullPivLu().solve(F2);

					value_type al0 = 0.0, al1 = X2[0], al2 = X2[1];

					tot_rhs_arr += al1*evol_correction[0] + al2*evol_correction[1];
			}

			if(rank == 0) 
			{
				testc_4 = clock();
				col_dc2_4 += testc_4 - testc_3;
				testc_3 = testc_4;
			}

			//Final step to calculate the updated f
			Lf_out_arr[sp_a_id] = col3_M_solver_arr[sp_a_id][ix]->solve(tot_rhs_arr);
			if(rank == 0) col_dc2_5 += clock() - testc_3;
		}
	}

	//Diagnostic routine 
	value_type diag_ratio = double(nstep)/double(col_diag_1d_period*500);
	int diag_ratio_int = max(1, int(diag_ratio));
	if (rank ==0 && t_step%col_diag_1d_period == 0 && iter_step == 0 && sub_step == 0)
	{
		if (ix == 0 && sub_step == 0)
		{
			ofstream fout;
			for(int sp_a_id = 0; sp_a_id < tot_species_num; sp_a_id++)
			{
				if(sp_kinetic_op_arr[sp_a_id] == 1)
				{
					string fname = "./col_diag_arr_" + to_string(sp_a_id) + ".txt";

					if (t_step == 0)
					{
						fout.open(fname);

							fout << "#1 : t_step, 2 : mass, 3 : V parallel, 4 : T_total, 5 : dbs_dvp_ov_bs, 6 : 3.0 + v_para*dbs_dvp_ov_bs, 7 : T_parallel, 8 : T_perp, 9 : total mom, 10 : total en, 11 : time (R/vT), 12 : real time (s)" << endl;
						fout.close();
					}

					fout.open(fname, ios_base::out | ios_base::app);
					if (fout.fail())
					{
						fout.open(fname);
						fout.close();
						fout.open(fname, ios_base::out | ios_base::app);
					}

					Vector tmp_avg = col_diag_cell_mat[sp_a_id][ix]*coeff_arr[sp_a_id];
					tmp_avg[1] /= tmp_avg[0];
					tmp_avg[2] *= Ms_arr[sp_a_id]/tmp_avg[0];
					tmp_avg[5] *= Ms_arr[sp_a_id]/tmp_avg[0];
					tmp_avg[6] *= Ms_arr[sp_a_id]/tmp_avg[0];


					fout.precision(16);
					fout << t_step << " ";
					for(int i = 0; i < col_diag_quantity_num; i++)
					{
						fout << tmp_avg[i] << " ";
					}
					fout << Ms_arr[sp_a_id]*tmp_avg[1]*tmp_avg[0] << " ";
					fout << tmp_avg[2]*tmp_avg[0] << " ";

					value_type vti0 = eq_reader->ph_const.get_property(Ph_const::vti0);
					value_type cf_x = eq_reader->ph_const.get_property(Ph_const::cf_x);
					
						fout << t_step*system_dt << " ";
						fout << t_step*system_dt*cf_x/vti0 << " ";
					fout << endl;

					fout.close();

					for(int sp_b_id = 0; sp_b_id < tot_species_num; sp_b_id++)
					{
						if(col3_sp_ab_ix_col_flag_arr[ix][sp_a_id][sp_b_id] == 1)
						{
							if (t_step%(col_diag_1d_period*diag_ratio_int) == 0)
								cout << "ion col : a, b, mom del, en del : " << sp_a_id << " " << sp_b_id << " " << del_quantity_arr[sp_a_id][sp_b_id][0] << " " << del_quantity_arr[sp_a_id][sp_b_id][1] << endl;

						}
					}
				}
			}
		}
	}

	if(rank == 0) col_dc5 += clock() - testc_5;
}

value_type Collision_dg::get_w_loc_mod(const value_type &ratio_in)
{

	value_type ratio_out, ratio_in_abs;

	value_type w_power_in_col = 1.0;

	ratio_out = abs(ratio_in);
	ratio_out = pow(abs(ratio_in), w_power_in_col);


	return ratio_out;
}

value_type Collision_dg::ie_mom_en_transfer_fn(int i1, int i2, int i3, int i4)
{
	return col3_ie_mom_en_transfer_arr[i1][i2][i3][i4];
}

void Collision_dg::ie_mom_en_transfer_write_fn(int i1, int i2, int i3, int i4, value_type loc_val)
{
	col3_ie_mom_en_transfer_arr[i1][i2][i3][i4] = loc_val;
}

/* Coulomb log functioins ------------------------------------------ */
/* Zi, Mi : ion charge and mass in unit of proton 
   ni, Ti : normalized density and temperature */
value_type Collision_dg::coulomb_log_ii(const value_type &Zi, const value_type &Mi, const value_type &ni_N, const value_type &Ti_N)
{
	value_type ni0 =eq_reader->ph_const.get_property(Ph_const::ni00);
	value_type Ti0 = eq_reader->ph_const.get_property(Ph_const::Ti0);
	value_type ni = ni_N*ni0;
	value_type Ti = Ti_N*Ti0;

	value_type cl_val = 23.0 - log((Zi*Zi/Ti)*sqrt(2.0*ni*Zi*Zi/Ti));
	cl_val *= col_mult_fac;
	return cl_val;
}

value_type Collision_dg::coulomb_log_ee(const value_type &ne_N, const value_type &Te_N)
{
	value_type ni0 =eq_reader->ph_const.get_property(Ph_const::ni00);
	value_type Ti0 = eq_reader->ph_const.get_property(Ph_const::Ti0);
	value_type ne = ne_N*ni0;
	value_type Te = Te_N*Ti0;

	value_type cl_val = 24.0 - log(sqrt(ne)/Te);
	cl_val *= col_mult_fac;
	return cl_val;
}

value_type Collision_dg::coulomb_log_Ii(const value_type &ZI, const value_type &MI, const value_type &nI_N, const value_type &TI_N, const value_type &Zi, const value_type &Mi, const value_type &ni_N, const value_type &Ti_N)
{
	value_type ni0 =eq_reader->ph_const.get_property(Ph_const::ni00);
	value_type Ti0 = eq_reader->ph_const.get_property(Ph_const::Ti0);
	value_type nI = nI_N*ni0;
	value_type TI = TI_N*Ti0;
	value_type ni = ni_N*ni0;
	value_type Ti = Ti_N*Ti0;

	value_type cl_val = 23.0 - log((Zi*ZI*(Mi+MI)/(Mi*TI+MI*Ti))*
			sqrt(ni*Zi*Zi/Ti + nI*ZI*ZI/TI));
	cl_val *= col_mult_fac;
	return cl_val;
}

/* basic ion collision time (sec) : NRL 37p */
value_type Collision_dg::taui_coll_time_cgs(const value_type &ni_N, const value_type &Ti_N, const value_type &Mi_ov_Mp, const value_type &Zi)
{
	value_type cl;

	cl = coulomb_log_ii(Zi, Mi_ov_Mp, ni_N, Ti_N);

	value_type ni0 =eq_reader->ph_const.get_property(Ph_const::ni00);
	value_type Ti0 = eq_reader->ph_const.get_property(Ph_const::Ti0);
	value_type ni = ni_N*ni0;
	value_type Ti = Ti_N*Ti0;

	return 2.0851e7*Ti*sqrt(Ti*Mi_ov_Mp)/(ni*cl);
}

value_type Collision_dg::taui_ab_coll_time_cgs(const int sp_a, const value_type &na_N, const value_type &Ta_N, const int sp_b, const value_type &nb_N, const value_type &Tb_N)
{
	value_type Ma_ov_Mp = Ms_ov_Mp_arr[sp_a];
	value_type Za = Zs_ov_e_arr[sp_a];
	value_type Mb_ov_Mp = Ms_ov_Mp_arr[sp_b];
	value_type Zb = Zs_ov_e_arr[sp_b];

	value_type cl;
	if(sp_a == 0)
	{
		cl = coulomb_log_ee(na_N, Ta_N);
	}
	else if (sp_b == 0)
	{
		cl = coulomb_log_ee(nb_N, Tb_N);
	}
	else
	{
		cl = coulomb_log_Ii(Za, Ma_ov_Mp, na_N, Ta_N, Zb, Mb_ov_Mp, nb_N, Tb_N);
	}

	value_type ni0 =eq_reader->ph_const.get_property(Ph_const::ni00);
	value_type Ti0 = eq_reader->ph_const.get_property(Ph_const::Ti0);
	value_type ni_b = nb_N*ni0;
	value_type Ti_a = Ta_N*Ti0;

	return 2.0851e7*Ti_a*sqrt(Ti_a*Ma_ov_Mp)/(Za*Za*Zb*Zb*ni_b*cl);
}

value_type Collision_dg::tau_ii_coll_time_norm(const value_type &ni_N, const value_type &Ti_N, const value_type &Mi_ov_Mp, const value_type &Zi_ov_e)
{
	value_type tau_cgs = taui_coll_time_cgs(ni_N, Ti_N, Mi_ov_Mp, Zi_ov_e)/sqrt(2);

	value_type vti0 = eq_reader->ph_const.get_property(Ph_const::vti0);
	value_type cf_x = eq_reader->ph_const.get_property(Ph_const::cf_x);
	return tau_cgs/(cf_x/vti0);
}

value_type Collision_dg::taui_ab_coll_time_norm(const int sp_a, const value_type &na_N, const value_type &Ta_N, const int sp_b, const value_type &nb_N, const value_type &Tb_N)
{
	value_type tau_cgs = taui_ab_coll_time_cgs(sp_a, na_N, Ta_N, sp_b, nb_N, Tb_N)/sqrt(2);

	value_type vti0 = eq_reader->ph_const.get_property(Ph_const::vti0);
	value_type cf_x = eq_reader->ph_const.get_property(Ph_const::cf_x);
	return tau_cgs/(cf_x/vti0);
}

/* Gamma_ab (cm^6/s^4) : (4*pi*q_a*q_b/m_a)^2*cln */
value_type Collision_dg::gamma_ab_cgs(const value_type &ni_N, const value_type &Ti_N, const value_type &Ms_ov_Mp, const value_type &Zi)
{
	value_type cl;

	cl = coulomb_log_ii(Zi, Ms_ov_Mp, ni_N, Ti_N);

	return 3.0044e12*(Zi*Zi*Zi*Zi/Ms_ov_Mp/Ms_ov_Mp)*cl;
}

value_type Collision_dg::gamma_ab_norm(const value_type &ni_N, const value_type &Ti_N, const value_type &Ms_ov_Mp, const value_type &Zi_ov_e)
{
	value_type gamma_cgs = gamma_ab_cgs(ni_N, Ti_N, Ms_ov_Mp, Zi_ov_e);

	value_type vti0 = eq_reader->ph_const.get_property(Ph_const::vti0);
	value_type cf_x = eq_reader->ph_const.get_property(Ph_const::cf_x);
	value_type ni0 =eq_reader->ph_const.get_property(Ph_const::ni00);

	return gamma_cgs/(vti0*vti0*vti0*vti0)*cf_x*ni0;
}

value_type Collision_dg::gamma_ab_gen_cgs(int sp_a, const value_type &na_N, const value_type &Ta_N, int sp_b, const value_type &nb_N, const value_type &Tb_N)
{
	value_type cl;
	value_type Ma_ov_Mp = species_data->mass_ov_mp(sp_a);
	value_type Za = species_data->charge_ov_e(sp_a);
	value_type Mb_ov_Mp = species_data->mass_ov_mp(sp_b);
	value_type Zb = species_data->charge_ov_e(sp_b);


	if(sp_a == 0)
	{
		cl = coulomb_log_ee(na_N, Ta_N);
	}
	else if (sp_b == 0)
	{
		cl = coulomb_log_ee(nb_N, Tb_N);
	}
	else
	{
		cl = coulomb_log_Ii(Za, Ma_ov_Mp, na_N, Ta_N, Zb, Mb_ov_Mp, nb_N, Tb_N);
	}

	return 3.0044e12*(Za*Za*Zb*Zb/Ma_ov_Mp/Ma_ov_Mp)*cl;
}

value_type Collision_dg::gamma_ab_gen_norm(int sp_a, const value_type &na_N, const value_type &Ta_N, int sp_b, const value_type &nb_N, const value_type &Tb_N)
{
	value_type gamma_cgs = gamma_ab_gen_cgs(sp_a, na_N, Ta_N, sp_b, nb_N, Tb_N);

	value_type vti0 = eq_reader->ph_const.get_property(Ph_const::vti0);
	value_type cf_x = eq_reader->ph_const.get_property(Ph_const::cf_x);
	value_type ni0 =eq_reader->ph_const.get_property(Ph_const::ni00);

	return gamma_cgs/(vti0*vti0*vti0*vti0)*cf_x*ni0;
}

//collision module speed diagnostic
void Collision_dg::col_time_out()
{
	if(rank == 0)
	{
		cout << "CPU time used in tot col steps : " << double(tot_dc1)/CLOCKS_PER_SEC << " s" << endl;
		cout << "CPU time used in col fhat, U, T calculation : " << double(tot_dc2)/CLOCKS_PER_SEC << " s" << endl;
		cout << "CPU time used in col volume terms calculation : " << double(tot_dc3)/CLOCKS_PER_SEC << " s" << endl;
		cout << "CPU time used in col edge terms calculation : " << double(tot_dc4)/CLOCKS_PER_SEC << " s" << endl;
		cout << "CPU time used in col consv mat  calculation : " << double(tot_dc5)/CLOCKS_PER_SEC << " s" << endl;
		cout << "CPU time used in col M solve : " << double(tot_dc6)/CLOCKS_PER_SEC << " s" << endl;
	}
}


