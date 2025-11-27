#include <gyrokinetic.h>

Gyrokinetic::Gyrokinetic(const Sokuri &sokuri, const MPI_Comm &comm) : comm(comm), Indexer(*sokuri.mesh, *sokuri.basis)
{
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &nproc);

	sp_tot_num = species_data->get_tot_sp_num();

	collision_dg = new Collision_dg(*mesh, *eq_reader, *config, *species_data, comm);

	for(int sp_i = 0; sp_i < sp_tot_num; sp_i++) 
	{	
		//setup module
		procs_mult_setup(sp_i);
	}

	collision_dg->initial_h0g0_setup(sp_dg_coeff);
}

//setup module for M matrix
void Gyrokinetic::procs_mult_setup(int sp_id)
{
	if (sp_id == 0)
	{
		//Configuration parameters
		dt = config->get_option<double>("sml_dt");
		nstep = config->get_option<int>("sml_nstep");
	}

	if (sp_kinetic_op[sp_id] == 1)
	{
		//Configuration miscellaneous parameters
		nx = nx_owned_arr[sp_id] = mesh_sp->own_size<ElementX>();
		nv = nv_arr[sp_id] = mesh_sp->size<ElementV>();
		num_elm = nx*nv;

		dof = dof_arr[sp_id] = basis_sp->get_dof();

		vector<value_type> tmp_M_mat_arr(dof*dof), tmp_E_mat_arr(dof*dof);
		vector<value_type> tmp_finit_arr(dof);

		vector<SparseMatrix_Triplet> tripletList_M_A_col;
		tripletList_M_A_col.reserve(nv*dof*dof);

		collision_dg->collision_dg_coeff_setup(*mesh_sp, *quadrature_sp, *basis_sp, *flux_sp, sp_id);

		//M matrix calculation 
		for(int ix = 0; ix < nx; ix++)
		{
			for(int iv = 0; iv < nv; iv++)
			{
				int k = ix*nv + iv;
				integration_MSE_mat_cal->ME_mat_cal_f_init(ElementX(ix), ElementV(iv), sp_id, tmp_M_mat_arr, tmp_E_mat_arr, tmp_finit_arr);

				for(int j = 0; j < dof; j++)
				{
					for(int i = 0; i < dof; i++)
					{
						tripletList_M_A_col.push_back(SparseMatrix_Triplet(iv*dof + i, iv*dof + j, tmp_M_mat_arr[j*dof + i]));
					}
				}

			}

			SparseMatrix M_A_col(dof*nv,dof*nv);
			M_A_col.setFromTriplets(tripletList_M_A_col.begin(), tripletList_M_A_col.end());
			collision_dg->col_consv_mat_setup(M_A_col, sp_id, ix);

			tripletList_M_A_col.clear();
		}
	}
}

//main loop module
void Gyrokinetic::procs_mult(void)
{
	for(int step = init_t_step; step < nstep; step++)
	{
		//lowest order part update
		collision_dg->h0g0fM_update(sp_dg_coeff, step, 0);

		//collision module
		collision_dg->RK_implicit_col(sp_dg_coeff, step);
	}
}




