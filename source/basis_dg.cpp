#include "basis_dg.h"

Basis_DG_quadratic_Seo::Basis_DG_quadratic_Seo()
{
	dofx = dofv = 6;
	dof = dofx*dofv;
	dof_S1_x = 4;
	dof_S1 = dof_S1_x*dofv; 
	dim = 4;

	fhat_v_basis_group_id_default = 1;
	x_basis_group_id = 1;
	v_basis_group_id = 1;
}

value_type Basis_DG_quadratic_Seo::operator()(const Point4 &p, int i, int diff) const
{
	assert(i >= 0 && i < dof);

	int r = i % dofx;
	int q = (i - r) / dofx;

	vector<value_type> xf(dofx), vf(dofv);

	switch(diff)
	{
		case 0:
			xf = {1.0 - p[0] - p[1], p[0], p[1], (1 - p[0] - p[1])*(1 - p[0] - p[1]), p[0]*p[0], p[1]*p[1]};
			vf = {1.0, p[2], p[3], p[2]*p[2], p[2]*p[3], p[3]*p[3]};
			break;

		case 1:
			xf = {-1.0, 1.0, 0.0, -2.0*(1 - p[0] - p[1]), 2.0*p[0], 0.0};
			vf = {1.0, p[2], p[3], p[2]*p[2], p[2]*p[3], p[3]*p[3]};
			break;

		case 2:
			xf = {-1.0, 0.0, 1.0,  -2.0*(1 - p[0] - p[1]), 0.0, 2.0*p[1]};
			vf = {1.0, p[2], p[3], p[2]*p[2], p[2]*p[3], p[3]*p[3]};
			break;

		case 3:
			xf = {1.0 - p[0] - p[1], p[0], p[1], (1 - p[0] - p[1])*(1- p[0] - p[1]), p[0]*p[0], p[1]*p[1]};
			vf = {0.0, 1.0, 0.0, 2.0*p[2], p[3], 0.0};
			break;

		case 4:
			xf = {1.0 - p[0] - p[1], p[0], p[1], (1 - p[0] - p[1])*(1- p[0] - p[1]), p[0]*p[0], p[1]*p[1]};
			vf = {0.0, 0.0, 1.0, 0.0, p[2], 2.0*p[3]};
			break;

		default:
			cerr << "Do you really need higher order derivatives?" << endl;
			exit(EXIT_FAILURE);
	}

	return xf[r]*vf[q];
}


value_type Basis_DG_quadratic_Seo::operator()(const Point4 &p, const Point4 &p_org, const ElementX &ti, const ElementV &tj, int i, int diff) const
{
	assert(i >= 0 && i < dof);

	int r = i % dofx;
	int q = (i - r) / dofx;

	vector<value_type> xf(dofx), vf(dofv);

	switch(diff)
	{
		case 0:
			xf = {1.0 - p[0] - p[1], p[0], p[1], (1 - p[0] - p[1])*(1 - p[0] - p[1]), p[0]*p[0], p[1]*p[1]};
			vf = {1.0, p[2], p[3], p[2]*p[2], p[2]*p[3], p[3]*p[3]};
			break;

		case 1:
			xf = {-1.0, 1.0, 0.0, -2.0*(1 - p[0] - p[1]), 2.0*p[0], 0.0};
			vf = {1.0, p[2], p[3], p[2]*p[2], p[2]*p[3], p[3]*p[3]};
			break;

		case 2:
			xf = {-1.0, 0.0, 1.0,  -2.0*(1 - p[0] - p[1]), 0.0, 2.0*p[1]};
			vf = {1.0, p[2], p[3], p[2]*p[2], p[2]*p[3], p[3]*p[3]};
			break;

		case 3:
			xf = {1.0 - p[0] - p[1], p[0], p[1], (1 - p[0] - p[1])*(1- p[0] - p[1]), p[0]*p[0], p[1]*p[1]};
			vf = {0.0, 1.0, 0.0, 2.0*p[2], p[3], 0.0};
			break;

		case 4:
			xf = {1.0 - p[0] - p[1], p[0], p[1], (1 - p[0] - p[1])*(1- p[0] - p[1]), p[0]*p[0], p[1]*p[1]};
			vf = {0.0, 0.0, 1.0, 0.0, p[2], 2.0*p[3]};
			break;

		default:
			cerr << "Do you really need higher order derivatives?" << endl;
			exit(EXIT_FAILURE);
	}

	return xf[r]*vf[q];
}

Vector Basis_DG_quadratic_Seo::ind_val(const Point4 &p, int i, int grad_op) const
{
	assert(i >= 0 && i < dof);
	int dim = 4;
	Vector d(dim);

	int r = i % dofx;
	int q = (i - r) / dofx;

	int X_dim = 2;
	vector<vector<value_type>> xf;
	xf.resize(X_dim);
	for(int j = 0; j < X_dim; j++)
		xf[j].resize(dofx);

	vector<value_type> vf(dofv), vf_p(dofv), vf_u(dofv);

	switch(grad_op)
	{
		case 0:
			xf[0] = {1.0 - p[0] - p[1], p[0], p[1], (1 - p[0] - p[1])*(1 - p[0] - p[1]), p[0]*p[0], p[1]*p[1]};
			xf[1] = {1.0 - p[0] - p[1], p[0], p[1], (1 - p[0] - p[1])*(1 - p[0] - p[1]), p[0]*p[0], p[1]*p[1]};
			vf = {1.0, p[2], p[3], p[2]*p[2], p[2]*p[3], p[3]*p[3]};
			vf_p = {1.0, p[2], 1.0, p[2]*p[2], p[2], 1.0};
			vf_u = {1.0, 1.0, p[3], 1.0, p[3], p[3]*p[3]};
			break;

		case 1:
			xf[0] = {-1.0, 1.0, 0.0, -2.0*(1 - p[0] - p[1]), 2.0*p[0], 0.0};
			xf[1] = {-1.0, 0.0, 1.0,  -2.0*(1 - p[0] - p[1]), 0.0, 2.0*p[1]};
			vf_p = {0.0, 1.0, 0.0, 2.0*p[2], 1.0, 0.0};
			vf_u = {0.0, 0.0, 1.0, 0.0, 1.0, 2.0*p[3]};
			break;

		default:
			cerr << "Unknown grdp_op in ind_val : " << grad_op << endl;
			exit(EXIT_FAILURE);
	}

	d[0] = xf[0][r];
	d[1] = xf[1][r];
	d[2] = vf_p[q];
	d[3] = vf_u[q];

	return d;
}

value_type Basis_DG_quadratic_Seo::ind_val_S1(const Point4 &p, int i) const
{
	assert(i >= 0 && i < dof_S1);

	int r = i % dof_S1_x;
	int q = (i - r) / dof_S1_x;

	vector<value_type> xf(dof_S1_x), vf(dofv);
	xf = {0.5*(1.0 - p[0]), 0.5*(1.0 + p[0]), 0.25*(1.0 - p[0])*(1.0 - p[0]), 0.25*(1.0 + p[0])*(1.0 + p[0])};
	vf = {1.0, p[2], p[3], p[2]*p[2], p[2]*p[3], p[3]*p[3]};

	return xf[r]*vf[q]; 
}

