#pragma once

template<class T>
Vector Basis::grad(const T &p, int i) const
{
	Vector d(dim);

	for(int j = 0; j < dim; j++)
		d[j] = operator()(p, i, j+1);

	return d;
}


template<class T>
Vector Basis::average_grad(const T &p, int i) const
{
	return 0.5*grad(p, i);
}


template<class T>
value_type Basis::jump(const T &p, int i) const
{
	return operator()(p, i, 0);
}


template<class T>
Vector Basis::grad(const T &p, int i, vector<value_type> &raw_basis) const
{
	Vector d(dim);

	for(int j = 0; j < dim; j++)
		d[j] = operator()(p, i, raw_basis, j+1);

	return d;
}

template<class T>
Vector Basis::grad(const T &p, const Point4 &p_org, const ElementX &ti, const ElementV &tj, int i) const
{
	Vector d(dim);

	for(int j = 0; j < dim; j++)
		d[j] = operator()(p, p_org, ti, tj, i, j+1);

	return d;
}
