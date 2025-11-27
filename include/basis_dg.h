#pragma once

#include "basis.h"

using namespace std;

class Basis_DG_quadratic_Seo : public Basis
{
	public:
		Basis_DG_quadratic_Seo();
		~Basis_DG_quadratic_Seo() = default;
		
		value_type operator()(const Point4 &, int, int) const;
		value_type operator()(const Point4 &, const Point4 &, const ElementX &, const ElementV &, int, int = 0) const;
		Vector ind_val(const Point4 &, int, int) const;
		value_type ind_val_S1(const Point4 &, int) const;

	protected:
};


