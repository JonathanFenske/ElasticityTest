#ifndef _BASIS_FUNS_TPP_
#define _BASIS_FUNS_TPP_

#include "basis_funs.h"

namespace BasisFun
{
  using namespace dealii;

  template <int dim>
  BasisQ1<dim>::BasisQ1(const BasisQ1<dim> &basis)
    : Function<dim>(dim)
    , index_basis(0)
    , coeff_matrix(basis.coeff_matrix)
  {}


  template <>
  BasisQ1<2>::BasisQ1(
    const typename Triangulation<2>::active_cell_iterator &cell)
    : Function<2>(2)
    , index_basis(0)
    , coeff_matrix(4, 4)
  {
    FullMatrix<double> point_matrix(4, 4);

    for (unsigned int i = 0; i < 4; ++i)
      {
        const Point<2> &p = cell->vertex(i);

        point_matrix(i, 0) = 1;
        point_matrix(i, 1) = p(0);
        point_matrix(i, 2) = p(1);
        point_matrix(i, 3) = p(0) * p(1);
      }

    // Columns of coeff_matrix are the coefficients of the polynomial
    coeff_matrix.invert(point_matrix);
  }


  template <>
  BasisQ1<3>::BasisQ1(
    const typename Triangulation<3>::active_cell_iterator &cell)
    : Function<3>(3)
    , index_basis(0)
    , coeff_matrix(8, 8)
  {
    FullMatrix<double> point_matrix(8, 8);

    for (unsigned int i = 0; i < 8; ++i)
      {
        const Point<3> &p = cell->vertex(i);

        point_matrix(i, 0) = 1;
        point_matrix(i, 1) = p(0);
        point_matrix(i, 2) = p(1);
        point_matrix(i, 3) = p(2);
        point_matrix(i, 4) = p(0) * p(1);
        point_matrix(i, 5) = p(1) * p(2);
        point_matrix(i, 6) = p(0) * p(2);
        point_matrix(i, 7) = p(0) * p(1) * p(2);
      }

    // Columns of coeff_matrix are the coefficients of the polynomial
    coeff_matrix.invert(point_matrix);
  }


  template <int dim>
  void
  BasisQ1<dim>::set_index(unsigned int index)
  {
    Assert(
      index <=
        static_cast<unsigned int>(GeometryInfo<dim>::vertices_per_cell * dim),
      ExcMessage(
        "Index of basis function must be smaller than dim * numer of corner points."));

    // This is integer division (returns only the integer part)
    index_basis = index / dim;

    // This is modulo division (returns only the rest)
    component_basis = index % GeometryInfo<dim>::vertices_per_cell;
  }


  template <>
  void
  BasisQ1<2>::vector_value(const Point<2> &p,
                           Vector<double> &vector_value) const
  {
    // set every component to zero
    vector_value = 0;

    double value = coeff_matrix(0, index_basis) +
                   coeff_matrix(1, index_basis) * p(0) +
                   coeff_matrix(2, index_basis) * p(1) +
                   coeff_matrix(3, index_basis) * p(0) * p(1);

    vector_value(component_basis) = value;
  }


  template <>
  void
  BasisQ1<3>::vector_value(const Point<3> &p,
                           Vector<double> &vector_value) const
  {
    // set every component to zero
    vector_value = 0;

    double value = coeff_matrix(0, index_basis) +
                   coeff_matrix(1, index_basis) * p(0) +
                   coeff_matrix(2, index_basis) * p(1) +
                   coeff_matrix(3, index_basis) * p(2) +
                   coeff_matrix(4, index_basis) * p(0) * p(1) +
                   coeff_matrix(5, index_basis) * p(1) * p(2) +
                   coeff_matrix(6, index_basis) * p(0) * p(2) +
                   coeff_matrix(7, index_basis) * p(0) * p(1) * p(2);

    vector_value(component_basis) = value;
  }


  template <>
  void
  BasisQ1<2>::vector_value_list(const std::vector<Point<2>> &points,
                                std::vector<Vector<double>> &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    for (unsigned int p = 0; p < points.size(); ++p)
      {
        vector_value(points[p], values[p]);
      } // end ++p
  }


  template <>
  void
  BasisQ1<3>::vector_value_list(const std::vector<Point<3>> &points,
                                std::vector<Vector<double>> &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    for (unsigned int p = 0; p < points.size(); ++p)
      {
        vector_value(points[p], values[p]);
      } // end ++p
  }
} // namespace BasisFun

#endif // _BASIS_FUNS_TPP_