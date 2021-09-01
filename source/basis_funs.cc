#include "basis_funs.h"
#include "basis_funs.tpp"

namespace BasisFun
{
  using namespace dealii;

  // template <>
  // BasisQ1Grad<2>::BasisQ1Grad(
  //   const typename Triangulation<2>::active_cell_iterator &cell)
  //   : Function<2>(2)
  //   , q_point(0)
  //   , coeff_matrix(4, 4)
  // {
  //   FullMatrix<double> point_matrix(4, 4);

  //   for (unsigned int i = 0; i < 4; ++i)
  //     {
  //       const Point<2> &p = cell->vertex(i);

  //       point_matrix(i, 0) = 1;
  //       point_matrix(i, 1) = p(0);
  //       point_matrix(i, 2) = p(1);
  //       point_matrix(i, 3) = p(0) * p(1);
  //     }

  //   // Columns of coeff_matrix are the coefficients of the polynomial
  //   coeff_matrix.invert(point_matrix);
  // }


  // template <>
  // BasisQ1Grad<3>::BasisQ1Grad(
  //   const typename Triangulation<3>::active_cell_iterator &cell)
  //   : Function<3>(3)
  //   , q_point(0)
  //   , coeff_matrix(8, 8)
  // {
  //   FullMatrix<double> point_matrix(8, 8);

  //   for (unsigned int i = 0; i < 8; ++i)
  //     {
  //       const Point<3> &p = cell->vertex(i);

  //       point_matrix(i, 0) = 1;
  //       point_matrix(i, 1) = p(0);
  //       point_matrix(i, 2) = p(1);
  //       point_matrix(i, 3) = p(2);
  //       point_matrix(i, 4) = p(0) * p(1);
  //       point_matrix(i, 5) = p(1) * p(2);
  //       point_matrix(i, 6) = p(0) * p(2);
  //       point_matrix(i, 7) = p(0) * p(1) * p(2);
  //     }

  //   // Columns of coeff_matrix are the coefficients of the polynomial
  //   coeff_matrix.invert(point_matrix);
  // }


  // template <>
  // void
  //   BasisQ1Grad<2>::vector_value(const Point<2> &p, Vector<double> &value) const
  // {
  //   value(0) =
  //     coeff_matrix(1, q_point) + coeff_matrix(3, q_point) * p(1);
  //   value(1) =
  //     coeff_matrix(2, q_point) + coeff_matrix(3, q_point) * p(0);
  // }


  // template <>
  // void
  //   BasisQ1Grad<3>::vector_value(const Point<3> &p, Vector<double> &value) const
  // {
  //   value(0) = coeff_matrix(1, q_point) +
  //              coeff_matrix(4, q_point) * p(1) +
  //              coeff_matrix(6, q_point) * p(2) +
  //              coeff_matrix(7, q_point) * p(1) * p(2);
  //   value(1) = coeff_matrix(2, q_point) +
  //              coeff_matrix(4, q_point) * p(0) +
  //              coeff_matrix(5, q_point) * p(2) +
  //              coeff_matrix(7, q_point) * p(0) * p(2);
  //   value(2) = coeff_matrix(3, q_point) +
  //              coeff_matrix(5, q_point) * p(1) +
  //              coeff_matrix(6, q_point) * p(0) +
  //              coeff_matrix(7, q_point) * p(0) * p(1);
  // }


  // template <>
  // void
  //   BasisQ1Grad<2>::vector_value_list(const std::vector<Point<2>> &points,
  //                                     std::vector<Vector<double>> &values) const
  // {
  //   Assert(points.size() == values.size(),
  //          ExcDimensionMismatch(points.size(), values.size()));

  //   for (unsigned int p = 0; p < points.size(); ++p)
  //     {
  //       vector_value(points[p], values[p]);
  //     } // end ++p
  // }


  // template <>
  // void
  //   BasisQ1Grad<3>::vector_value_list(const std::vector<Point<3>> &points,
  //                                     std::vector<Vector<double>> &values) const
  // {
  //   Assert(points.size() == values.size(),
  //          ExcDimensionMismatch(points.size(), values.size()));

  //   for (unsigned int p = 0; p < points.size(); ++p)
  //     {
  //       vector_value(points[p], values[p]);
  //     } // end ++p
  // }


  // template <>
  // void
  //   BasisQ1Grad<2>::tensor_value_list(const std::vector<Point<2>> &points,
  //                                     std::vector<Tensor<1, 2>> &  values) const
  // {
  //   Assert(points.size() == values.size(),
  //          ExcDimensionMismatch(points.size(), values.size()));

  //   Vector<double> value_tmp(2);
  //   for (unsigned int p = 0; p < points.size(); ++p)
  //     {
  //       value_tmp = 0;
  //       vector_value(points[p], value_tmp);

  //       values[p][0] = value_tmp(0);
  //       values[p][1] = value_tmp(1);
  //     } // end ++p
  // }


  // template <>
  // void
  //   BasisQ1Grad<3>::tensor_value_list(const std::vector<Point<3>> &points,
  //                                     std::vector<Tensor<1, 3>> &  values) const
  // {
  //   Assert(points.size() == values.size(),
  //          ExcDimensionMismatch(points.size(), values.size()));

  //   Vector<double> value_tmp(3);
  //   for (unsigned int p = 0; p < points.size(); ++p)
  //     {
  //       value_tmp = 0;
  //       vector_value(points[p], value_tmp);

  //       values[p][0] = value_tmp(0);
  //       values[p][1] = value_tmp(1);
  //       values[p][2] = value_tmp(2);
  //     } // end ++p
  // }

  template class BasisQ1Grad<2>;
  template class BasisQ1Grad<3>;
}