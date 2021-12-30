#ifndef _INCLUDE_FINE_MS_SOLUTION_FUNCTION_TPP_
#define _INCLUDE_FINE_MS_SOLUTION_FUNCTION_TPP_

#include "fine_ms_solution_function.h"

namespace Elasticity
{
  using namespace dealii;

  // template <int dim>
  // FineMsSolution<dim>
  // FineMsSolution<dim>::operator=(const FineMsSolution<dim> &other)
  // {
  //   triangulation.copy_triangulation(other.triangulation);
  //   cell_basis_map = other.cell_basis_map;
  //   return *this;
  // }

  template <int dim>
  void
  FineMsSolution<dim>::vector_value(const Point<dim> &p,
                                    Vector<double>   &vector_value) const
  {
    typename Triangulation<dim>::active_cell_iterator cell =
      GridTools::find_active_cell_around_point(triangulation, p);

    (function_map.at(cell->id()))->vector_value(p, vector_value);
  }

  template <int dim>
  void
  FineMsSolution<dim>::vector_gradient(
    const Point<dim>            &p,
    std::vector<Tensor<1, dim>> &gradients) const
  {
    typename Triangulation<dim>::active_cell_iterator cell =
      GridTools::find_active_cell_around_point(triangulation, p);

    (function_map.at(cell->id()))->vector_gradient(p, gradients);
  }
} // namespace Elasticity

#endif // _INCLUDE_FINE_MS_SOLUTION_FUNCTION_TPP_