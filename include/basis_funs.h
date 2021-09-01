#ifndef _BASIS_FUNS_h_
#define _BASIS_FUNS_h_

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

namespace BasisFun
{
  using namespace dealii;

  /*!
   * @class BasisQ1Grad
   *
   * @brief Gradient of \f$Q_1\f$ basis on given cell
   *
   * Class implements curls of vectorial Nedelec basis for a given
   * quadrilateral.
   *
   * @note The gradients of \f$H(\mathrm{grad})\f$ conforming functions are in \f$H(\mathrm{curl})\f$.
   * So we need covariant transforms to guarantee \f$H(\mathrm{curl})\f$
   * conformity.
   */
  template <int dim>
  class BasisQ1Grad : public Function<dim>
  {
  public:
    BasisQ1Grad() = delete;

    /*!
     * Constructor.
     *
     * @param cell
     */
    BasisQ1Grad(const typename Triangulation<dim>::active_cell_iterator &cell);

    /*!
     * Copy constructor.
     */
    BasisQ1Grad(const BasisQ1Grad<dim> &);

    /*!
     * Set the q_point of the basis function to be evaluated.
     *
     * @param q_point
     */
    void
      set_q_point(unsigned int q_point);

    /*!
     * Compute the gradient of a Q1 function.
     *
     * @param p
     * @param value
     */
    virtual void
      vector_value(const Point<dim> &p, Vector<double> &value) const override;

    /*!
     * Compute the gradient of a Q1 function for a list of points.
     *
     * @param[in] points
     * @param[out] values
     */
    virtual void
      vector_value_list(const std::vector<Point<dim>> &points,
                        std::vector<Vector<double>> &  values) const override;

    /*!
     * Compute the gradient of a Q1 function for a list of points. Return
     * tensors.
     *
     * @param[in] points
     * @param[out] values
     */
    void
      tensor_value_list(const std::vector<Point<dim>> &points,
                        std::vector<Tensor<1, dim>> &  values) const;

  private:
    /*!
     * q_point of current basis function to be evaluated.
     */
    unsigned int q_point;

    /*!
     * Matrix columns hold coefficients of basis functions.
     */
    FullMatrix<double> coeff_matrix;
  };

  // declare specializations
  template <>
  BasisQ1Grad<2>::BasisQ1Grad(
    const typename Triangulation<2>::active_cell_iterator &cell);

  template <>
  BasisQ1Grad<3>::BasisQ1Grad(
    const typename Triangulation<3>::active_cell_iterator &cell);

  template <>
  void
    BasisQ1Grad<2>::vector_value(const Point<2> &p,
                                 Vector<double> &value) const;

  template <>
  void
    BasisQ1Grad<3>::vector_value(const Point<3> &p,
                                 Vector<double> &value) const;

  template <>
  void
    BasisQ1Grad<2>::vector_value_list(
      const std::vector<Point<2>> &points,
      std::vector<Vector<double>> &values) const;

  template <>
  void
    BasisQ1Grad<3>::vector_value_list(
      const std::vector<Point<3>> &points,
      std::vector<Vector<double>> &values) const;


  template <>
  void
    BasisQ1Grad<2>::tensor_value_list(const std::vector<Point<2>> &points,
                                      std::vector<Tensor<1, 2>> &values) const;

  template <>
  void
    BasisQ1Grad<3>::tensor_value_list(const std::vector<Point<3>> &points,
                                      std::vector<Tensor<1, 3>> &values) const;

  // exernal template instantiations
  extern template class BasisQ1Grad<2>;
  extern template class BasisQ1Grad<3>;
}

#endif // _BASIS_FUNS_h_