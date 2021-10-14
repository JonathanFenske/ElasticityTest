#ifndef _INCLUDE_POSTPROCESSING_H_
#define _INCLUDE_POSTPROCESSING_H_

#include <deal.II/base/utilities.h>

#include <deal.II/numerics/data_postprocessor.h>

#include "process_parameter_file.h"

/**
 * @file postprocessing.h
 *
 * @brief Postprocessing
 *
 * @see https://www.dealii.org/9.2.0/doxygen/deal.II/classDataPostprocessor.html
 */

namespace Elasticity
{
  using namespace dealii;

  /****************************************************************************/
  /* Postprocessing */

  // (The uncommented parts are my try to combine stress and strain
  // into a single class.)
  /**
   * @brief Class that enables the output of the linearized strain tensor.
   *
   * @tparam dim Space dimension
   *
   * This class enables the output of the linearized strain tensor.
   *
   * Based on the documentation for the deal.ii class
   * <a href="https://www.dealii.org/9.2.0/doxygen/deal.II/classData
   * PostprocessorTensor.html">DataPostprocessorTensor</a>.
   */
  template <int dim>
  class StrainPostprocessor : public DataPostprocessor<dim>
  {
  public:
    /**
     * @brief Construct a new StrainPostprocessor object for ElaBasis.
     *
     * @param basis_index Index of the basis function
     */
    StrainPostprocessor(unsigned int basis_index);

    /**
     * @brief Construct a new StrainPostprocessor object.
     */
    StrainPostprocessor();

    /**
     * @brief Copy constructor for StrainPostprocessor objects
     *
     * @param other
     */
    StrainPostprocessor(const StrainPostprocessor<dim> &other) = default;

    /**
     * @brief Evaluate the vector field and
     *        construct the linearized strain tensor.
     *
     * @param input_data Input data
     * @param computed_quantities Vector which will be overridden with the
     *                            linearized strain tensor.
     */
    virtual void
    evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &input_data,
      std::vector<Vector<double>> &computed_quantities) const override;

    /**
     * @brief Specify how the computed contities can be interpreted.
     *
     * @return std::vector<
     * DataComponentInterpretation::DataComponentInterpretation>
     */
    virtual std::vector<
      DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation() const override;

    /**
     * @brief Specify the names of the components of the computed data
     *
     * @return std::vector<std::string>
     */
    virtual std::vector<std::string>
    get_names() const override;

    /**
     * @brief Specify which update flags are needed
     *        for the shape functions in
     *        get_data_component_interpretation()
     *
     * @return UpdateFlags
     */
    virtual UpdateFlags
    get_needed_update_flags() const override;

  private:
    /**
     * String that contains the index of the local shape
     * shape function in ElaBasis.
     */
    std::string basis_str;
  };


  /**
   * @brief Class that enables the output of the linearized stress tensor.
   *
   * @tparam dim Space dimension
   */
  template <int dim>
  class StressPostprocessor : public DataPostprocessor<dim>
  {
  public:
    /**
     * @brief Construct a new (empty) StressPostprocessor object.
     */
    StressPostprocessor();

    /**
     * @brief Construct a new StressPostprocessor object for ElaBasis.
     *
     * @param basis_index Index of the basis function
     * @param global_parameters Parameters that many classes need
     */
    StressPostprocessor(unsigned int                 basis_index,
                        const GlobalParameters<dim> &global_parameters);

    /**
     * @brief Construct a new Stress Postprocessor object
     *
     * @param global_parameters Parameters that many classes need
     */
    StressPostprocessor(const GlobalParameters<dim> &global_parameters);

    /**
     * @brief Copy constructor for StressPostprocessor objects
     *
     * @param other Other StressPostprocessor
     */
    StressPostprocessor(const StressPostprocessor<dim> &other) = default;

    /**
     * @brief Evaluate the vector field and
     *        construct the linearized stress tensor.
     *
     * @param input_data Input data
     * @param computed_quantities Vector which will be overridden with the
     *                            linearized stress tensor.
     */
    virtual void
    evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &input_data,
      std::vector<Vector<double>> &computed_quantities) const override;

    /**
     * @brief Specify how the computed contities can be interpreted.
     *
     * @return std::vector<
     * DataComponentInterpretation::DataComponentInterpretation>
     */
    virtual std::vector<
      DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation() const override;

    /**
     * @brief Specify the names of the components of the computed data
     *
     * @return std::vector<std::string>
     */
    virtual std::vector<std::string>
    get_names() const override;

    /**
     * @brief Specify which update flags are needed
     *        for the shape functions in
     *        get_data_component_interpretation()
     *
     * @return UpdateFlags
     */
    virtual UpdateFlags
    get_needed_update_flags() const override;

  private:
    /**
     * String that contains the index of the local shape
     * shape function in ElaBasis.
     */
    std::string           basis_str;
    GlobalParameters<dim> parameters;
  };
} // namespace Elasticity

#endif // _INCLUDE_POSTPROCESSING_H_