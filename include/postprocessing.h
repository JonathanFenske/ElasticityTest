#ifndef _INCLUDE_POSTPROCESSING_H_
#define _INCLUDE_POSTPROCESSING_H_

#include <deal.II/base/utilities.h>

#include <deal.II/numerics/data_postprocessor.h>

#include "process_parameter_file.h"

namespace Elasticity
{
  using namespace dealii;

  /****************************************************************************/
  /* Postprocessing */

  // Class that enables the output of the linearized strain tensor.
  // Taken from the documentation for the deal.ii class DataPostprocessorTensor.
  // (The uncommented parts are my try to combine stress and strain
  // into a single class.)
  template <int dim>
  class StrainPostprocessor : public DataPostprocessor<dim>
  {
  public:
    StrainPostprocessor(unsigned int basis_index);

    StrainPostprocessor();

    StrainPostprocessor(const StrainPostprocessor<dim> &other);

    virtual void
    evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &input_data,
      std::vector<Vector<double>> &computed_quantities) const override;

    virtual std::vector<
      DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation() const override;

    virtual std::vector<std::string>
    get_names() const override;

    virtual UpdateFlags
    get_needed_update_flags() const override;

  private:
    std::string basis_str;
  };


  // Class that enables the output of the linearized stress tensor.
  // Taken from the documentation for the deal.ii class DataPostprocessorTensor.
  template <int dim>
  class StressPostprocessor : public DataPostprocessor<dim>
  {
  public:
    StressPostprocessor();

    StressPostprocessor(unsigned int                 basis_index,
                        const GlobalParameters<dim> &global_parameters);

    StressPostprocessor(const GlobalParameters<dim> &global_parameters);

    StressPostprocessor(const StressPostprocessor<dim> &other);

    // void
    // operator=(StressPostprocessor<dim> &other);

    virtual void
    evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &input_data,
      std::vector<Vector<double>> &computed_quantities) const override;

    virtual std::vector<
      DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation() const override;

    virtual std::vector<std::string>
    get_names() const override;

    virtual UpdateFlags
    get_needed_update_flags() const override;

  private:
    std::string           basis_str;
    GlobalParameters<dim> parameters;
  };
} // namespace Elasticity

#endif // _INCLUDE_POSTPROCESSING_H_