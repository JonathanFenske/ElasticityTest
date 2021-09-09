#ifndef _INCLUDE_POSTPROCESSING_H_
#define _INCLUDE_POSTPROCESSING_H_

#include <deal.II/numerics/data_postprocessor.h>

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
  };


  // Class that enables the output of the linearized stress tensor.
  // Taken from the documentation for the deal.ii class DataPostprocessorTensor.
  template <int dim>
  class StressPostprocessor : public DataPostprocessor<dim>
  {
  public:
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
  };
}

#endif // _INCLUDE_POSTPROCESSING_H_