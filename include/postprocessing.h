#ifndef _INCLUDE_POSTPROCESSING_H_
#define _INCLUDE_POSTPROCESSING_H_

#include <deal.II/base/utilities.h>

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
    StrainPostprocessor(unsigned int basis_index) :
      basis_str("_" + Utilities::int_to_string(basis_index,2))
      {}

    StrainPostprocessor() :
      basis_str("")
      {}

    StrainPostprocessor(const StrainPostprocessor<dim> &other) :
      basis_str(other.basis_str)
      {}

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
    StressPostprocessor(unsigned int basis_index) :
      basis_str("_" + Utilities::int_to_string(basis_index,2))
      {}

    StressPostprocessor() :
      basis_str("")
      {}

    StressPostprocessor(const StressPostprocessor<dim> &other) :
      basis_str(other.basis_str)
      {}
    
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
}

#endif // _INCLUDE_POSTPROCESSING_H_