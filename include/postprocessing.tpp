#ifndef _INCLUDE_POSTPROCESSING_TPP_
#define _INCLUDE_POSTPROCESSING_TPP_

#include "postprocessing.h"

namespace Elasticity
{
  using namespace dealii;

  /****************************************************************************/
  /* Postprocessing */

  template <int dim>
  StrainPostprocessor<dim>::StrainPostprocessor(unsigned int basis_index)
    : basis_str("_" + Utilities::int_to_string(basis_index, 2))
  {}


  template <int dim>
  StrainPostprocessor<dim>::StrainPostprocessor()
    : basis_str("")
  {}


  template <int dim>
  void
  StrainPostprocessor<dim>::evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &input_data,
    std::vector<Vector<double>>                &computed_quantities) const
  {
    AssertDimension(input_data.solution_gradients.size(),
                    computed_quantities.size());
    // mu<dim> mu;
    // lambda<dim> lambda;
    // std::vector<double> mu_values(input_data.evaluation_points.size()),
    // lambda_values(input_data.evaluation_points.size());
    // mu.value_list(input_data.evaluation_points, mu_values);
    // lambda.value_list(input_data.evaluation_points, lambda_values);
    for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {
        AssertDimension(computed_quantities[p].size(),
                        (Tensor<2, dim>::n_independent_components));
        for (unsigned int d = 0; d < dim; ++d)
          {
            for (unsigned int e = 0; e < dim; ++e)
              {
                computed_quantities[p]
                                   [Tensor<2, dim>::component_to_unrolled_index(
                                     TableIndices<2>(d, e))] =
                                     (input_data.solution_gradients[p][d][e] +
                                      input_data.solution_gradients[p][e][d]) /
                                     2;
                // if (d == e)
                // {
                //   computed_quantities[p][
                //       Tensor<2,dim>::component_to_unrolled_index(
                //         TableIndices<2>(dim+d,dim+e))]
                //     = (lambda_values[p] + 2*mu_values[p])
                //       * input_data.solution_gradients[p][d][e];
                //     for (unsigned int i=1; i<dim; ++i)
                //     {
                //       computed_quantities[p][
                //         Tensor<2,dim>::component_to_unrolled_index(
                //           TableIndices<2>(dim+d,dim+e))]
                //         += lambda_values[p]
                //           * input_data.solution_gradients[p][(d+i)%dim][(e+i)%dim];
                //     }
                // }
                // else
                // {
                //   //std::cout << p << " : " << dim+d << "," << dim+e <<
                //   std::endl; computed_quantities[p][
                //       Tensor<2,dim>::component_to_unrolled_index(
                //         TableIndices<2>(dim+d,dim+e))]
                //     = mu_values[p]
                //       * (input_data.solution_gradients[p][d][e]
                //         + input_data.solution_gradients[p][e][d]);
                // }
              }
          }
      }
  }


  // function that specifies how many entries of the strain
  // belong to the tensor representing the strain
  template <int dim>
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  StrainPostprocessor<dim>::get_data_component_interpretation() const
  {
    return std::vector<
      DataComponentInterpretation::DataComponentInterpretation>(
      dim * dim, DataComponentInterpretation::component_is_part_of_tensor);
  }


  template <int dim>
  std::vector<std::string>
  StrainPostprocessor<dim>::get_names() const
  {
    // std::vector<std::string> names(dim*dim, "strain");
    // std::vector<std::string> stress(dim*dim, "stress");
    // names.insert(names.end(), stress.begin(), stress.end());
    // return names;
    return std::vector<std::string>(dim * dim, "strain" + basis_str);
  }


  template <int dim>
  UpdateFlags
  StrainPostprocessor<dim>::get_needed_update_flags() const
  {
    return update_gradients;
  }


  template <int dim>
  StressPostprocessor<dim>::StressPostprocessor()
    : basis_str("")
  {}


  template <int dim>
  StressPostprocessor<dim>::StressPostprocessor(
    unsigned int              basis_index,
    const ElaParameters<dim> &ela_parameters)
    : basis_str("_" + Utilities::int_to_string(basis_index, 2))
    , parameters(ela_parameters)
  {}


  template <int dim>
  StressPostprocessor<dim>::StressPostprocessor(
    const ElaParameters<dim> &ela_parameters)
    : basis_str("")
    , parameters(ela_parameters)
  {}


  // function that computes the linearized stress (Hooke's law)
  template <int dim>
  void
  StressPostprocessor<dim>::evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &input_data,
    std::vector<Vector<double>>                &computed_quantities) const
  {
    AssertDimension(input_data.solution_gradients.size(),
                    computed_quantities.size());
    std::vector<double> mu_values(input_data.evaluation_points.size()),
      lambda_values(input_data.evaluation_points.size());
    parameters.mu->value_list(input_data.evaluation_points, mu_values);
    parameters.lambda->value_list(input_data.evaluation_points, lambda_values);
    for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
      {
        AssertDimension(computed_quantities[p].size(),
                        (Tensor<2, dim>::n_independent_components));
        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            {
              if (d == e)
                {
                  computed_quantities
                    [p][Tensor<2, dim>::component_to_unrolled_index(
                      TableIndices<2>(d, e))] =
                      (lambda_values[p] + 2 * mu_values[p]) *
                      input_data.solution_gradients[p][d][e];
                  for (unsigned int i = 1; i < dim; ++i)
                    computed_quantities
                      [p][Tensor<2, dim>::component_to_unrolled_index(
                        TableIndices<2>(d, e))] +=
                      lambda_values[p] *
                      input_data
                        .solution_gradients[p][(d + i) % dim][(e + i) % dim];
                }
              else
                {
                  computed_quantities
                    [p][Tensor<2, dim>::component_to_unrolled_index(
                      TableIndices<2>(d, e))] =
                      mu_values[p] * (input_data.solution_gradients[p][d][e] +
                                      input_data.solution_gradients[p][e][d]);
                }
            }
      }
  }


  // function that specifies how many entries of the stress
  // belong to the tensor representing the stress
  template <int dim>
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  StressPostprocessor<dim>::get_data_component_interpretation() const
  {
    return std::vector<
      DataComponentInterpretation::DataComponentInterpretation>(
      dim * dim, DataComponentInterpretation::component_is_part_of_tensor);
  }


  template <int dim>
  std::vector<std::string>
  StressPostprocessor<dim>::get_names() const
  {
    return std::vector<std::string>(dim * dim, "stress" + basis_str);
  }


  template <int dim>
  UpdateFlags
  StressPostprocessor<dim>::get_needed_update_flags() const
  {
    return update_gradients | update_quadrature_points;
  }
} // namespace Elasticity

#endif // _INCLUDE_POSTPROCESSING_TPP_