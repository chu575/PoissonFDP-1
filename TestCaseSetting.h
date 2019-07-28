#ifndef _TestCaseSetting_H_                 
#define _TestCaseSetting_H_                 
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/sparse_direct.h>                       
#include <deal.II/numerics/vector_tools.h>            
#include <fstream>
#include <iostream>
#include <cmath>                                      


namespace MMSassignment 
{  
  using namespace dealii;

  // Exact solution of this test case. In MMS, the assignment of the exact solution is the first priority 
  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution () : Function<dim>(dim) {}

    virtual double value (const Point<dim> &p,
                          const unsigned int component = 0) const;
    
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;
  };
  template class ExactSolution<2>;


  // Right hand side of U equation
  template <int dim>
  class RightHandSide_u : public Function<dim>                       
  {
    public:
      RightHandSide_u () : Function<dim>(dim) {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;
  };
  template class RightHandSide_u<2>;


  // Right hand side of Lambda equation
  template <int dim>
  class RightHandSide_l : public Function<dim>                          
  {
    public:
      RightHandSide_l () : Function<dim>(dim) {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;
  };
  template class RightHandSide_l<2>;


  // define range for error computation, value function is for L2 norm of error, vector value function is for H1 semi-norm of error
  template <int dim>
  class ErrorComputationRange : public Function<dim>
  {
  public:
    ErrorComputationRange () : Function<dim>(dim){}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
 
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;   
  };
  template class ErrorComputationRange<2>;
}// end of namespace MMSassignment
  

namespace IBimplementation
{
  using namespace dealii;
 
  template <int dim>
  class IBlocation : public Function<dim>
  {
    public:
      IBlocation () : Function<dim>(dim) {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;
  }; 
  template class IBlocation<2>;     
}// end of namespace IB implementation

#endif          
