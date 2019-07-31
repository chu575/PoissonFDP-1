#include "TestCaseSetting.h"

namespace MMSassignment
{
  // This function describes the exact solution in this test case. It is valid in the needed domain only
  template <int dim>
  double ExactSolution<dim>::value (const Point<dim> &p,
                                     const unsigned int /*component*/) const
  {
    double value = 0.0;
    const Point<dim> center (0.0,0.0);
    double dis = center.distance(p);
    double g1 = 6;
    double ri = 0.6;

    if (dis < ri)    
      value = g1*(p[0]*p[0]+p[1]*p[1])/(ri*ri);
    
    return value;
  }
 

  // gradient equals to g1*(2*x)/(ri*ri) in x-direction and g1*(2*y)/(ri*ri) in y-direction. It is valid in the needed domain only
  template <int dim>
  Tensor<1,dim> ExactSolution<dim>::gradient (const Point<dim>   &p,
                                               const unsigned int) const
  { 
    Tensor<1,dim> exact_solution_gradient;
  
    const Point<dim> center (0.0,0.0);
    double dis = center.distance(p);
    double g1 = 6;
    double ri = 0.6;
 
    if (dis < ri)
    {
      exact_solution_gradient[0] = g1*(2*p[0])/(ri*ri); 
      exact_solution_gradient[1] = g1*(2*p[1])/(ri*ri);
    }
            
    return exact_solution_gradient;
  } 


  // This function is set for the R.H.S. of the u-equation, which is f_tilda in the dissertation.
  // However, here, I plan to let the weight factor, IBlocation to decide its value
  template <int dim>
  double RightHandSide_u<dim>::value (const Point<dim>  &p,
				     const unsigned int /*component*/) const
  {
    double value;  //////
    const Point<dim> center (0.0,0.0);
    double dis = center.distance(p);
    double ri = 0.6;    // radius of the needed domain
    double g1 = 6.0;    // constraint on the immersed boundary

    value= -4*g1/(ri*ri);
     
    return value;            //////
  } 


  // This class is set for the R.H.S. of the Lambda-equation, which is g_tilda in the dissertation
  // However, here, I plan to let the weight factor, IBlocation to decide its value
  template <int dim>
  double RightHandSide_l<dim>::value (const Point<dim>  &p,
                                    const unsigned int /*component*/) const
  {
    double value = 6.0;      //////                                                                                                                           
    return value;
  }

 
  template <int dim>
  double ErrorComputationRange<dim>::value (const Point<dim>  &p,
                                    const unsigned int component) const
  {
    int value = 0;
    const Point<dim> center (0.0,0.0);
    double ri = 0.6;
    double dis = center.distance(p);                                                                                                          
        
    if (component == 0)                                                                          
      if (dis < ri) 
        value = 1;
     
    return value;
  }


  template <int dim>                                                                       
  void ErrorComputationRange<dim>::vector_value (const Point<dim> &p,                                 
                                         Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = ErrorComputationRange<dim>::value (p, c);
  }  
}// end of namespace MMSassignment











