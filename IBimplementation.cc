#include "TestCaseSetting.h"

namespace IBimplementation
{
  // In approach 1, immersed boundary condition is not only valid on the immersed boundary, 
  // but also valid in the area outside the needed domain
  template <int dim>
  double IBlocation_approach1<dim>::value (const Point<dim>  &p,
                                  const unsigned int /*component*/) const
  {
    int value = 0;
    const Point<dim> center (0.0,0.0);
    double dis = center.distance(p);
    double ri = 0.6;                     
   
    if (dis >= ri) 
      value = 1;           // Immersed boundary condition is valid when IBCweightfactor is 1

    return value;
  }   

  template <int dim>
  double smoothfunction_approach1<dim>::value (const Point<dim>  &p,
				               const unsigned int component) const
  {
    double value = 1.0;
            
    return value;
  } 


  // In approach 1, immersed boundary condition is valid in a confined region. In addition, a smooth function is needed to distribute 
  // the contribution from the immersed boundary condition
  template <int dim>
  double IBlocation_approach2<dim>::value (const Point<dim>  &p,
                                           const unsigned int /*component*/) const
  {
    int value = 0;
    const Point<dim> center (0.0,0.0);
    double dis = center.distance(p);
    double ri = 0.6;                     
    double half_width = 0.0, mesh_size = 0.0, temp = 0.0;
    temp = pow(2,refinement_number);
    mesh_size = 2.0/temp;
    half_width = width_const * mesh_size;

    if (dis >= ri-half_width && dis <= ri+half_width) 
      value = 1;
     
    return value;
  }   

  
  template <int dim>
  double smoothfunction_approach2_const<dim>::value (const Point<dim>  &p,
				                     const unsigned int component) const
  {
    double value = 0.0;
    double total_width = 0.0, mesh_size = 0.0, temp = 0.0;
    temp = pow(2,refinement_number);
    mesh_size = 2.0/temp;
    total_width = 2 * width_const * mesh_size;
    value = 1.0 / total_width;                          
  
    return value; 
  } 

  
  template <int dim>
  double smoothfunction_approach2_trian<dim>::value (const Point<dim>  &p,
				                     const unsigned int component) const
  {
    const Point<dim> center (0.0,0.0);
    double dis = center.distance(p);
    double ri = 0.6;
    
    double value = 0.0;
    double half_width = 0.0, mesh_size = 0.0, temp = 0.0;
    temp = pow(2,refinement_number);
    mesh_size = 2.0/temp;
    half_width = width_const * mesh_size;
    value = (1.0-fabs(dis-ri)/half_width)*(1.0/half_width);                          

    return value; 
  } 


  template <int dim>
  double smoothfunction_approach2_gauss<dim>::value (const Point<dim>  &p,
				                     const unsigned int component) const
  {
    const Point<dim> center (0.0,0.0);
    double dis = center.distance(p);
    double ri = 0.6;
    
    double value = 0.0;
    double half_width = 0.0, mesh_size = 0.0, temp = 0.0;
    temp = pow(2,refinement_number);
    mesh_size = 2.0/temp;
    half_width = width_const * mesh_size;

    double d = 1./3.;            // d = 1/3 means boundary region gamma_bar contains 6 standard deviations                    
     
    double standard_deviation = 0.0;                                                                                                          
    standard_deviation = d * half_width;

    // The following computation can be checked in any data regarding Gaussian function
    double pi = 3.1415926535897;
    double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    tmp1 = dis-ri;
    tmp2 = -pow(tmp1,2.0);
    tmp3 = 2 * pow(standard_deviation,2);
    tmp4 = tmp2/tmp3;
    tmp5 = exp(tmp4);
    tmp6 = standard_deviation * pow(2*pi,0.5);
    value = tmp5/tmp6;    

    return value; 
  } 
}// end of namespace IBimplementation











