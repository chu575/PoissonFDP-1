#include "TestCaseSetting.h"

namespace IBimplementation
{
  // In approach 1, immersed boundary condition is not only valid on the immersed boundary, 
  // but also valid in the area outside the needed domain
  template <int dim>
  double IBlocation<dim>::value (const Point<dim>  &p,
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
}// end of namespace IBimplementation











