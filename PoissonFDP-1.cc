/* This code is to solve the designed Poisson problem using the fictitious domain method. 
 * The approach used in this implementation in Approach1 which extends the contraint to be valid not only on the original domain (i.e.,   
 * immersed boundary), but also valid in the area outtside the needed domain. Detailed descriptions can be seen in Chapter 2 of my PhD   
 * dissertation 
 * In 2019 July, I decided to modified the original code. The objective of this modification is to make the code simple and more organized. */


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
#include <deal.II/lac/sparse_direct.h>           //for SparseUMFPACK             
#include <deal.II/numerics/vector_tools.h>       // for error computation         
#include <fstream>
#include <iostream>
#include <cmath>                                 // for absolute value and exponential function computation 
#include "TestCaseSetting.h"                     
                  

namespace PoissonFictitiousDomainProblem
{
  using namespace dealii;
  using namespace MMSassignment;                
  using namespace IBimplementation;             

  template <int dim>
  class PFDP
  {
    public:
      PFDP (const unsigned int degree, const unsigned int refinement_number);  
      void run ();

    private:
      void setup_dofs ();
      void assemble_system ();
      void compute_errors () const;                                          
      void solve ();
      void output_results () const;

      const unsigned int   degree;
      const unsigned int   refinement_number;   

      Triangulation<dim>   triangulation;
      FESystem<dim>        fe;
      DoFHandler<dim>      dof_handler;

      BlockSparsityPattern      sparsity_pattern;
      BlockSparseMatrix<double> system_matrix;

      BlockVector<double>       solution;
      BlockVector<double>       system_rhs;
  };
  

  template <int dim>
  PFDP<dim>::PFDP (const unsigned int degree, const unsigned int refinement_number)
    :
    degree (degree),
    fe (FE_Q<dim>(degree), 1,                                                                      
	FE_Q<dim>(degree), 1),                                                                     
    dof_handler (triangulation),
    refinement_number (refinement_number)         
  {}


  template <int dim>
  void PFDP<dim>::setup_dofs ()
  {
    dof_handler.distribute_dofs (fe);

    DoFRenumbering::component_wise (dof_handler);

    std::vector<unsigned int> dofs_per_component (dim);
    DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);
    const unsigned int n_u = dofs_per_component[0],
                       n_l = dofs_per_component[1];

    std::cout << "Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl
              << "Total number of cells: "
              << triangulation.n_cells()
              << std::endl
              << "Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << " (" << n_u << '+' << n_l << ')'
              << std::endl;

    const unsigned int
    n_couplings = dof_handler.max_couplings_between_dofs();

    sparsity_pattern.reinit (2,2);
    sparsity_pattern.block(0,0).reinit (n_u, n_u, n_couplings);
    sparsity_pattern.block(1,0).reinit (n_l, n_u, n_couplings);
    sparsity_pattern.block(0,1).reinit (n_u, n_l, n_couplings);
    sparsity_pattern.block(1,1).reinit (n_l, n_l, n_couplings);
    sparsity_pattern.collect_sizes();

    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    sparsity_pattern.compress();

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (2);
    solution.block(0).reinit (n_u);
    solution.block(1).reinit (n_l);
    solution.collect_sizes ();

    system_rhs.reinit (2);
    system_rhs.block(0).reinit (n_u);
    system_rhs.block(1).reinit (n_l);
    system_rhs.collect_sizes ();
  }


  template <int dim>
  void PFDP<dim>::assemble_system ()
  {
    QGauss<dim>   quadrature_formula(degree+2);    
      
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    | update_gradients |
                             update_quadrature_points  | update_JxW_values);
    
    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();
    
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);
   
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    const RightHandSide_u<dim>          right_hand_side_u;                                           
    std::vector<double> rhs_values_u (n_q_points);                         

    const RightHandSide_l<dim>          right_hand_side_l;    
    std::vector<double> rhs_values_l (n_q_points);                         
                                            
    const IBlocation<dim> ib_location;                      
    std::vector<double> ib_location_weight (n_q_points);  
           
    const FEValuesExtractors::Scalar u (0);
    const FEValuesExtractors::Scalar lamda (1);                           

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      right_hand_side_u.value_list (fe_values.get_quadrature_points(), rhs_values_u);                            
      right_hand_side_l.value_list (fe_values.get_quadrature_points(), rhs_values_l);                                       
      ib_location.value_list (fe_values.get_quadrature_points(), ib_location_weight);
        
      for (unsigned int q=0; q<n_q_points; ++q)	 
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const double                phi_i_u      = fe_values[u].value (i, q);
          const Tensor<1,dim>         grad_phi_i_u = fe_values[u].gradient (i, q);
          const double                phi_i_l      = fe_values[lamda].value (i, q); 

          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const double             phi_j_u      = fe_values[u].value (j, q);
            const Tensor<1,dim>      grad_phi_j_u = fe_values[u].gradient (j, q);
            const double             phi_j_l      = fe_values[lamda].value (j, q);

            local_matrix(i,j) += (  grad_phi_i_u * grad_phi_j_u
                 		  - ib_location_weight[q] * phi_i_u * phi_j_l     
				  + ib_location_weight[q] * phi_i_l * phi_j_u                
                                  + (1-ib_location_weight[q]) * phi_i_l * phi_j_l                                   
                                 ) * fe_values.JxW(q);                                       
        }

          local_rhs(i) += (  phi_i_u * rhs_values_u[q]   
	                   + ib_location_weight[q] * phi_i_l * rhs_values_l[q]
		           + (1-ib_location_weight[q]) * phi_i_l * 0.0                                            
                          ) * fe_values.JxW(q);
        }
       
        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               local_matrix(i,j));

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += local_rhs(i);
    }   
  }


  // The solver used here is UMFPACK. This function comes from the tutorial step29
  template <int dim>              
  void PFDP<dim>::solve ()
  {
    SparseDirectUMFPACK  A_direct;
    A_direct.initialize(system_matrix);

    A_direct.vmult (solution, system_rhs);
  }


  template <int dim>
  void PFDP<dim>::compute_errors () const
  {
    ExactSolution<dim> exact_solution;
   
    Vector<double> cellwise_errors (triangulation.n_active_cells());

    QTrapez<1>     q_trapez;
    QIterated<dim> quadrature (q_trapez, degree+2);   

    ErrorComputationRange<dim> ecr;                          
       
    VectorTools::integrate_difference (dof_handler, 
                                       solution,
                                       exact_solution,
                                       cellwise_errors, 
                                       quadrature,
                                       VectorTools::Linfty_norm,
                                       &ecr
                                      );
    const double u_linfty_error = cellwise_errors.linfty_norm();
    std::cout << "Errors: ||e_u||_Linf = " << u_linfty_error << std::endl; 

    VectorTools::integrate_difference (dof_handler, 
                                       solution,
                                       exact_solution,
                                       cellwise_errors, 
                                       quadrature,
                                       VectorTools::L1_norm,
                                       &ecr
                                      );
    const double u_l1_error = cellwise_errors.l1_norm();
    std::cout << "Errors: ||e_u||_L1 = " << u_l1_error << std::endl;
    
    
    VectorTools::integrate_difference (dof_handler, 
                                       solution,
                                       exact_solution,
                                       cellwise_errors, 
                                       quadrature,
                                       VectorTools::L2_norm,
                                       &ecr
                                      );
    const double u_l2_error = cellwise_errors.l2_norm();
    std::cout << "Errors: ||e_u||_L2 = " << u_l2_error << std::endl;

    VectorTools::integrate_difference (dof_handler,
                                       solution,
                                       exact_solution,
                                       cellwise_errors,
                                       quadrature,
                                       VectorTools::H1_seminorm,
                                       &ecr
                                      );
    const double u_H1_error = cellwise_errors.l2_norm();
    std::cout << "Errors: ||e_u||_H1 = " << u_H1_error << std::endl;
  }
 

  // This is based on the function "output_results" in step-21
  template <int dim>
  void PFDP<dim>::output_results () const
  {
    std::vector<std::string> solution_names;

    switch (dim)
    {
      case 2:
        solution_names.push_back ("u");
        solution_names.push_back ("l");
        break;

      case 3:
        solution_names.push_back ("u");
        solution_names.push_back ("l");
        break;

      default:
        Assert (false, ExcNotImplemented());
    }

    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, solution_names);

    data_out.build_patches (degree);                              // if we choose degree+1, the outer surface might not be very smooth

    std::ofstream output ("solution.vtk");
    data_out.write_vtk (output);

  }


  template <int dim>
  void PFDP<dim>::run ()
  {
    GridGenerator::hyper_cube (triangulation, -1, 1);
    triangulation.refine_global (refinement_number);                                                               

    setup_dofs ();
    
    assemble_system ();
    
    solve ();
    
    compute_errors ();
    
    output_results ();
  }
}


int main ()
{
  using namespace dealii;
  using namespace PoissonFictitiousDomainProblem;    

  deallog.depth_console (0);

  unsigned int refinement_number = 0;                       
  std::cout << " mesh refinement number= " << std::endl;   
  std::cin >> refinement_number;          
  unsigned int degree = 0;                  
  std::cout << " degree for Finite Element class= " << std::endl;   // degree 1 represents Q1 element
  std::cin >> degree;          

  std::cout << "start to test this code" << std::endl;

  PFDP<2> pfdp(degree, refinement_number);  
  pfdp.run ();
  
  return 0;
}
