// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================
//
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    giovanni.filomeno@tum.de   $
//   Last modified by:    $Co-Author: giovanni.filomeno@tum.de   $
//   Date:                $Date:                      Decem 2016 $
//   Revision:            $Revision:                         0.0 $
//
// ==============================================================================

#ifndef CAD_MAPPER_H
#define CAD_MAPPER_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../../kratos/includes/define.h"
#include "../../kratos/processes/process.h"
#include "../../kratos/includes/node.h"
#include "../../kratos/includes/element.h"
#include "../../kratos/includes/model_part.h"
#include "../../kratos/includes/kratos_flags.h"
#include "../../kratos/spatial_containers/spatial_containers.h"
#include "../../kratos/utilities/binbased_fast_point_locator.h"
#include "linear_solvers/linear_solver.h"
#include "patch.h"
#include "control_point.hpp"
#include "shape_optimization_application.h"
// ==============================================================================

namespace Kratos
{
class CADMapper
{
  public:
    ///@name Type Definitions
    ///@{

    // ==========================================================================
    // Type definitions for linear algebra including sparse systems
    // ==========================================================================
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef typename SparseSpaceType::MatrixType SparseMatrixType;
    typedef typename SparseSpaceType::VectorType VectorType;
    typedef std::vector<double> DoubleVector;
    typedef std::vector<int> IntVector;
    typedef std::vector<ControlPoint> ControlPointVector;
	typedef std::vector<Patch> PatchVector;
	typedef Node<3> NodeType;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
	

	// For python extraction
    typedef boost::python::extract<double> takeDouble;
    typedef boost::python::extract<int> takeInt;
    typedef boost::python::extract<bool> takeBool;
	
	// for tree search
	typedef std::vector<double> DistanceVector;
    typedef std::vector<double>::iterator DistanceIterator;
	    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /// Pointer definition of CADMapper
    KRATOS_CLASS_POINTER_DEFINITION(CADMapper);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CADMapper(ModelPart& fe_model_part, boost::python::dict cad_model, unsigned int domain_size, LinearSolverType::Pointer linear_solver)
	: mr_fe_model_part(fe_model_part),
	  mr_cad_model(cad_model),
	  m_domain_size(domain_size),
	  m_linear_solver(linear_solver)
    {
		// Initialize CAD model
		initialize_patches();

		// Initialize FEM model
		initialize_fe_model_for_mapping();
    }

    /// Destructor.
    virtual ~CADMapper()
    {
    }

    // --------------------------------------------------------------------------
    void initialize_patches()
    {
		// Counter for control points which provides a global Id for each control point
		m_n_control_points = 0;

		// loop over faces
		for (int i = 0; i < boost::python::len(mr_cad_model["faces"]); i++)
		{
			// each face is a patch: mr_cad_model["faces"][i]

			// 1. initialize "ingredients" ==================================
			DoubleVector knot_vector_u;
			DoubleVector knot_vector_v;
			int p;
			int q;
			ControlPointVector control_points;

			// 2. create "ingredients" ======================================
			// fill in knot_vector_u
			for (int u_idx = 0;
			u_idx < boost::python::len(mr_cad_model["faces"][i]["surface"][0]["knot_vectors"][0]);
			u_idx++)
			{
			// read knot
			double knot = takeDouble(mr_cad_model["faces"][i]["surface"][0]["knot_vectors"][0][u_idx]);
			knot_vector_u.push_back(knot);
			}

			//    	    		for( auto entris : knot_vector_u)
			//    	    		{
			//    	    			std::cout << "my test: " << entris << std::endl;
			//    	    		}
			// fill in knot_vector_v

			for (int v_idx = 0; v_idx < boost::python::len(mr_cad_model["faces"][i]["surface"][0]["knot_vectors"][1]); v_idx++)
			{
			// read knot
			double knot = takeDouble(mr_cad_model["faces"][i]["surface"][0]["knot_vectors"][1][v_idx]);
			knot_vector_v.push_back(knot);
			}

			//    	    		for( auto entris : knot_vector_v)
			//    	    		{
			//    	    			std::cout << "my test: " << entris << std::endl;
			//    	    		}
			// get p and q
			p = takeInt(mr_cad_model["faces"][i]["surface"][0]["degrees"][0]);
			q = takeInt(mr_cad_model["faces"][i]["surface"][0]["degrees"][1]);

			// fill in control_points

			// Control points in each patch get a global and local ID as well as a mapping matrix id
			// global Id: Unique Id for each control point (given by m_n_control_points)
			// local Id: given by input .json file
			// mapping matrix id: specifies position in global mapping matrix
			for (int cp_idx = 0; cp_idx < boost::python::len(mr_cad_model["faces"][i]["surface"][0]["control_points"]); cp_idx++)
			{
				//    	    			std::cout << cp_idx << std::endl;

				unsigned int local_id = takeInt(mr_cad_model["faces"][i]["surface"][0]["control_points"][cp_idx][0]);
				unsigned int global_id = ++m_n_control_points;
				unsigned int mapping_matrix_id = global_id-1;
				double x = takeDouble(mr_cad_model["faces"][i]["surface"][0]["control_points"][cp_idx][1][0]);
				double y = takeDouble(mr_cad_model["faces"][i]["surface"][0]["control_points"][cp_idx][1][1]);
				double z = takeDouble(mr_cad_model["faces"][i]["surface"][0]["control_points"][cp_idx][1][2]);
				double w = takeDouble(mr_cad_model["faces"][i]["surface"][0]["control_points"][cp_idx][1][3]);

				// read knot
				ControlPoint myControlPoint(x, y, z, w, local_id, global_id, mapping_matrix_id);
				control_points.push_back(myControlPoint);
			}

			// 3. create patch  =============================================
			Patch patch(knot_vector_u, knot_vector_v, p, q, control_points);
			m_patches.push_back(patch);
		}
    }

	// --------------------------------------------------------------------------
	void initialize_fe_model_for_mapping()
	{
		// Assign each FE node a uniqe Id in the mapping matrix (according to iterator over nodes)
		m_n_fem_points = 0;
		for (ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i != mr_fe_model_part.NodesEnd(); ++node_i)
			node_i->SetValue(MAPPING_MATRIX_ID,m_n_fem_points++);
	}

    // --------------------------------------------------------------------------
    void compute_mapping_matrix(const unsigned int u_resolution, const unsigned int v_resolution)
    {
		std::cout << "\n> Starting computation of mapping matrix..." << std::endl;
		boost::timer function_timer;

		// Initialize matrices to zero matrices of correct size
		m_mapping_matrix_CAD_CAD.resize(m_n_control_points,m_n_control_points);
		m_mapping_matrix_CAD_FEM.resize(m_n_control_points,m_n_fem_points); 
		m_mapping_matrix_CAD_CAD.clear();
		m_mapping_matrix_CAD_FEM.clear(); 		

		// 1st step: Create for each patch coarse cloud of CAD points in x-y space for neighbor search later 

		// Each correspondingly created point is stored in a list.
		// Further lists in the same order are created to store the respective u & v parameter as well as the patch id
		// As list iterator we use a counter for the number of CAD nodes
		unsigned int cad_node_counter = 0;
		NodeVector list_of_cad_nodes;
		DoubleVector list_of_us_of_cad_nodes;
		DoubleVector list_of_vs_of_cad_nodes;
		IntVector list_of_patch_ids_of_cad_nodes;

		//Loop over all patches / faces
		for (int patch_id = 0; patch_id < boost::python::len(mr_cad_model["faces"]); patch_id++)
		{

			std::cout << "\n> Processing Patch " << patch_id << std::endl;

			DoubleVector& knot_vec_u_i = m_patches[patch_id].get_knot_vector_u();
			DoubleVector& knot_vec_v_i = m_patches[patch_id].get_knot_vector_v();
			unsigned int knot_vector_u_dimension = knot_vec_u_i.size();
			unsigned int knot_vector_v_dimension = knot_vec_v_i.size();

			double u_min = knot_vec_u_i[0];
			double u_max = knot_vec_u_i[knot_vector_u_dimension-1];
			double v_min = knot_vec_v_i[0];
			double v_max = knot_vec_v_i[knot_vector_v_dimension-1];

			double delta_u = (u_max-u_min) / u_resolution;
			double delta_v = (v_max-v_min) / v_resolution;

			// Loop over all u & v according to specified resolution
			for(unsigned int i=0; i<=u_resolution; i++)
			{
				// current u-value
				double u_i = u_min + i*delta_u;

				for(unsigned int j=0; j<=v_resolution; j++)
				{
					// current v-value
					double v_j = v_min + j*delta_v;

					// compute unique point in CAD-model for given u&v
					++cad_node_counter;					
					Point<3> cad_point_coordinates;
					m_patches[patch_id].evaluate_surface_point(cad_point_coordinates, u_i, v_j);

					// Add id to point --> node. Add node to list of CAD nodes
					NodeType::Pointer new_cad_node = Node < 3 > ::Pointer(new Node<3>(cad_node_counter, cad_point_coordinates));
					list_of_cad_nodes.push_back(new_cad_node);

					// Store for cad node the corresponding cad information in separate vectors
					list_of_us_of_cad_nodes.push_back(u_i);
					list_of_vs_of_cad_nodes.push_back(v_j);
					list_of_patch_ids_of_cad_nodes.push_back(patch_id);
				}
			}
		}

		// 2nd step: Construct KD-Tree with all cad nodes
		std::cout << "\n> Starting construction of search-tree..." << std::endl;
		boost::timer timer;
        typedef Bucket< 3, NodeType, NodeVector, NodeType::Pointer, NodeIterator, DistanceIterator > BucketType;
        typedef Tree< KDTreePartition<BucketType> > tree;
        int bucket_size = 20;
        tree nodes_tree(list_of_cad_nodes.begin(), list_of_cad_nodes.end(), bucket_size);
		std::cout << "> Time needed for constructing search-tree: " << timer.elapsed() << " s" << std::endl;

		// 3rd step: Compute mapping matrix

		// Loop over all integration points of fe-model-part and find corresponding closest neighbors of cad-model
		// We assume the surface in the fe-model to be mapped is described by conditions (e.g. ShapeOptimizationConditions)
        for (ModelPart::ConditionsContainerType::iterator cond_i = mr_fe_model_part.ConditionsBegin(); cond_i != mr_fe_model_part.ConditionsEnd(); ++cond_i)
        {
        	// Get geometry information of current condition
			Condition::GeometryType& geom_i = cond_i->GetGeometry();
			unsigned int n_fem_nodes = geom_i.size();

			// Get and store mapping matrix ids of nodes of current condition
			Vector mapping_matrix_ids_fem = ZeroVector(n_fem_nodes);
			for(unsigned int i=0; i<n_fem_nodes;i++)
				mapping_matrix_ids_fem[i] = geom_i[i].GetValue(MAPPING_MATRIX_ID);

			// Set integration method and evaluate shape functions of FE model accordingly
        	const Condition::GeometryType::IntegrationMethod integration_method = GeometryData::GI_GAUSS_2;
        	const Condition::GeometryType::IntegrationPointsArrayType& integration_points = geom_i.IntegrationPoints(integration_method);
			const unsigned int number_of_integration_points = integration_points.size();
        	const Matrix& N_container = geom_i.ShapeFunctionsValues(integration_method);

        	// std::cout << "########################################" << std::endl;
        	// std::cout << "cond_i->Id() = " << cond_i->Id() << std::endl;

        	for ( unsigned int PointNumber = 0; PointNumber < number_of_integration_points; PointNumber++ )
        	{
				// std::cout << "------------------------------" << std::endl;
        		// std::cout << "PointNumber = " << PointNumber << std::endl;

        		// Compute global coordinates of current integration point and get corresponding weight
        		NodeType::CoordinatesArrayType ip_coordinates = geom_i.GlobalCoordinates(ip_coordinates, integration_points[PointNumber].Coordinates());
        		NodeType::Pointer gauss_point_i = Node < 3 > ::Pointer(new Node<3>(PointNumber, ip_coordinates ));
				double integration_weight = integration_points[PointNumber].Weight();

        		// Search nearest cad neighbor of current integration point
        		NodeType resulting_nearest_point;
        		NodeType::Pointer nearest_point = nodes_tree.SearchNearestPoint( *gauss_point_i );

				// Recover CAD information of nearest point
				double u_of_nearest_point = list_of_us_of_cad_nodes[nearest_point->Id()-1];
				double v_of_nearest_point = list_of_vs_of_cad_nodes[nearest_point->Id()-1];
				int patch_id_of_nearest_point = list_of_patch_ids_of_cad_nodes[nearest_point->Id()-1];

				// Perform Newton-Raphson detailed search



				// To be developed by Massimo and Giovanni...

				// Some useful functions

				// //Computing inverse
				// Matrix H; // Hessian obtained from patch.h
				// Matrix InvH;
				// double det_H;
    			// MathUtils<double>::InvertMatrix( H, InvH, det_H );

				// // Multiplication
				// Vector some_vector;
				// Vector result = prod( invH, some_vector );

				// // Get NURBS shape functions derivatives
				// matrix<double> dR, matrix<double> ddR
				// m_patches[patch_id_of_nearest_point].deriv_Nurbsbasisfunc(0, 0, u_of_nearest_point, v_of_nearest_point, dR, ddR);


				// Prepare FEM side for computation of mapping matrix

        		// Get FEM-shape-function-value for current integration point
        		Vector N_FEM_GPi = row( N_container, PointNumber);
				
				// Prepare CAD side for computation of mapping matrix

				// Get CAD-shape-function-value for all control points affecting the nearest cad point
				// Additionally get the corresponding ids of control points in the mapping matrix
				matrix<double> R_CAD_Pi;
				matrix<unsigned int> mapping_matrix_ids_cad;
				m_patches[patch_id_of_nearest_point].eval_Nurbsbasis(0, 0, u_of_nearest_point, v_of_nearest_point, R_CAD_Pi, mapping_matrix_ids_cad);

        		// Some output
        		// // // KRATOS_WATCH(integration_points[PointNumber]);
        		KRATOS_WATCH(*gauss_point_i);
        		KRATOS_WATCH(*nearest_point);
				KRATOS_WATCH(u_of_nearest_point);
				KRATOS_WATCH(v_of_nearest_point);
				// KRATOS_WATCH(R_CAD_Pi);
				// KRATOS_WATCH(mapping_matrix_ids_cad);
        		// KRATOS_WATCH(N_FEM_GPi);
				// KRATOS_WATCH(mapping_matrix_ids_fem);

				// Assemble mapping matrix
        		
				// First we assemble CAD-FEM matrix
				for(unsigned int i=0; i<mapping_matrix_ids_cad.size1();i++)
				{
					for(unsigned int j=0; j<mapping_matrix_ids_cad.size2();j++)
					{

						double R_id = mapping_matrix_ids_cad(i,j);
						double R = R_CAD_Pi(i,j);

						for (unsigned int k=0; k<n_fem_nodes;k++)
						{
							double N_id = mapping_matrix_ids_fem[k];
							double N = N_FEM_GPi[k];

							m_mapping_matrix_CAD_FEM(R_id,N_id) += integration_weight * R * N;
						}
					}
				}

				// Then we assemble CAD-CAD matrix
				for(unsigned int i=0; i<mapping_matrix_ids_cad.size1();i++)
				{
					for(unsigned int j=0; j<mapping_matrix_ids_cad.size2();j++)
					{
						double R_row_id = mapping_matrix_ids_cad(i,j);
						double R_row = R_CAD_Pi(i,j);

						for(unsigned int k=0; k<mapping_matrix_ids_cad.size1();k++)
						{
							for(unsigned int l=0; l<mapping_matrix_ids_cad.size2();l++)
							{
								double R_coll_id = mapping_matrix_ids_cad(k,l);
								double R_coll = R_CAD_Pi(k,l);

								m_mapping_matrix_CAD_CAD(R_row_id,R_coll_id) += integration_weight * R_row * R_coll;
							}
						}
					}
				}
        	}
        }
		std::cout << "\n> Finished computation of mapping matrix in " << function_timer.elapsed() << " s." << std::endl;

		// Write CAD nodes into a file (only for verification)
		std::cout << "\n> Writing cad nodes..." << std::endl;
		std::ofstream temp_file("cad_node_coordinates.txt");
		for(NodeVector::iterator it =  list_of_cad_nodes.begin(); it!=list_of_cad_nodes.end(); it++)
		{
			NodeType::Pointer cad_node_i = *it;
			temp_file << cad_node_i->X() << " " << cad_node_i->Y() << " " << cad_node_i->Z() << std::endl;
		}
		temp_file.close();
    }

    // --------------------------------------------------------------------------
    void map_to_cad_space()
    {
		// Initialize vectors needed later
		Vector dx = ZeroVector(m_n_fem_points);
		Vector dy = ZeroVector(m_n_fem_points);
		Vector dz = ZeroVector(m_n_fem_points);
		Vector mapping_rhs_vector_x = ZeroVector(m_n_control_points);
		Vector mapping_rhs_vector_y = ZeroVector(m_n_control_points);
		Vector mapping_rhs_vector_z = ZeroVector(m_n_control_points);
		Vector dsx = ZeroVector(m_n_control_points);
		Vector dsy = ZeroVector(m_n_control_points);
		Vector dsz = ZeroVector(m_n_control_points);

		// Prepare RHS vector of mapping system of equation
		for (ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i != mr_fe_model_part.NodesEnd(); ++node_i)
		{
			unsigned int mapping_id = node_i->GetValue(MAPPING_MATRIX_ID);

			dx[mapping_id] = node_i->GetValue(SHAPE_UPDATE_X);
			dy[mapping_id] = node_i->GetValue(SHAPE_UPDATE_Y);
			dz[mapping_id] = node_i->GetValue(SHAPE_UPDATE_Z);
		}
		SparseSpaceType::Mult(m_mapping_matrix_CAD_FEM,dx,mapping_rhs_vector_x);
		SparseSpaceType::Mult(m_mapping_matrix_CAD_FEM,dy,mapping_rhs_vector_y);
		SparseSpaceType::Mult(m_mapping_matrix_CAD_FEM,dz,mapping_rhs_vector_z);


		// Solve linear systems to obtain mapped quantities each in X,Y,Z-direction separately
		// Note that an alternative would be to solve a big block structured matrix at once
		m_linear_solver->Solve(m_mapping_matrix_CAD_CAD, dsx, mapping_rhs_vector_x);
		m_linear_solver->Solve(m_mapping_matrix_CAD_CAD, dsy, mapping_rhs_vector_y);
		m_linear_solver->Solve(m_mapping_matrix_CAD_CAD, dsz, mapping_rhs_vector_z);

		std::cout << "\n> WARNING!!!!!!!!!!!!!!!!!!!: Think about solution precision. Why does the solution only show six digits and why there is a -0?....\n" << std::endl;

		// Sort solution (displacement of control points) into cad data set
		for(unsigned int i = 0; i<m_patches.size();i++)
		{
			ControlPointVector& control_points = m_patches[i].get_control_points();

			for(unsigned int j = 0; j<control_points.size();j++)
			{
				ControlPoint& cp_j = control_points[j];
				unsigned int cp_mapping_matrix_id = cp_j.getMappingMatrixId();

				cp_j.setdX( dsx[cp_mapping_matrix_id] );
				cp_j.setdY( dsy[cp_mapping_matrix_id] );
				cp_j.setdZ( dsz[cp_mapping_matrix_id] );
			}
		}

		// Test solution
		Vector rhs_test_z = ZeroVector(m_n_control_points);
		SparseSpaceType::Mult(m_mapping_matrix_CAD_CAD,dsz,rhs_test_z);
		KRATOS_WATCH(dz);
		KRATOS_WATCH(dsz);
		// KRATOS_WATCH(m_mapping_matrix_CAD_FEM);
		// KRATOS_WATCH(m_mapping_matrix_CAD_CAD);
		KRATOS_WATCH(mapping_rhs_vector_z);
		KRATOS_WATCH(rhs_test_z);
	}	

    // --------------------------------------------------------------------------
    void output_control_point_displacements()
    {
		// Outputs the displacements of the control points in a format that may be read by Gid

		std::cout << "\n> Starting to write displacement of control points nodes..." << std::endl;
		std::ofstream output_file("control_point_displacements.post.res");

		output_file << "Rhino Post Results File 1.0" << std::endl;
		output_file << "Result \"Displacement\" \"Load Case\" 0 Vector OnNodes" << std::endl;
		output_file << "Values" << std::endl;

		for(unsigned int i = 0; i<m_patches.size();i++)
		{
			ControlPointVector& control_points = m_patches[i].get_control_points();

			for(unsigned int j = 0; j<control_points.size();j++)
			{
				ControlPoint& cp_j = control_points[j];
				output_file << cp_j.getGlobalId() << " " << cp_j.getdX() << " " << cp_j.getdY() << " " << cp_j.getdZ() << std::endl;
			}
		}

		output_file << "End Values" << std::endl;

		output_file.close();
	}

    // --------------------------------------------------------------------------
    boost::python::list compute_point(unsigned int patch_id, double u, double v)
    {
		boost::python::list point;
		Point<3> cad_point_coordinates;
		m_patches[patch_id].evaluate_surface_point(cad_point_coordinates, u, v);

		point.append(cad_point_coordinates[0]);
		point.append(cad_point_coordinates[1]);
		point.append(cad_point_coordinates[2]);

		return point;
    }

    // ==============================================================================

    /// Turn back information as a string.
    virtual std::string Info() const
    {
		return "CADMapper";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
		rOStream << "CADMapper";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

  private:
    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ModelPart &mr_fe_model_part;
    boost::python::dict mr_cad_model;
	unsigned int m_domain_size;

    // ==============================================================================
    // General working arrays
    // ==============================================================================
    PatchVector m_patches;
	unsigned int m_n_control_points;
	unsigned int m_n_fem_points;
    SparseMatrixType m_mapping_matrix_CAD_CAD;
	SparseMatrixType m_mapping_matrix_CAD_FEM;

	// ==============================================================================
    // Solver and strategies
    // ==============================================================================
	LinearSolverType::Pointer m_linear_solver;

    /// Assignment operator.
    //      CADMapper& operator=(CADMapper const& rOther);

    /// Copy constructor.
    //      CADMapper(CADMapper const& rOther);

}; // Class CADMapper
} // namespace Kratos.

#endif // CAD_MAPPER_H
