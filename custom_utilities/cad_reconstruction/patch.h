// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 [Released on march 05, 2007].

 Copyright [c] 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  [the
 "Software"], to  deal in  the Software without  restriction, including
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

#ifndef PATCH_H
#define PATCH_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <numeric>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "control_point.hpp"
// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

 */

class Patch
{
public:
	// ==========================================================================
	// Type definitions for linear algebra including sparse systems
	// ==========================================================================
	typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	typedef typename SparseSpaceType::MatrixType SparseMatrixType;
	typedef typename SparseSpaceType::VectorType VectorType;
	typedef std::vector<int> IntVector;
	typedef std::vector<double> DoubleVector;
	typedef boost::python::extract<double> takeDouble;
	typedef boost::python::extract<int> takeInt;
	typedef boost::python::extract<bool> takeBool;
	typedef std::vector<ControlPoint> ControlPointVector;

	/// Pointer definition of Patch
	//    KRATOS_CLASS_POINTER_DEFINITION[Patch];

	/// Default constructor.
	Patch(DoubleVector &m_knot_vector_u,
			DoubleVector &m_knot_vector_v,
			int p, int q,
			ControlPointVector &m_control_points) : m_knot_vector_u(m_knot_vector_u),
					m_knot_vector_v(m_knot_vector_v),
					m_p(p), m_q(q),
					m_control_points(m_control_points)
	{
		m_n_u = m_knot_vector_u.size() - m_p - 1;
		m_n_v = m_knot_vector_v.size() - m_q - 1;

		if (m_control_points.size() != m_n_u * m_n_v)
		{
			std::cout << "Invalid Patch" << std::endl;
		}
	}

	/// Destructor.
	virtual ~Patch()
	{
	}










	// ###################################################################################################
	// Algorithm of Massimo and Giovanni
	// ###################################################################################################
	// --------------------------------------------------------------------------
	// int getCPIndex(const int i, const int j)
	// {
	// 	return j * m_n_u + i;
	// }

	// // --------------------------------------------------------------------------
	// double R(const double &u, const double &v, const int &i, const int &j)
	// {
	// 	//    	for(int k=5; k>= 0 ; k--)
	// 	//    	{
	// 	//    	for(int i=0; i<k; i++)
	// 	//    	{
	// 	//    		std::cout << "i = " << i << std::endl;
	// 	//    		std::cout << N( m_knot_vector_v, 0, i, 5-k) << std::endl;
	// 	//    		std::cout << N( m_knot_vector_v, 2.5, i, 5-k) << std::endl;
	// 	//    		std::cout << N( m_knot_vector_v, 5.0, i, 5-k) << std::endl << std::endl;
	// 	//    	}
	// 	//    	}

	// 	// max index of knot vector (from 0 to m-1)

	// 	double lo = 0.0;
	// 	double n_tmp = 0.0;

	// 	//        std::cout << "I'm in R" << std::endl;

	// 	DoubleVector Nii;
	// 	DoubleVector Njj;

	// 	for (unsigned int jj = 0; jj < m_n_v; jj++)
	// 	{
	// 		auto Nj = N(m_knot_vector_v, v, jj, m_q);
	// 		Njj.push_back(Nj);
	// 	}

	// 	//        if( v == m_knot_vector_v.back())
	// 	//        {
	// 	//        	Njj[ m_n_v - 1] = 1.0;
	// 	//        }

	// 	for (unsigned int ii = 0; ii < m_n_u; ii++)
	// 	{
	// 		double N_ii = N(m_knot_vector_u, u, ii, m_p);
	// 		Nii.push_back(N_ii);
	// 	}

	// 	//        if( u == m_knot_vector_u.back())
	// 	//        {
	// 	//           	Nii[ m_n_u - 1] = 1.0;
	// 	//        }

	// 	for (unsigned int ii = 0; ii < m_n_u; ii++)
	// 	{
	// 		for (unsigned int jj = 0; jj < m_n_v; jj++)
	// 		{
	// 			n_tmp = Nii[ii] * Njj[jj];
	// 			int index = getCPIndex(ii, jj);
	// 			lo = lo + n_tmp * m_control_points[index].getWeight();
	// 			//        		lo = lo + n_tmp;
	// 		}
	// 	}

	// 	//        std::cout << "I'm in R" << std::endl;

	// 	if (fabs(lo) < TOL)
	// 	{
	// 		return 0;
	// 	}

	// 	int index = getCPIndex(i, j);
	// 	double up = N(m_knot_vector_u, u, i, m_p) * N(m_knot_vector_v, v, j, m_q) * m_control_points[index].getWeight();
	// 	//        double up = N(m_knot_vector_u, u, i, p) * N(m_knot_vector_v, v, j, q);
	// 	//        std:: cout << "I am in R " << up/lo << std::endl;
	// 	return up / lo;
	// }

	// // --------------------------------------------------------------------------
	// // Returns the R without weight in u/v direction

	// // U = knot vector
	// // t = parameter in the parameter space
	// // p = degree
	// // i = number of shape function

	// double N(DoubleVector &U, const double &t, const int &i, const int &p_degree)
	// {
	// 	// max index of knot vector (from 0 to m)
	// 	int m = U.size();

	// 	//		int lastZeroDegree = m - 1;

	// 	if (p_degree == 0 && t == U.back() && i == m - m_p - 2)
	// 	{
	// 		return 1.0;
	// 	}

	// 	//		std::cout << std::endl << "BEGIN" << std::endl;
	// 	//
	// 	//		std::cout << "N" << i << "," << p_degree << std::endl;

	// 	// check if i is in interval 0 <= i <= n. (with n = m-p-1)
	// 	// and if t is in the interval of the knot vector
	// 	if ((i < 0) || (i >= (m - m_p - 1)) || t < U.front() || t > U.back())
	// 	{

	// 		return 0.0;
	// 	}

	// 	// treat case of p_degree == 0
	// 	if (p_degree == 0)
	// 	{
	// 		// check if t is in interval from t_i <= t < t_(i+1)
	// 		if ((U[i] <= t) && (t < U[i + 1]))
	// 		{

	// 			return 1.0;
	// 		}
	// 		else
	// 		{

	// 			return 0.0;
	// 		}
	// 	}
	// 	else
	// 	{
	// 		double e1, e2, up, lo;

	// 		// compute division term 1
	// 		up = t - U[i];
	// 		lo = U[i + p_degree] - U[i];

	// 		if (fabs(lo) < TOL)
	// 		{
	// 			e1 = 0.0;
	// 		}
	// 		else
	// 		{
	// 			e1 = up / lo;
	// 		}

	// 		//			std::cout<<"Up=" <<up << " LO=" << lo <<std::endl;
	// 		// compute division term 2
	// 		up = U[i + p_degree + 1] - t;
	// 		lo = U[i + p_degree + 1] - U[i + 1];

	// 		if (fabs(lo) < TOL)
	// 		{
	// 			e2 = 0.0;
	// 		}
	// 		else
	// 		{
	// 			e2 = up / lo;
	// 		}
	// 		//			std::cout<<"Up=" <<up << " LO=" << lo <<std::endl;
	// 		//			std::cout <<"N" << i << "," << p_degree << "-->e1 " << e1 << " e2 " << e2 << std::endl;
	// 		// compute basis function recursively
	// 		double sn = 0.0;

	// 		if (fabs(e1) > TOL)
	// 		{
	// 			//				std::cout << "calling" << std::endl;
	// 			sn += e1 * N(U, t, i, p_degree - 1);
	// 		}

	// 		if (fabs(e2) > TOL)
	// 		{
	// 			//				std::cout << "calling" << std::endl;
	// 			sn += e2 * N(U, t, i + 1, p_degree - 1);
	// 		}

	// 		return sn;
	// 	}
	// }

	// // --------------------------------------------------------------------------
	// // computes the cartesian coordinates of the point for given parameters
	// void S(Point<3>& rResultingPoint, double u, double v)
	// {
	// 	rResultingPoint[0] = 0;
	// 	rResultingPoint[1] = 0;
	// 	rResultingPoint[2] = 0;

	// 	for (unsigned int i = 0; i < m_n_u; i++)
	// 	{
	// 		for (unsigned int j = 0; j < m_n_v; j++)
	// 		{
	// 			int index = getCPIndex(i, j);
	// 			double RR = R(u, v, i, j);
	// 			rResultingPoint[0] += RR * m_control_points[index].getX();
	// 			rResultingPoint[1] += RR * m_control_points[index].getY();
	// 			rResultingPoint[2] += RR * m_control_points[index].getZ();
	// 		}
	// 	}
	// }
	// ###################################################################################################
	// Eof algorithm of Massimo and Giovanni
	// ###################################################################################################











	//  #####################################################################################
	// #######################################################################################
	///
	///  \details    returns the cartesian coordinates (global) for a specific point 
	///              located on the NURBS surface S(u=fixed and v=fixed)
	///              Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
	///              Algorithm A4.3
	///
	/// ======================================================================================
	///  \param[in]  rSurfacePoint   	evaluated point
	///  \param[in]  _uoi    			local parameter in u-direction
	///  \param[in]  _voi    			local parameter in v-direction
	///
	/// ======================================================================================
	///  \author     Daniel Baumgärtner (12/2016)
	//
	//########################################################################################
	void evaluate_surface_point(Point<3>& rSurfacePoint, double u, double v)
	{
		rSurfacePoint[0] = 0;
		rSurfacePoint[1] = 0;
		rSurfacePoint[2] = 0;

		int span_u=find_Knot_Span(m_knot_vector_u,u,m_p,m_n_u);
		int span_v=find_Knot_Span(m_knot_vector_v,v,m_q,m_n_v);

		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				matrix<double> R;
				eval_Nurbsbasis(span_u, span_v, u, v, R);

				rSurfacePoint[0] += R(b,c) * m_control_points[control_point_index].getX();
				rSurfacePoint[1] += R(b,c) * m_control_points[control_point_index].getY();
				rSurfacePoint[2] += R(b,c) * m_control_points[control_point_index].getZ();
			}
		}
	}

	//  #####################################################################################
	// #######################################################################################
	//
	///  \details    returns the basis functions of NURBS basis function w.r.t. u,v
	///              span_u, span_v are the knot span indices. if unknown, insert 0!
	///
	/// ======================================================================================
	///  \param[in]  span_u     knotspan index in u-direction
	///  \param[in]  span_v     knotspan index in v-direction
	///  \param[in]  _u         local parameter in u-direction
	///  \param[in]  _v         local parameter in v-direction
	///  \param[out] R         basis func
	///
	/// ======================================================================================
	///  \author     Daniel Baumgärtner (12/2016)
	//
	//########################################################################################
	void eval_Nurbsbasis(int span_u, int span_v, double _u, double _v, matrix<double>& R)
	{
		if(span_u==0) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==0) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);

		matrix<double> N_matrix;
		matrix<double> M_matrix;
		Vector N;
		Vector M;

		R.resize( m_p+1,m_q+1 );
		noalias(R) = zero_matrix<double>( m_p+1 ,m_q+1 );

		// Evaluate basis functions with derivatives
		eval_nonzero_basis_function_with_derivatives(N_matrix, m_knot_vector_u, _u, span_u, m_p, 0);
		eval_nonzero_basis_function_with_derivatives(M_matrix, m_knot_vector_v, _v, span_v, m_q, 0);
		N = row( N_matrix, 0);
		M = row( M_matrix, 0);

		// eval_Nonzero_Basis_Fct(N,m_knot_vector_u,_u,span_u,m_p);
		// eval_Nonzero_Basis_Fct(M,m_knot_vector_v,_v,span_v,m_q);

		double sum = 0.0;

		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				// Evaluate basis function
				R(b,c) = N(b)*M(c)*m_control_points[control_point_index].getWeight();
				sum +=R(b,c);
			}
		}

		// divide by sum only required in terms of rational basis functions
		//if (fabs(sum-weight)> cepsilon) //Breitenberger 18.06.2014
		double inv_sum = 1/sum;
		// divide through by sum
		for(int c=0;c<=m_q;c++)
			for(int b=0;b<=m_p;b++)
				R(b,c) = inv_sum*R(b,c);
	}

	//  #####################################################################################
	// #######################################################################################
	//
	///  \details    returns the basis functions of NURBS basis function w.r.t. u,v
	///              span_u, span_v are the knot span indices. if unknown, insert 0!
	///				 Additionally returns mapping matrix ids of control points involved
	///
	/// ======================================================================================
	///  \param[in]  span_u     			knotspan index in u-direction
	///  \param[in]  span_v     			knotspan index in v-direction
	///  \param[in]  _u         			local parameter in u-direction
	///  \param[in]  _v         			local parameter in v-direction
	///  \param[out] R         				basis func
	///  \param[out] cp_mapping_matrix_ids  lists mapping matrix ids of control point corresponding to basis func
	///
	/// ======================================================================================
	///  \author     Daniel Baumgärtner (12/2016)
	//
	//########################################################################################
	void eval_Nurbsbasis(int span_u, int span_v, double _u, double _v, matrix<double>& R, matrix<unsigned int>& mapping_matrix_ids)
	{
		if(span_u==0) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==0) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);

		matrix<double> N_matrix;
		matrix<double> M_matrix;
		Vector N;
		Vector M;

		R.resize( m_p+1,m_q+1 );
		mapping_matrix_ids.resize( m_p+1,m_q+1 );
		noalias(R) = zero_matrix<double>( m_p+1 ,m_q+1 );
		noalias(mapping_matrix_ids) = zero_matrix<unsigned int>( m_p+1 ,m_q+1 );

		// Evaluate basis functions with derivatives
		eval_nonzero_basis_function_with_derivatives(N_matrix, m_knot_vector_u, _u, span_u, m_p, 0);
		eval_nonzero_basis_function_with_derivatives(M_matrix, m_knot_vector_v, _v, span_v, m_q, 0);
		N = row( N_matrix, 0);
		M = row( M_matrix, 0);

		// eval_Nonzero_Basis_Fct(N,m_knot_vector_u,_u,span_u,m_p);
		// eval_Nonzero_Basis_Fct(M,m_knot_vector_v,_v,span_v,m_q);

		double sum = 0.0;

		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				// Store control point Id in corresponding matrix
				mapping_matrix_ids(b,c) = m_control_points[control_point_index].getMappingMatrixId();

				// Evaluate basis function
				R(b,c) = N(b)*M(c)*m_control_points[control_point_index].getWeight();
				sum +=R(b,c);
			}
		}

		// divide by sum only required in terms of rational basis functions
		//if (fabs(sum-weight)> cepsilon) //Breitenberger 18.06.2014
		double inv_sum = 1/sum;
		// divide through by sum
		for(int c=0;c<=m_q;c++)
			for(int b=0;b<=m_p;b++)
				R(b,c) = inv_sum*R(b,c);
	}	

	// // #######################################################################################
	// ///
	// ///  \details    returns the first and second derivative of NURBS basis function w.r.t. u,v
	// ///              span_u,span_v are the knot span indices. if unknown, insert 0!
	// ///
	// /// ======================================================================================
	// ///  \param[in]  span_u     knotspan index in u-direction
	// ///  \param[in]  _u     	local parameter in u-direction
	// ///  \param[in]  span_v 	knotspan index in v-direction
	// ///  \param[in]  _v         local parameter in v-direction
	// ///  \param[out] _dR        1st derivatives
	// ///  \param[out] _ddR       2nd derivatives
	// ///
	// /// ======================================================================================
	// ///  \author     from M.Breitenberger in Carat (12/2009)
	// //
	// //########################################################################################
	
	 void deriv_Nurbsbasisfunc(int span_u, int span_v, double _u, double _v, matrix<double>& _dR, matrix<double>& _ddR)
	 {
	 	if(span_u==0) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
	 	if(span_v==0) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);
	
	 	int ne = (m_p+1)*(m_q+1); // Control Points per element
	 	matrix<double> N;              // Basisfunc at _u
	 	matrix<double> M;              // Basisfunc at _v
	 	eval_nonzero_basis_function_with_derivatives(N, m_knot_vector_u, _u, span_u, m_p, 2);
	 	eval_nonzero_basis_function_with_derivatives(M, m_knot_vector_v, _v, span_v, m_q, 2);
		
	 	vector<double> r(ne);
	 	r.clear();
	 	_dR.resize(ne,2);
	 	_ddR.resize(ne,3);
	 	double sum = 0.0;
	 	Vector dsum = ZeroVector(2);
	 	Vector ddsum = ZeroVector(3);
	 	int k=0;
	 	double weight;
		
	 	for(int c=0;c<=m_q;c++)
	 	{
	 		for(int b=0;b<=m_p;b++)
	 		{
	 			// the control point vector is filled up by first going over u, then over v
	 			int ui = span_u-m_p+b;
	 			int vi = span_v-m_q+c;
	 			int control_point_index =vi*m_n_u + ui;

	 			// Evaluate basis function
	 			weight = m_control_points[control_point_index].getWeight();

	 			r[k] = N(0,b)*M(0,c)*weight;
	 			sum +=r[k];
	 			//First derivatives
	 			_dR(k,0) = N(1,b)*M(0,c)*weight;
	 			dsum[0] +=_dR(k,0);
	 			_dR(k,1) = N(0,b)*M(1,c)*weight;
	 			dsum(1) +=_dR(k,1);
	 			//Second derivatives  1-du^2, 2-dv^2, 3-dudv
	 			_ddR(k,0) = N(2,b)*M(0,c)*weight;
	 			ddsum(0)  = ddsum(0) + _ddR(k,0);
	 			_ddR(k,1) = N(0,b)*M(2,c)*weight;
	 			ddsum(1)  = ddsum(1) + _ddR(k,1);
	 			_ddR(k,2) = N(1,b)*M(1,c)*weight;
	 			ddsum(2)  = ddsum(2) + _ddR(k,2);
	 			k++;
	 		}
	 	}
	 	double sum_2 = pow(sum,2);
	 	double sum_3 = pow(sum,3);
	 	// divide through by sum
	 	for(int k=0;k<ne;k++)
	 	{
		
	 		_ddR(k,0) = _ddR(k,0)/sum - 2.0*_dR(k,0)*dsum[0]/sum_2
	 				-r[k]*ddsum[0]/sum_2 + 2.0*r[k]*dsum[0]*dsum[0]/sum_3;
	 		_ddR(k,1) = _ddR(k,1)/sum - 2.0*_dR(k,1)*dsum[1]/sum_2
	 				-r[k]*ddsum[1]/sum_2 + 2.0*r[k]*dsum[1]*dsum[1]/sum_3;
	 		_ddR(k,2) = _ddR(k,2)/sum - _dR(k,0)*dsum[1]/sum_2 - _dR(k,1)*dsum[0]/sum_2
	 				-r[k]*ddsum[2]/sum_2 + 2.0*r[k]*dsum[0]*dsum[1]/sum_3;
	 		_dR(k,0) = _dR(k,0)/sum - r[k]*dsum[0]/sum_2;
	 		_dR(k,1) = _dR(k,1)/sum - r[k]*dsum[1]/sum_2;
	 	}
	 }

    //  #####################################################################################
    // #######################################################################################
	///
    ///   \details   evaluates all non-zero B-Spline basis functions (recursive formula)
    ///              based on a specified u and its corresponding knot span
    ///              Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
    ///              Algorithm A2.2
    ///
    /// ======================================================================================
    ///   \param[in]   _NBasisFct vector of non-zero B-Spline basis functions
    ///                           N = [Ni-p,p(u),...,Ni,p(u)]
    ///   \param[in]  _knotVec    knot vector in _par direction
    ///   \param[in]  _par        local parameter (curve -> 1D, isoline -> 2D)
    ///   \param[in]  _span       corresponding knot span
    ///   \param[in]  _pDeg       polynominal degree in _par direction
    ///
    /// ======================================================================================
	///  \author     Daniel Baumgärtner (12/2016)
    //
    //########################################################################################
	void eval_Nonzero_Basis_Fct(Vector& _NBasisFct,DoubleVector _knotVec,double _par,int _span,int _pDeg)
	{
		Vector left;
		Vector right;
		double saved;
		double temp;

		left.resize(_pDeg+1);
		right.resize(_pDeg+1);
		_NBasisFct.resize(_pDeg+1);
		_NBasisFct(0) = 1.00;

		for(int i=1;i<=_pDeg;i++)
		{
			left(i) = _par-_knotVec[_span+1-i];
			right(i) = _knotVec[_span+i]-_par;

			saved = 0.00;
			for(int j=0;j<i;j++)
			{
				temp = _NBasisFct(j)/(right[j+1]+left(i-j));
				_NBasisFct(j) = saved+right(j+1)*temp;
				saved = left(i-j)*temp;
			}
			_NBasisFct(i) = saved;
		}
	}

	//  #####################################################################################
	// #######################################################################################
	///
	///   \details   evaluates all non-zero B-Spline basis functions (recursive formula)
	///              and its n-th derivative based on a specified u and its corresponding
	///              knot span
	///              Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
	///              Algorithm A2.3
	///
	/// ======================================================================================
	///   \param[in] _dNBasisFct   matrix including the kth derivatives of the non-zero
	///                            B-spline basis fuctions (reals)
	///                            dNdu = |Ni-p,p(u)          ,...,  Ni,p(u)         |
	///                                   |dNdu_i-p,p(u)      ,...,  dNdu_i,p(u)     |
	///                                   |  :                  :       :            |
	///                                   |(dNdu_i-p,p(u))^k  ,...,  (dNdu_i,p(u))^k |
	///   \param[in]  _knotVec    knot vector in _par direction
	///   \param[in]  _par        local parameter (curve -> 1D, isoline -> 2D) (knot of interest)
	///   \param[in]  _span       corresponding knot span
	///   \param[in]  _pDeg       polynominal degree in _par direction
	///
	/// ======================================================================================
	///  \author     Daniel Baumgärtner (12/2016)
	//
	//########################################################################################
	void eval_nonzero_basis_function_with_derivatives(matrix<double>& _dNBasisFct,  DoubleVector _knotVec,double _par, int _span,int _pDeg,int _kth)
	{
	matrix<double> ndu;
	matrix<double> a;
	Vector left;
	Vector right;
	int s1,s2,j1,j2,jk,pk,jj,ll;
	double d,saved,temp;
	
	_dNBasisFct.resize(_kth+1,_pDeg+1);
	ndu.resize(_pDeg+1,_pDeg+1);
	a.resize(2,_pDeg+1);
	left.resize(_pDeg+1);
	right.resize(_pDeg+1);
	ndu(0,0) = 1.00;
	
	for(int i=1;i<=_pDeg;i++){
	left(i) = _par-_knotVec[_span+1-i];
	right(i) = _knotVec[_span+i]-_par;
	saved = 0.00;
	for(int j=0;j<i;j++){
		ndu(i,j) = right(j+1)+left(i-j);
		temp = ndu(j,i-1)/ndu(i,j);
		ndu(j,i) = saved+right(j+1)*temp;
		saved = left(i-j)*temp;
	}
	ndu(i,i) = saved;
	}
	
	for(int i=0;i<=_pDeg;i++){
	_dNBasisFct(0,i) = ndu(i,_pDeg);
	}
	for(int j=0;j<=_pDeg;j++){
	s1 = 0;
	s2 = 1;
	a(0,0) = 1.00;
	for(int k=1;k<=_kth;k++){
		d = 0.00;
		jk = j-k;
		pk = _pDeg-k;
		if(j >= k){
		a(s2,0) = a(s1,0)/ndu(pk+1,jk);
		d = a(s2,0)*ndu(jk,pk);
		}
		if(jk >= -1){
		j1 = 1;
		}
		else{
		j1 = -jk;
		}
		if(j-1 <= pk){
		j2 = k-1;
		}
		else{
		j2 = _pDeg-j;
		}
		for(int l=j1;l<=j2;l++){
		a(s2,l) = (a(s1,l)-a(s1,l-1))/ndu(pk+1,jk+l);
		d += a(s2,l)*ndu(jk+l,pk);
		}
		if(j<=pk){
		a(s2,k) = -a(s1,k-1)/ndu(pk+1,j);
		d += a(s2,k)*ndu(j,pk);
		}
		_dNBasisFct(k,j) = d;
		ll = s1;
		s1 = s2;
		s2 = ll;
	}
	}
	jj = _pDeg;
	for(int k=1;k<=_kth;k++){
	for(int l=0;l<=_pDeg;l++){
		_dNBasisFct(k,l) *= jj;
	}
	jj *= (_pDeg-k);
	}
	}
	

	//  #####################################################################################
	// #######################################################################################
	//#
	//#                  ++++++++++++++++++++++++++++++++++++
	//#                  +++  NurbsBasis::find_Knot_Span  +++
	//#                  ++++++++++++++++++++++++++++++++++++
	//#
	///   \details   returns the corresponding knot span based on a specified u
	///              Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
	///              Algorithm A2.1
	///
	/// ======================================================================================
	///   \return     index of the first knot in the knot span in _par direction
	///   \param[in]  _knotVec  knot vector in _par direction
	///   \param[in]  _par      local parameter (curve -> 1D, isoline -> 2D)
	///   \param[in]  _pDeg     polynominal degree in _par direction
	///   \param[in]  _nCtrl    number of control points in _par direction
	///
	/// ======================================================================================
	///  \author     Daniel Baumgärtner (12/2016)
	//
	//########################################################################################
	int find_Knot_Span(DoubleVector _knotVec, double _par, int _pDeg, int _nCtrl)
	{

		if(_par<=_knotVec[0])
			_par = _knotVec[0]+1.0e-13;

		if(_par>=_knotVec[_knotVec.size()-1])
			_par = _knotVec[_knotVec.size()-1]-1.0e-13;

		int low = _pDeg;
		int high = _nCtrl+1;
		int span=(low+high)/2;

		while (_par < _knotVec[span] || _par >= _knotVec[span+1])
		{
			if (_par < _knotVec[span])
			{
				high = span;
			}
			else
			{
				low = span;
			}
			span = (low+high)/2;
		}

		return span;
	}

	// #######################################################################################
	///
	///  \details    returns global ids of control points which are affected by given u,v.
	///				 Ids are returned as matrix ordered in the same way as the corresponding
	///				 basis function. I.e. shape function value of control point specified by 
	/// 			 global_cp_ids(i,j) corresponds to R(i,j).
	///
	/// ======================================================================================
	///  \param[in]  span_u     	knotspan index in u-direction
	///  \param[in]  span_v     	knotspan index in v-direction	
	///  \param[in]  _u    			local parameter in u-direction
	///  \param[in]  _v    			local parameter in v-direction
	///  \param[out] global_cp_ids 	matrix storing global control point ids
	///
	/// ======================================================================================
	///  \author     Daniel Baumgärtner (12/2016)
	//
	//########################################################################################
	void get_global_control_point_ids(int span_u, int span_v, double _u, double _v, matrix<int>& global_cp_ids)
	{

		if(span_u==0) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==0) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);

		global_cp_ids.resize( m_p+1,m_q+1 );
		noalias(global_cp_ids) = zero_matrix<int>( m_p+1 ,m_q+1 );

		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				// Store control point Id in corresponding matrix
				global_cp_ids(b,c) = m_control_points[control_point_index].getGlobalId();
			}
		}
	}	

	// #######################################################################################
	///
	///  \details    returns local ids of control points which are affected by given u,v.
	///				 Ids are returned as matrix ordered in the same way as the corresponding
	///				 basis function. I.e. shape function value of control point specified by 
	/// 			 local_cp_ids(i,j) corresponds to R(i,j).
	///
	/// ======================================================================================
	///  \param[in]  span_u    	 	knotspan index in u-direction
	///  \param[in]  span_v     	knotspan index in v-direction	
	///  \param[in]  _u    			local parameter in u-direction
	///  \param[in]  _v    			local parameter in v-direction
	///  \param[out] local_cp_ids 	matrix storing global control point ids
	///
	/// ======================================================================================
	///  \author     Daniel Baumgärtner (12/2016)
	//
	//########################################################################################
	void get_local_control_point_ids(int span_u, int span_v, double _u, double _v, matrix<int>& local_cp_ids)
	{

		if(span_u==0) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==0) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);

		local_cp_ids.resize( m_p+1,m_q+1 );
		noalias(local_cp_ids) = zero_matrix<double>( m_p+1 ,m_q+1 );

		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				// Store control point Id in corresponding matrix
				local_cp_ids(b,c) = m_control_points[control_point_index].getLocalId();
			}
		}
	}		
	//-----------------------------------------------------------------------------------------
	// #######################################################################################
	///
	///  \details    Returns a DoubleVector (a.k.a. std::vector<double> given by c*M*b
	///
	/// ======================================================================================
	///  \param[in]  QminP   	 	Vector distance between Q and P
	///  \param[in]  H     			Hessian reference matrix, defined as matrix<double>
	///  \param[in]  Gradient  		Gradient reference vector, defined as std::vector<double> (a.k.a. DoubleVector)
    ///  \param[in]  u,v  			u and v value of the parameter space
	///
	/// ======================================================================================
	///  \author     Giovanni Filomeno (1/2017)
	//
	//########################################################################################
	void evaluate_Hessian_and_Gradient(DoubleVector QminP, matrix<double>& H, DoubleVector& Gradient , double u, double v)
	{
		matrix<double> dR;
		matrix<double> ddR;
		deriv_Nurbsbasisfunc(0,0, u, v, dR,ddR);
		DoubleVector dQdu(3,0.0);
		DoubleVector dQdv(3,0.0);
		DoubleVector dQdudu(3,0.0);
		DoubleVector dQdvdv(3,0.0);
		DoubleVector dQdudv(3,0.0);
		matrix<int> local_cp_ids;

		get_local_control_point_ids(0, 0, u, v, local_cp_ids);
		int k=0;
		for( int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				dQdu[0] += dR(k,0) * m_control_points[local_cp_ids(b,c)].getX();
				dQdu[1] += dR(k,0) * m_control_points[local_cp_ids(b,c)].getY();
				dQdu[2] += dR(k,0) * m_control_points[local_cp_ids(b,c)].getZ();

				dQdv[0] += dR(k,1) * m_control_points[local_cp_ids(b,c)].getX();
				dQdv[1] += dR(k,1) * m_control_points[local_cp_ids(b,c)].getY();
				dQdv[2] += dR(k,1) * m_control_points[local_cp_ids(b,c)].getZ();

				dQdudu[0] += ddR(k,0) * m_control_points[local_cp_ids(b,c)].getX();
				dQdudu[1] += ddR(k,0) * m_control_points[local_cp_ids(b,c)].getY();
				dQdudu[2] += ddR(k,0) * m_control_points[local_cp_ids(b,c)].getZ();

				dQdvdv[0] += ddR(k,1) * m_control_points[local_cp_ids(b,c)].getX();
				dQdvdv[1] += ddR(k,1) * m_control_points[local_cp_ids(b,c)].getY();
				dQdvdv[2] += ddR(k,1) * m_control_points[local_cp_ids(b,c)].getZ();

				dQdudv[0] += ddR(k,2) * m_control_points[local_cp_ids(b,c)].getX();
				dQdudv[1] += ddR(k,2) * m_control_points[local_cp_ids(b,c)].getY();
				dQdudv[2] += ddR(k,2) * m_control_points[local_cp_ids(b,c)].getZ();

				k++;
			}
		}

		H(0,0) = std::inner_product( dQdudu.begin(), dQdudu.end(),  QminP.begin(), 0);
		H(0,1) = std::inner_product( dQdudv.begin(), dQdudv.end(),  QminP.begin(), 0);
		H(1,0) = std::inner_product( dQdudv.begin(), dQdudv.end(),  QminP.begin(), 0);
		H(1,1) = std::inner_product( dQdvdv.begin(), dQdvdv.end(),  QminP.begin(), 0);


		Gradient[0] = std::inner_product( dQdu.begin(), dQdu.end(),  QminP.begin(), 0);
		Gradient[1] = std::inner_product( dQdv.begin(), dQdv.end(),  QminP.begin(), 0);
//		for (int c=0;c<=m_q;c++)
//		{
//			for (int b=0;b<=m_p;b++)
//			{
//				// the control point vector is filled up by first going over u, then over v
//				int ui = find_Knot_Span(m_knot_vector_u,u,m_p,m_n_u)-m_p+b;
//				int vi = find_Knot_Span(m_knot_vector_v,v,m_q,m_n_v)-m_q+c;
//				int control_point_index =vi*m_n_u + ui;
//
//				// Store control point Id in corresponding matrix
//				local_cp_ids(b,c) = m_control_points[control_point_index].getLocalId();
//			}
//		}

	}
	// --------------------------------------------------------------------------
	DoubleVector& get_knot_vector_u()
	{
		return m_knot_vector_u;
	}

	// --------------------------------------------------------------------------
	DoubleVector& get_knot_vector_v()
	{
		return m_knot_vector_v;
	}

	ControlPointVector& get_control_points()
	{
		return m_control_points;
	}

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "Patch";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "Patch";
	}

	// ==============================================================================
	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}

	int get_m_p()
	{
		return m_p;
	}

	int get_m_q()
	{
		return m_q;
	}


private:
	// ==============================================================================
	// Initialized by class constructor
	// ==============================================================================
	DoubleVector m_knot_vector_u;
	DoubleVector m_knot_vector_v;
	int m_p;
	int m_q;
	ControlPointVector m_control_points;
	unsigned int m_n_u; // number of control points in u-direction
	unsigned int m_n_v; // number of control points in v-direction
	double m_epsilon = 1e-13;

	// double TOL = 10e-4;

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      Patch& operator=[Patch const& rOther];

	/// Copy constructor.
	//      Patch[Patch const& rOther];

}; // Class Patch

} // namespace Kratos.

#endif // PATCH_H
