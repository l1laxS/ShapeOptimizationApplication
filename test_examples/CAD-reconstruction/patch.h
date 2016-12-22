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
#include <Python.h>
#include <cmath>
#include <math.h>

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
#include "shape_optimization_application.h"
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
    typedef std::vector<double> DoubleVector;
    typedef boost::python::extract<double> takeDouble;
    typedef boost::python::extract<int> takeInt;
    typedef boost::python::extract<bool> takeBool;
    typedef std::vector< std::vector<double> > DoubleMatrix;
    typedef std::vector<ControlPoint> controlPointVcr;
    typedef std::vector<int> IntVector;

    /// Pointer definition of Patch
//    KRATOS_CLASS_POINTER_DEFINITION[Patch];

    /// Default constructor.
    Patch( DoubleVector& knot_vector_u,
    	   DoubleVector& knot_vector_v,
		  int p, int q,
		  controlPointVcr& control_points):
			  knot_vector_u(knot_vector_u),
			  knot_vector_v(knot_vector_v),
			  p(p), q(q),
			  control_points(control_points)
    {
    	m_u = knot_vector_u.size();
    	m_v = knot_vector_v.size();

    	n_u = m_u - p -1;
    	n_v = m_v - q -1;

//    	std::cout << m_u << m_v << n_u << n_v;

    	size_t nu = n_u;
    	size_t nv = n_v;

    	if( control_points.size() != nu*nv )
    	{
    		std::cout << "Invalid Patch" << std::endl;
    	}

    }

    /// Destructor.
    virtual ~Patch()
    {
    }

    // ==============================================================================
//    // returns the value of the NURBS function given the parameters
//    void eval_Nurbsbasis(int i, int j, double u, double v, DoubleMatrix& R)
//    {
//    	//
//    	if(i==0)
//    	{
//    		i=find_Knot_Span(knot_vector_u,u,p,n_v); // da dove sbucano fuori?
//    	}
//    	if(j==0)
//    	{
//    		j=find_Knot_Span(knot_vector_v,v,q,n_u); //da dove vengono pt2?
//    	}
//
//    	  //matrix<cfloat> N;
//    	  //matrix<cfloat> M;
//    	DoubleVector N;
//    	DoubleVector M;
//    	  //eval_Derivative_NonzeroBasis_Fct[N,knot_vector_u,_u,_i,p,2];
//    	  //eval_Derivative_NonzeroBasis_Fct[M,knot_vector_v,_v,_j,q,2];
//    	R.resize[p+1][q+1];
//    	eval_Nonzero_Basis_Fct(N,knot_vector_u,u,i,p);
//    	eval_Nonzero_Basis_Fct(M,knot_vector_v,v,j,q);
//
//    	double sum = 0.0;
//    	double weight;
//
//    	for (int c=0;c <= q; c++)
//    	{
//    	  for (int b = 0; b <= p; b++)
//    	  {
//    	    weight = Ctrl_Pt_Net[i-p+b][j-q+c]->get_Weight();
//    	    R[b][c] = N[b]*M[c]*weight;
//    	    sum +=R[b][c];
//    	  }
//    	}
//
//    	  // divide by sum only required in terms of rational basis functions
//    	  //if [fabs[sum-weight]> cepsilon] //Breitenberger 18.06.2014
//    	double inv_sum = 1/sum;
//    	   // divide through by sum
//    	for(int c=0; c<=q; c++)
//    	{
//    	 	for(int b = 0 ; b <= p; b++)
//    	  	{
//    	   		R[b][c] = inv_sum*R[b][c];
//    	   	}
//    	}
//
//    	//hint for computing X from given parameter
//    	double x=0;
//    	double y=0;
//    	double z=0;
//
//		for(int l = 0; l <= q; l++)
//		{
//			for(int a = 0; a <=q; a++ )
//			{
//				x = x + R[a][l] * Ctrl_Pt_Net[ i - p + a ][ j - q + l ]->getX();
//				y = y + R[a][l] * Ctrl_Pt_Net[ i - p + a ][ j - q + l ]->getY();
//				z = z + R[a][l] * Ctrl_Pt_Net[ i - p + a ][ j - q + l ]->getZ();
//			}
//		}
//
//    };

//    // returns the value of the NURBS derivative given the parameters
//    void eval_Derivative_NonzeroBasis_Fct(DoubleMatrix dNBasisFct,
//    		DoubleVector knotVec, double par, int span, int pDeg, int kth)
//    {
//    	DoubleMatrix dnFct;
//    	DoubleMatrix ndu;
//    	DoubleMatrix a;
//    	DoubleVector left;
//    	DoubleVector right;
//    	int s1,s2,j1,j2,jk,pk,jj,ll;
//    	double d,saved,temp;
//
////    	  dnFct.resize[_kth+1,_pDeg+1];
////    	  ndu.resize[_pDeg+1,_pDeg+1];
////    	  a.resize[2,_pDeg+1];
////    	  left.resize[_pDeg+1];
////    	  right.resize[_pDeg+1];
//    	ndu[0][0] = 1.00;
//
//    	for(int i=1; i <= pDeg; i++)
//    	{
//    	  left[i]= par-knotVec[span+1-i];
//    	  right[i]= knotVec[span+i]-par;
//    	  saved = 0.00;
//    	  for(int j = 0; j < i; j++)
//    	  {
//    	    ndu[i][j] = right[j+1]+left[i-j];
//    	    temp = ndu[j][i-1]/ndu[i][j];
//    	    ndu[j][i] = saved+right[j+1]*temp;
//    	    saved = left[i-j]*temp;
//    	  }
//    	  ndu[i][i] = saved;
//    	  }
//
//    	  for(int i=0;i <= pDeg;i++)
//    	  {
//    		  dnFct[0][i] = ndu[i][pDeg];
//    	  }
//    	  for(int j=0;j <= pDeg;j++)
//    	  {
//    		  s1 = 0;
//    		  s2 = 1;
//    		  a[0][0] = 1.00;
//    		  for(int k=1 ;k <= kth;k++)
//    		  {
//    			  d = 0.00;
//    			  jk = j-k;
//    			  pk = pDeg-k;
//    			  if(j >= k)
//    			  {
//    				  a[s2][0] = a[s1][0]/ndu[pk+1][jk];
//    				  d = a[s2][0]*ndu[jk][pk];
//    			  }
//    			  if(jk >= -1)
//    			  {
//    				  j1 = 1;
//    			  }
//    			  else
//    			  {
//    				  j1 = -jk;
//    			  }
//    			  if(j-1 <= pk)
//    			  {
//    				  j2 = k-1;
//    			  }
//    			  else
//    			  {
//    				  j2 = pDeg-j;
//    			  }
//    			  for(int l=j1;l<=j2;l++)
//    			  {
//    				  a[s2][l] = (a[s1][l]-a[s1][l-1])/ndu[pk+1][jk+l];
//    				  d += a[s2][l]*ndu[jk+l][pk];
//    			  }
//    			  if(j<=pk)
//    			  {
//    				  a[s2][k] = -a[s1][k-1]/ndu[pk+1][j];
//    				  d += a[s2][k]*ndu[j][pk];
//    			  }
//    			  dnFct[k][j] = d;
//    		ll = s1;
//    	    s1 = s2;
//    	    s2 = ll;
//    	  }
//    	  }
//    	  jj = pDeg;
//    	  for(int k=1;k<=kth;k++)
//    	  {
//    		  for(int l = 0;l <= pDeg;l++)
//    		  {
//    			  dnFct[k][l] *= jj;
//    		  }
//    	  jj *= (pDeg-k);
//    	  }
//    	  dNBasisFct = dnFct;
//    }

    // returns the coordinates of the point given the parameters
    void S(double& x, double& y, double& z, double u, double v)
    {
    	x = 0;
    	y = 0;
    	z = 0;

    	for(int i = 0; i < n_u; i++)
    	{
    		for(int j = 0; j < n_v; j++)
    		{
    			int index = getCPIndex( i, j);
    			double RR = R( u, v, i, j);
    			x = x + RR * control_points[ index ].getX();
    			y = y + RR * control_points[ index ].getY();
    			z = z + RR * control_points[ index ].getZ();
    		}
    	}
    }

    int getCPIndex(const int i,const int j)
    {
    	return j * n_u + i;
    }


    double R(const double& u, const double& v,const int& i, const int& j)
    {
//    	for(int k=5; k>= 0 ; k--)
//    	{
//    	for(int i=0; i<k; i++)
//    	{
//    		std::cout << "i = " << i << std::endl;
//    		std::cout << N( knot_vector_v, 0, i, 5-k) << std::endl;
//    		std::cout << N( knot_vector_v, 2.5, i, 5-k) << std::endl;
//    		std::cout << N( knot_vector_v, 5.0, i, 5-k) << std::endl << std::endl;
//    	}
//    	}

        // max index of knot vector (from 0 to m-1)

        double lo = 0.0;
        double n_tmp = 0.0;

//        std::cout << "I'm in R" << std::endl;

        DoubleVector Nii;
        DoubleVector Njj;

        for(int jj=0; jj < n_v; jj++)
        {
         	 auto Nj = N(knot_vector_v, v, jj, q);
         	 Njj.push_back( Nj );
        }

//        if( v == knot_vector_v.back())
//        {
//        	Njj[ n_v - 1] = 1.0;
//        }

        for(int ii=0; ii< n_u; ii++)
        {
        	double N_ii = N(knot_vector_u, u, ii, p);
        	Nii.push_back( N_ii );
        }

//        if( u == knot_vector_u.back())
//        {
//           	Nii[ n_u - 1] = 1.0;
//        }



        for(int ii=0; ii < n_u; ii++)
        {
        	for(int jj=0; jj < n_v; jj++)
        	{
        		n_tmp = Nii[ii] * Njj[jj];
        		int index  = getCPIndex( ii, jj);
        		lo = lo + n_tmp * control_points[ index ].getWeight();
//        		lo = lo + n_tmp;
        	}
        }

//        std::cout << "I'm in R" << std::endl;

        if(fabs(lo)<TOL) {
            return 0;
        }

        int index = getCPIndex( i, j);
        double up = N(knot_vector_u, u, i, p) * N(knot_vector_v, v, j, q) * control_points[ index ].getWeight();
//        double up = N(knot_vector_u, u, i, p) * N(knot_vector_v, v, j, q);
//        std:: cout << "I am in R " << up/lo << std::endl;
        return up/lo;



    }
    // Returns the R without weight in u/v direction


    // U = knot vector
    // t = parameter in the parameter space
    // p = degree
    // i = number of shape function

	double N(DoubleVector& U,const double& t,const int& i,const int& p_degree)
	{
		// max index of knot vector (from 0 to m)
		int m = U.size();

//		int lastZeroDegree = m - 1;

		if( p_degree == 0 && t == U.back() && i == m-p-2)
		{
			return 1.0;
		}

//		std::cout << std::endl << "BEGIN" << std::endl;
//
//		std::cout << "N" << i << "," << p_degree << std::endl;

		// check if i is in interval 0 <= i <= n. (with n = m-p-1)
		// and if t is in the interval of the knot vector
		if ((i < 0) || (i >= (m - p - 1)) || t < U.front() || t > U.back()) {

			return 0.0;
		}

		// treat case of p_degree == 0
		if (p_degree == 0) {
			// check if t is in interval from t_i <= t < t_(i+1)
			if ((U[i] <= t) && (t < U[i + 1])) {

				return 1.0;
			} else {

				return 0.0;
			}

		} else {
			double e1, e2, up, lo;

			// compute division term 1
			up = t - U[i];
			lo = U[i + p_degree] - U[i];

			if (fabs(lo) < TOL) {
				e1 = 0.0;
			} else {
				e1 = up / lo;
			}

//			std::cout<<"Up=" <<up << " LO=" << lo <<std::endl;
			// compute division term 2
			up = U[i + p_degree + 1] - t;
			lo = U[i + p_degree + 1] - U[i + 1];

			if (fabs(lo) < TOL) {
				e2 = 0.0;
			} else {
				e2 = up / lo;
			}
//			std::cout<<"Up=" <<up << " LO=" << lo <<std::endl;
//			std::cout <<"N" << i << "," << p_degree << "-->e1 " << e1 << " e2 " << e2 << std::endl;
			// compute basis function recursively
			double sn = 0.0;

			if (fabs(e1) > TOL) {
//				std::cout << "calling" << std::endl;
				sn += e1 * N(U, t, i, p_degree - 1);
			}

			if (fabs(e2) > TOL) {
//				std::cout << "calling" << std::endl;
				sn += e2 * N(U, t, i + 1, p_degree - 1);
			}


			return sn;
		}
	}
    // ==============================================================================
    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Patch";
    }

    // ==============================================================================
    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Patch";
    }

    // ==============================================================================
    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    int find_Knot_Span(DoubleVector knotVec, double par,int pDeg, int nCtrl) const
    {

        if(par<=knotVec[0])
        {
        	par = knotVec[0]+1.0e-14; // what is cepsilon?? maybe 1.0e-14
        }

        if(par >= knotVec[knotVec.size()-1])
        {
        	par = knotVec[knotVec.size()-1]-1.0e-14;
        }

//      float check;
//      check = fabs(par - knotVec[nCtrl+1] );

      int low = pDeg;
      int high = nCtrl+1;
      int span= ( low + high )/2;

      while (par < knotVec[span] || par >= knotVec[span+1])
      {
        if (par < knotVec[span])
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

private:
    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    DoubleVector knot_vector_u;
    DoubleVector knot_vector_v;
    int p;
    int q;
    controlPointVcr control_points;
    // Number of control points in u & v
    int n_u;
    int n_v;
    int m_u;
    int m_v;
    double TOL = 10e-4;

    // ==============================================================================
    // General working arrays
    // ==============================================================================
    /// Assignment operator.
//      Patch& operator=[Patch const& rOther];

    /// Copy constructor.
//      Patch[Patch const& rOther];


}; // Class Patch

}  // namespace Kratos.

#endif // PATCH_H
