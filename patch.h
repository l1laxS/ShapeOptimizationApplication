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
    }

    /// Destructor.
    virtual ~Patch()
    {
    }

    // ==============================================================================
    // returns the value of the NURBS function given the parameters
//    void eval_Nurbsbasis(int i, int j, double u, double v, DoubleMatrix& R)
//    {
//    	//
//    	if(i==0)
//    	{
//    		i=find_Knot_Span(knot_vector_u,u,p,N_Ctrl); // da dove sbucano fuori?
//
//    		if(j==0)
//    		{
//    			j=find_Knot_Span(knot_vector_v,v,q,M_Ctrl); //da dove vengono pt2?
//    		}
//
//    	  //matrix<cfloat> N;
//    	  //matrix<cfloat> M;
//    	  DoubleVector N;
//    	  DoubleVector M;
//    	  //eval_Derivative_NonzeroBasis_Fct[N,knot_vector_u,_u,_i,p,2];
//    	  //eval_Derivative_NonzeroBasis_Fct[M,knot_vector_v,_v,_j,q,2];
//    	  R.resize[p+1][q+1];
//    	  eval_Nonzero_Basis_Fct(N,knot_vector_u,u,i,p);
//    	  eval_Nonzero_Basis_Fct(M,knot_vector_v,v,j,q);
//
//    	  double sum = 0.0;
//    	  double weight;
//
//    	  for (int c=0;c <= q; c++)
//    	  {
//    	    for (int b = 0; b <= p; b++)
//    	    {
//    	      weight = Ctrl_Pt_Net[i-p+b][j-q+c]->get_Weight();
//    	      R[b][c] = N[b]*M[c]*weight;
//    	      sum +=R[b][c];
//    	    }
//    	  }
//
//    	  // divide by sum only required in terms of rational basis functions
//    	  //if [fabs[sum-weight]> cepsilon] //Breitenberger 18.06.2014
//    	    double inv_sum = 1/sum;
//    	    // divide through by sum
//    	    for(int c=0; c<=q; c++)
//    	    {
//    	    	for(int b = 0 ; b <= p; b++)
//    	    	{
//    	    		R[b][c] = inv_sum*R[b][c];
//    	    	}
//    	    }
//
//    	//hint for computing X from given parameter
//    	double x=0
//    	double y=0
//    	double z=0
//
//    	for (int c = 0; c <= q; c++ )
//    	  {
//    	    for (int b=0; b <= p; b++)
//    	    {
//    	    	x = x + R[b][c]* Ctrl_Pt_Net[i-p+b][j-q+c]->x;
//    	    	y = y + R[b][c] * Ctrl_Pt_Net[i-p+b][j-q+c]->y;
//    	    	z += R[b][c] * Ctrl_Pt_Net[i-p+b][j-q+c]->z;
//    	    }
//    	  }
//
//    	}
//    };

    // returns the value of the NURBS derivative given the parameters
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
//    	ndu[0,0] = 1.00;
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
//    	  for(int j=0;j<=pDeg;j++)
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
//    				  a[s2][l] = [a[s1][l]-a[s1][l-1]]/ndu[pk+1][jk+l];
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
//    	  jj *= [pDeg-k];
//    	  }
//    	  dNBasisFct = dnFct;
//    }

    // returns the coordinates of the point given the parameters
    void S(double& x, double& y, double& z, double u, double v)
    {

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

    int find_Knot_Span (DoubleVector knotVec, double par,int pDeg, int nCtrl) const
    {

        if(par<=knotVec[0])
        {
        	par = knotVec[0]+1.0e-14; // what is cepsilon?? maybe 1.0e-14
        }

        if(par >= knotVec[knotVec.size()-1])
        {
        	par = knotVec[knotVec.size()-1]-1.0e-14;
        }

      float check;
      check = fabs(par - knotVec[nCtrl+1] );

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
