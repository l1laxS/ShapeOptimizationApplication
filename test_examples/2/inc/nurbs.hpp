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
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//   Last modified by:    $Co-Author: daniel.baumgaertner@tum.de $
//   Date:                $Date:                      March 2016 $
//   Revision:            $Revision:                         0.0 $
//
// ==============================================================================

#ifndef NURBS_HPP_
#define NURBS_HPP_

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <iomanip>      // for std::setprecision
#include "/home/giovanni/workspace/kratos/applications/ShapeOptimizationApplication/test_examples/CAD_reconstruction/supportTypes/supportTypes.hpp"
// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
//#include <boost/python.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
//#include "../../kratos/includes/define.h"
//#include "../../kratos/processes/process.h"
//#include "../../kratos/includes/node.h"
//#include "../../kratos/includes/element.h"
//#include "../../kratos/includes/model_part.h"
//#include "../../kratos/includes/kratos_flags.h"
//#include "../../kratos/utilities/normal_calculation_utils.h"
//#include "shape_optimization_application.h"

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
#define TOL 0.00001
class Nurbs
{
public:
    ///@name Type Definitions
    ///@{

    // ==========================================================================
    // Type definitions for better reading later
    // ==========================================================================
    typedef std::vector<double> DoubleVector;

    /// Pointer definition of Nurbs
//    KRATOS_CLASS_POINTER_DEFINITION(Nurbs);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Nurbs(int grade)
    {
    	degree=grade;
    }

    /// Destructor.
    virtual ~Nurbs()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================

    // Set X coordinates of control points
    void setCPX( const DoubleVector& cPX)
    {
        myCPX = cPX;
    }

    // --------------------------------------------------------------------------


    // Set Y coordinates of control points
    void setCPY( const DoubleVector& cPY )
    {
       myCPY = cPY;
    }

    // Set Z coordinates of control points
    void setCPZ( const DoubleVector& cPZ )
    {
       myCPZ = cPZ;
    }

    // ==============================================================================

    // Set weight of control point
    void setCPW(const DoubleVector& cPW)
    {
    	myCPW = cPW;
    }

    // Set knot vector
    void setKnotVector(const DoubleVector& knot)
    {
    	U.clear();
    	DoubleVector::const_iterator iter;

    	// entries is a smart iterator
    	for(iter = knot.begin(); iter!=knot.end(); )
    	{
    		U.push_back(*iter); //push back in my knot
    	}
    }


    // check if a knot vector is valid,
    //valid means that p+1 are equal (at begining and end)
    // but no each other
    bool checkKnotVector()
    {
    	// m = n+p+1
    	// m number of knots
    	// p degree
    	// n number of control points
    	// n+p+1 = m
    	// My knot is empty
    	if( U.empty() == 0)
    	{
    		return false;
    	}

    	int m = U.size();
    	int n = myCPW.size();

    	if(n+degree+1 == m)
    	{
    		return true;
    	}else
    	{
    		return false;
    	}
    }

    double N(const double& t, const int& i, const int& p)
    {
        // max index of knot vector (from 0 to m)
        int m = U.size() - 1;

        // check if i is in interval 0 <= i <= n. (with n = m-p-1)
        // and if t is in the interval of the knot vector
        if((i < 0) || (i > (m-p-1)) || t<U.front() || t>U.back())
        {
        	return 0.0;
        }

        // treat case of p == 0
        if(p == 0) {
            // check if t is in interval from t_i <= t < t_(i+1)
            if((U[i] <= t) && (t < U[i+1]))
            {
            	// If yes, return 1
            	return 1.0;
            }
            else
            {
            	// If I am out of interval, return zero
            	return 0.0;
            }

        }
        else
        {
            double e1,e2, up,lo;

            // compute division term 1
            up = t - U[i];
            lo = U[i+p] - U[i];

            if(fabs(lo)<TOL)
            {
            	// According to the definition 0/0 = 0
            	e1 = 0.0;
            }
            else
            {
            	e1 = up/lo;
            }

            // compute division term 2
            up = U[i+p+1] - t;
            lo = U[i+p+1] - U[i+1];

            if(fabs(lo)<TOL)
            {
            	e2 = 0.0;
            }
            else
            {
            	e2 = up/lo;
            }

            // compute basis function recursively
            double sn = 0.0;
            if(fabs(e1)>TOL)
            {
            	sn += e1 * N(t,i  ,p-1); // recursively function
            }
            if(fabs(e2)>TOL)
            {
            	sn += e2 * N(t,i+1,p-1);
            }

            return sn;
        }
    }

    double dN(const double& t, const int& i, const int& p)
    {
    	double up = p;
    	double lo = U[i+p] - U[i];
    	double e1;
    	double e2;
    	double dN = 0;

    	if(fabs(lo) < TOL)
    	{
    		e1 = 0;
    	}else
    	{
    		e1 = up / lo;
    	}

    	lo = U[i+p+1] - U[i+1];

    	if(fabs(lo) < TOL)
    	{
    		e2 = 0;
    	}else
    	{
    		e2 = up/lo;
    	}

    	if( fabs(e1) > TOL)
    	{
    		dN = dN + e1*N(t, i, p-1);
    	}

    	if ( fabs(e2) > TOL)
    	{
    		dN = dN + e2*N(t, i+1, p-1);
    	}

    	return dN;
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Nurbs";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Nurbs";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }



    ///@}

protected:


private:

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================

    // Polynomial degree
    int degree;

    // knot vector
    DoubleVector U;

    // Characteristics of control points
    DoubleVector myCPX;
    DoubleVector myCPY;
    DoubleVector myCPZ;
    DoubleVector myCPW;

    /// Assignment operator.
//      Nurbs& operator=(Nurbs const& rOther);

    /// Copy constructor.
//      Nurbs(Nurbs const& rOther);

}; // Class Nurbs

}  // namespace Kratos.

#endif // NURBS_HPP_
