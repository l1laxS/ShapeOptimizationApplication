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

#ifndef CONTROL_POINT_H
#define CONTROL_POINT_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <Python.h>

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
// ==============================================================================

namespace Kratos
{
class ControlPoint
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
    typedef std::vector<controlPoint> controlPointVcr;
    typedef std::vector<int> IntVector;

    /// Pointer definition of ControlPoint
    KRATOS_CLASS_POINTER_DEFINITION(ControlPoint);

    /// Default constructor.
    ControlPoint( int id,
    			  double x,
    			  double y,
				  double z,
				  double w):
				id(id),
    			x(x),
				y(y),
				z(z),
				w(w)
    {
    }

    /// Destructor.
    virtual ~ControlPoint()
    {
    }

    // ==============================================================================
    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ControlPoint";
    }

    // ==============================================================================
    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ControlPoint";
    }

    // ==============================================================================
    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    int id;
    double x;
    double y;
    double z;
    double w;

    // ==============================================================================
    // General working arrays
    // ==============================================================================
    /// Assignment operator.
//      ControlPoint& operator=(ControlPoint const& rOther);

    /// Copy constructor.
//      ControlPoint(ControlPoint const& rOther);


}; // Class ControlPoint

}  // namespace Kratos.

#endif // CONTROL_POINT_H
