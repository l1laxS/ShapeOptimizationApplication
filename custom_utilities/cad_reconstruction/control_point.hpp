// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Giovanni Filomeno
                     giovanni.filomeno@tum.de
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

#ifndef CONTROL_POINT_H_
#define CONTROL_POINT_H_

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

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

    /// Pointer definition of ControlPoint
    KRATOS_CLASS_POINTER_DEFINITION(ControlPoint);

    /// Default constructor.
    ControlPoint( double X, double Y, double Z, double W, unsigned int local_id, unsigned int global_id, unsigned int mapping_matrix_id)
    {
    	m_coordinates.push_back( X );
    	m_coordinates.push_back( Y );
    	m_coordinates.push_back( Z );
        m_displacements.push_back( 0.0 );
    	m_displacements.push_back( 0.0 );
    	m_displacements.push_back( 0.0 );
		m_w = W;
    	m_local_id = local_id;
        m_global_id = global_id;
        m_mapping_matrix_id = mapping_matrix_id;
    }

    double getX()
    {
    	return m_coordinates[0];
    }

    double getY()
    {
     	return m_coordinates[1];
    }

    double getZ()
    {
    	return m_coordinates[2];
    }

    double getdX()
    {
    	return m_displacements[0];
    }

    double getdY()
    {
     	return m_displacements[1];
    }

    double getdZ()
    {
    	return m_displacements[2];
    }  

    void setdX(double dX)
    {
        m_displacements[0] = dX;
    }  

    void setdY(double dY)
    {
        m_displacements[1] = dY;
    }  

    void setdZ(double dZ)
    {
        m_displacements[2] = dZ;
    }          

    double getWeight()
    {
    	return m_w;
    }

    int getLocalId()
    {
    	return m_local_id;
    }

    int getGlobalId()
    {
    	return m_global_id;
    }

    int getMappingMatrixId()
    {
    	return m_mapping_matrix_id;
    }

    /// Destructor.
    virtual ~ControlPoint()
    {
    }
    // ==============================================================================


//    evaluat

    // ==============================================================================


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ControlPoint";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ControlPoint";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


private:
    // ==============================================================================
    // General working arrays
    // ==============================================================================
    DoubleVector m_coordinates;
    DoubleVector m_displacements;
    double m_w;
    unsigned int m_local_id;
    unsigned int m_global_id;
    unsigned int m_mapping_matrix_id;


}; // Class ControlPoint

}  // namespace Kratos.

#endif // CONTROL_POINT_H_
