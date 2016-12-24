//// ==============================================================================
///*
// KratosShapeOptimizationApplication
// A library based on:
// Kratos
// A General Purpose Software for Multi-Physics Finite Element Analysis
// (Released on march 05, 2007).
//
// Copyright (c) 2016: Daniel Baumgaertner
//                     daniel.baumgaertner@tum.de
//                     Chair of Structural Analysis
//                     Technische Universitaet Muenchen
//                     Arcisstrasse 21 80333 Munich, Germany
//
// Permission is hereby granted, free  of charge, to any person obtaining
// a  copy  of this  software  and  associated  documentation files  (the
// "Software"), to  deal in  the Software without  restriction, including
// without limitation  the rights to  use, copy, modify,  merge, publish,
// distribute,  sublicense and/or  sell copies  of the  Software,  and to
// permit persons to whom the Software  is furnished to do so, subject to
// the following condition:
//
// Distribution of this code for  any  commercial purpose  is permissible
// ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.
//
// The  above  copyright  notice  and  this permission  notice  shall  be
// included in all copies or substantial portions of the Software.
//
// THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
// EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
// CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
// TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//*/
////==============================================================================
////
////   Project Name:        KratosShape                            $
////   Created by:          $Author:    giovanni.filomeno@tum.de   $
////   Last modified by:    $Co-Author: giovanni.filomeno@tum.de   $
////   Date:                $Date:                      Decem 2016 $
////   Revision:            $Revision:                         0.0 $
////
//// ==============================================================================
//
//#ifndef IBRA_MAPPER_H
//#define IBRA_MAPPER_H
//
//// ------------------------------------------------------------------------------
//// System includes
//// ------------------------------------------------------------------------------
//#include <iostream>
//#include <string>
//#include <algorithm>
//#include <Python.h>
//
//// ------------------------------------------------------------------------------
//// External includes
//// ------------------------------------------------------------------------------
//#include <boost/python.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>
//
//// ------------------------------------------------------------------------------
//// Project includes
//// ------------------------------------------------------------------------------
//#include "../../kratos/includes/define.h"
//#include "../../kratos/processes/process.h"
//#include "../../kratos/includes/node.h"
//#include "../../kratos/includes/element.h"
//#include "../../kratos/includes/model_part.h"
//#include "../../kratos/includes/kratos_flags.h"
//#include "shape_optimization_application.h"
//#include "/home/giovanni/workspace/kratos/applications/ShapeOptimizationApplication/custom_utilities/controlPoint.hpp"
//#include "/home/giovanni/workspace/kratos/applications/ShapeOptimizationApplication/test_examples/CAD_reconstruction/inc/nurbs.hpp"
//// ==============================================================================
//
//namespace Kratos
//{
//
/////@name Kratos Globals
/////@{
//
/////@}
/////@name Type Definitions
/////@{
//
//
/////@}
/////@name  Enum's
/////@{
//
/////@}
/////@name  Functions
/////@{
//
/////@}
/////@name Kratos Classes
/////@{
//
///// Short class definition.
///** Detail class definition.
//
//*/
//
//class IBRAMapper
//{
//public:
//    ///@name Type Definitions
//    ///@{
//
//
//    // ==========================================================================
//    // Type definitions for linear algebra including sparse systems
//    // ==========================================================================
//    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//    typedef typename SparseSpaceType::MatrixType SparseMatrixType;
//    typedef typename SparseSpaceType::VectorType VectorType;
//    typedef std::vector<double> DoubleVector;
//    typedef boost::python::extract<double> takeDouble;
//    typedef boost::python::extract<int> takeInt;
//    typedef boost::python::extract<bool> takeBool;
//    typedef std::vector<controlPoint> controlPointVcr;
//    typedef std::vector<int> IntVector;
////    static PyObject *;
//
//    /// Pointer definition of IBRAMapper
//    KRATOS_CLASS_POINTER_DEFINITION(IBRAMapper);
//
//    ///@}
//    ///@name Life Cycle
//    ///@{
//
//    /// Default constructor.
//    IBRAMapper( ModelPart& model_part,
//    		boost::python::dict self)
//		: mr_model_part(model_part),
//          myPythondict(self)
//    {
//        // Initialize filter matrix
//        //m_mapping_matrix.resize(m_number_of_design_variables*3,m_number_of_design_variables);
//
//    }
//
//    /// Destructor.
//    virtual ~IBRAMapper()
//    {
//    }
//
//
//    ///@}
//    ///@name Operators
//    ///@{
//
//
//    ///@}
//    ///@name Operations
//    ///@{
//
////    void compute_mapping_matrix()
////    {
////    	DoubleVector knotOne = {0,0,0,5,5,5};
////    	DoubleVector knotTwo = { 0,0,0,5,5,5};
////
////    	int grade = 2;
////
////    	Nurbs *myTestForNurbs = new Nurbs (grade);
////    	Nurbs *mySecondTest = new Nurbs (grade);
////    	myTestForNurbs->setKnotVector( knotOne );
////    	mySecondTest->setKnotVector( knotTwo );
////
////    	DoubleVector XCP;
////    	XCP.push_back(0);
////    	XCP.push_back(2.5);
////    	XCP.push_back(5);
////    	XCP.push_back(0);
////    	XCP.push_back(2.5);
////    	XCP.push_back(5);
////		XCP.push_back(0);
////		XCP.push_back(2.5);
////		XCP.push_back(5);
////    	DoubleVector YCP = { 0,0,0,2.5,2.5,2.5,5,5,5 };
////    	DoubleVector ZCP ={ 0,0,0,0,0,0,0,0,0};
////    	DoubleVector WCP = {1,1,1,1,1,1,1,1,1};
////    	myTestForNurbs->setCPX( XCP );
////    	myTestForNurbs->setCPY( YCP );
////    	myTestForNurbs->setCPZ( ZCP );
////    	myTestForNurbs->setCPW( WCP );
////
////    	double t = 0;
////    	double x=0;
////    	int i;
////
////    	for(int j=0; j < 9; j++)
////    	{
////    		if(XCP[j] == 0)
////    		{
////    			i = 0;
////    		}
////    		if(XCP[j] == 2.5)
////    		{
////    			i = 1;
////    		}
////    		if(XCP[j] == 5)
////    		{
////    		 	i = 2;
////    		}
////    		x = x + myTestForNurbs->N(t,i,2);
////    	}
////
////    	std::cout << "Result" << x << std::endl;
////
////    }
//    // ==============================================================================
//    void compute_mapping_matrix()
//    {
//        KRATOS_TRY;
//
////        VectorType knot_vector = ZeroVector(3);
//
//        double test;
//
//        DoubleVector coordinates;
//        controlPointVcr myControlPoint;
//        IntVector degree;
//        DoubleVector knotVector;
//
//
//
//
//        for( int k = 0; k < 9; k++)
//        {
//
//        	std::cout << "checkmark 100" << std::endl;
//
//        	int list_length = boost::python::len(myPythondict["faces"]);
//
//
//
//
//
//        	std::cout << "list_length = " << list_length << std::endl;
//
//
//
//
//
//
//
//
//        	takeInt ID(myPythondict["faces"][0]["surface"][0]["control_points"][k][0]);
//
//        	for(int i = 0; i < 3; i++)
//        	{
//        		std::cout << "Porco dio nel secondo for" << std::endl;
//        	  	double myCoord = takeDouble (myPythondict["faces"][0]["surface"][0]["control_points"][0][k][1][i]);
//        	   	coordinates.push_back( myCoord );
//        	}
//
//        	double weight = takeDouble (myPythondict["faces"][0]["surface"][0]["control_points"][0][k][1][3]);
//
//        	controlPoint controlPointFromJson( coordinates, weight, ID);
//        	myControlPoint.push_back( controlPointFromJson );
//
//        }
//
//        for(int k = 0; k < 2; k++)
//        {
//        	int deg = takeInt (myPythondict["faces"][0]["surface"][0]["degree"][k]);
//        	degree.push_back( deg );
//        }
//
////        bool is_rational = takeBool (myPythondict["faces"][0]["surface"][0]["is_rational"]);
////        bool is_trimmed = takeBool (myPythondict["faces"][0]["surface"][0]["is_trimmed"]);
//
//        for(int i = 0; i < 2; i++)
//        {
//        	for(int k = 0; k < 6; k++)
//        	{
//        		double knot = takeDouble (myPythondict["faces"][0]["surface"][0]["knot_vectors"][i][k]);
//        		knotVector.push_back( knot );
//        		std::cout << "Knot " << knot << std::endl;
//        	}
//        }
//
//
////        double test = boost::python::extract<double> (myPythondict["faces"]);
////
////        std::cout<<"Test1 " << test << std::endl;
//
//
////        int i=0;
////        while (i < myPythondict["knot_vector"].attr("length") )
////       {
////        	boost::python::extract<double> test2 (myPythondict["knot_vector"][i]);
////        	std::cout <<"vector test " << test2 << std::endl;
////        	i++;
////       }
//
////        boost::python::extract<const double *> test2 (myPythondict["knot_vector"]);
////
////        std::cout<<"vector test " << test2 << std::endl;
//
//
//
//
////        string test2;
////
////        boost::python::extract<const char*> test2(myPythondict["loop"]);
////
////        std::cout<<"Test2 " << test2 << std::endl;
////
////        boost::python::extract<double> test3(myPythondict["knot_vector"][0]);
////        std::cout << "Test3 " << test3 << std::endl;
////
////        boost::python::extract<Vec2&> test4(myPythondict["knot_vector"]);
////        std::cout<<"Test4 " << test4 << std::endl;
////        KRATOS_WATCH(test4);
//
////        KRATOS_WATCH(knot_vector);
//
//        KRATOS_CATCH("");
//    }
//
////    evaluat
//
//    // ==============================================================================
//
//    ///@}
//    ///@name Access
//    ///@{
//
//
//    ///@}
//    ///@name Inquiry
//    ///@{
//
//
//    ///@}
//    ///@name Input and output
//    ///@{
//
//    /// Turn back information as a string.
//    virtual std::string Info() const
//    {
//        return "IBRAMapper";
//    }
//
//    /// Print information about this object.
//    virtual void PrintInfo(std::ostream& rOStream) const
//    {
//        rOStream << "IBRAMapper";
//    }
//
//    /// Print object's data.
//    virtual void PrintData(std::ostream& rOStream) const
//    {
//    }
//
//
//    ///@}
//    ///@name Friends
//    ///@{
//
//
//    ///@}
//
//protected:
//    ///@name Protected static Member Variables
//    ///@{
//
//
//    ///@}
//    ///@name Protected member Variables
//    ///@{
//
//
//    ///@}
//    ///@name Protected Operators
//    ///@{
//
//
//    ///@}
//    ///@name Protected Operations
//    ///@{
//
//
//    ///@}
//    ///@name Protected  Access
//    ///@{
//
//
//    ///@}
//    ///@name Protected Inquiry
//    ///@{
//
//
//    ///@}
//    ///@name Protected LifeCycle
//    ///@{
//
//
//    ///@}
//
//private:
//    ///@name Static Member Variables
//    ///@{
//
//
//    ///@}
//    ///@name Member Variables
//    ///@{
//
//    // ==============================================================================
//    // Initialized by class constructor
//    // ==============================================================================
//    ModelPart& mr_model_part;
//    boost::python::dict myPythondict;
//
//    // ==============================================================================
//    // General working arrays
//    // ==============================================================================
//    SparseMatrixType m_mapping_matrix;
//
//    ///@}
//    ///@name Private Operators
//    ///@{
//
//
//    ///@}
//    ///@name Private Operations
//    ///@{
//
//
//    ///@}
//    ///@name Private  Access
//    ///@{
//
//
//    ///@}
//    ///@name Private Inquiry
//    ///@{
//
//
//    ///@}
//    ///@name Un accessible methods
//    ///@{
//
//    /// Assignment operator.
////      IBRAMapper& operator=(IBRAMapper const& rOther);
//
//    /// Copy constructor.
////      IBRAMapper(IBRAMapper const& rOther);
//
//
//    ///@}
//
//}; // Class IBRAMapper
//
/////@}
//
/////@name Type Definitions
/////@{
//
//
/////@}
/////@name Input and output
/////@{
//
/////@}
//
//
//}  // namespace Kratos.
//
//#endif // IBRA_MAPPER_H
