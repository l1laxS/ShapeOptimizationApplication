from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *

# Read the FE model
fe_model_part = ModelPart("name_of_empty_mdpa")
model_part_io = ModelPartIO("FEM_mesh")
model_part_io.ReadModelPart(fe_model_part)


my_dict = {'faces': 1, 'loop':'Outer', 'knot_vector':[1,2,3],
           'some_bool': True};

CAD_Reconstructor = IBRAMapper(fe_model_part, my_dict)
CAD_Reconstructor.compute_mapping_matrix()
