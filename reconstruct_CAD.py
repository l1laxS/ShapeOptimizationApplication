# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------
# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Read model part
fe_model_part = ModelPart("name_of_empty_mdpa")
model_part_io = ModelPartIO("FEM_mesh")
model_part_io.ReadModelPart(fe_model_part)

for node in fe_model_part.Nodes:
	print(node)






