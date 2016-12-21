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

import json
from pprint import pprint

#########################################################################

###############################################################################

# Read the FE model
fe_model_part = ModelPart("name_of_empty_mdpa")
model_part_io = ModelPartIO("FEM_mesh")
model_part_io.ReadModelPart(fe_model_part)

# Read CAD model
# path = "simple_example_geometry.json"
path = "hock_geometry.json"
face = None
with open(path) as data_file:
    data = json.load(data_file)
    pprint(data)

mapper = IBRAMapper(fe_model_part,data)

mapper.initialize_patches()

#mapper.Surface()


patch_id = 1
u = 2.5
v = 2.5

point = mapper.compute_point(patch_id, u, v)
print(point)

# mapper.map_to_CAD_space()




