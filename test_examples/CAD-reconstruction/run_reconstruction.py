# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import json as json

# ======================================================================================================================================
# Parameters
# ======================================================================================================================================

fem_input_filename = "FEM_mesh_simple_example_triangles"
# "FEM_mesh_simple_example" 
# "FEM_mesh_simple_example_triangles" 
# "FEM_mesh_hook" !!!!!!!!!!!!!! FEM shape function evaluation not working for this triangular elements
cad_input_filename = "simple_example_geometry.json" 
# "simple_example_geometry.json"
# "simple_example_geometry_pq3x3"
# "simple_example_geometry_pq3x3_ele4x4"
# "hock_geometry.json" 

# ======================================================================================================================================
# Reconstruction part
# ======================================================================================================================================

# Read the FE model
fe_model_part = ModelPart("name_of_empty_mdpa")
model_part_io = ModelPartIO(fem_input_filename)
model_part_io.ReadModelPart(fe_model_part)

# Read CAD model
cad_model = False
with open(cad_input_filename) as cad_data:
    cad_model = json.load(cad_data)

# Create CAD-mapper
linear_solver = SuperLUSolver()
# DiagPrecond = DiagonalPreconditioner()
# linear_solver =  BICGSTABSolver(1e-9, 300, DiagPrecond)
# linear_solver = AMGCLSolver(AMGCLSmoother.GAUSS_SEIDEL, AMGCLIterativeSolverType.BICGSTAB, 1e-9, 300, 2, 10)
domain_size = 3
mapper = CADMapper(fe_model_part,cad_model,domain_size,linear_solver)

# ======================================================================================================================================
# Only for flat plate model
# ======================================================================================================================================

# Compute mapping matrix
mapper.compute_mapping_matrix(1000,1000)

# Set shape update to map

# Constant displacement
for node in fe_model_part.Nodes:
    shape_update = Vector(3)
    shape_update[0] = 0.0
    shape_update[1] = 0.0
    shape_update[2] = 1.0
    node.SetValue(SHAPE_UPDATE,shape_update)

# # Move middle node
# shape_update = Vector(3)
# shape_update[0] = 0.0
# shape_update[1] = 0.0
# shape_update[2] = 1.0
# fe_model_part.Nodes[4].SetValue(SHAPE_UPDATE,shape_update)

# # Linear boundary update
# for node in fe_model_part.Nodes:
#     shape_update = Vector(3)
#     shape_update[0] = 0.0
#     shape_update[1] = 0.0
#     shape_update[2] = node.Y
#     node.SetValue(SHAPE_UPDATE,shape_update)

# Perform mapping
mapper.map_to_cad_space()

# Output control point update in gid-format 
mapper.output_control_point_displacements() 

# ======================================================================================================================================
# Only for circular shell model
# ======================================================================================================================================

# ======================================================================================================================================
# Only for hook model
# ======================================================================================================================================

