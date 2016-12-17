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

###########################################################################
class Face:
    def __init__(self, json_dict):
        # attributes
        self.boundary_loops = None
        self.brep_id = None
        self.surface = None
        self.swapped_surface_normal = None
        
        self.__dict__ = json_dict
        self.fix_boundary_loops()
        self.fix_surface()

    def fix_boundary_loops(self):
        self.boundary_loops = [BoundaryLoop(json_dict)
                               for json_dict in self.boundary_loops]
        
    def fix_surface(self):
        json_dict, = self.surface
        self.surface = Surface(json_dict)

class Surface:
    def __init__(self, json_dict):
        self.__dict__ = json_dict
    
class BoundaryLoop:
    def __init__(self, json_dict):
        self.boundary_edges = None
        self.brep_id = None
        self.loop_type = None

        self.__dict__ = json_dict
        self.fix_boundary_edges()

    def fix_boundary_edges(self):
        self.boundary_edges = [BoundaryEdge(json_dict)
                               for json_dict in self.boundary_edges]

class BoundaryEdge:
    def __init__(self, json_dict):
        # attributes
        self.boundary_vertices = None
        self.brep_id = None
        self.parameter_curve = None
        
        self.__dict__ = json_dict
        self.fix_boundary_vertices()
        self.fix_parameter_curve()

    def fix_boundary_vertices(self):
        self.boundary_vertices = [BoundaryVertex(json_dict)
                                  for json_dict in self.boundary_vertices]
    def fix_parameter_curve(self):
        self.parameter_curve = ParameterCurve(self.parameter_curve)

class BoundaryVertex:
    def __init__(self, json_dict):
        # attributes
        self.parameter_value = None
        
        self.__dict__ = json_dict

class ParameterCurve:
    def __init__(self, json_dict):
        # attributes
        self.control_points = None
        self.degrees = None
        self.is_rational = None
        self.u_vec = None
        
        self.__dict__ = json_dict
        self.fix_control_points()
        
    def fix_control_points(self):
        self.control_points = [ControlPoint(parsed_list)
                               for parsed_list in self.control_points]            

class ControlPoint:
    def __init__(self, parsed_list):
        self.id = parsed_list[0]    # Control Point id
        self.x = parsed_list[1]     # x-coordinate
        self.y = parsed_list[2]     # y-coordinate
        self.z = parsed_list[3]     # z-coordinate
        self.w = parsed_list[4]     # weight

###############################################################################

# Read the FE model
fe_model_part = ModelPart("name_of_empty_mdpa")
model_part_io = ModelPartIO("FEM_mesh")
model_part_io.ReadModelPart(fe_model_part)

# Read CAD model
path = "/home/giovanni/workspace/kratos/applications/ShapeOptimizationApplication/test_examples/CAD_reconstruction/simple example_geometry.json"
face = None
with open(path) as data_file:
    data = json.load(data_file)
    #pprint(data)
    json_face, = data['faces']

    face = Face(json_face)

mapper = IBRAMapper(fe_model_part,face.__dict__)

mapper.compute_mapping_matrix()

# mapper.map_to_CAD_space()




