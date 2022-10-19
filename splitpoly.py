# -*- coding: utf-8 -*-

# * Copyright (c) 2009-2018. Authors: see NOTICE file.
# *
# * Licensed under the Apache License, Version 2.0 (the "License");
# * you may not use this file except in compliance with the License.
# * You may obtain a copy of the License at
# *
# *      http://www.apache.org/licenses/LICENSE-2.0
# *
# * Unless required by applicable law or agreed to in writing, software
# * distributed under the License is distributed on an "AS IS" BASIS,
# * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# * See the License for the specific language governing permissions and
# * limitations under the License.

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import logging
import sys
from argparse import ArgumentParser

import os
import numpy as np
from shapely.geometry import shape, box, Polygon, Point, MultiPolygon, LineString
from shapely import wkt
from shapely.ops import split
import geopandas

from glob import glob
from tifffile import imread

from cytomine import Cytomine, CytomineJob
from cytomine.models import Property, Annotation, AnnotationTerm, AnnotationCollection, Job, TermCollection
# from cytomine.models.ontology import Ontology, OntologyCollection, Term, RelationTerm, TermCollection
# from cytomine.models.property import Tag, TagCollection, PropertyCollection
# from cytomine.utilities.software import parse_domain_list, str2bool, setup_classify, stringify


from PIL import Image
import matplotlib.pyplot as plt
import time
import cv2
import math

from shapely import affinity
from shapely.geometry.multipolygon import MultiPolygon
from scipy.spatial import Voronoi



__author__ = "WSH Munirah W Ahmad <wshmunirah@gmail.com>"
__version__ = "1.0.0"
# Date created: 25 August 2021

def _quadrat_cut_geometry(geometry, quadrat_width, min_num=3):
    """
    Split a Polygon or MultiPolygon up into sub-polygons of a specified size.

    Parameters
    ----------
    geometry : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the geometry to split up into smaller sub-polygons
    quadrat_width : numeric
        the linear width of the quadrats with which to cut up the geometry (in
        the units the geometry is in)
    min_num : int
        the minimum number of linear quadrat lines (e.g., min_num=3 would
        produce a quadrat grid of 4 squares)

    Returns
    -------
    geometry : shapely.geometry.MultiPolygon
    """
    # create n evenly spaced points between the min and max x and y bounds
    west, south, east, north = geometry.bounds
    x_num = math.ceil((east - west) / quadrat_width) + 1
    y_num = math.ceil((north - south) / quadrat_width) + 1
    x_points = np.linspace(west, east, num=max(x_num, min_num))
    y_points = np.linspace(south, north, num=max(y_num, min_num))

    # create a quadrat grid of lines at each of the evenly spaced points
    vertical_lines = [LineString([(x, y_points[0]), (x, y_points[-1])]) for x in x_points]
    horizont_lines = [LineString([(x_points[0], y), (x_points[-1], y)]) for y in y_points]
    lines = vertical_lines + horizont_lines

    # recursively split the geometry by each quadrat line
    for line in lines:
        geometry = MultiPolygon(split(geometry, line))

    return geometry

#============================================

def main(argv):
    with CytomineJob.from_cli(argv) as conn:
    # with Cytomine(argv) as conn:
        print(conn.parameters)
        id_user=conn.parameters.cytomine_id_user

        conn.job.update(status=Job.RUNNING, progress=0, statusComment="Initialization...")
        base_path = "{}".format(os.getenv("HOME")) # Mandatory for Singularity
        working_path = os.path.join(base_path,str(conn.job.id))


        terms = TermCollection().fetch_with_filter("project", conn.parameters.cytomine_id_project)
        conn.job.update(status=Job.RUNNING, progress=1, statusComment="Terms collected...")
        print(terms)
        for term in terms:
            print("ID: {} | Name: {}".format(
                term.id,
                term.name
            )) 

        id_project=conn.parameters.cytomine_id_project
        id_image = conn.parameters.cytomine_id_images
        id_user=conn.parameters.cytomine_id_user
        
        
        id_term = conn.parameters.cytomine_id_roi_term
        id_term_poly = conn.parameters.cytomine_id_roipoly_term
        
        poly_sides = conn.parameters.cytomine_poly_sides

        print('parameters:',id_project, id_image, id_term, id_term_poly)


        roi_annotations = AnnotationCollection()
        roi_annotations.project = id_project
        roi_annotations.image = id_image
        roi_annotations.term = id_term
        roi_annotations.showWKT = True
        roi_annotations.showMeta = True
        roi_annotations.showGIS = True
        roi_annotations.showTerm = True
        if id_user:
            roi_annotations.user = id_user
        roi_annotations.fetch()
        print(roi_annotations)

        conn.job.update(status=Job.RUNNING, progress=10, statusComment="Running splitpoly on ROI-WSI...")

        for i, roi in enumerate(roi_annotations):
                #Get Cytomine ROI coordinates for remapping to whole-slide
                #Cytomine cartesian coordinate system, (0,0) is bottom left corner                
                print("----------------------------Cells------------------------------")
                roi_geometry = wkt.loads(roi.location)
                # print("ROI Geometry from Shapely: {}".format(roi_geometry))
                print("ROI Bounds")
                print(roi_geometry.bounds)
                minx=roi_geometry.bounds[0]
                miny=roi_geometry.bounds[3]
                #Dump ROI image into local PNG file
                # roi_path=os.path.join(working_path,str(roi_annotations.project)+'/'+str(roi_annotations.image)+'/'+str(roi.id))
                roi_path=os.path.join(working_path,str(roi_annotations.project)+'/'+str(roi_annotations.image)+'/')
                print(roi_path)
                roi_png_filename=os.path.join(roi_path+str(roi.id)+'.png')
                # conn.job.update(status=Job.RUNNING, progress=20, statusComment=roi_png_filename)
                print("roi_png_filename: %s" %roi_png_filename)
                output = _quadrat_cut_geometry(roi_geometry, quadrat_width=poly_sides, min_num=1)  
                print(output)
#                 cytomine_annotations = AnnotationCollection()
                
                annotations = AnnotationCollection()
                for annotation_poly in output:                  
                    annotations.append(Annotation(
                        location=annotation_poly.wkt,
                        id_terms=[id_term_poly],
                        id_project=id_project,
                        id_image=id_image
                    ))
                annotations.save()
                print(".",end = '',flush=True)

 
        conn.job.update(status=Job.TERMINATED, progress=100, statusComment="Finished.")

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])

                  






