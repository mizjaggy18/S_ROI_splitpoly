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
import logging
import shutil

import os
import numpy as np
from shapely.geometry import shape, box, Polygon, Point, MultiPolygon, LineString
from shapely import wkt
from shapely.ops import split
import geopandas

from glob import glob
from tifffile import imread

import cytomine
from cytomine import Cytomine, CytomineJob
from cytomine.models import Property, Annotation, AnnotationTerm, AnnotationCollection, Job, JobData, TermCollection, ImageInstanceCollection


from PIL import Image
import matplotlib.pyplot as plt
import time
import cv2
import math

from shapely import affinity
from shapely.geometry.multipolygon import MultiPolygon
from scipy.spatial import Voronoi



__author__ = "WSH Munirah W Ahmad <wshmunirah@gmail.com>"
__version__ = "2.1.0"
# Date created: 20 Oct 2022

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

def save_annotations (output,id_image,id_project,id_term_poly):
    annotations = AnnotationCollection()
    # print("Annotation collections: ", annotations)
    for i, annotation_poly in enumerate(output): 
        annotations.append(Annotation(location=annotation_poly.wkt,
                                    id_image=id_image,
                                    id_project=id_project,
                                    id_terms=[id_term_poly]))
        print(".",end = '',flush=True)
    annotations.save()

def run(cyto_job, parameters):
    logging.info("----- ROI Split-poly v%s -----", __version__)
    logging.info("Entering run(cyto_job=%s, parameters=%s)", cyto_job, parameters)

    job = cyto_job.job
    project = cyto_job.project
    id_user = parameters.cytomine_id_user

    job.update(status=Job.RUNNING, progress=10, statusComment="Initialization...")

    terms = TermCollection().fetch_with_filter("project", parameters.cytomine_id_project)    
    print(terms)
    for term in terms:
        print("ID: {} | Name: {}".format(
            term.id,
            term.name
        )) 
    job.update(status=Job.RUNNING, progress=20, statusComment="Terms collected...")
    
    images = ImageInstanceCollection().fetch_with_filter("project", project.id)    
    list_imgs = []
    if parameters.cytomine_id_images == 'all':
        for image in images:
            list_imgs.append(int(image.id))
    else:
        list_imgs = [int(id_img) for id_img in parameters.cytomine_id_images.split(',')]
        print('Images: ', list_imgs)    
    job.update(status=Job.RUNNING, progress=30, statusComment="Images gathered...")
         
    id_project = parameters.cytomine_id_project
    id_user = parameters.cytomine_id_user   
    id_term = parameters.cytomine_id_roi_term
    id_term_poly = parameters.cytomine_id_roipoly_term
    poly_sides = parameters.cytomine_poly_sides
    
    working_path = os.path.join("tmp", str(job.id))
    
    if not os.path.exists(working_path):
        logging.info("Creating working directory: %s", working_path)
        os.makedirs(working_path)
    try:

        for id_image in list_imgs:
            print('Parameters (id_project, id_image, id_term, id_term_poly):',id_project, id_image, id_term, id_term_poly)

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

            job.update(status=Job.RUNNING, progress=40, statusComment="Running splitpoly on ROI-WSI...")

            for i, roi in enumerate(roi_annotations):
                    #Get Cytomine ROI coordinates for remapping to whole-slide
                    #Cytomine cartesian coordinate system, (0,0) is bottom left corner                
                    print("----------------------------Cells------------------------------")
                    roi_geometry = wkt.loads(roi.location)
                    # print("ROI Geometry from Shapely: {}".format(roi_geometry))
                    print("ROI Bounds")
                    print(roi_geometry.bounds)
                    minx=roi_geometry.bounds[0]
                    maxx=roi_geometry.bounds[2]
                    miny=roi_geometry.bounds[1]
                    maxy=roi_geometry.bounds[3]
                    widthROI=maxx-minx
                    heightROI=maxy-miny
                    #Dump ROI image into local PNG file
                    # roi_path=os.path.join(working_path,str(roi_annotations.project)+'/'+str(roi_annotations.image)+'/'+str(roi.id))
                    roi_path=os.path.join(working_path,str(roi_annotations.project)+'/'+str(roi_annotations.image)+'/')
                    print(roi_path)
                    roi_png_filename=os.path.join(roi_path+str(roi.id)+'.png')
                    print("roi_png_filename: %s" %roi_png_filename)
                    
                    if poly_sides <= 512: 
                        if widthROI > 10000 or heightROI > 10000:
                            print("Large ROI")
                            largepoly_sides=4096
                            job.update(status=Job.RUNNING, progress=50, statusComment="Split ROI-WSI into 4096px sides...")
                            
                            outputlarge = _quadrat_cut_geometry(roi_geometry, quadrat_width=largepoly_sides, min_num=1)
                            outputlarge_poly = len(list(outputlarge))
                            print("Output large polygons: ", outputlarge_poly)

                            for i, roi_large in enumerate(outputlarge):
                                print("Current ROI4096:", str(i), "out of: ", str(outputlarge_poly))
                                job.update(status=Job.RUNNING, progress=60, statusComment="Split ROI-WSI into user-defined sides...")

                                output = _quadrat_cut_geometry(roi_large, quadrat_width=poly_sides, min_num=1)
                                output_poly = list(output)
                                print("Output polygons: ", len(output_poly))                                
                                save_annotations (output,id_image,id_project,id_term_poly)

                    else:
                        output = _quadrat_cut_geometry(roi_geometry, quadrat_width=poly_sides, min_num=1)
                        output_poly = list(output)
                        print("Output polygons: ", len(output_poly))
                        save_annotations (output,id_image,id_project,id_term_poly)
                    
#                     output = _quadrat_cut_geometry(roi_geometry, quadrat_width=poly_sides, min_num=1)  
#                     output_poly = list(output)
#                     print("Output polygons: ", len(output_poly))
                    
#                     annotations = AnnotationCollection()
#                     print("Annotation collections: ", annotations)
#                     for annotation_poly in output: 
#                         annotation=Annotation(location=annotation_poly.wkt,
#                                    id_image=id_image,
#                                    id_project=id_project).save()
#                         AnnotationTerm(annotation.id, id_term_poly).save()
# #                         annotations.append(Annotation(location=annotation_poly.wkt,
# #                                                       id_image=id_image,
# #                                                       id_project=id_project,
# #                                                       id_terms=[id_term_poly]))
# #                         print(".",end = '',flush=True)
# #                     annotations.save()
                            
                                           
    finally:
        job.update(progress=100, statusComment="Run complete.")
        shutil.rmtree(working_path, ignore_errors=True)
        logging.debug("Leaving run()")
        
if __name__ == "__main__":
    logging.debug("Command: %s", sys.argv)

    with cytomine.CytomineJob.from_cli(sys.argv) as cyto_job:
        run(cyto_job, cyto_job.parameters)

                  






