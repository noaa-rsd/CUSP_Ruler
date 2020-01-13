CUSP Ruler
----------

CUSP Ruler is a NOAA Remote Sensing Division ArcPro tool that calculates the length of Continually Updated Shoreline Product (CUSP) shoreline relative to a specified reference shoreline dataset.

GUI
===

- *Reference Shoreline*:  feature class containing reference shoreline
- *CUSP*:  feature class containing CUSP
- *Output Directory*:  output directory where results table and products will be generated

.. image:: ./support_files/gui.png

Basic Workflow
==============

The CUSP Ruler workflow has four main steps:

1. Generalize the CUSP datashet
2. Buffer the generalized CUSP datashet
3. Clip reference shoreline with CUSP buffer
4. Measure geodesic lengths of clipped reference shoreline

1.  Generalize the CUSP datashet
++++++++++++++++++++++++++++++++
The CUSP dataset is generalized to computationaly simplify the proceeding steps.  The Geopandas simplify algorithm is used, with a tolerance of 0.002 (degrees) and perserve_topology=False.

2. Buffer the generalized CUSP datashet
+++++++++++++++++++++++++++++++++++++++
The CUSP data are buffered as way to determine which reference shoreline segments the CUSP corresponds to.  The reference shoreline segments intersecting the CUSP buffer are considered to correspsond to the associated CUSP.  The simplified CUSP are buffered using the shapely buffer method, with each line segment transformed to a local azimuthal equidistand projection.



3. Clip reference shoreline with CUSP buffer
++++++++++++++++++++++++++++++++++++++++++++
The reference shoreline are clipped with the CUSP buffer using the Geopandas intersecion method.

4. Measure geodesic lengths of clipped reference shoreline
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The geodesic lengths of the clipped reference shoreline segment are calcualted using the inverse method of the cartopy package.

Output
======

-

.. image:: ./support_files/CUSP_Overview_Map.png
