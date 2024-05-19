from math import log10

import numpy as np
import openmc
import to_cubit_journal



#Arbitrarily rotated plane
if True:
    box = openmc.model.RectangularParallelepiped(*3*[-10, 10])
    p = openmc.XPlane(x0=5)
    cell = openmc.Cell(region=-box & +p)
    geometry = openmc.Geometry([cell])
    geometry.export_to_xml()
    to_cubit_journal.to_cubit_journal( geometry, world=[30,30,30], filename="box_pane.jou" )

if True:
    box = openmc.model.RectangularParallelepiped(*3*[-10, 10])
    p = openmc.XPlane(x0=5)
    cell = openmc.Cell(region=-box & +p)
    cell.region = -box & -p
    geometry = openmc.Geometry([cell])
    geometry.export_to_xml()
    to_cubit_journal.to_cubit_journal( geometry, world=[30,30,30], filename="box_pane2.jou"  )


if True:
    box = openmc.model.RectangularParallelepiped(*3*[-10, 10])
    p = openmc.XPlane(x0=5)
    new_p = p.rotate((10, 20, 30)) # in degrees
    cell = openmc.Cell(region=-box & +p)
    cell.region = -box & -new_p
    geometry = openmc.Geometry([cell])
    geometry.export_to_xml()
    to_cubit_journal.to_cubit_journal( geometry, world=[60,60,60], filename="box_rotated_pane.jou"  )


#Arbitrarily rotated cylinder

if True:
    cyl = openmc.ZCylinder(r=5)
    print( cyl )
    cell = openmc.Cell(region=-cyl)
    cell.region = -cyl
    geometry = openmc.Geometry([cell])
    geometry.export_to_xml()
    to_cubit_journal.to_cubit_journal( geometry, world=[30,30,30], filename="cylinder.jou"  )

if True:
    cyl = openmc.ZCylinder(r=5)
    cell = openmc.Cell(region=-cyl)
    new_cyl = cyl.rotate((30, 20, 30)) 
    print( new_cyl )
    cell.region = -new_cyl
    geometry = openmc.Geometry([cell])
    geometry.export_to_xml()
    to_cubit_journal.to_cubit_journal( geometry, world=[30,30,30] , filename="rotated_cylinder.jou" )

# Translation Example

if True:
    cyl = openmc.ZCylinder(r=5)
    cell = openmc.Cell(region=-cyl)
    new_cyl = cyl.rotate((30, 20, 30)) 
    new_new_cyl = new_cyl.translate((3, 0, 0))
    print( new_new_cyl )
    cell.region = -new_new_cyl
    geometry = openmc.Geometry([cell])
    geometry.export_to_xml()
    to_cubit_journal.to_cubit_journal( geometry, world=[30,30,30] , filename="rotated_translated_cylinder.jou" )
                                

#Cone example

if True:
    cone = openmc.XCone(r2=1.0)
    start = openmc.XPlane(x0=1.0)
    stop = openmc.XPlane(x0=3.0)
    cone_cell = openmc.Cell(region=+start & -stop & -cone)
    geometry = openmc.Geometry([cone_cell])
    geometry.export_to_xml( path="cone.xml" )
    to_cubit_journal.to_cubit_journal( geometry, world=[30,30,30], filename="cone.jou"  )


if True:
    cone = openmc.XCone(r2=1.0)
    start = openmc.XPlane(x0=1.0)
    start.x0 = 0
    stop = openmc.XPlane(x0=3.0)
    cone_cell = openmc.Cell(region=+start & -stop & -cone)
    geometry = openmc.Geometry([cone_cell])
    geometry.export_to_xml()
    to_cubit_journal.to_cubit_journal( geometry, world=[30,30,30], filename="cone2.jou"  )


if True:
    stop = openmc.XPlane(x0=3.0)
    cone = openmc.model.XConeOneSided(r2=1.0)
    cone_cell = openmc.Cell(region=-cone)
    cone_cell.region &= -stop
    geometry = openmc.Geometry([cone_cell])
    geometry.export_to_xml()
    to_cubit_journal.to_cubit_journal( geometry, world=[30,30,30], filename="cone3.jou"  )

#Complement example
if True:
    top = openmc.ZPlane(z0=10)
    bottom = openmc.ZPlane(z0=-10)
    cyl1 = openmc.ZCylinder(r=5)
    cyl2 = openmc.ZCylinder(r=2)
    inside_cyl = -cyl1 & + cyl2 & -top & +bottom
    cell = openmc.Cell(region=inside_cyl)

    geometry = openmc.Geometry([cell])
    geometry.export_to_xml()
    to_cubit_journal.to_cubit_journal( geometry, world=[30,30,30], filename="complement.jou" )

# same plot using the complement of this region
if True:
    top = openmc.ZPlane(z0=10)
    bottom = openmc.ZPlane(z0=-10)
    cyl1 = openmc.ZCylinder(r=5)
    cyl2 = openmc.ZCylinder(r=2)
    inside_cyl = -cyl1 & + cyl2 & -top & +bottom
    cell = openmc.Cell(region=inside_cyl)
    cell.region = ~cell.region

    geometry = openmc.Geometry([cell])
    geometry.export_to_xml()
    to_cubit_journal.to_cubit_journal( geometry, world=[30,30,20], filename="complement2.jou"  )
