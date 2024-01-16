#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import sys
import re
from pathlib import Path
import numpy as np

import gmsh
from gmshairfoil2d.airfoil_func import (NACA_4_digit_geom, get_airfoil_points,
                                        get_all_available_airfoil_names)
from gmshairfoil2d.geometry_def import (Airfoil, Circle, PlaneSurface, AirfoilBoundaryCondition,
                                        Rectangle, MeshExtrusion, BoundaryCondition, AirfoilStructuredRegion)

def main():
    # Instantiate the parser
    parser = argparse.ArgumentParser(
        description="Optional argument description",
        usage=argparse.SUPPRESS,
        formatter_class=lambda prog: argparse.HelpFormatter(
            prog, max_help_position=80, width=99
        ),
    )

    parser.add_argument(
        "--list",
        action="store_true",
        help="Display all airfoil available in the database : https://m-selig.ae.illinois.edu/ads/coord_database.html",
    )

    parser.add_argument(
        "--naca",
        type=str,
        metavar="4DIGITS",
        nargs="?",
        help="NACA airfoil 4 digit (default 0012)",
    )

    parser.add_argument(
        "--airfoil",
        type=str,
        metavar="NAME",
        nargs="?",
        help="Name of an airfoil profile in the database (database available with the --list argument)",
    )

    parser.add_argument(
        "--blayer",
        type=bool,
        metavar="NAME",
        nargs="?",
        help="Set it to true for activating boundary layer prisms",
    )


    parser.add_argument(
        "--blayer_thickness",
        type=float,
        metavar="NAME",
        nargs="?",
        help="Thickness of the boundary layer",
    )

    parser.add_argument(
        "--blayer_size",
        type=float,
        metavar="NAME",
        nargs="?",
        help="Wall-normal size of the first boudnary layer cell",
    )


    parser.add_argument(
        "--blayer_ratio",
        type=float,
        metavar="NAME",
        nargs="?",
        help="Cell expansion ratio within the boundary layer",
    )

    parser.add_argument(
        "--aoa",
        type=float,
        nargs="?",
        help="Angle of attack [deg] (default: 0 [deg])",
        default=0.0,
    )

    parser.add_argument(
        "--farfield",
        type=float,
        metavar="RADIUS",
        nargs="?",
        default=10,
        help="Create a circular farfield mesh of given radius [m] (default 10m)",
    )
    parser.add_argument(
        "--box",
        type=str,
        metavar="LENGTHxWIDTH",
        nargs="?",
        help="Create a box mesh of dimensions [length]x[height] [m]",
    )
    parser.add_argument(
        "--airfoil_mesh_size",
        type=float,
        metavar="SIZE",
        nargs="?",
        default=0.01,
        help="Mesh size of the airfoil countour [m]  (default 0.01m)",
    )

    parser.add_argument(
        "--ext_mesh_size",
        type=float,
        metavar="SIZE",
        nargs="?",
        default=0.2,
        help="Mesh size of the external domain [m] (default 0.2m)",
    )

    parser.add_argument(
        "--format",
        type=str,
        nargs="?",
        default="su2",
        help="format of the mesh file, e.g: msh, vtk, wrl, stl, mesh, cgns, su2, dat (default su2)",
    )

    parser.add_argument(
        "--output",
        type=str,
        metavar="PATH",
        nargs="?",
        default=".",
        help="output path for the mesh file (default : current dir)",
    )

    parser.add_argument(
        "--ui", action="store_true", help="Open GMSH user interface to see the mesh"
    )
    
    parser.add_argument(
        "--quad", action="store_true", help="Create quadrangular mesh based on Frontal-Delaunay for Quads (experimental) and the Simple Full-Quad recombination algorithm"
    )

    parser.add_argument(
        "--hopr", action="store_true", help="Create quadrangular mesh based on Frontal-Delaunay for Quads (experimental) and the Simple Full-Quad recombination algorithm"
    )
    
    parser.add_argument(
        "--extrusion",
        type=str,
        metavar="VALUE/NELEMENTS",
        nargs="?",
        help="Extrusion of the 2D mesh in z by [width] [m]",
    )

    parser.add_argument(
        "--le_refinement",
        type=float,
        metavar="LENGTH",
        nargs="?",
        help="Add a refinement in the LE area and match the airfoil mesh size at [length]",
    )
    
    parser.add_argument(
        "--bc",
        type=str,
        metavar="NAME:START..END-SIDE",
        nargs="*",
        help="Create a boundary condition named [name] on the airfoil ranging from the x coord [start] to [end] [%%c] (start < end) on the suction side (side = SS) or the pressure side (side = PS).",
    )

    parser.add_argument(
        "--high_order",
        type=int,
        metavar="ORDER",
        nargs="?",
        help="Add a refinement in the LE area and match the airfoil mesh size at [length]",
    )

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    if args.list:
        get_all_available_airfoil_names()
        sys.exit()

    # Airfoil choice
    cloud_points = None
    if args.naca:
        airfoil_name = args.naca
        cloud_points = NACA_4_digit_geom(airfoil_name)

    if args.airfoil:
        airfoil_name = args.airfoil
        cloud_points = get_airfoil_points(airfoil_name)

    if cloud_points is None:
        print("\nNo airfoil profile specified, exiting")
        print("You must use --naca or --airfoil\n")
        parser.print_help()
        sys.exit()

    blayer = False
    blayer_ratio = 5e-5
    blayer_size = 5e-5
    blayer_thickness = 0.01 
    if args.blayer:
        print('Activating boundary layer prisms')
        blayer = True
        if args.blayer_size:
            blayer_size = args.blayer_size
        if args.blayer_ratio:
            blayer_ratio = args.blayer_ratio
        if args.blayer_thickness:
            blayer_thickness = args.blayer_thickness
        print('        Size:      '+str(blayer_size))
        print('        Ratio:     '+str(blayer_ratio))
        print('        Thickness: '+str(blayer_thickness))


    # Angle of attack
    aoa = -args.aoa * (math.pi / 180)

    # Generate Geometry
    gmsh.initialize()
    gmsh.option.setNumber("General.NumThreads", 32)

    # Boundary conditions
    airfoilBCs = []
    if args.bc:
        delimiters = [':','..', '-']
        for bc in args.bc:
            if not re.match(r".*:\d*\.\d*\.\.\d*\.\d*-(SS|PS)$", bc):
                raise SyntaxError("Boundary conditions were provided with wrong syntax. For more information use --help")
            name, start, end, side = re.split('|'.join(map(re.escape, delimiters)), bc)
            
            # Assign boundary condition values to dictionary
            if float(start) < 1 and float(end) < 1 and float(start) < float(end):
                airfoilBCs.append(AirfoilBoundaryCondition(name, float(start), float(end), side))
            else:
                raise ValueError("Given values are outside of the allowed specification. Make sure they are between 0 and 1 and start < end.")
                
        for bc in airfoilBCs:
            print(f"{bc.name}: {bc.start}..{bc.end} - {bc.side}")
    

    airfoil_mesh_size = args.airfoil_mesh_size
    extMeshSize = args.ext_mesh_size
    if args.blayer_size:
        blayer_size = args.blayer_size
    if args.blayer_thickness:
        blayer_thickness = args.blayer_thickness
    if args.quad:
        airfoil_mesh_size = 2*airfoil_mesh_size
        extMeshSize = 2*extMeshSize
        blayer_size = 2*blayer_size
        gmsh.option.setNumber("Mesh.Algorithm", 8)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1) 
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        if args.extrusion:
            gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)
        else:
            gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)


    # Airfoil
    airfoil = Airfoil(cloud_points, airfoil_mesh_size, airfoilBCs)
    airfoil.rotation(aoa, (0.5, 0, 0), (0, 0, 1))
    airfoil.generate()
    if args.blayer:
        airfoil.setBoundaryLayer(blayer,blayer_size,blayer_ratio,blayer_thickness)
    
    # External domain
    if args.box:
        length, width = [float(value) for value in args.box.split("x")]
        ext_domain = Rectangle(0.5, 0, 0, length, width, mesh_size=extMeshSize)
    else:
        ext_domain = Circle(0.5, 0, 0, radius=args.farfield, mesh_size=extMeshSize)
    
    offsetTrigger = True
    surface_domain = None
    offset = None
    # Generate domain
    if offsetTrigger:
        offset = AirfoilStructuredRegion(airfoil, 0.05, 10, 0.1, 20)
        for planeSurface in offset.planeSurfaces:
            planeSurface.define_bc()

        surface_domain = PlaneSurface([ext_domain, offset])
    else:
        surface_domain = PlaneSurface([ext_domain, airfoil])

    # Synchronize and generate BC marker
    gmsh.model.occ.synchronize()
    ext_domain.define_bc()
    airfoil.define_bc()
    surface_domain.define_bc()
    
    # Set high order methods
    if args.high_order:
        meshOrder = args.high_order
        print(f"Polynomial degree: {meshOrder}")
        gmsh.option.setNumber("Mesh.ElementOrder", meshOrder)    
        # gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)
    else:
        meshOrder = 1

    # Mesh settings
    airfoil.setTransfinite()
    if offsetTrigger:
        offset.setTransfinite(blayer_size, blayer_thickness, meshOrder)
    
    # 3D extrusion and meshing
    if args.extrusion:
        delimiters = ['/']
        if re.match(r"\d*\.?\d*\/\d*$", args.extrusion):
            extrusionValue, nElements = re.split('|'.join(map(re.escape, delimiters)), args.extrusion)
        elif re.match(r"\d*\.?\d*", args.extrusion):
            extrusionValue = args.extrusion.strip()
            nElements = 1
        else:
            raise SyntaxError("Extrusion parameters were provided with wrong syntax. For more information use --help")
    
        extrusionValue = float(extrusionValue)
        nElements = int(nElements)

        if args.quad and divmod(nElements, 2)[1] != 0:
            raise Exception("Quad mesh uses subdivision algorithms to ensure only quads. Hence, the number of extruded elements have to be a multiple of 2!")
        elif args.quad:
            nElements = nElements/2

        if offsetTrigger:
            planeSurfaces = offset.planeSurfaces + [surface_domain]
            extrusion = MeshExtrusion(planeSurfaces, extrusionValue, nElements)
        else:
            extrusion = MeshExtrusion([surface_domain], extrusionValue, nElements)
        extrusion.define_bc()
        BoundaryCondition.generatePhysicalGroups(2,3)
        
        gmsh.model.mesh.generate(3)
    else:
        BoundaryCondition.generatePhysicalGroups(1,2)
        gmsh.model.mesh.generate(2)

    # Open user interface of GMSH
    if args.ui:
        gmsh.fltk.run()

    # Mesh file name and output
    mesh_path = Path(args.output, f"mesh_airfoil_{airfoil_name}.{args.format}")
    gmsh.write(str(mesh_path))
    gmsh.finalize()


if __name__ == "__main__":
    main()
