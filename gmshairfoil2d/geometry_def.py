"""
This script contain the definition of geometrical objects needed to build the geometry.
"""

from operator import attrgetter
import gmsh
import numpy as np
import math


class Point:
    """
    A class to represent the point geometrical object of gmsh

    ...

    Attributes
    ----------
    x : float
        position in x
    y : float
        position in y
    z : float
        position in z
    mesh_size : float
        If mesh_size is > 0, add a meshing constraint
            at that point
    """

    def __init__(self, x, y, z, mesh_size=None):

        self.x = x
        self.y = y
        self.z = z

        self.mesh_size = mesh_size
        self.dim = 0

        # create the gmsh object and store the tag of the geometric object
        if mesh_size != None:
            self.tag = gmsh.model.occ.addPoint(self.x, self.y, self.z, self.mesh_size)
        else:
            self.tag = gmsh.model.occ.addPoint(self.x, self.y, self.z)

    def rotation(self, angle, origin, axis):
        """
        Methode to rotate the object Point
        ...

        Parameters
        ----------
        angle : float
            angle of rotation in rad
        origin : tuple
            tuple of point (x,y,z) which is the origin of the rotation
        axis : tuple
            tuple of point (x,y,z) which represent the axis of rotation
        """
        gmsh.model.occ.rotate(
            [(self.dim, self.tag)],
            *origin,
            *axis,
            angle,
        )

    def translation(self, vector):
        """
        Methode to translate the object Point
        ...

        Parameters
        ----------
        direction : tuple
            tuple of point (x,y,z) which represent the direction of the translation
        """
        gmsh.model.occ.translate([(self.dim, self.tag)], *vector)


class Line:
    """
    A class to represent the Line geometrical object of gmsh

    ...

    Attributes
    ----------
    start_point : Point
        first point of the line
    end_point : Point
        second point of the line
    """

    def __init__(self, start_point, end_point):
        self.start_point = start_point
        self.end_point = end_point

        self.dim = 1

        # create the gmsh object and store the tag of the geometric object
        self.tag = gmsh.model.occ.addLine(self.start_point.tag, self.end_point.tag)

    def rotation(self, angle, origin, axis):
        """
        Methode to rotate the object Line
        ...

        Parameters
        ----------
        angle : float
            angle of rotation in rad
        origin : tuple
            tuple of point (x,y,z) which is the origin of the rotation
        axis : tuple
            tuple of point (x,y,z) which represent the axis of rotation
        """
        gmsh.model.occ.rotate(
            [(self.dim, self.tag)],
            *origin,
            *axis,
            angle,
        )

    def translation(self, vector):
        """
        Methode to translate the object Line
        ...

        Parameters
        ----------
        direction : tuple
            tuple of point (x,y,z) which represent the direction of the translation
        """
        gmsh.model.occ.translate([(self.dim, self.tag)], *vector)


class Spline:
    """
    A class to represent the Spine geometrical object of gmsh

    ...

    Attributes
    ----------
    points_list : list(Point)
        list of Point object forming the Spline
    """

    def __init__(self, point_list):
        self.point_list = point_list

        # generate the Lines tag list to folow
        self.tag_list = [point.tag for point in self.point_list]
        self.dim = 1
        # create the gmsh object and store the tag of the geometric object
        self.tag = gmsh.model.occ.addSpline(self.tag_list)

    def rotation(self, angle, origin, axis):
        """
        Methode to rotate the object Spline

        Rotate the spline itself (curve, starpoint,endpoint), then rotate the indermediate points
        ...

        Parameters
        ----------
        angle : float
            angle of rotation in rad
        origin : tuple
            tuple of point (x,y,z) which is the origin of the rotation
        axis : tuple
            tuple of point (x,y,z) which represent the axis of rotation
        """
        gmsh.model.occ.rotate(
            [(self.dim, self.tag)],
            *origin,
            *axis,
            angle,
        )

        [
            interm_point.rotation(angle, origin, axis)
            for interm_point in self.point_list[1:-1]
        ]

    def translation(self, vector):
        """
        Methode to translate the object Line

        Translate the spline itself (curve, starpoint,endpoint), then translate the indermediate points
        ...

        Parameters
        ----------
        direction : tuple
            tuple of point (x,y,z) which represent the direction of the translation
        """
        gmsh.model.occ.translate([(self.dim, self.tag)], *vector)
        [interm_point.translation(vector) for interm_point in self.point_list[1:-1]]


class CurveLoop:
    """
    A class to represent the CurveLoop geometrical object of gmsh
    Curveloop object are an addition entity of the existing line that forms it
    Curveloop must be created when the geometry is in its final layout

    ...

    Attributes
    ----------
    line_list : list(Line)
        List of Line object, in the order of the wanted CurveLoop and closed
    """

    def __init__(self, line_list):

        self.line_list = line_list
        self.dim = 1

        # generate the Lines tag list to folow
        self.tag_list = [line.tag for line in self.line_list]
        # create the gmsh object and store the tag of the geometric object
        self.tag = gmsh.model.occ.addCurveLoop(self.tag_list)


class Circle:
    """
    A class to represent a Circle geometrical object, composed of many arcCircle object of gmsh

    ...

    Attributes
    ----------
    xc : float
        position of the center in x
    yc : float
        position of the center in y
    z : float
        position in z
    radius : float
        radius of the circle
    mesh_size : float
        determine the mesh resolution and how many segment the
        resulting circle will be composed of
    """

    def __init__(self, xc, yc, zc, radius, mesh_size):
        # Position of the disk center
        self.xc = xc
        self.yc = yc
        self.zc = zc

        self.radius = radius
        self.mesh_size = mesh_size
        self.dim = 1

        # create a structured arcCricle to merge in one curveloop
        self.distribution = math.floor((np.pi * 2 * self.radius) / self.mesh_size)
        self.arcCircle_list = [
            gmsh.model.occ.addCircle(
                self.xc,
                self.yc,
                self.zc,
                self.radius,
                angle1=2 * np.pi / self.distribution * i,
                angle2=2 * np.pi / self.distribution * (1 + i),
            )
            for i in range(0, self.distribution)
        ]
        # Remove the duplicated points generated by the arcCircle
        gmsh.model.occ.synchronize()
        gmsh.model.occ.removeAllDuplicates()

    def close_loop(self):
        """
        Method to form a close loop with the current geometrical object

        Returns
        -------
        _ : int
            return the tag of the CurveLoop object
        """
        return gmsh.model.occ.addCurveLoop(self.arcCircle_list)

    def define_bc(self):
        """
        Method that define the marker of the circle
        for the boundary condition
        -------
        """

        self.bc = gmsh.model.addPhysicalGroup(self.dim, self.arcCircle_list)
        self.physical_name = gmsh.model.setPhysicalName(self.dim, self.bc, "farfield")

    def rotation(self, angle, origin, axis):
        """
        Methode to rotate the object Circle
        ...

        Parameters
        ----------
        angle : float
            angle of rotation in rad
        origin : tuple
            tuple of point (x,y,z) which is the origin of the rotation
        axis : tuple
            tuple of point (x,y,z) which represent the axis of rotation
        """
        [
            gmsh.model.occ.rotate(
                [(self.dim, arccircle)],
                *origin,
                *axis,
                angle,
            )
            for arccircle in self.arcCircle_list
        ]

    def translation(self, vector):
        """
        Methode to translate the object Circle
        ...

        Parameters
        ----------
        direction : tuple
            tuple of point (x,y,z) which represent the direction of the translation
        """
        [
            gmsh.model.occ.translate([(self.dim, arccircle)], *vector)
            for arccircle in self.arcCircle_list
        ]


class Rectangle:
    """
    A class to represent a rectangle geometrical object, composed of 4 Lines object of gmsh

    ...

    Attributes
    ----------
    xc : float
        position of the center in x
    yc : float
        position of the center in y
    z : float
        position in z
    dx: float
        length of the rectangle along the x direction
    dy: float
        length of the rectangle along the y direction
    mesh_size : float
        attribute given for the class Point
    """

    def __init__(self, xc, yc, z, dx, dy, mesh_size):

        self.xc = xc
        self.yc = yc
        self.z = z

        self.dx = dx
        self.dy = dy

        self.mesh_size = mesh_size
        self.dim = 1

        # Generate the 4 corners of the rectangle
        self.points = [
            Point(self.xc - self.dx / 2, self.yc - self.dy / 2, z, self.mesh_size),
            Point(self.xc + self.dx / 2, self.yc - self.dy / 2, z, self.mesh_size),
            Point(self.xc + self.dx / 2, self.yc + self.dy / 2, z, self.mesh_size),
            Point(self.xc - self.dx / 2, self.yc + self.dy / 2, z, self.mesh_size),
        ]

        # Generate the 4 lines of the rectangle
        self.lines = [
            Line(self.points[0], self.points[1]),
            Line(self.points[1], self.points[2]),
            Line(self.points[2], self.points[3]),
            Line(self.points[3], self.points[0]),
        ]

    def close_loop(self):
        """
        Method to form a close loop with the current geometrical object

        Returns
        -------
        _ : int
            return the tag of the CurveLoop object
        """
        return CurveLoop(self.lines).tag

    def define_bc(self):
        """
        Method that define the different markers of the rectangle for the boundary condition
        self.lines[0] => wall_bot
        self.lines[1] => outlet
        self.lines[2] => wall_top
        self.lines[3] => inlet
        -------
        """

        self.bc_in = gmsh.model.addPhysicalGroup(self.dim, [self.lines[3].tag], tag=-1)
        gmsh.model.setPhysicalName(self.dim, self.bc_in, "inlet")

        self.bc_out = gmsh.model.addPhysicalGroup(self.dim, [self.lines[1].tag])
        gmsh.model.setPhysicalName(self.dim, self.bc_out, "outlet")

        self.bc_wall = gmsh.model.addPhysicalGroup(
            self.dim, [self.lines[0].tag, self.lines[2].tag]
        )
        gmsh.model.setPhysicalName(self.dim, self.bc_wall, "wall")

        self.bc = [self.bc_in, self.bc_out, self.bc_wall]

    def rotation(self, angle, origin, axis):
        """
        Methode to rotate the object Rectangle
        ...

        Parameters
        ----------
        angle : float
            angle of rotation in rad
        origin : tuple
            tuple of point (x,y,z) which is the origin of the rotation
        axis : tuple
            tuple of point (x,y,z) which represent the axis of rotation
        """
        [line.rotation(angle, origin, axis) for line in self.lines]

    def translation(self, vector):
        """
        Methode to translate the object Rectangle
        ...

        Parameters
        ----------
        direction : tuple
            tuple of point (x,y,z) which represent the direction of the translation
        """
        [line.translation(vector) for line in self.lines]


class Airfoil:
    """
    A class to represent and airfoil as a CurveLoop object formed with lines

    ...

    Attributes
    ----------
    point_cloud : list(list(float))
        List of points forming the airfoil in the order,
        each point is a list containing in the order
        its postion x,y,z
    mesh_size : float
        attribute given for the class Point,Note that a mesh size larger
        than the resolution given by the cloud of points
        will not be taken into account
    name : str
        name of the marker that will be associated to the airfoil
        boundary condition
    """

    def __init__(self, point_cloud, mesh_size, name="airfoil"):

        self.name = name
        self.dim = 1
        # Generate Points object from the point_cloud
        self.points = [
            Point(point_cord[0], point_cord[1], point_cord[2], mesh_size)
            for point_cord in point_cloud
        ]

    def gen_skin(self):
        """
        Method to generate the line forming the foil, Only call this function when the points
        of the airfoil are in their final position
        -------
        """
        self.lines = [
            Line(self.points[i], self.points[i + 1])
            for i in range(-1, len(self.points) - 1)
        ]
        self.lines_tag = [line.tag for line in self.lines]

    def close_loop(self):
        """
        Method to form a close loop with the current geometrical object

        Returns
        -------
        _ : int
            return the tag of the CurveLoop object
        """
        return CurveLoop(self.lines).tag

    def define_bc(self):
        """
        Method that define the marker of the airfoil for the boundary condition
        -------
        """

        self.bc = gmsh.model.addPhysicalGroup(self.dim, self.lines_tag)
        gmsh.model.setPhysicalName(self.dim, self.bc, self.name)

    def rotation(self, angle, origin, axis):
        """
        Methode to rotate the object CurveLoop
        ...

        Parameters
        ----------
        angle : float
            angle of rotation in rad
        origin : tuple
            tuple of point (x,y,z) which is the origin of the rotation
        axis : tuple
            tuple of point (x,y,z) which represent the axis of rotation
        """
        [point.rotation(angle, origin, axis) for point in self.points]

    def translation(self, vector):
        """
        Methode to translate the object CurveLoop
        ...

        Parameters
        ----------
        direction : tuple
            tuple of point (x,y,z) which represent the direction of the translation
        """
        [point.translation(vector) for point in self.points]


class AirfoilSpline:
    """
    A class to represent and airfoil as a CurveLoop object formed with Splines
    ...

    Attributes
    ----------
    point_cloud : list(list(float))
        List of points forming the airfoil in the order,
        each point is a list containing in the order
        its postion x,y,z
    mesh_size : float
        attribute given for the class Point,Note that a mesh size larger
        than the resolution given by the cloud of points
        will not be taken into account
    name : str
        name of the marker that will be associated to the airfoil
        boundary condition
    """

    def __init__(self, point_cloud, mesh_size, bc_data = None, name="airfoil"):
        
        self.dim = 1
        self.mesh_size = mesh_size

        # Generate Points object from the point_cloud
        self.points = [
            Point(point_cord[0], point_cord[1], point_cord[2], mesh_size)
            for point_cord in point_cloud
        ]
        # Find leading and trailing edge location
        # in space
        self.le = min(self.points, key=attrgetter("x"))
        self.te = max(self.points, key=attrgetter("x"))

        # in the list of point
        self.le_indx = self.points.index(self.le)
        self.te_indx = self.points.index(self.te)
        
        # generate separat point sets for upper and lower side
        self.upper_splines = dict(pointSet = [self.gen_upper_pointSet()], name = [name], spline =[])
        self.lower_splines = dict(pointSet = [self.gen_lower_pointSet()], name = [name], spline =[])
        
        # split point sets
        if bc_data:            
            for i in range(len(bc_data["start"])):
                splittingPoints = self.gen_bcSplittingPoints(bc_data['start'][i], bc_data['end'][i], bc_data['side'][i])
                if bc_data["side"][i] == "SS":
                    self.upper_splines = self.split_pointSet(splittingPoints, self.upper_splines, bc_data["name"][i])
                else:
                    self.lower_splines = self.split_pointSet(splittingPoints, self.lower_splines, bc_data["name"][i])
        
        # set class variable self.all_splines containing all airfoil splines
        self.update_all_splines()

    def update_all_splines(self):
        self.all_splines = {key: [] for key in self.upper_splines.keys()}
        for key, value in self.upper_splines.items():
            self.all_splines[key].extend(value)
        for key, value in self.lower_splines.items():
            self.all_splines[key].extend(value)

    def gen_skin(self, *args):
        """
        Method to generate the two splines forming the foil, Only call this function when the points
        of the airfoil are in their final position
        -------
        """
        # Create the Splines depending on the le and te location in point_cloud
        for i in range(len(self.upper_splines["pointSet"])):
            self.upper_splines["spline"].insert(i,Spline(self.upper_splines["pointSet"][i]))
        for i in range(len(self.lower_splines["pointSet"])):
            self.lower_splines["spline"].insert(i,Spline(self.lower_splines["pointSet"][i]))
        
        self.update_all_splines()
        
        for i in range(len(self.all_splines[next(iter(self.all_splines.keys()))])):
                    print(f"{self.all_splines['name'][i]}:")
                    for point in self.all_splines['pointSet'][i]:
                        print(f"x = {point.x}, y = {point.y}, z = {point.z}")
        
        if args[0]:            
            f = gmsh.model.mesh.field.add('BoundaryLayer')
            gmsh.model.mesh.field.setNumbers(f, 'CurvesList', [spline.tag for spline in self.all_splines['spline']])
            gmsh.model.mesh.field.setNumber(f, 'Size', args[1])
            gmsh.model.mesh.field.setNumber(f, 'Ratio', args[2])
            gmsh.model.mesh.field.setNumber(f, 'Quads', 1)
            gmsh.model.mesh.field.setNumber(f, 'Thickness', args[3])
            gmsh.option.setNumber('Mesh.BoundaryLayerFanElements', 7)
            gmsh.model.mesh.field.setNumbers(f, 'FanPointsList', [self.te.tag])
            gmsh.model.mesh.field.setAsBoundaryLayer(f)

        return self.upper_splines, self.lower_splines
        # form the curvedloop

    def gen_upper_pointSet(self):
        if self.le_indx < self.te_indx:
            upper_pointSet = self.points[self.le_indx : self.te_indx + 1]
        else:
            upper_pointSet = self.points[self.le_indx :] + self.points[: (self.te_indx + 1)]
        
        upper_pointSet = sorted(upper_pointSet, key=lambda obj: obj.x)
        return upper_pointSet

    def gen_lower_pointSet(self):
        if self.le_indx < self.te_indx:
            lower_pointSet = self.points[self.te_indx :] + self.points[: (self.le_indx) + 1]
        else:
            lower_pointSet = self.points[self.te_indx : self.le_indx + 1]
        
        lower_pointSet = sorted(lower_pointSet, key=lambda obj: obj.x)
        return lower_pointSet
    
    def gen_bcSplittingPoints(self, startPoint_xCoord, endPoint_xCoord, airfoilSide):
        # Get coordinates for splitting points:
        if airfoilSide == "SS":
            pointSet = self.gen_upper_pointSet()
        elif airfoilSide == "PS":
            pointSet = self.gen_lower_pointSet()
        
        startPoint_yCoord = None
        endPoint_yCoord = None
        
        for point in pointSet:
            if pointSet.index(point) < len(pointSet) - 1:
                nextPoint = pointSet[pointSet.index(point) + 1]                
                if point.x <= startPoint_xCoord and startPoint_xCoord < nextPoint.x:
                    startPoint_yCoord = point.y + (nextPoint.y - point.y)/(nextPoint.x - point.x)*(startPoint_xCoord - point.x)
                    print(f"{point.x} <= x = {startPoint_xCoord} < {nextPoint.x}  ->  y = {startPoint_yCoord}")
                elif point.x <= endPoint_xCoord and endPoint_xCoord < nextPoint.x:
                    endPoint_yCoord = point.y + (nextPoint.y - point.y)/(nextPoint.x - point.x)*(endPoint_xCoord - point.x)
                    print(f"{point.x} <= x = {endPoint_xCoord} < {nextPoint.x}  ->  y = {endPoint_yCoord}")
               
        startPoint = Point(startPoint_xCoord, startPoint_yCoord, 0, self.mesh_size)
        endPoint = Point(endPoint_xCoord, endPoint_yCoord, 0, self.mesh_size)
            
        return startPoint, endPoint
    
    def split_pointSet(self, splittingPoints, pointSets, bc_name):
        pointSet_to_split = None
        index_pointSet_to_split = None
        
        for pointSet in pointSets["pointSet"]:
            if pointSet[0].x < splittingPoints[0].x and splittingPoints[1].x < pointSet[-1].x:
                pointSet_to_split = pointSet
                index_pointSet_to_split = pointSets["pointSet"].index(pointSet)
                break
        
        if pointSet_to_split == None or index_pointSet_to_split == None:
            raise ValueError(f"Error: boundary conditions are overlapping  ->  {bc_name} can not be created!")
            return pointSets
                
        pointSets["pointSet"][index_pointSet_to_split:index_pointSet_to_split + 1] = ([[splittingPoints[0]],list(splittingPoints),[splittingPoints[1]]])
        pointSets["name"][index_pointSet_to_split:index_pointSet_to_split + 1] = (["airfoil", bc_name, "airfoil"])

        for point in pointSet:
            if point.x < splittingPoints[0].x:
                pointSets["pointSet"][index_pointSet_to_split].append(point)
            elif splittingPoints[0].x < point.x and point.x < splittingPoints[1].x:
                pointSets["pointSet"][index_pointSet_to_split + 1].append(point)
            else:
                pointSets["pointSet"][index_pointSet_to_split + 2].append(point)
        
        for i in range(len(pointSets["pointSet"])):
            pointSets["pointSet"][i] = sorted(pointSets["pointSet"][i], key=lambda obj: obj.x)
        
        return pointSets
    
    def close_loop(self):
        """
        Method to form a close loop with the current geometrical object

        Returns
        -------
        _ : int
            return the tag of the CurveLoop object
        """                
        return CurveLoop(self.all_splines["spline"]).tag

    def define_bc(self):
        """
        Method that define the marker of the airfoil for the boundary condition
        -------
        """
        # Find all different boundaries
        
        self.bcs = {"name": [], "spline tags": [], "physical group": []}
        for i in range(len(self.all_splines[next(iter(self.all_splines.keys()))])):
            if self.all_splines["name"][i] in self.bcs["name"]:
                self.bcs["spline tags"][self.bcs["name"].index(self.all_splines["name"][i])].append(self.all_splines["spline"][i].tag)
            else:
                self.bcs["name"].append(self.all_splines["name"][i])
                self.bcs["spline tags"].append([self.all_splines["spline"][i].tag])
                
        for i in range(len(self.bcs[next(iter(self.bcs.keys()))])):
            self.bcs["physical group"].append(gmsh.model.addPhysicalGroup(self.dim, self.bcs["spline tags"][i]))
            gmsh.model.setPhysicalName(self.dim, self.bcs["physical group"][i], self.bcs["name"][i])

    def rotation(self, angle, origin, axis):
        """
        Methode to rotate the object AirfoilSpline
        ...

        Parameters
        ----------
        angle : float
            angle of rotation in rad
        origin : tuple
            tuple of point (x,y,z) which is the origin of the rotation
        axis : tuple
            tuple of point (x,y,z) which represent the axis of rotation
        """       
        rot_pointSet = []
        for pointSet in self.all_splines["pointSet"]:
            for point in pointSet:
                if point not in rot_pointSet:
                    rot_pointSet.append(point)
                    
        for point in rot_pointSet:
            print(f"x = {point.x}, y = {point.y}, z = {point.z}")
                    
        [point.rotation(angle, origin, axis) for point in rot_pointSet]

    def translation(self, vector):
        """
        Methode to translate the object AirfoilSpline
        ...

        Parameters
        ----------
        direction : tuple
            tuple of point (x,y,z) which represent the direction of the translation
        """
        trans_pointSet = []
        for pointSet in self.all_splines["pointSet"]:
            for point in pointSet:
                if point not in trans_pointSet:
                    trans_pointSet.append(point)
                    
        for point in trans_pointSet:
            print(f"x = {point.x}, y = {point.y}, z = {point.z}")
                    
        [point.translation(vector) for point in trans_pointSet]


class PlaneSurface:
    """
    A class to represent the PlaneSurface geometrical object of gmsh


    ...

    Attributes
    ----------
    geom_objects : list(geom_object)
        List of geometrical object able to form closedloop,
        First the object will be closed in ClosedLoop
        the first curve loop defines the exterior contour; additional curve loop
        define holes in the surface domaine

    """

    def __init__(self, geom_objects):

        self.geom_objects = geom_objects
        # close_loop() will form a close loop object and return its tag
        self.tag_list = [geom_object.close_loop() for geom_object in self.geom_objects]

        self.dim = 2

        # create the gmsh object and store the tag of the geometric object
        self.tag = gmsh.model.occ.addPlaneSurface(self.tag_list)

    def define_bc(self):
        """
        Method that define the domain marker of the surface
        -------
        """
        self.ps = gmsh.model.addPhysicalGroup(self.dim, [self.tag])
        gmsh.model.setPhysicalName(self.dim, self.ps, "fluid")


class MeshExtrusion:
    """
    A class to represent the PlaneSurface geometrical object of gmsh


    ...

    Attributes
    ----------
    geom_objects : list(geom_object)
        List of geometrical object able to form closedloop,
        First the object will be closed in ClosedLoop
        the first curve loop defines the exterior contour; additional curve loop
        define holes in the surface domaine

    """

    def __init__(self, plane_surface, extrusion_value):
        
        self.dim = 3
        self.extrusion_value = extrusion_value
        self.plane_surface = plane_surface

        ov = gmsh.model.occ.extrude([(plane_surface.dim,plane_surface.tag)], 0, 0, extrusion_value, [1],recombine = True)
        gmsh.model.occ.synchronize()
        
        self.bc_volume = dict(name = "fluid", tags = [ov[1][1]], physical_group = None)
        self.bc_surface = dict(name = ["side_z-","side_z+"], tags = [[plane_surface.tag],[ov[0][1]]], physical_group = [])
        
        for bc in gmsh.model.getPhysicalGroups(1):
            bc_name = gmsh.model.getPhysicalName(bc[0],bc[1])
            bc_entities = gmsh.model.getEntitiesForPhysicalGroup(bc[0],bc[1])
            
            if bc_name not in self.bc_surface['name']:
                self.bc_surface['name'].append(bc_name)
                self.bc_surface['tags'].append([])
                for extruded_surface_tag in ov[2:len(ov)]:
                    if gmsh.model.getBoundary([extruded_surface_tag], combined=True, oriented=False, recursive=False)[0][1] in bc_entities:
                        self.bc_surface['tags'][self.bc_surface['name'].index(bc_name)].append(extruded_surface_tag[1])

    def define_bc(self):
        """
        Method that define the domain marker of the surface
        -------
        """
        gmsh.model.removePhysicalGroups([(self.plane_surface.dim,self.plane_surface.ps)])
        
        self.bc_volume['physical_group'] = gmsh.model.addPhysicalGroup(self.dim, self.bc_volume['tags'])        
        gmsh.model.setPhysicalName(self.dim, self.bc_volume['physical_group'], self.bc_volume['name'])
        
        for i in range(len(self.bc_surface['name'])):
            self.bc_surface['physical_group'].append(gmsh.model.addPhysicalGroup(self.plane_surface.dim, self.bc_surface['tags'][i]))
            gmsh.model.setPhysicalName(self.plane_surface.dim, self.bc_surface['physical_group'][i], self.bc_surface['name'][i])