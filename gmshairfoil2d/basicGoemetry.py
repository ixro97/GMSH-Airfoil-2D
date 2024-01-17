from operator import attrgetter
import gmsh
import numpy as np
import math
from gmshairfoil2d.Utils import InstanceTracker

class Point(metaclass=InstanceTracker):
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

    def __init__(self, x, y, z, tag=None):

        self.x = x
        self.y = y
        self.z = z

        self.dim = 0

        self.tag = tag

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

        # Define needed vectors
        posVector = np.array([self.x, self.y, self.z])
        axisVector = np.array([axis[0], axis[1], axis[2]])
        axisVector_normalized = axisVector/np.linalg.norm(axisVector)
        originVector = np.array([origin[0], origin[1], origin[2]])

        # Perform the rotation using the rotation matrix (Rodrigues' rotation formula)
        posVector_rotated = originVector + ((posVector - originVector)*math.cos(angle)                                     \
                            + np.cross(axisVector_normalized, posVector - originVector)*math.sin(angle)                    \
                            + axisVector_normalized*np.dot(axisVector_normalized, (posVector - originVector))*(1 - math.cos(angle)))
        
        # Assign the rotated vector components to the point.<coordinates>
        self.x = round(posVector_rotated[0], 12)
        self.y = round(posVector_rotated[1], 12)
        self.z = round(posVector_rotated[2], 12)

        # Rotate the point in the gmsh model
        if self.tag != None:
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

    def generate(self):
        # create the gmsh object and store the tag of the geometric object
        if self.tag == None:
            self.tag = gmsh.model.occ.addPoint(self.x, self.y, self.z)

    def remove(self):
        """
        Methode to translate the object Point
        ...

        Parameters
        ----------
        direction : tuple
            tuple of point (x,y,z) which represent the direction of the translation
        """
        self.tag = None
        gmsh.model.occ.synchronize()
        gmsh.model.occ.remove([(self.dim, self.tag)])
        gmsh.model.occ.synchronize()

    @staticmethod
    def interpolate(x, lastpoint, nextPoint, type):
        if type == "linear2D":
            y = lastpoint.y + (nextPoint.y - lastpoint.y)/(nextPoint.x - lastpoint.x)*(x - lastpoint.x)
            print(f"{lastpoint.x} < x = {x} < {nextPoint.x}  ->  y = {y}")
            resultingPoint = Point(x, y, 0)
        
        return resultingPoint

    @staticmethod
    def getInstance(x, y, z):
        for point in Point.instances:
            if point.x == x and point.y == y and point.z == z:
                return point
        else:
            return None

class Line(metaclass=InstanceTracker):
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
        self.points = [start_point, end_point]

        self.dim = 1
        self.length = math.sqrt((self.points[1].x - self.points[0].x)**2 + (self.points[1].y - self.points[0].y)**2 + (self.points[1].z - self.points[0].z)**2)

        # create the gmsh object and store the tag of the geometric object
        self.tag = None

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
        if self.tag != None:
            gmsh.model.occ.rotate(
                [(self.dim, self.tag)],
                *origin,
                *axis,
                angle,
            )
        else:
            [point.rotation(angle, origin, axis) for point in self.points]

    def translation(self, vector):
        """
        Methode to translate the object Line
        ...

        Parameters
        ----------
        direction : tuple
            tuple of point (x,y,z) which represent the direction of the translation
        """
        if self.tag != None:
            gmsh.model.occ.translate([(self.dim, self.tag)], *vector)
        else:
            [point.translation(vector) for point in self.points]
    
    def generate(self):
        self.tag = gmsh.model.occ.addLine(self.points[0].tag, self.points[1].tag)

    def setTransfinite(self, startMeshSize = None, endMeshSize = None, nPnts = None, growthRate = None):
        if startMeshSize != None and endMeshSize != None:
            growthRate = (self.length - startMeshSize)/(self.length - endMeshSize)

            self.startMeshSize = startMeshSize
            self.endMeshSize = endMeshSize

            if growthRate == 1:
                self.nTransfinitePoints = round(self.length/startMeshSize - 1)
            else:
                self.nTransfinitePoints = round(math.log(endMeshSize/startMeshSize, growthRate) + 1)

        elif nPnts != None and growthRate != None:
            if growthRate == 1:
                self.startMeshSize = self.length / (nPnts + 1)
            else:
                self.startMeshSize = self.length * (growthRate - 1)/(growthRate**(nPnts + 1) - 1)
            self.endMeshSize = self.startMeshSize * growthRate**(nPnts)
            self.nTransfinitePoints = nPnts

        elif endMeshSize != None and nPnts != None:
            # Estimator for progression factor:
            estimatedLength = 0
            growthRate = 1
            toleranceLen = 0.01
            direction = 0
            exponent = 0
            while math.fabs(self.length - estimatedLength) > self.length*toleranceLen:
                if growthRate == 1:
                    estimatedLength = endMeshSize * (nPnts + 1)
                else:                
                    estimatedLength = endMeshSize * (growthRate**(nPnts + 1) - 1)/(growthRate**(nPnts + 1) - growthRate**(nPnts))

                if self.length < estimatedLength and direction != 1:
                    exponent += 1
                    direction = 1
                elif self.length > estimatedLength and direction != -1:
                    exponent += 1
                    direction = -1
                elif self.length == estimatedLength: break
                
                growthRate += direction*10**(-exponent)

            self.startMeshSize = endMeshSize/growthRate**nPnts
            self.endMeshSize = endMeshSize
            self.nTransfinitePoints = nPnts
            
        elif startMeshSize != None and nPnts != None:
            # Estimator for progression factor:
            estimatedLength = 0
            growthRate = 1
            toleranceLen = 0.01
            direction = 0
            exponent = 0
            while math.fabs(self.length - estimatedLength) > self.length*toleranceLen:
                if growthRate == 1:
                    estimatedLength = startMeshSize * (nPnts + 1)
                else:                
                    estimatedLength = startMeshSize * (growthRate**(nPnts + 1) - 1)/(growthRate - 1)

                if self.length > estimatedLength and direction != 1:
                    exponent += 1
                    direction = 1
                elif self.length < estimatedLength and direction != -1:
                    exponent += 1
                    direction = -1
                elif self.length == estimatedLength: break
                
                growthRate += direction*10**(-exponent)

            self.startMeshSize = startMeshSize
            self.endMeshSize = startMeshSize*growthRate**nPnts
            self.nTransfinitePoints = nPnts

        gmsh.model.mesh.setTransfiniteCurve(self.tag, self.nTransfinitePoints, "Progression", growthRate)

class Spline(metaclass=InstanceTracker):
    """
    A class to represent the Spine geometrical object of gmsh

    ...

    Attributes
    ----------
    points_list : list(Point)
        list of Point object forming the Spline
    """

    def __init__(self, points):
        self.points = points
        self.length = 0

        # calcutlate spline length 
        for i in range(len(self.points) - 1):
            currentPoint = self.points[i]
            nextPoint = self.points[i + 1]

            self.length += math.sqrt((nextPoint.x - currentPoint.x)**2 + (nextPoint.y - currentPoint.y)**2 + (nextPoint.z - currentPoint.z)**2)

        # generate the Lines tag list to follow
        self.dim = 1
        # create the gmsh object and store the tag of the geometric object
        self.tag = None

        self.nTransfinitePoints = None
        self.startMeshSize = None
        self.endMeshSize = None

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
        if self.tag != None:
            [
                point.rotation(angle, origin, axis)
                for point in self.points[1:-1]
            ]
            gmsh.model.occ.rotate(
                [(self.dim, self.tag)],
                *origin,
                *axis,
                angle,
            )
        else:
            [point.rotation(angle, origin, axis) for point in self.points]
        
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
        [interm_point.translation(vector) for interm_point in self.points[1:-1]]

    def generate(self):
        for point in self.points:
            point.generate()
        if self.tag == None:
            self.tag = gmsh.model.occ.addSpline([point.tag for point in self.points])

    def setTransfinite(self, startMeshSize = None, endMeshSize = None, nPnts = None, growthRate = None, meshType = "Progression"):
        if startMeshSize != None and endMeshSize != None:
            growthRate = (self.length - startMeshSize)/(self.length - endMeshSize)

            self.startMeshSize = startMeshSize
            self.endMeshSize = endMeshSize
            
            print(self.length,self.startMeshSize,self.endMeshSize,growthRate)
            if growthRate == 1:
                self.nTransfinitePoints = round(self.length/startMeshSize + 1)
            else:
                self.nTransfinitePoints = round(math.log(endMeshSize/startMeshSize, growthRate) + 1) + 1
            self.growthRate = growthRate

        elif nPnts != None and growthRate != None:
            if growthRate == 1:
                self.startMeshSize = self.length / (nPnts + 1)
            else:
                self.startMeshSize = self.length * (growthRate - 1)/(growthRate**(nPnts + 1) - 1)
            self.endMeshSize = self.startMeshSize * growthRate**(nPnts)
            self.nTransfinitePoints = nPnts
            self.growthRate = growthRate

        elif startMeshSize != None and growthRate != None:
            self.startMeshSize = startMeshSize
            self.endMeshSize = self.length - (self.length - self.startMeshSize)/growthRate
            self.nTransfinitePoints = round(math.log(self.endMeshSize/self.startMeshSize, growthRate) + 1) + 1
            self.growthRate = growthRate

        elif endMeshSize != None and nPnts != None:
            # Estimator for progression factor:
            estimatedLength = 0
            growthRate = 1
            toleranceLen = 0.001
            direction = 0
            exponent = 0
            while math.fabs(self.length - estimatedLength) > self.length*toleranceLen:
                if growthRate == 1:
                    estimatedLength = endMeshSize * (nPnts + 1)
                else:                
                    estimatedLength = endMeshSize * (growthRate**(nPnts + 1) - 1)/(growthRate**(nPnts + 1) - growthRate**(nPnts))

                if self.length < estimatedLength and direction != 1:
                    exponent += 1
                    direction = 1
                elif self.length > estimatedLength and direction != -1:
                    exponent += 1
                    direction = -1
                elif self.length == estimatedLength: break
                
                growthRate += direction*10**(-exponent)

            self.startMeshSize = endMeshSize/growthRate**nPnts
            self.endMeshSize = endMeshSize
            self.nTransfinitePoints = nPnts
            self.growthRate = growthRate

        elif startMeshSize != None and nPnts != None:
            # Estimator for progression factor:
            estimatedLength = 0
            growthRate = 1
            toleranceLen = 0.001
            direction = 0
            exponent = 0
            while math.fabs(self.length - estimatedLength) > self.length*toleranceLen:
                if growthRate == 1:
                    estimatedLength = startMeshSize * (nPnts + 1)
                else:                
                    estimatedLength = startMeshSize * (growthRate**(nPnts + 1) - 1)/(growthRate - 1)

                if self.length > estimatedLength and direction != 1:
                    exponent += 1
                    direction = 1
                elif self.length < estimatedLength and direction != -1:
                    exponent += 1
                    direction = -1
                elif self.length == estimatedLength: break
                
                growthRate += direction*10**(-exponent)

            self.startMeshSize = startMeshSize
            self.endMeshSize = startMeshSize*growthRate**nPnts
            self.nTransfinitePoints = nPnts
            self.growthRate = growthRate

        gmsh.model.mesh.setTransfiniteCurve(self.tag, self.nTransfinitePoints, meshType, growthRate)

    @staticmethod
    def getInstance(points):
        for spline in Spline.instances:
            if spline.points == points:
                return spline
        else:
            return None

class AirfoilSpline(Spline):
    def __init__(self, points, side):
        super().__init__(points)    

        self.offsetSpline = None
        self.side = side
        self.controlSplines = None
        self.bc = None
        
class ExtensionSpline(Spline):
    def __init__(self, points):
        super().__init__(points)    

        self.offsetSplines = None
        self.controlSplines = None

class CircleArc:
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

    def __init__(self, xc, yc, zc, radius, startAngle, endAngle):
        # Position of the disk center
        self.xc = xc
        self.yc = yc
        self.zc = zc

        self.radius = radius
        self.length = (endAngle - startAngle)*radius
        
        self.startAngle = startAngle
        self.endAngle = endAngle
        
        # Calculate points
        self.center = Point.getInstance(self.xc, self.yc, self.zc)
        if self.center is None:
            self.center = Point(self.xc, self.yc, self.zc)
            
        self.start = Point.getInstance(self.xc + self.radius*math.cos(startAngle), self.yc + self.radius*math.sin(startAngle), self.zc)
        if self.start is None:
            self.start = Point(self.xc + self.radius*math.cos(startAngle), self.yc + self.radius*math.sin(startAngle), self.zc)
        
        self.end = Point.getInstance(self.xc + self.radius*math.cos(endAngle), self.yc + self.radius*math.sin(endAngle), self.zc)
        if self.end is None:
            self.end = Point(self.xc + self.radius*math.cos(endAngle), self.yc + self.radius*math.sin(endAngle), self.zc)    
        
        self.dim = 1
        self.tag = None
    
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
        gmsh.model.occ.rotate([(self.dim, self.dim)], *origin, *axis, angle)
        
    def translation(self, vector):
        """
        Methode to translate the object Circle
        ...

        Parameters
        ----------
        direction : tuple
            tuple of point (x,y,z) which represent the direction of the translation
        """

        gmsh.model.occ.translate([(self.dim, self.tag)], *vector)

    def generate(self):
        if self.endAngle - self.startAngle < math.pi:
            self.start.generate()
            self.center.generate()
            self.end.generate()
            self.tag = gmsh.model.occ.addCircleArc(self.start.tag, self.center.tag, self.end.tag)
        else:
            self.tag = gmsh.model.occ.addCircle(self.xc, self.yc, self.zc, self.radius, angle1=self.startAngle, angle2=self.endAngle)
            gmsh.model.occ.synchronize()
            if len(gmsh.model.getBoundary([(self.dim, self.tag)])) != 0:
                self.start.tag = gmsh.model.getBoundary([(self.dim, self.tag)])[0][1]
                self.end.tag = gmsh.model.getBoundary([(self.dim, self.tag)])[1][1]
            
    def setTransfinite(self, meshSize):
        nPnts = math.floor(self.length / meshSize)
        
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.setTransfiniteCurve(self.tag, nPnts)

class Circle:
    def __init__(self, xc, yc, zc, radius):
        self.xc = xc
        self.yc = yc
        self.zc = zc

        self.radius = radius
        
        self.arcList = [CircleArc(xc, yc, zc, radius, startAngle=0., endAngle=2*math.pi)]
    
    def split(self, angle1, angle2):
        #Check if the Circle is already splitted
        if len(self.arcList) != 1:
            raise SyntaxError("Circle is already splitted")
        
        for arc in self.arcList:
            if arc.startAngle <= angle1 and angle1 <= arc.endAngle and arc.startAngle <= angle2 and angle2 <= arc.endAngle:
                # Resize base circle arc to the correct angles
                if angle1 < angle2:
                    arc.startAngle = angle1
                    arc.endAngle = angle2
                elif angle1 > angle2:
                    arc.startAngle = angle2
                    arc.endAngle = angle1
                
                arc.start.x = arc.xc + arc.radius*math.cos(arc.startAngle)
                arc.start.y = arc.yc + arc.radius*math.sin(arc.startAngle)
                if Point.getInstance(arc.start.x, arc.start.y, arc.start.z):
                    arc.start = Point.getInstance(arc.start.x, arc.start.y, arc.start.z)
                
                arc.end.x = arc.xc + arc.radius*math.cos(arc.endAngle)
                arc.end.y = arc.yc + arc.radius*math.sin(arc.endAngle)
                if Point.getInstance(arc.end.x, arc.end.y, arc.end.z):
                    arc.start = Point.getInstance(arc.end.x, arc.end.y, arc.end.z)
                
                # Create arc
                self.arcList.append(CircleArc(self.xc, self.yc, self.zc, self.radius, arc.endAngle, arc.startAngle))
        
    def generate(self):
        for arc in self.arcList:
            arc.generate()

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

        self.bc_in = BoundaryCondition(self.dim, "inlet", [self.lines[3].tag])
        self.bc_out = BoundaryCondition(self.dim, "outlet", [self.lines[1].tag])
        self.bc_wall = BoundaryCondition(self.dim, "wall", [self.lines[0].tag, self.lines[2].tag])

        self.bcs = [self.bc_in, self.bc_out, self.bc_wall]

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