"""
This script contain the definition of geometrical objects needed to build the geometry.
"""

from operator import attrgetter
import gmsh
import numpy as np
import math
from gmshairfoil2d.Utils import InstanceTracker
from gmshairfoil2d.basicGoemetry import (Point, Line, Spline, AirfoilSpline, 
                                         ExtensionSpline, Circle, CircleArc,
                                         Rectangle)

class ExternalDomain:
    def __init__(self) -> None:
        pass
    def close_loop(self):
        """
        Method to form a close loop with the current geometrical object

        Returns
        -------
        _ : int
            return the tag of the CurveLoop object
        """
        return gmsh.model.occ.addCurveLoop([self.tag])

    def define_bc(self):
        """
        Method that define the marker of the circle
        for the boundary condition
        -------
        """

        self.bc = BoundaryCondition(self.dim, "farfield", [self.tag])


class SplineInterpolation:
    @staticmethod 
    def generateEndPoint(startPoint, tangVector, transitionLength, transitionOffset = None):
        # Define or calculate the position vectors
        if(transitionOffset is None):
            if tangVector[0]/tangVector[1] > 0:
                endPoint_x = startPoint.x + transitionLength
                endPoint_y = startPoint.y + transitionLength*(-tangVector[0]/tangVector[1] + np.sqrt(1 + (tangVector[0]/tangVector[1])**2))
                endPoint_z = 0
            elif tangVector[0]/tangVector[1] < 0:
                endPoint_x = startPoint.x + transitionLength
                endPoint_y = startPoint.y + transitionLength*(-tangVector[0]/tangVector[1] - np.sqrt(1 + (tangVector[0]/tangVector[1])**2))
                endPoint_z = 0
        else:
            endPoint_x = startPoint.x + transitionLength
            endPoint_y = startPoint.y + transitionOffset
            endPoint_z = 0

        endPoint = Point(endPoint_x, endPoint_y, endPoint_z)

        print(f"Created end point for spline transition: x = {endPoint.x}, y = {endPoint.y}, z = {endPoint.z}")

        return endPoint
    
    @staticmethod
    def generateControlPoint(startPoint, tangVector, endPoint):
        controlPoint_x = startPoint.x + tangVector[0]/tangVector[1]*(endPoint.y - startPoint.y)
        controlPoint_y = endPoint.y
        controlPoint_z = 0

        controlPoint = Point(controlPoint_x, controlPoint_y, controlPoint_z)

        print(f"Created control point for spline transition: x = {controlPoint.x}, y = {controlPoint.y}, z = {controlPoint.z}")

        return controlPoint
 
    @staticmethod
    def quadBezierInterpolation(startPoint, controlPoint, endPoint):
        startVector = np.array([startPoint.x, startPoint.y, startPoint.z])
        controlVector = np.array([controlPoint.x, controlPoint.y, controlPoint.z])
        endVector = np.array([endPoint.x, endPoint.y, endPoint.z])

        interpolationPoints = [startPoint]

        # Performing spline interpolation on <n_interpolation> points
        length = np.linalg.norm(endVector - startVector)
        n_interpolation = round(20*length)
        
        for n in range(n_interpolation - 1):
            # Kurvenparameter
            t = (n+1)*1/n_interpolation
            
            # Quadratic Bezier
            interpolationPoint = (1 - t)**2*startVector + 2*(1 - t)*t*controlVector + t**2*endVector
            interpolationPoints.append(Point(interpolationPoint[0], interpolationPoint[1], interpolationPoint[2]))
        
        interpolationPoints.append(endPoint)

        return interpolationPoints

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

class Airfoil:
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

    def __init__(self, point_cloud, bcs):
        
        self.dim = 1
        self.bcs = bcs
        self.refinementLengthLE = 0.2

        # Generate Points object from the point_cloud
        self.points = [
            Point(point_cord[0], point_cord[1], point_cord[2])
            for point_cord in point_cloud
        ]

        self.pointLE = min(self.points, key=attrgetter("x"))
        self.pointTE = max(self.points, key=attrgetter("x"))

        self.idxLE = self.points.index(self.pointLE)
        self.idxTE = self.points.index(self.pointTE)
        
        self.splines = []
        self.generateSplines()

        self.offsetPoints = None
        self.offsetSplines = None

        self.extensionSplines = None

    def addBoundaryPoint(self, x, side):
        resultingPoint = None

        # Get points for given side
        if side == "SS":
            points = self.points[self.idxLE : self.idxTE + 1]
            points = sorted(points, key=lambda obj: obj.x)
        elif side == "PS":
            points = self.points[self.idxTE :] + self.points[: (self.idxLE) + 1]
            points = sorted(points, key=lambda obj: obj.x)
        else:
            raise Exception(f"Side {side} is not known")
        
        for point in points:
            if point.x == x:
                resultingPoint = point
            elif points.index(point) < len(points) - 1:
                nextPoint = points[points.index(point) + 1]
                if point.x < x and x < nextPoint.x:         
                    resultingPoint = Point.interpolate(x, point, nextPoint, "linear2D")

        if resultingPoint not in self.points and side == "SS":
            self.points.insert(self.idxLE, resultingPoint)
            self.points[self.idxLE:self.idxTE + 1] = sorted(self.points[self.idxLE:self.idxTE + 1], key=lambda obj: obj.x)
        elif resultingPoint not in self.points and side == "PS":
            self.points.append(resultingPoint)
            self.points[self.idxTE:] = sorted(self.points[self.idxTE:], key=lambda obj: obj.x, reverse=True)

        self.idxLE = self.points.index(self.pointLE)
        self.idxTE = self.points.index(self.pointTE)

        return resultingPoint                                   

    def generateSplines(self):
        # Add splitting points
        splittingPoints = [self.pointLE, self.pointTE]

        for bc in self.bcs:
            startPoint = self.addBoundaryPoint(bc.start, bc.side)
            endPoint = self.addBoundaryPoint(bc.end, bc.side)

            if startPoint not in splittingPoints and bc.side == "SS":
                splittingPoints.insert(0, startPoint)
            elif startPoint not in splittingPoints and bc.side == "PS":
                splittingPoints.append(startPoint)
            if endPoint not in splittingPoints and bc.side == "SS":
                splittingPoints.insert(0, endPoint)
            elif endPoint not in splittingPoints and bc.side == "PS":
                splittingPoints.append(endPoint)
        
        if self.refinementLengthLE != None:
            ssPoint = self.addBoundaryPoint(self.refinementLengthLE, "SS")
            psPoint = self.addBoundaryPoint(self.refinementLengthLE, "PS")

            if ssPoint not in splittingPoints:
                splittingPoints.insert(0, ssPoint)
            if psPoint not in splittingPoints:
                splittingPoints.append(psPoint)

        idxTE = splittingPoints.index(self.pointTE)
        splittingPoints[0:idxTE + 1] = sorted(splittingPoints[0:idxTE + 1], key=lambda obj: obj.x)
        splittingPoints[idxTE:] = sorted(splittingPoints[idxTE:], key=lambda obj: -obj.x)

        for point in splittingPoints:
            print(self.points.index(point))

        for idxPoint, point in enumerate(splittingPoints):
            startPoint = point
            idxStartPoint = self.points.index(startPoint)
            if idxPoint == len(splittingPoints) - 1:
                endPoint = splittingPoints[0]
                idxEndPoint = self.points.index(endPoint)
                print(f"\n{idxStartPoint} / {idxEndPoint}")

                splinePoints = self.points[idxStartPoint:]
                splinePoints.append(self.points[0])
                for point in splinePoints:
                    print(f"x = {point.x}, y = {point.y}, z = {point.z} -> Tag: {point.tag}")
                side = "PS"
                splinePoints.reverse()
                self.splines.append(AirfoilSpline(splinePoints, "PS"))
            else:
                endPoint = splittingPoints[idxPoint + 1]
                idxEndPoint = self.points.index(endPoint)
                
                print(f"\n{idxStartPoint} / {idxEndPoint} -> {self.idxTE}")

                splinePoints = self.points[idxStartPoint:idxEndPoint + 1]
                
                for point in splinePoints:
                    print(f"x = {point.x}, y = {point.y}, z = {point.z} -> Tag: {point.tag}")

                if idxEndPoint <= self.idxTE:
                    side = "SS"
                    self.splines.append(AirfoilSpline(splinePoints, side))
                else:
                    side = "PS"
                    splinePoints.reverse()
                    self.splines.append(AirfoilSpline(splinePoints, side))

        airfoilBC = AirfoilBoundaryCondition("airfoil", None, None, None)
        for spline in self.splines:
            for bc in self.bcs:
                if spline.points[0].x == bc.start and spline.points[-1].x == bc.end:
                    print(f"{spline.points[0].x} == {bc.start} and {spline.points[-1].x} == {bc.end}")
                    spline.bc = bc
                    break
            else:
                spline.bc = airfoilBC

        for spline in self.splines:
            print(spline.bc.name)
    
    def close_loop(self):
        """
        Method to form a close loop with the current geometrical object

        Returns
        -------
        _ : int
            return the tag of the CurveLoop object
        """                
        return CurveLoop(self.splines).tag

    def define_bc(self):
        """
        Method that define the marker of the airfoil for the boundary condition
        -------
        """
        # Find all different boundaries
        for spline in self.splines:
            self.bcs = []
            bc_name = spline.bc.name
            bc_entities = [spline.tag]

            for bc in BoundaryCondition.instances:
                if bc.name == bc_name and bc.dim == self.dim:
                    self.bcs.append(bc)
                    bc.tag_list.extend(bc_entities)
                    break
            else:
                self.bcs.append(BoundaryCondition(self.dim, bc_name, bc_entities))

    def setBoundaryLayer(self, *args):
        f = gmsh.model.mesh.field.add('BoundaryLayer')
        gmsh.model.mesh.field.setNumbers(f, 'CurvesList', [spline.tag for spline in self.splines])
        gmsh.model.mesh.field.setNumber(f, 'Size', args[0])
        gmsh.model.mesh.field.setNumber(f, 'Ratio', args[1])
        gmsh.model.mesh.field.setNumber(f, 'Quads', 1)
        gmsh.model.mesh.field.setNumber(f, 'Thickness', args[2])
        gmsh.option.setNumber('Mesh.BoundaryLayerFanElements', 7)
        gmsh.model.mesh.field.setNumbers(f, 'FanPointsList', [self.pointTE.tag])
        gmsh.model.mesh.field.setAsBoundaryLayer(f)

    def rotation(self, angle, origin, axis):
        """
        Methode to rotate the object Airfoil
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
        self.aoa = angle
        # [spline.rotation(angle, origin, axis) for spline in self.splines]
        [point.rotation(angle, origin, axis) for point in self.points]

    def translation(self, vector):
        """
        Methode to translate the object Airfoil
        ...

        Parameters
        ----------
        direction : tuple
            tuple of point (x,y,z) which represent the direction of the translation
        """
        [spline.translation(vector) for spline in self.splines]     

    def offset(self, offsetValue):
        self.offsetValue = offsetValue
        self.offsetPoints = []
        self.offsetSplines = []

        vector_z = np.array([0,0,1])
        
        for idxPoint, point in enumerate(self.points):
            if point == self.pointLE:
                lastPoint = self.points[-1]
                nextPoint = self.points[idxPoint + 1]
            elif idxPoint == len(self.points) - 1:
                lastPoint = self.points[idxPoint - 1]
                nextPoint = self.points[0]
            else:
                lastPoint = self.points[idxPoint - 1]
                nextPoint = self.points[idxPoint + 1]

            vector_lastPoint = np.array([lastPoint.x - point.x, lastPoint.y - point.y, 0])
            vector_lastPoint_normalized = vector_lastPoint/np.linalg.norm(vector_lastPoint)
            vector_normal_lastPoint = np.cross(vector_lastPoint_normalized, vector_z)

            vector_nextPoint = np.array([nextPoint.x - point.x, nextPoint.y - point.y, 0])
            vector_nextPoint_normalized = vector_nextPoint/np.linalg.norm(vector_nextPoint)
            vector_normal_nextPoint = np.cross(vector_z, vector_nextPoint_normalized)

            vector_normal = (vector_normal_lastPoint + vector_normal_nextPoint)/np.linalg.norm(vector_normal_lastPoint + vector_normal_nextPoint)

            if point == self.pointTE:
                vector_offsetPoint = np.array([point.x, point.y, point.z]) + offsetValue*vector_normal_lastPoint
                offsetPoint = Point(vector_offsetPoint[0], vector_offsetPoint[1], 0)
                self.offsetPoints.append(offsetPoint)

                vector_offsetPoint = np.array([point.x, point.y, point.z]) + offsetValue*vector_normal_nextPoint
                offsetPoint = Point(vector_offsetPoint[0], vector_offsetPoint[1], 0)
                self.offsetPoints.append(offsetPoint)
            else:
                vector_normal = (vector_normal_lastPoint + vector_normal_nextPoint)/np.linalg.norm(vector_normal_lastPoint + vector_normal_nextPoint)
                vector_offsetPoint = np.array([point.x, point.y, point.z]) + offsetValue*vector_normal
                offsetPoint = Point(vector_offsetPoint[0], vector_offsetPoint[1], 0)
                self.offsetPoints.append(offsetPoint)

        idxEndPoint = self.points.index(self.pointLE)
        for idxSpline, spline in enumerate(self.splines):
            if self.pointTE in spline.points and spline.side == "PS":
                idxStartPoint = idxEndPoint + 1
            else:
                idxStartPoint = idxEndPoint
            idxEndPoint = idxStartPoint + len(spline.points) - 1
        
            if idxSpline == len(self.splines) - 1:
                points = self.offsetPoints[idxStartPoint:idxEndPoint + 1] + [self.offsetPoints[self.points.index(self.pointLE)]]
            else:
                points = self.offsetPoints[idxStartPoint:idxEndPoint + 1]

            if spline.side == "PS":
                points.reverse()

            spline.offsetSpline = Spline(points)
            
            self.offsetSplines.append(spline.offsetSpline)                

        for spline in self.splines:
            print(f"First Point: {spline.points[0].tag} / Last Point: {spline.points[-1].tag}")
        for spline in self.offsetSplines:
            print(f"First Point: {spline.points[0].tag} / Last Point: {spline.points[-1].tag}")

    def extendTE(self, extensionLength, transitionLength, extension_offset = None):
        self.extensionSplines = []
        print("TE extension")

        lastPoint_SS = self.points[self.idxTE - 1]
        lastPoint_PS = self.points[self.idxTE + 1]
 
        tangVector_SS = np.array([self.pointTE.x - lastPoint_SS.x, self.pointTE.y - lastPoint_SS.y, 0])
        tangVector_SS_normalized = tangVector_SS/np.linalg.norm(tangVector_SS)

        tangVector_te_PS = np.array([self.pointTE.x - lastPoint_PS.x, self.pointTE.y - lastPoint_PS.y, 0])
        tangVector_te_PS_normalized = tangVector_te_PS/np.linalg.norm(tangVector_te_PS)

        tangVector_normalized = (tangVector_SS_normalized + tangVector_te_PS_normalized)/np.linalg.norm(tangVector_SS_normalized + tangVector_te_PS_normalized)

        # Build te extension spline
        endPoint = SplineInterpolation.generateEndPoint(self.pointTE, tangVector_normalized, transitionLength, extension_offset)
        controlPoint = SplineInterpolation.generateControlPoint(self.pointTE, tangVector_normalized, endPoint)
        interpolationPoints = SplineInterpolation.quadBezierInterpolation(self.pointTE, controlPoint, endPoint)

        extensionSpline = ExtensionSpline(interpolationPoints)
        self.extensionSplines.append(extensionSpline)
            
        # Build point set for upper extension line:
        startPoint = endPoint
        endPoint = Point(self.pointTE.x + extensionLength, startPoint.y, startPoint.z)

        extensionSpline = ExtensionSpline([startPoint, endPoint])
        self.extensionSplines.append(extensionSpline)

    def generate(self):
        for point in self.points:
            point.generate()
        for spline in self.splines:
            spline.generate()                

    def setTransfinite(self, meshSize):
        self.meshSize = meshSize
        
        refinementSize = 1/4*self.meshSize
        for spline in self.splines:
            if self.pointLE in spline.points and self.refinementLengthLE != None:
                spline.setTransfinite(startMeshSize = refinementSize, endMeshSize = self.meshSize)
            else:
                spline.setTransfinite(startMeshSize = self.meshSize, endMeshSize = self.meshSize)

        if self.extensionSplines:
            self.extensionSplines[0].setTransfinite(startMeshSize = self.meshSize, endMeshSize = self.meshSize)
            self.extensionSplines[1].setTransfinite(startMeshSize = self.extensionSplines[0].endMeshSize, endMeshSize = self.extensionSplines[1].length/10*2)

class AirfoilStructuredRegion:
    def __init__(self, airfoil, offsetValue, extensionLength, transitionLength, conicalWakeAngle):
        airfoil.offset(offsetValue)
        airfoil.extendTE(extensionLength, transitionLength)
        for spline in airfoil.splines:
            spline.offsetSpline.generate()
        for spline in airfoil.extensionSplines:
            spline.generate()
            
        gmsh.model.occ.synchronize()
        gmsh.fltk.run()
        self.offsetValue = airfoil.offsetValue
        self.airfoil = airfoil

        self.points = []
        self.splines = []
        self.wakeCurves = airfoil.extensionSplines
        self.controlCurves = []
        self.planeSurfaces = []
        
        self.createWakeRegion(conicalWakeAngle)
        self.createPlaneSurfaces()

        self.upper_connection_lines = []
        self.lower_connection_lines = []

    def createWakeRegion(self, conicalWakeAngle):
        for spline in self.airfoil.extensionSplines:
            spline.offsetSplines = []

        # Lines
        angle = conicalWakeAngle * (math.pi / 180)
        direction = 1
        refStartPoint = self.airfoil.extensionSplines[-1].points[0]
        refEndPoint = self.airfoil.extensionSplines[-1].points[1]

        for idx in range(2):
            print(idx)
            startPoint = self.airfoil.offsetPoints[self.airfoil.idxTE + idx]
            endPoint = Point(refStartPoint.x, refStartPoint.y + direction*self.offsetValue, refStartPoint.z)
            lastPoint = self.airfoil.offsetPoints[self.airfoil.idxTE + idx - direction]

            tangVector = np.array([startPoint.x - lastPoint.x, startPoint.y - lastPoint.y, 0])
            tangVector = tangVector/np.linalg.norm(tangVector)

            controlPoint = SplineInterpolation.generateControlPoint(startPoint, tangVector, endPoint)
            interpolationPoints = SplineInterpolation.quadBezierInterpolation(startPoint, controlPoint, endPoint)
            spline = Spline(interpolationPoints)
            self.airfoil.extensionSplines[0].offsetSplines.append(spline)
            spline.generate()
            
            startPoint = endPoint
            endPoint = Point(refEndPoint.x, refEndPoint.y + direction*(self.offsetValue + (refEndPoint.x - refStartPoint.x)*math.tan(angle)), refEndPoint.z)
            spline = Spline([startPoint, endPoint])
            self.airfoil.extensionSplines[1].offsetSplines.append(spline)
            spline.generate()
            
            direction = -1
            
    def createPlaneSurfaces(self):
        for spline in self.airfoil.splines:
            spline.controlSplines = []
            offsetSpline = spline.offsetSpline
            
            firstControlSpline = Spline.get([spline.points[0], offsetSpline.points[0]])
            if not firstControlSpline:
                firstControlSpline = Spline([spline.points[0], offsetSpline.points[0]])
                firstControlSpline.generate()
            spline.controlSplines.append(firstControlSpline)

            lastControlSpline = Spline.get([spline.points[-1], offsetSpline.points[-1]])
            if not lastControlSpline:
                lastControlSpline = Spline([spline.points[-1], offsetSpline.points[-1]])
                lastControlSpline.generate()
            spline.controlSplines.append(lastControlSpline)
            
            curveList = [spline, spline.controlSplines[0], offsetSpline, spline.controlSplines[1]]
            self.planeSurfaces.append(PlaneSurface([CurveLoop(curveList)]))

        for spline in self.airfoil.extensionSplines:
            spline.controlSplines = []

            for offsetSpline in spline.offsetSplines:            
                firstControlSpline = Spline.get([spline.points[0], offsetSpline.points[0]])
                if not firstControlSpline:
                    firstControlSpline = Spline([spline.points[0], offsetSpline.points[0]])
                    firstControlSpline.generate()
                spline.controlSplines.append(firstControlSpline)
            
                lastControlSpline = Spline.get([spline.points[-1], offsetSpline.points[-1]])
                if not lastControlSpline:
                    lastControlSpline = Spline([spline.points[-1], offsetSpline.points[-1]])
                    lastControlSpline.generate()
                spline.controlSplines.append(lastControlSpline)
            
                curveList = [spline, spline.controlSplines[-2], offsetSpline, spline.controlSplines[-1]]
                print([curve.tag for curve in curveList])
                self.planeSurfaces.append(PlaneSurface([CurveLoop(curveList)]))

    def close_loop(self):
        """
        Method to form a close loop with the current geometrical object

        Returns
        -------
        _ : int
            return the tag of the CurveLoop object
        """                           
        outer_loop_curve_list = []
        list = []
        for airfoilSpline in self.airfoil.splines:
            if airfoilSpline.side == "SS":
                list.append(airfoilSpline.offsetSpline)
        outer_loop_curve_list.extend(list)

        list = []
        for airfoilSpline in self.airfoil.splines:
            if airfoilSpline.side == "PS":
                list.append(airfoilSpline.offsetSpline)
        list.reverse()
        outer_loop_curve_list.extend(list)
        
        outer_loop_curve_list.extend([extensionSpline.offsetSplines[0] for extensionSpline in self.airfoil.extensionSplines])
        outer_loop_curve_list.extend([extensionSpline.offsetSplines[1] for extensionSpline in self.airfoil.extensionSplines])
        outer_loop_curve_list.append(self.airfoil.extensionSplines[-1].controlSplines[-3])
        outer_loop_curve_list.append(self.airfoil.extensionSplines[-1].controlSplines[-1])
        
        return CurveLoop(outer_loop_curve_list).tag

    def setTransfinite(self, deltaYWall, boundaryLayerThickness, N):
        # Face normal transfinite mesh settings
        growthRate = ((boundaryLayerThickness - deltaYWall/2)/(boundaryLayerThickness - 10*deltaYWall/2))**(1.5*N)
        print(f"Boundary layer parameter:   y+Wall = {deltaYWall/2}, delta = {boundaryLayerThickness}, growth rate = {growthRate}")

        deltaYWall = deltaYWall*(N+1)

        for airfoilSpline in self.airfoil.splines:
            for controlSpline in airfoilSpline.controlSplines:
                if controlSpline.nTransfinitePoints:
                    continue
                else:
                    controlSpline.setTransfinite(startMeshSize=deltaYWall, growthRate=growthRate)

        for airfoilSpline in self.airfoil.extensionSplines:
            for controlSpline in airfoilSpline.controlSplines:
                if controlSpline.nTransfinitePoints:
                    continue
                else:
                    controlSpline.setTransfinite(startMeshSize=deltaYWall, growthRate=growthRate)

        self.airfoil.extensionSplines[-1].controlSplines[-1].setTransfinite(nPnts = self.airfoil.extensionSplines[-1].controlSplines[-2].nTransfinitePoints, growthRate = self.airfoil.extensionSplines[-1].controlSplines[-2].growthRate)
        self.airfoil.extensionSplines[-1].controlSplines[-3].setTransfinite(nPnts = self.airfoil.extensionSplines[-1].controlSplines[-4].nTransfinitePoints, growthRate = self.airfoil.extensionSplines[-1].controlSplines[-2].growthRate)

        # Calculate the tangential mesh distribution
        if self.airfoil.refinementLengthLE != None:
            for airfoilSpline in self.airfoil.splines:
                if self.airfoil.pointLE not in airfoilSpline.points:
                    airfoilSpline.offsetSpline.setTransfinite(growthRate = 1, nPnts = airfoilSpline.nTransfinitePoints)

            for idxAirfoilSpline, airfoilSpline in enumerate(self.airfoil.splines):
                if self.airfoil.pointLE in airfoilSpline.points:
                    if airfoilSpline.side == "SS":
                        nextSpline = self.airfoil.splines[idxAirfoilSpline + 1]
                    elif airfoilSpline.side == "PS":
                        nextSpline = self.airfoil.splines[idxAirfoilSpline - 1]
                    
                    airfoilSpline.offsetSpline.setTransfinite(endMeshSize = nextSpline.offsetSpline.endMeshSize, nPnts = airfoilSpline.nTransfinitePoints)
        else:
            for airfoilSpline in self.airfoil.splines:
                airfoilSpline.offsetSpline.setTransfinite(growthRate = 1, nPnts = airfoilSpline.nTransfinitePoints)

        idxTESpline = None
        for spline in self.airfoil.splines:
            if self.airfoil.pointTE in spline.points and spline.side == "SS":
                idxTESpline = self.airfoil.splines.index(spline)

        for offsetSpline in self.airfoil.extensionSplines[0].offsetSplines:
            lastSpline = self.airfoil.splines[idxTESpline]
            offsetSpline.setTransfinite(startMeshSize = lastSpline.offsetSpline.startMeshSize, nPnts = self.airfoil.extensionSplines[0].nTransfinitePoints)
            idxTESpline += 1

        for idxOffsetSpline, offsetSpline in enumerate(self.airfoil.extensionSplines[-1].offsetSplines):
            lastSpline = self.airfoil.extensionSplines[0]
            offsetSpline.setTransfinite(startMeshSize = self.airfoil.extensionSplines[-1].startMeshSize, nPnts = self.airfoil.extensionSplines[-1].nTransfinitePoints)

        # Set and recombine transfinite surfaces
        for surface in self.planeSurfaces:
            gmsh.model.mesh.setTransfiniteSurface(surface.tag)
            gmsh.model.mesh.setRecombine(2, surface.tag)
            # gmsh.model.mesh.setSmoothing(2, surface.tag, 20)

        print("All transfinte conditions set")

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
        self.dim = 2
        
        if all(isinstance(geom_object, CurveLoop) for geom_object in geom_objects):
            self.tag_list = [geom_object.tag for geom_object in self.geom_objects]
            
            # create the gmsh object and store the tag of the geometric object
            self.tag = gmsh.model.occ.addPlaneSurface(self.tag_list)
        else:
            # close_loop() will form a close loop object and return its tag
            self.tag_list = [geom_object.close_loop() for geom_object in self.geom_objects]

            # create the gmsh object and store the tag of the geometric object
            self.tag = gmsh.model.occ.addPlaneSurface(self.tag_list)

    def define_bc(self):
        """
        Method that define the domain marker of the surface
        -------
        """
        bc_name = "fluid"
        bc_entities = [self.tag]

        for bc in BoundaryCondition.instances:
            if bc.name == bc_name and bc.dim == self.dim:
                self.bc = bc 
                bc.tag_list.append(self.tag)
                return
        
        self.bc = BoundaryCondition(self.dim, bc_name, bc_entities)

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
    def __init__(self, plane_surfaces, extrusion_value, numElements):
        
        self.dim = 3
        self.extrusion_value = extrusion_value
        self.plane_surfaces = plane_surfaces
        self.numElements = numElements

        self.tagList2D = [entity[1] for entity in gmsh.model.occ.getEntities(1)]
        dimtags = []
        for planeSurface in plane_surfaces:
            dimtags.append((planeSurface.dim,planeSurface.tag))

        self.extrude_dimtags = gmsh.model.occ.extrude(dimtags, 0, 0, extrusion_value, [self.numElements], [1], recombine = True)
        gmsh.model.occ.synchronize()
        self.bcs =[]

    def define_bc(self):
        """
        Method that define the domain marker of the surface
        -------
        """
        # Rename 2D fluid domain:
        bc_name, bc_dim = "fluid", 2
        currentBC = next((instance for instance in BoundaryCondition.instances if instance.name == bc_name and instance.dim == bc_dim), None)
        if currentBC != None:
            currentBC.name = "side_z-"

        for extrude_dimtag in self.extrude_dimtags:
            if extrude_dimtag[0] == 3:
                bc_name, bc_dim = "fluid", 3
                bc_entities = [extrude_dimtag[1]]
            elif gmsh.model.getBoundary([extrude_dimtag], combined=True, oriented=False, recursive=False)[0][1] not in self.tagList2D:
                bc_name, bc_dim = "side_z+", 2
                bc_entities = [extrude_dimtag[1]]
            else:
                for bc in BoundaryCondition.instances:
                    if gmsh.model.getBoundary([extrude_dimtag], combined=True, oriented=False, recursive=False)[0][1] in bc.tag_list and bc.dim == 1:
                        bc_name, bc_dim = bc.name, 2
                        bc_entities = [extrude_dimtag[1]]
                        break
                else:
                    continue
                    
            currentBC = next((instance for instance in BoundaryCondition.instances if instance.name == bc_name and instance.dim == bc_dim), None)
            if currentBC != None:
                currentBC.tag_list.extend(bc_entities)
            else:
                currentBC = BoundaryCondition(bc_dim, bc_name, bc_entities)
            
        for bc in BoundaryCondition.instances:
            print(f"{bc.name}, dim = {bc.dim}, tags = {bc.tag_list}")

class AirfoilBoundaryCondition():
    def __init__(self, name, start, end, side):
        self.name = name
        self.start = start
        self.end = end
        self.side = side

class BoundaryCondition(metaclass=InstanceTracker):
    def __init__(self, dim, name, tag_list):
        self.dim = dim
        self.name = name

        self.tag_list = tag_list
        self.tag = None

    @staticmethod
    def generatePhysicalGroups(*dims):
        for bc in BoundaryCondition.instances:
            if not dims:
                print(f"Created BoundaryCondition: {bc.name}")
                bc.tag = gmsh.model.addPhysicalGroup(bc.dim, bc.tag_list, name = bc.name)
                continue
            for dim in dims:
                match dim:
                    case 1:
                        if dim == bc.dim:
                            print(f"Created BoundaryCondition: {bc.name}")
                            bc.tag = gmsh.model.addPhysicalGroup(bc.dim, bc.tag_list, name = bc.name)
                    case 2:
                        if dim == bc.dim:
                            print(f"Created BoundaryCondition: {bc.name}")
                            bc.tag = gmsh.model.addPhysicalGroup(bc.dim, bc.tag_list, name = bc.name)
                    case 3:
                        if dim == bc.dim:
                            print(f"Created BoundaryCondition: {bc.name}")
                            bc.tag = gmsh.model.addPhysicalGroup(bc.dim, bc.tag_list, name = bc.name)
                    case _:
                        raise ValueError(f"{dim} is no allowed value for dimension!")

    @staticmethod
    def exists(dim, name):
        for bc in BoundaryCondition.instances:
            if bc.name == name and bc.dim == dim:
                return True
        else:
            return False
    
    def get(dim, name):
        for bc in BoundaryCondition.instances:
            if bc.name == name and bc.dim == dim:
                return bc
        else:
            raise ValueError(f"No boundary condition {name} of dimension {dim} found!")