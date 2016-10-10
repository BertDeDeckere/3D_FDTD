import os
os.environ['QT_API'] = 'pyqt'
# os.environ['ETS_TOOLKIT'] = 'qt4'

from pyface.qt import QtGui, QtCore

from traits.api import HasTraits, Instance, on_trait_change, Int, Dict
from traitsui.api import View, Item
from tvtk.pyface.api import Scene
from tvtk.api import tvtk
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor
from mayavi import mlab
from multiprocessing import cpu_count
from time import time

import matplotlib			# Qt compatibility
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# from matplotlib.patches import Circle, Wedge, Polygon
import matplotlib.patches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from numpy import pi, sin, cos, sqrt, amin, amax, ones, array, zeros, concatenate, arange, int32, arccos, matrix, cumsum, dot, roll, real, imag, \
	reshape, linspace, angle
from scipy.sparse.linalg import eigsh, eigs
from scipy import sparse
import scipy.io as io
from math import ceil
from triangle import triangulate
from Field import Field

class Cutplane():
	def __init__(self, orientation):
		self.orientation = orientation

	def computeNormal(self, P1, P2):		# P1 and P2 are vectors on a plane
		self.n = [P1[1]*P2[2]-P1[2]*P2[1], \
				P1[2]*P2[0]-P1[0]*P2[2], \
				P1[0]*P2[1]-P1[1]*P2[0]]

	def defineLineSegments(self, bottom, top):
		self.computeNormal([bottom[2][0]-bottom[1][0], bottom[2][1]-bottom[1][1], bottom[2][2]-bottom[1][2]], \
						[bottom[1][0]-bottom[0][0], bottom[1][1]-bottom[0][1], bottom[1][2]-bottom[0][2]])
		self.segmentsBottom = []
		self.segmentsTop = []
		if self.n.index(max(self.n)) == 0:	# Drop x-coordinate
			for k in range(len(bottom)):
				self.segmentsBottom.append([bottom[k][1], bottom[k][2]])
				self.segmentsTop.append([top[k][1], top[k][2]])
		elif self.n.index(max(self.n)) == 1:
			for k in range(len(bottom)):
				self.segmentsBottom.append([bottom[k][0], bottom[k][2]])
				self.segmentsTop.append([top[k][0], top[k][2]])
		elif self.n.index(max(self.n)) == 2:
			for k in range(len(bottom)):
				self.segmentsBottom.append([bottom[k][0], bottom[k][1]])
				self.segmentsTop.append([top[k][0], top[k][1]])

	def pointInPolygon(self, P, points):
		i = 0
		j = len(points)-1
		c = False
		while i < len(points):
			if ((points[i][1] > P[1]) != (points[j][1] > P[1])) and \
				(P[0] < (points[j][0] - points[i][0])*(P[1] - points[i][1])/(points[j][1] - points[i][1]) + points[i][0]):
				c = not(c)
			j = i
			i += 1
		return c

	def lineIntersection(self, P0, P1, Q0, Q1):
		A = P1[0]-P0[0]
		B = Q0[0]-Q1[0]
		C = P1[1]-P0[1]
		D = Q0[1]-Q1[1]
		E = Q0[0]-P0[0]
		F = Q0[1]-P0[1]

		det = A*D-B*C
		if det != 0:
			return ((D*E-B*F)/det, (A*F-E*C)/det)
		return (-1, -1)

	def lineInPolygon(self, P0, P1):
		if self.n.index(max(self.n)) == 0:	# Drop x-coordinate
			midPoint = [(P0[1]+P1[1])/2.0, (P0[2]+P1[2])/2.0]
			if self.pointInPolygon(midPoint, self.segmentsBottom) or self.pointInPolygon(midPoint, self.segmentsTop):
				Q = self.segmentsTop
				length = len(self.segmentsTop)
				for k in range(length):
					u,v = self.lineIntersection([P0[1], P0[2]], [P1[1], P1[2]], [Q[k][0], Q[k][1]], [Q[(k+1)%length][0], Q[(k+1)%length][1]])
					if u > 1E-10 and u < 1-1E-10 and v > 1E-10 and v < 1-1E-10:
						return False

				Q = self.segmentsBottom
				for k in range(length):
					u,v = self.lineIntersection([P0[1], P0[2]], [P1[1], P1[2]], [Q[k][0], Q[k][1]], [Q[(k+1)%length][0], Q[(k+1)%length][1]])
					if u > 1E-10 and u < 1-1E-10 and v > 1E-10 and v < 1-1E-10:
						return False
			else:
				return False
		elif self.n.index(max(self.n)) == 1:
			midPoint = [(P0[0]+P1[0])/2.0, (P0[2]+P1[2])/2.0]
			if self.pointInPolygon(midPoint, self.segmentsBottom) or self.pointInPolygon(midPoint, self.segmentsTop):
				Q = self.segmentsTop
				length = len(self.segmentsTop)
				for k in range(length):
					u,v = self.lineIntersection([P0[0], P0[2]], [P1[0], P1[2]], [Q[k][0], Q[k][1]], [Q[(k+1)%length][0], Q[(k+1)%length][1]])
					if u > 1E-10 and u < 1-1E-10 and v > 1E-10 and v < 1-1E-10:
						return False

				Q = self.segmentsBottom
				for k in range(length):
					u,v = self.lineIntersection([P0[0], P0[2]], [P1[0], P1[2]], [Q[k][0], Q[k][1]], [Q[(k+1)%length][0], Q[(k+1)%length][1]])
					if u > 1E-10 and u < 1-1E-10 and v > 1E-10 and v < 1-1E-10:
						return False
			else:
				return False

		elif self.n.index(max(self.n)) == 2:
			midPoint = [(P0[0]+P1[0])/2.0, (P0[1]+P1[1])/2.0]
			if self.pointInPolygon(midPoint, self.segmentsBottom) or self.pointInPolygon(midPoint, self.segmentsTop):
				Q = self.segmentsTop
				length = len(self.segmentsTop)
				for k in range(length):
					u,v = self.lineIntersection([P0[0], P0[1]], [P1[0], P1[1]], [Q[k][0], Q[k][1]], [Q[(k+1)%length][0], Q[(k+1)%length][1]])
					if u > 1E-10 and u < 1-1E-10 and v > 1E-10 and v < 1-1E-10:
						return False

				Q = self.segmentsBottom
				for k in range(length):
					u,v = self.lineIntersection([P0[0], P0[1]], [P1[0], P1[1]], [Q[k][0], Q[k][1]], [Q[(k+1)%length][0], Q[(k+1)%length][1]])
					if u > 1E-10 and u < 1-1E-10 and v > 1E-10 and v < 1-1E-10:
						return False
			else:
				return False

		return True

	def eq(self, a, b, eps=1E-8):
		return abs( a[0] - b[0] ) <= eps*abs(a[0]) and abs( a[1] - b[1] ) <= eps*abs(a[1])

	def computeIntersection(self, pos, bottom, top):		# P = points, faces contains triplets
		length = len(bottom)
		points = concatenate([bottom,top])
		faces = [[k, (k+1)%length, (k+1)%length + length, k+length] for k in range(length)]

		segments = []
		if self.orientation == 'xy':
			z = pos
			for face in faces:
				P0 = points[face[0]]
				P1 = points[face[1]]
				P2 = points[face[2]]
				P3 = points[face[3]]

				# Faces should have non-zero area, otherwise errors occur, because 3 points are addes in the intersection
				temp = []
				if z > min(P0[2], P1[2]) and z < max(P0[2], P1[2]):		# line segment between P0 and P1
					if P1[2] != P0[2]:
						u = (z-P0[2])/(P1[2] - P0[2])
						if u >= 0 and u <= 1:
							y = u*(P1[1] - P0[1]) + P0[1]
							x = u*(P1[0] - P0[0]) + P0[0]
							temp.append([x,y])

				if z >= min(P1[2], P2[2]) and z <= max(P1[2], P2[2]):		# line segment between P1 and P2
					if P2[2] != P1[2]:
						u = (z-P1[2])/(P2[2] - P1[2])
						if u >= 0 and u <= 1:
							y = u*(P2[1] - P1[1]) + P1[1]
							x = u*(P2[0] - P1[0]) + P1[0]
							temp.append([x,y])

				if z > min(P2[2], P3[2]) and z < max(P2[2], P3[2]):		# line segment between P1 and P2
					if P3[2] != P2[2]:
						u = (z-P2[2])/(P3[2] - P2[2])
						if u >= 0 and u <= 1:
							y = u*(P3[1] - P2[1]) + P2[1]
							x = u*(P3[0] - P2[0]) + P2[0]
							temp.append([x,y])

				if z >= min(P3[2], P0[2]) and z <= max(P3[2], P0[2]):		# line segment between P2 and P0
					if P0[2] != P3[2]:
						u = (z-P3[2])/(P0[2] - P3[2])
						if u >= 0 and u <= 1:
							y = u*(P0[1] - P3[1]) + P3[1]
							x = u*(P0[0] - P3[0]) + P3[0]
							temp.append([x,y])

				if len(temp) >= 2:
					segments.append(temp)

		elif self.orientation == 'yz':
			x = pos
			for face in faces:
				P0 = points[face[0]]
				P1 = points[face[1]]
				P2 = points[face[2]]
				P3 = points[face[3]]

				# Faces should have non-zero area, otherwise errors occur, because 3 points are addes in the intersection
				temp = []
				if x > min(P0[0], P1[0]) and x < max(P0[0], P1[0]):		# line segment between P0 and P1
					if P1[0] != P0[0]:
						u = (x-P0[0])/(P1[0] - P0[0])
						if u >= 0 and u <= 1:
							y = u*(P1[1] - P0[1]) + P0[1]
							z = u*(P1[2] - P0[2]) + P0[2]
							temp.append([y,z])

				if x >= min(P1[0], P2[0]) and x <= max(P1[0], P2[0]):		# line segment between P1 and P2
					if P2[0] != P1[0]:
						u = (x-P1[0])/(P2[0] - P1[0])
						if u >= 0 and u <= 1:
							y = u*(P2[1] - P1[1]) + P1[1]
							z = u*(P2[2] - P1[2]) + P1[2]
							temp.append([y,z])

				if x > min(P2[0], P3[0]) and x < max(P2[0], P3[0]):		# line segment between P1 and P2
					if P3[0] != P2[0]:
						u = (x-P2[0])/(P3[0] - P2[0])
						if u >= 0 and u <= 1:
							y = u*(P3[1] - P2[1]) + P2[1]
							z = u*(P3[2] - P2[2]) + P2[2]
							temp.append([y,z])

				if x >= min(P3[0], P0[0]) and x <= max(P3[0], P0[0]):		# line segment between P2 and P0
					if P0[0] != P3[0]:
						u = (x-P3[0])/(P0[0] - P3[0])
						if u >= 0 and u <= 1:
							y = u*(P0[1] - P3[1]) + P3[1]
							z = u*(P0[2] - P3[2]) + P3[2]
							temp.append([y,z])

				if len(temp) >= 2:
					segments.append(temp)

		elif self.orientation == 'zx':
			y = pos
			for face in faces:
				P0 = points[face[0]]
				P1 = points[face[1]]
				P2 = points[face[2]]
				P3 = points[face[3]]

				# Faces should have non-zero area, otherwise errors occur, because 3 points are addes in the intersection
				temp = []
				if y > min(P0[1], P1[1]) and y < max(P0[1], P1[1]):		# line segment between P0 and P1
					if P1[1] != P0[1]:
						u = (y-P0[1])/(P1[1] - P0[1])
						if u >= 0 and u <= 1:
							x = u*(P1[0] - P0[0]) + P0[0]
							z = u*(P1[2] - P0[2]) + P0[2]
							temp.append([x,z])

				if y >= min(P1[1], P2[1]) and y <= max(P1[1], P2[1]):		# line segment between P1 and P2
					if P2[1] != P1[1]:
						u = (y-P1[1])/(P2[1] - P1[1])
						if u >= 0 and u <= 1:
							x = u*(P2[0] - P1[0]) + P1[0]
							z = u*(P2[2] - P1[2]) + P1[2]
							temp.append([x,z])

				if y > min(P2[1], P3[1]) and y < max(P2[1], P3[1]):		# distinguish between < and <=, otherwise you end up with 3 points when the plane is at the edge
					if P3[1] != P2[1]:
						u = (y-P2[1])/(P3[1] - P2[1])
						if u >= 0 and u <= 1:
							x = u*(P3[0] - P2[0]) + P2[0]
							z = u*(P3[2] - P2[2]) + P2[2]
							temp.append([x,z])

				if y >= min(P3[1], P0[1]) and y <= max(P3[1], P0[1]):		# line segment between P2 and P0
					if P0[1] != P3[1]:
						u = (y-P3[1])/(P0[1] - P3[1])
						if u >= 0 and u <= 1:
							x = u*(P0[0] - P3[0]) + P3[0]
							z = u*(P0[2] - P3[2]) + P3[2]
							temp.append([x,z])

				if len(temp) >= 2:
					segments.append(temp)	

		# Connect the line segments
		polygonSet = []
		if len(segments) > 1:
			if self.eq(segments[0][1], segments[1][0]) or self.eq(segments[0][1], segments[1][1]):
				polygon = [segments[0][0], segments[0][1]]
			else:
				polygon = [segments[0][1], segments[0][0]]
			for k in range(1, len(segments)):
				if self.eq(polygon[-1], segments[k][0]):
					polygon.append(segments[k][1])
				elif self.eq(polygon[-1], segments[k][1]):
					polygon.append(segments[k][0])
				else:
					polygonSet.append(polygon)
					if k < len(segments)-1:
						if self.eq(segments[k][1], segments[k+1][0]) or self.eq(segments[k][1], segments[k+1][1]):
							polygon = [segments[k][0], segments[k][1]]
						else:
							polygon = [segments[k][1], segments[k][0]]
					else:
						polygon = [segments[k][1], segments[k][0]]

			polygonSet.append(polygon)

		# Look for any missing connections on top and bottom
		result = []
		index = -1
		self.defineLineSegments(bottom, top)
		filt = [True for _ in range(len(polygonSet))]
		for k in range(len(polygonSet)):
			if filt[k] == True:
				result.append(polygonSet[k])
				filt[k] = False
				index += 1
				for m in range(k+1, len(polygonSet)+1):
					if filt[m%len(polygonSet)] == True and k != m%len(polygonSet):
						if self.orientation == 'xy':
							P = [result[index][0][0], result[index][0][1], pos]
							Q1 = [polygonSet[m%len(polygonSet)][0][0], polygonSet[m%len(polygonSet)][0][1], pos]
							Q2 = [polygonSet[m%len(polygonSet)][-1][0], polygonSet[m%len(polygonSet)][-1][1], pos]
						elif self.orientation == 'yz':
							P = [pos, result[index][0][0], result[index][0][1]]
							Q1 = [pos, polygonSet[m%len(polygonSet)][0][0], polygonSet[m%len(polygonSet)][0][1]]
							Q2 = [pos, polygonSet[m%len(polygonSet)][-1][0], polygonSet[m%len(polygonSet)][-1][1]]
						elif self.orientation == 'zx':
							P = [result[index][0][0], pos, result[index][0][1]]
							Q1 = [polygonSet[m%len(polygonSet)][0][0], pos, polygonSet[m%len(polygonSet)][0][1]]
							Q2 = [polygonSet[m%len(polygonSet)][-1][0], pos, polygonSet[m%len(polygonSet)][-1][1]]
						if self.lineInPolygon(P, Q1):
							filt[m%len(polygonSet)] = False
							length = len(polygonSet[m%len(polygonSet)])
							for t in range(length):
								result[index].insert(0, polygonSet[m%len(polygonSet)][t])
						elif self.lineInPolygon(P, Q2):
							filt[m%len(polygonSet)] = False
							length = len(polygonSet[m%len(polygonSet)])
							for t in range(length):
								result[index].insert(0, polygonSet[m%len(polygonSet)][length-t-1])
		return result

	def rayCasting(self, polygons, height):		# Casts a horizontal ray at "height" and returns the intersections with polygons
		intersections = []
		for polygon in polygons:
			for k in range(len(polygon)):
				P0 = polygon[k]
				P1 = polygon[(k+1)%len(polygon)]

				if P1[1] != P0[1]:
					u = (height-P0[1])/(P1[1]-P0[1])
					if u >= 0 and u <= 1:
						intersections.append(P0[0] + u*(P1[0]-P0[0]))

		return sorted(intersections)

class Plotter2D():
	def __init__(self):
		self.figure2D = Figure(dpi=72, facecolor=(1,1,1), edgecolor=(0,0,0))
		self.axes2D = self.figure2D.add_subplot(1, 1, 1)
		self.axes2D.set_aspect('equal', adjustable='box')
		self.widget = FigureCanvas(self.figure2D)

	def area(self, points):		# Does not work for self intersecting polygons
		x = [points[k][0] for k in range(len(points))]
		y = [points[k][1] for k in range(len(points))]
		return 0.5*abs(dot(x,roll(y,1))-dot(y,roll(x,1)))

	def zeroArea(self, points):		# A line also gives False, which is incorrect, but in this case it's the required return value
		for k in range(1, len(points)):
			if points[k][0] != points[0][0] or points[k][1] != points[0][1]:
				return False
		return True

	def updatePlot(self, plane):
		self.axes2D.clear()

		if plane.type == 'Workplane':
			geometries = plane.geometries
			grid = plane.grid
			patches = []
			for geometry in geometries:
				if geometry.enable == True:
					if geometry.type == 'Polygon 2D':
						if geometry.material == None:
							color = (1, 0.5, 0)
						else:
							color = (geometry.material.color[0]/255.0, geometry.material.color[1]/255.0, geometry.material.color[2]/255.0)

						if len(geometry.P) > 2:		# At least 2 points are needed to make a polygon
							polygon = matplotlib.patches.Polygon(array(geometry.P)/grid.sizeScalefactor(), fc=color, alpha=geometry.opacity/100.0)
						if self.zeroArea(geometry.P) == False:	
							patches.append(polygon)

					if geometry.type == 'Circle':
						if geometry.material == None:
							color = (1, 0.5, 0)
						else:
							color = (geometry.material.color[0]/255.0, geometry.material.color[1]/255.0, geometry.material.color[2]/255.0)

						circle = matplotlib.patches.Circle((geometry.center[0]/grid.sizeScalefactor(), geometry.center[1]/grid.sizeScalefactor()), \
							radius=geometry.radius/grid.sizeScalefactor(), fc=color, alpha=geometry.opacity/100.0)
						if geometry.radius > 0:
							patches.append(circle)

			if len(patches) > 0:
				p = PatchCollection(patches, match_original=True)
				# p.set_array(array(colors))
				self.axes2D.set_xlabel("u ["+grid.sizeUnit+"]")
				self.axes2D.set_ylabel("v ["+grid.sizeUnit+"]")
				self.axes2D.add_collection(p)
				self.axes2D.autoscale()
				self.figure2D.canvas.draw()

		elif plane.type == 'Port Plane':
			grid = plane.project.grid
			patches = []
			planeType = ''
			workplanes = (x for x in plane.project.geometries if x.type == 'Workplane')
			for workplane in workplanes:
				extrusions = (x for x in workplane.geometries if x.type == 'Extrusion')
				for extrusion in extrusions:
					for geometry in workplane.geometries:
						top = []
						bottom = []
						if geometry.type == 'Circle':
							resolution = 100
							center = geometry.center
							radius = geometry.radius
							dalpha = 2*pi/resolution

							points = [[radius*cos(k*dalpha)+center[0], radius*sin(k*dalpha)+center[1]] for k in range(resolution)]

						if geometry.type == 'Polygon 2D':
							points = geometry.P

						if geometry.type != 'Extrusion' and len(points) > 0:
							length = len(points)

							bottom = workplane.uvwToxyz([[points[k][0], points[k][1], 0] for k in range(length)])
							top = workplane.uvwToxyz([[points[k][0], points[k][1], extrusion.height] for k in range(length)])

							cutplane = Cutplane(orientation=plane.plane)
							intersection = cutplane.computeIntersection(plane.pos, bottom, top)

							if geometry.material == None:
								color = (1, 0.5, 0)
							else:
								color = (geometry.material.color[0]/255.0, geometry.material.color[1]/255.0, geometry.material.color[2]/255.0)

							for pointset in intersection:
								if len(pointset) > 1 and self.zeroArea(pointset) == False:
									polygon = matplotlib.patches.Polygon(array(pointset)/grid.sizeScalefactor(), fc=color, alpha=geometry.opacity/100.0)
									patches.append(polygon)

			for portPlane in plane.project.excitations:
				if portPlane.type == 'Port Plane':
					planeType = portPlane.plane
					for port in portPlane.sources:
						if port.type == 'Waveguide Port' and len(port.v) == 0:
							points = [[port.pos[0][0], port.pos[1][0]], [port.pos[0][0], port.pos[1][1]], \
								[port.pos[0][1], port.pos[1][1]], [port.pos[0][1], port.pos[1][0]]]
							polygon = matplotlib.patches.Polygon(array(points)/grid.sizeScalefactor(), fc=(0.11, 0.56, 1), alpha=0.5)
							if self.zeroArea(points) == False:
								patches.append(polygon)
						else:
							# Ex = port.v[0:(port.sizeX)*(port.sizeY-1), port.mode]
							# Ex = Ex.reshape(port.sizeX, port.sizeY-1)
							# Ey = port.v[port.sizeX*(port.sizeY-1):port.sizeX*(port.sizeY-1) + (port.sizeX-1)*port.sizeY, port.mode]
							# Ey = Ey.reshape(port.sizeX-1, port.sizeY)
							# print port.v.shape, (port.sizeX)*(port.sizeY-1), port.sizeX*(port.sizeY-1) + (port.sizeX-1)*port.sizeY

							data = []
							if port.field == 'Ex':
								data = port.Ex[port.mode]
							elif port.field == 'Ey':
								data = port.Ey[port.mode]
							elif port.field == 'Ez':
								data = port.Ez[port.mode]

							# print array(data).shape
							self.axes2D.imshow(abs(data.T), extent=[port.pos[0][0]/grid.sizeScalefactor(), port.pos[0][1]/grid.sizeScalefactor(), \
								port.pos[1][0]/grid.sizeScalefactor(), port.pos[1][1]/grid.sizeScalefactor()], origin='lower', cmap=cm.hot)

			if planeType != '':
				self.axes2D.set_xlabel("u ["+grid.sizeUnit+"]")
				self.axes2D.set_ylabel("v ["+grid.sizeUnit+"]")
			if planeType == 'xy':
				self.axes2D.set_xlim([0, grid.sizeX/grid.sizeScalefactor()])
				self.axes2D.set_ylim([0, grid.sizeY/grid.sizeScalefactor()])
			elif planeType == 'yz':
				self.axes2D.set_xlim([0, grid.sizeY/grid.sizeScalefactor()])
				self.axes2D.set_ylim([0, grid.sizeZ/grid.sizeScalefactor()])
			elif planeType == 'zx':
				self.axes2D.set_xlim([0, grid.sizeX/grid.sizeScalefactor()])
				self.axes2D.set_ylim([0, grid.sizeZ/grid.sizeScalefactor()])

			if len(patches) > 0:
				# print len(patches), patches[0], patches[0].get_xy()
				p = PatchCollection(patches, match_original=True)
				self.axes2D.add_collection(p)
				# self.axes2D.autoscale()
			self.figure2D.canvas.draw()

class Visualization(HasTraits):
	scene = Instance(MlabSceneModel, ())
	
	def __init__(self):
		self.objects = []
		# self.update = False

	# def rectilinearGrid(self, grid, data):
	# 	r = tvtk.RectilinearGrid()
	# 	r.point_data.scalars = data.ravel()
	# 	r.point_data.scalars.name = 'scalars'
	# 	r.dimensions = data.shape
	# 	r.x_coordinates = array(cumsum(grid.dx))
	# 	r.y_coordinates = array(cumsum(grid.dy))
	# 	r.z_coordinates = array(cumsum(grid.dz))
	# 	return r

	@on_trait_change('scene.activated')
	def update_plot(self, project=None):
		fig = mlab.gcf()
		fig.scene.disable_render = True
		if project != True:
			if project.refreshPlot == True:
				mlab.clf()

			for g in project.geometries:
				if g.enable == True and g.type == 'Workplane':
					vertices = []
					color = []
					opacity = []
					for item in g.geometries:		# First collect all vertices and their connectivity
						if item.type == 'Polygon 2D':
							l = len(item.P)
							vertices.append({'vertices':item.P, 'segments': array([[i%l,(i+1)%l] for i in range(l)], dtype=int32)})

							if item.material != None:
								color.append((item.material.color[0]/255.0, item.material.color[1]/255.0, item.material.color[2]/255.0))
							else:
								color.append((1, 0.5, 0))

							opacity.append(item.opacity/99.0)

						if item.type == 'Circle':
							resolution = 50		# Resolution for the approximation of the circle as polygon
							dalpha = 2.0*pi/resolution
							P = [[item.radius*cos(k*dalpha)+item.center[0], item.radius*sin(k*dalpha)+item.center[1]] for k in range(resolution)]
							vertices.append({'vertices':P, 'segments': array([[i%resolution,(i+1)%resolution] for i in range(resolution)], dtype=int32)})

							if item.material != None:
								color.append((item.material.color[0]/255.0, item.material.color[1]/255.0, item.material.color[2]/255.0))
							else:
								color.append((1, 0.5, 0))

							opacity.append(item.opacity/99.0)

					for item in g.geometries:		# Then do the post-processing
						if item.type == 'Extrusion':
							for v, c, o in zip(vertices, color, opacity):
								cndt = triangulate(v, 'pq')	# Constrained delaunay triangulation, alternative is ear clipping
								triangles = cndt['triangles']

								points1 = item.workplane.uvwToxyz([[cndt['vertices'][k][0], cndt['vertices'][k][1], 0] for k in range(len(cndt['vertices']))])
								points2 = item.workplane.uvwToxyz([[cndt['vertices'][k][0], cndt['vertices'][k][1], item.height] for k in range(len(cndt['vertices']))])

								length = len(v['vertices'])
								points3 = item.workplane.uvwToxyz([[v['vertices'][k][0], v['vertices'][k][1], 0] for k in range(length)])
								points3 += item.workplane.uvwToxyz([[v['vertices'][k][0], v['vertices'][k][1], item.height] for k in range(length)])

								a = [[k, (k+1)%length, k+length] for k in range(length)]
								b = [[(k+1)%length, k+length, (k+1)%length + length] for k in range(length)]
								triangles3 = [val for pair in zip(a, b) for val in pair]	# interleave a and b to get rid of rendering artefacts

								mesh1 = tvtk.PolyData(points=points1, polys=triangles)
								mesh2 = tvtk.PolyData(points=points2, polys=triangles)
								mesh3 = tvtk.PolyData(points=points3, polys=triangles3)
								surf1 = mlab.pipeline.surface(mesh1, opacity=o)
								surf2 = mlab.pipeline.surface(mesh2, opacity=o)
								surf3 = mlab.pipeline.surface(mesh3, opacity=o)
								mlab.pipeline.surface(surf1, color=c, opacity=o)
								mlab.pipeline.surface(surf2, color=c, opacity=o)
								mlab.pipeline.surface(surf3, color=c, opacity=o)

				if g.enable == True and g.type == 'Block':
					if g.material != None:
						color = (g.material.color[0]/255.0, g.material.color[1]/255.0, g.material.color[2]/255.0)
					else:
						color = (1, 0.5, 0)

					points = array([[g.x0, g.y0, g.z0], [g.x0, g.y0, g.z1], [g.x0, g.y1, g.z0], [g.x0, g.y1, g.z1], 
									[g.x1, g.y0, g.z0], [g.x1, g.y0, g.z1], [g.x1, g.y1, g.z0], [g.x1, g.y1, g.z1]])
					triangles = array([[0,1,2], [1,2,3], [4,5,7], [4,6,7], [0,2,4], [2,4,6], 
										[1,3,7], [1,5,7], [0,1,5], [0,4,5], [2,3,7], [2,6,7]])

					mesh = tvtk.PolyData(points=points, polys=triangles)
					surf = mlab.pipeline.surface(mesh, opacity=g.opacity/99.0)
					mlab.pipeline.surface(surf, color=color, opacity=g.opacity/99.0)
					# mlab.outline()

			for excitation in project.excitations:
				if excitation.type == 'Port Plane':
					if excitation.plane == 'xy':
						points = [[0, 0, excitation.pos], [0, project.grid.sizeY, excitation.pos], [project.grid.sizeX, 0, excitation.pos], \
							[project.grid.sizeX, project.grid.sizeY, excitation.pos]]
						triangles = [[0, 1, 2], [1, 2, 3]]
						mesh = tvtk.PolyData(points=points, polys=triangles)
						surf = mlab.pipeline.surface(mesh, opacity=0.5)
						mlab.pipeline.surface(surf, color=(0.1, 0.1, 0.1), opacity=0.2)
					elif excitation.plane == 'yz':
						points = [[excitation.pos, 0, 0], [excitation.pos, project.grid.sizeY, 0], [excitation.pos, 0, project.grid.sizeZ], \
							[excitation.pos, project.grid.sizeY, project.grid.sizeZ]]
						triangles = [[0, 1, 2], [1, 2, 3]]
						mesh = tvtk.PolyData(points=points, polys=triangles)
						surf = mlab.pipeline.surface(mesh, opacity=0.5)
						mlab.pipeline.surface(surf, color=(0.1, 0.1, 0.1), opacity=0.2)
					elif excitation.plane == 'zx':
						points = [[0, excitation.pos, 0], [0, excitation.pos, project.grid.sizeZ], [project.grid.sizeX, excitation.pos, 0], \
							[project.grid.sizeX, excitation.pos, project.grid.sizeZ]]
						triangles = [[0, 1, 2], [1, 2, 3]]
						mesh = tvtk.PolyData(points=points, polys=triangles)
						surf = mlab.pipeline.surface(mesh, opacity=0.5)
						mlab.pipeline.surface(surf, color=(0.1, 0.1, 0.1), opacity=0.2)

					for port in excitation.sources:
						if port.type == 'Waveguide Port':
							if excitation.plane == 'xy':
								points = [[port.pos[0][0], port.pos[1][0], excitation.pos], [port.pos[0][0], port.pos[1][1], excitation.pos], \
								[port.pos[0][1], port.pos[1][0], excitation.pos], [port.pos[0][1], port.pos[1][1], excitation.pos]]
								triangles = [[0, 1, 2], [1, 2, 3]]
								mesh = tvtk.PolyData(points=points, polys=triangles)
								surf = mlab.pipeline.surface(mesh, opacity=0.5)
								mlab.pipeline.surface(surf, color=(0.11, 0.56, 1), opacity=0.5)
							elif excitation.plane == 'yz':
								points = [[excitation.pos, port.pos[0][0], port.pos[1][0]], [excitation.pos, port.pos[0][1], port.pos[1][0]], \
								[excitation.pos, port.pos[0][0], port.pos[1][1]], [excitation.pos, port.pos[0][1], port.pos[1][1]]]
								mesh = tvtk.PolyData(points=points, polys=triangles)
								surf = mlab.pipeline.surface(mesh, opacity=0.5)
								mlab.pipeline.surface(surf, color=(0.11, 0.56, 1), opacity=0.5)
							elif excitation.plane == 'zx':
								points = [[port.pos[0][0], excitation.pos, port.pos[1][0]], [port.pos[0][0], excitation.pos, port.pos[1][1]], \
									[port.pos[0][1], excitation.pos, port.pos[1][0]], [port.pos[0][1], excitation.pos, port.pos[1][1]]]
								triangles = [[0, 1, 2], [1, 2, 3]]
								mesh = tvtk.PolyData(points=points, polys=triangles)
								surf = mlab.pipeline.surface(mesh, opacity=0.5)
								mlab.pipeline.surface(surf, color=(0.11, 0.56, 1), opacity=0.5)


			if project.result.showFields == True:
				n = project.result.timestep
				for p in project.plots:
					Ex = project.solver.Ex[n,:,:,:]
					Ey = project.solver.Ey[n,:,:,:]
					Ez = project.solver.Ez[n,:,:,:]
					Hx = project.solver.Hx[n,:,:,:]
					Hy = project.solver.Hy[n,:,:,:]
					Hz = project.solver.Hz[n,:,:,:]

					if p.type == 'Scalar Cutplane':
						data = eval(p.data)
						if p.plotted == False or project.refreshPlot == True:
							# dataset = self.rectilinearGrid(project.grid, data)

							# p.plot = mlab.pipeline.scalar_cut_plane(dataset, plane_orientation=p.orientation, 
							# 	vmin=p.minimum, vmax=p.maximum, colormap='seismic')
							p.plot = mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(data), plane_orientation=p.orientation, 
								slice_index=project.grid.cellsX/2, vmin=p.minimum, vmax=p.maximum, colormap='seismic')
							p.plotted = True
						else:
							p.plot.mlab_source.scalars = data

					elif p.type == 'Scalar Volume':
						data = eval(p.data)
						if p.plotted == False or project.refreshPlot == True:
							p.plot = mlab.pipeline.volume(mlab.pipeline.scalar_field(data), vmin=p.minimum, vmax=p.maximum)
							p.plotted = True
						else:
							p.plot.mlab_source.scalars = data

					elif p.type == 'Vector Volume':
						x = eval(p.xComponent)
						y = eval(p.yComponent)
						z = eval(p.zComponent)

						if p.plotted == False or project.refreshPlot == True:
							p.plot = mlab.quiver3d(x,y,z, mask_points=p.maskPoints, scale_factor=1, \
								vmin=p.minimum, vmax=p.maximum, mode='cone', scale_mode='scalar', resolution=20)
							p.plotted = True
						else:
							p.plot.mlab_source.set(u=x, v=y, w=z)

					elif p.type == 'Vector Cutplane':
						x = eval(p.xComponent)
						y = eval(p.yComponent)
						z = eval(p.zComponent)

						if p.plotted == False or project.refreshPlot == True:
							p.plot = mlab.pipeline.vector_cut_plane(mlab.pipeline.vector_field(x, y, z), mask_points=p.maskPoints, \
									scale_factor=1, vmin=p.minimum, vmax=p.maximum, mode='cone', scale_mode='scalar', resolution=20, \
									plane_orientation=p.orientation)
							p.plotted = True
						else:
							p.plot.mlab_source.set(u=x, v=y, w=z)

					elif p.type == 'Contour':
						data = eval(p.data)
						if p.plotted == False or project.refreshPlot == True:
							p.plot = mlab.contour3d(data, contours=p.contours, vmin=p.minimum, vmax=p.maximum, opacity=p.opacity, colormap='seismic')
							p.plotted = True
						else:
							p.plot.mlab_source.scalars = data

					# mlab.axes(extent=[0, project.grid.cellsX, 0, project.grid.cellsY, 0, project.grid.cellsZ])

			project.refreshPlot = False

			mlab.colorbar(orientation='vertical')
			if len(fig.scene.actor_list) > 0:
				outline = mlab.outline(line_width=1)
				outline.bounds = (0, project.grid.sizeX, 0, project.grid.sizeY, 0, project.grid.sizeZ)
		fig.scene.disable_render = False

	# the layout of the dialog screated
	view = View(Item('scene', editor=SceneEditor(),	# scene_class = Scene or MayaviScene
			height=250, width=300, show_label=False), resizable=True) # We need this to resize with the parent widget

class MayaviQWidget(QtGui.QWidget):
	def __init__(self, parent=None):
		QtGui.QWidget.__init__(self, parent)
		layout = QtGui.QVBoxLayout(self)
		layout.setContentsMargins(0,0,0,0)
		layout.setSpacing(0)
		self.visualization = Visualization()
		# self.visualization.configure_traits()

		# The edit_traits call will generate the widget to embed.
		self.ui = self.visualization.edit_traits(parent=self, kind='subpanel').control
		layout.addWidget(self.ui)
		self.ui.setParent(self)

class Grid():
	def __init__(self):
		self.fmax = 1E9 			# All values are stored with their SI unit
		self.sizeX = 1
		self.sizeY = 1 
		self.sizeZ = 1

		self.eps0 = 8.854187817E-12
		self.mu0 = 1.2566370614E-6
		self.c = 1/sqrt(self.eps0*self.mu0)
		self.courant = 0.99

		wavelength = self.c/self.fmax
		self.cellsX = (int)(ceil(self.sizeX/(wavelength/10)))
		self.cellsY = (int)(ceil(self.sizeY/(wavelength/10)))
		self.cellsZ = (int)(ceil(self.sizeZ/(wavelength/10)))

		self.dx = (wavelength/10)*ones(self.cellsX+1, dtype=float)
		self.dy = (wavelength/10)*ones(self.cellsY+1, dtype=float)
		self.dz = (wavelength/10)*ones(self.cellsZ+1, dtype=float)

		self.distanceX = cumsum(self.dx)
		self.distanceY = cumsum(self.dy)
		self.distanceZ = cumsum(self.dz)

		self.dt = self.courant/(self.c*sqrt(1.0/amin(self.dx)**2 + 1.0/amin(self.dy)**2 + 1.0/amin(self.dz)**2))

		self.timesteps = 50
		self.time = self.timesteps*self.dt

		self.PMLlayers = 10

		self.freqUnit = 'GHz'
		self.sizeUnit = 'mm'

		self.freqScalingDictionary = {'THz':1E12, 'GHz':1E9, 'MHz':1E6, 'kHz':1E3, 'Hz':1, 'mHz':1E-3}
		self.timeScalingDictionary = {'ps':1E-12, 'ns':1E-9, u"\u00B5s":1E-6, 'ms':1E-3, 's':1, 'ks':1E3}
		self.freqTimeConversion = {'THz':'ps', 'GHz':'ns', 'MHz':u"\u00B5s", 'kHz':'ms', 'Hz':'s', 'mHz':'ks'}
		self.sizeScalingDictionary = {'km':1E3, 'm':1, 'mm':1E-3, u"\u00B5m":1E-6, 'nm':1E-9}

	def computeDt(self):
		self.dt = self.courant/(self.c*sqrt(1.0/amin(self.dx)**2 + 1.0/amin(self.dy)**2 + 1.0/amin(self.dz)**2))
		return self.dt

	def computeCourant(self):
		self.courant = self.dt*(self.c*sqrt(1.0/amin(self.dx)**2 + 1.0/amin(self.dy)**2 + 1.0/amin(self.dz)**2))
		return self.courant

	def freqScalefactor(self):
		return self.freqScalingDictionary[self.freqUnit]

	def sizeScalefactor(self):
		return self.sizeScalingDictionary[self.sizeUnit]

	def timeScaleFactor(self):
		return self.timeScalingDictionary[self.freqTimeConversion[self.freqUnit]]

class GridWidget():
	def __init__(self, project):
		self.grid = project.grid
		self.project = project
		self.hold = False
		self.createWidget()

	def projectNameChanged(self, name):
		self.project.name = name
		self.project.item.setText(name)

	def freqUnitChanged(self, value):
		self.grid.freqUnit = self.freqUnitComboBox.itemText(value)

	def sizeUnitChanged(self, value):
		self.grid.sizeUnit = self.sizeUnitComboBox.itemText(value)

	def sizeXGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.sizeX = value*self.grid.sizeScalefactor()
			self.grid.cellsX = (int)(self.grid.sizeX/amax(self.grid.dx))
			self.cellsX.setValue(self.grid.cellsX)
			self.hold = False

	def sizeYGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.sizeY = value*self.grid.sizeScalefactor()
			self.grid.cellsY = (int)(self.grid.sizeY/amax(self.grid.dy))
			self.cellsY.setValue(self.grid.cellsY)
			self.hold = False

	def sizeZGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.sizeZ = value*self.grid.sizeScalefactor()
			self.grid.cellsZ = (int)(self.grid.sizeZ/amax(self.grid.dz))
			self.cellsZ.setValue(self.grid.cellsZ)
			self.hold = False

	def fmaxGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.fmax = value*self.grid.freqScalefactor()
			wavelength = self.grid.c/self.grid.fmax

			self.grid.cellsX = (int)(ceil(self.grid.sizeX/(wavelength/10)))
			self.grid.cellsY = (int)(ceil(self.grid.sizeY/(wavelength/10)))
			self.grid.cellsZ = (int)(ceil(self.grid.sizeZ/(wavelength/10)))

			self.cellsX.setValue(self.grid.cellsX)
			self.cellsY.setValue(self.grid.cellsY)
			self.cellsZ.setValue(self.grid.cellsZ)

			self.grid.dx = (wavelength/10)*ones(self.grid.cellsX+1)
			self.grid.dy = (wavelength/10)*ones(self.grid.cellsY+1)
			self.grid.dz = (wavelength/10)*ones(self.grid.cellsZ+1)

			self.dx.setValue(wavelength/10/self.grid.sizeScalefactor())
			self.dy.setValue(wavelength/10/self.grid.sizeScalefactor())
			self.dz.setValue(wavelength/10/self.grid.sizeScalefactor())

			self.grid.distanceX = cumsum(self.grid.dx)
			self.grid.distanceY = cumsum(self.grid.dy)
			self.grid.distanceZ = cumsum(self.grid.dz)

			self.grid.computeDt()
			self.dt.setValue(self.grid.dt/self.grid.timeScaleFactor())

			self.grid.time = self.grid.timesteps*self.grid.dt
			self.time.setValue(self.grid.time/self.grid.timeScaleFactor())
			self.hold = False

	def courantGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.courant = value
			self.grid.computeDt()
			self.dt.setValue(self.grid.dt/self.grid.timeScaleFactor())
			self.hold = False

	def timestepsGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.timesteps = (int)(value)
			self.grid.time = self.grid.timesteps*self.grid.dt
			self.time.setValue(self.grid.time/self.grid.timeScaleFactor())
			self.hold = False

	def PMLlayersGridChanged(self, value):
		self.grid.PMLlayers = value

	def cellsXGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.cellsX = (int)(value)
			self.grid.sizeX = self.grid.cellsX*amax(self.grid.dx)
			self.sizeX.setValue(self.grid.sizeX/self.grid.sizeScalefactor())
			self.hold = False

	def cellsYGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.cellsY = (int)(value)
			self.grid.sizeY = self.grid.cellsY*amax(self.grid.dy)
			self.sizeY.setValue(self.grid.sizeY/self.grid.sizeScalefactor())
			self.hold = False

	def cellsZGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.cellsZ = (int)(value)
			self.grid.sizeZ = self.grid.cellsZ*amax(self.grid.dz)
			self.sizeZ.setValue(self.grid.sizeZ/self.grid.sizeScalefactor())
			self.hold = False

	def dxGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.dx = value*ones(self.grid.cellsX+1)*self.grid.sizeScalefactor()
			self.grid.distanceX = cumsum(self.grid.dx)
			self.grid.computeCourant()
			self.courant.setValue(self.grid.courant)
			self.hold = False

	def dyGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.dy = value*ones(self.grid.cellsY+1)*self.grid.sizeScalefactor()
			self.grid.distanceY = cumsum(self.grid.dy)
			self.grid.computeCourant()
			self.courant.setValue(self.grid.courant)
			self.hold = False

	def dzGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.dz = value*ones(self.grid.cellsZ+1)*self.grid.sizeScalefactor()
			self.grid.distanceZ = cumsum(self.grid.dz)
			self.grid.computeCourant()
			self.courant.setValue(self.grid.courant)
			self.hold = False

	def dtGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.dt = value*self.grid.timeScaleFactor()
			self.grid.computeCourant()
			self.courant.setValue(self.grid.courant)
			self.grid.timesteps = self.grid.time/self.grid.dt 
			self.timesteps.setValue(self.grid.timesteps)
			self.hold = False

	def timeGridChanged(self, value):
		if self.hold == False:
			self.hold = True
			self.grid.time = value*self.grid.timeScaleFactor()
			self.grid.timesteps = self.grid.time/self.grid.dt
			self.timesteps.setValue(self.grid.timesteps)
			self.hold = False

	def createWidget(self):
		self.tab = QtGui.QWidget()
		self.sizeX = QtGui.QDoubleSpinBox()
		self.sizeY = QtGui.QDoubleSpinBox()
		self.sizeZ = QtGui.QDoubleSpinBox()
		self.fmax = QtGui.QDoubleSpinBox()
		self.courant = QtGui.QDoubleSpinBox()
		self.timesteps = QtGui.QDoubleSpinBox()
		self.PMLlayers = QtGui.QDoubleSpinBox()
		self.cellsX = QtGui.QDoubleSpinBox()
		self.cellsY = QtGui.QDoubleSpinBox()
		self.cellsZ = QtGui.QDoubleSpinBox()
		self.dx = QtGui.QDoubleSpinBox()
		self.dy = QtGui.QDoubleSpinBox()
		self.dz = QtGui.QDoubleSpinBox()
		self.dt = QtGui.QDoubleSpinBox()
		self.time = QtGui.QDoubleSpinBox()

		maximum = 1e9
		minimum = 0
		self.sizeX.setRange(minimum, maximum)
		self.sizeY.setRange(minimum, maximum)
		self.sizeZ.setRange(minimum, maximum)
		self.fmax.setRange(0, 100E9)
		self.courant.setRange(minimum, maximum)
		self.timesteps.setRange(0, 1E9)
		self.PMLlayers.setRange(0, 30)
		self.cellsX.setRange(minimum, maximum)
		self.cellsY.setRange(minimum, maximum)
		self.cellsZ.setRange(minimum, maximum)
		self.dx.setRange(minimum, 1E9)
		self.dy.setRange(minimum, 1E9)
		self.dz.setRange(minimum, 1E9)
		self.dt.setRange(minimum, 1E9)
		self.time.setRange(minimum, maximum)

		self.courant.setDecimals(3)
		self.PMLlayers.setDecimals(0)
		self.cellsX.setDecimals(0)
		self.cellsY.setDecimals(0)
		self.cellsZ.setDecimals(0)
		self.timesteps.setDecimals(0)
		self.dx.setDecimals(4)
		self.dy.setDecimals(4)
		self.dz.setDecimals(4)
		self.dt.setDecimals(4)

		self.sizeX.setValue(self.grid.sizeX/self.grid.sizeScalefactor())
		self.sizeY.setValue(self.grid.sizeY/self.grid.sizeScalefactor())
		self.sizeZ.setValue(self.grid.sizeZ/self.grid.sizeScalefactor())
		self.fmax.setValue(self.grid.fmax/self.grid.freqScalefactor())
		self.courant.setValue(self.grid.courant)
		self.timesteps.setValue(self.grid.timesteps)
		self.PMLlayers.setValue(self.grid.PMLlayers)
		self.cellsX.setValue(self.grid.cellsX)
		self.cellsY.setValue(self.grid.cellsY)
		self.cellsZ.setValue(self.grid.cellsZ)
		self.dx.setValue(amax(self.grid.dx)/self.grid.sizeScalefactor())
		self.dy.setValue(amax(self.grid.dy)/self.grid.sizeScalefactor())
		self.dz.setValue(amax(self.grid.dz)/self.grid.sizeScalefactor())
		self.dt.setValue(self.grid.dt/self.grid.timeScaleFactor())
		self.time.setValue(self.grid.time/self.grid.timeScaleFactor())

		self.sizeX.valueChanged.connect(self.sizeXGridChanged)
		self.sizeY.valueChanged.connect(self.sizeYGridChanged)
		self.sizeZ.valueChanged.connect(self.sizeZGridChanged)
		self.fmax.valueChanged.connect(self.fmaxGridChanged)
		self.courant.valueChanged.connect(self.courantGridChanged)
		self.timesteps.valueChanged.connect(self.timestepsGridChanged)
		self.PMLlayers.valueChanged.connect(self.PMLlayersGridChanged)
		self.cellsX.valueChanged.connect(self.cellsXGridChanged)
		self.cellsY.valueChanged.connect(self.cellsYGridChanged)
		self.cellsZ.valueChanged.connect(self.cellsZGridChanged)
		self.dx.valueChanged.connect(self.dxGridChanged)
		self.dy.valueChanged.connect(self.dyGridChanged)
		self.dz.valueChanged.connect(self.dzGridChanged)
		self.dt.valueChanged.connect(self.dtGridChanged)
		self.time.valueChanged.connect(self.timeGridChanged)

		self.sizeXLabel = QtGui.QLabel("Size x ["+self.grid.sizeUnit+"]:")
		self.sizeYLabel = QtGui.QLabel("Size y ["+self.grid.sizeUnit+"]:")
		self.sizeZLabel = QtGui.QLabel("Size z ["+self.grid.sizeUnit+"]:")
		self.fmaxLabel = QtGui.QLabel(u"f\u2098 ["+self.grid.freqUnit+"]:")
		self.courantLabel = QtGui.QLabel("Courant number:")
		self.timestepsLabel = QtGui.QLabel(u"n\u00B0 timesteps:")
		self.PMLLabel = QtGui.QLabel(u"n\u00B0 PML layers:")

		self.cellsXLabel = QtGui.QLabel("Cells x:")
		self.cellsYLabel = QtGui.QLabel("Cells y:")
		self.cellsZLabel = QtGui.QLabel("Cells z:")
		self.deltaXLabel = QtGui.QLabel(u"\u0394x ["+self.grid.sizeUnit+"]:")
		self.deltaYLabel = QtGui.QLabel(u"\u0394y ["+self.grid.sizeUnit+"]:")
		self.deltaZLabel = QtGui.QLabel(u"\u0394z ["+self.grid.sizeUnit+"]:")
		self.deltaTLabel = QtGui.QLabel(u"\u0394t ["+self.grid.freqTimeConversion[self.grid.freqUnit]+"]:")
		self.timeLabel = QtGui.QLabel("Time ["+self.grid.freqTimeConversion[self.grid.freqUnit]+"]:")

		self.freqUnitLabel = QtGui.QLabel("Frequency unit:")
		self.sizeUnitLabel = QtGui.QLabel("size unit:")
		self.freqUnitComboBox = QtGui.QComboBox()
		self.sizeUnitComboBox = QtGui.QComboBox()

		freqUnits = ('THz', 'GHz', 'MHz', 'kHz', 'Hz', 'mHz')
		self.freqUnitComboBox.addItems(freqUnits)
		for k in range(len(freqUnits)):
			if freqUnits[k] == self.grid.freqUnit:
				self.freqUnitComboBox.setCurrentIndex(k)
				break

		sizeUnits = ('km', 'm', 'mm', u"\u00B5m", 'nm')
		self.sizeUnitComboBox.addItems(sizeUnits)
		for k in range(len(sizeUnits)):
			if sizeUnits[k] == self.grid.sizeUnit:
				self.sizeUnitComboBox.setCurrentIndex(k)
				break

		self.freqUnitComboBox.currentIndexChanged.connect(self.freqUnitChanged)
		self.sizeUnitComboBox.currentIndexChanged.connect(self.sizeUnitChanged)

		self.projectNameLabel = QtGui.QLabel("Project name:")
		self.projectNameEntry = QtGui.QLineEdit(self.project.name)
		self.projectNameEntry.textChanged.connect(self.projectNameChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)
		horizontalLine = QtGui.QFrame()
		horizontalLine.setFrameStyle(QtGui.QFrame.HLine)
		horizontalLine.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)

		horizontalLine2 = QtGui.QFrame()
		horizontalLine2.setFrameStyle(QtGui.QFrame.HLine)
		horizontalLine2.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(2, 1)
		layout.setColumnStretch(4, 1)
		layout.addWidget(self.projectNameLabel, 0, 1, QtCore.Qt.AlignRight)
		layout.addWidget(self.projectNameEntry, 0, 2, 1, 3)

		layout.addWidget(horizontalLine,1,0,1,5)

		layout.addWidget(self.sizeXLabel, 2, 3, QtCore.Qt.AlignRight)
		layout.addWidget(self.sizeX, 2, 4)

		layout.addWidget(self.sizeYLabel, 3, 3, QtCore.Qt.AlignRight)
		layout.addWidget(self.sizeY, 3, 4)

		layout.addWidget(self.sizeZLabel, 4, 3, QtCore.Qt.AlignRight)
		layout.addWidget(self.sizeZ, 4, 4)

		layout.addWidget(self.fmaxLabel, 5, 3, QtCore.Qt.AlignRight)
		layout.addWidget(self.fmax, 5, 4)

		layout.addWidget(self.courantLabel, 6, 3, QtCore.Qt.AlignRight)
		layout.addWidget(self.courant, 6, 4)

		layout.addWidget(self.timestepsLabel, 7, 3, QtCore.Qt.AlignRight)
		layout.addWidget(self.timesteps, 7, 4)

		layout.addWidget(self.PMLLabel, 8, 3, QtCore.Qt.AlignRight)
		layout.addWidget(self.PMLlayers, 8, 4)

		layout.addWidget(self.cellsXLabel, 2, 1, QtCore.Qt.AlignRight)
		layout.addWidget(self.cellsX, 2, 2)

		layout.addWidget(self.cellsYLabel, 3, 1, QtCore.Qt.AlignRight)
		layout.addWidget(self.cellsY, 3, 2)

		layout.addWidget(self.cellsZLabel, 4, 1, QtCore.Qt.AlignRight)
		layout.addWidget(self.cellsZ, 4, 2)

		layout.addWidget(self.deltaXLabel, 5, 1, QtCore.Qt.AlignRight)
		layout.addWidget(self.dx, 5, 2)

		layout.addWidget(self.deltaYLabel, 6, 1, QtCore.Qt.AlignRight)
		layout.addWidget(self.dy, 6, 2)

		layout.addWidget(self.deltaZLabel, 7, 1, QtCore.Qt.AlignRight)
		layout.addWidget(self.dz, 7, 2)

		layout.addWidget(self.deltaTLabel, 8, 1, QtCore.Qt.AlignRight)
		layout.addWidget(self.dt, 8, 2)

		layout.addWidget(self.timeLabel, 9, 1, QtCore.Qt.AlignRight)
		layout.addWidget(self.time, 9, 2)

		layout.addWidget(horizontalLine2,10,0,1,5)

		layout.addWidget(self.freqUnitLabel, 11, 1, QtCore.Qt.AlignRight)
		layout.addWidget(self.freqUnitComboBox, 11, 2)

		layout.addWidget(self.sizeUnitLabel, 11, 3, QtCore.Qt.AlignRight)
		layout.addWidget(self.sizeUnitComboBox, 11, 4)

		layout.addItem(spacer, 12, 2)
		self.tab.setLayout(layout)

class Sensor():
	def __init__(self, P, name='Sensor'):
		self.P = P 			# Point
		self.name = name

class Block():
	def __init__(self, P1=(0,0,0), P2=(0,0,0), name='Block'):
		self.type = 'Block'
		self.name = name
		self.enable = True
		self.opacity = 99
		self.material = None

		self.x0 = min(P1[0], P2[0])		# Allow for oblique cubes in the future
		self.x1 = max(P1[0], P2[0])

		self.y0 = min(P1[1], P2[1])
		self.y1 = max(P1[1], P2[1])

		self.z0 = min(P1[2], P2[2])
		self.z1 = max(P1[2], P2[2])

	def changeMaterial(material):
		self.material = material

	def pointInGeometry(P):
		if P[0] >= self.x0 and P[0] <= self.x1 and \
			P[1] >= self.y0 and P[1] <= self.y1 and \
			P[2] >= self.z0 and P[2] <= self.z1:
			return True
		else:
			return False

class BlockWidget():
	def __init__(self, project, mayavi, view, model, block=None):
		self.project = project
		self.mayavi = mayavi
		self.view = view
		self.model = model
		self.block = block

		if block==None:
			self.addBlock()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def addBlockEnable(self, state):
		self.project.refreshPlot = True
		if state:
			self.block.enable = True
		else:
			self.block.enable = False

		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def opacityBlockChanged(self, position):
		self.project.refreshPlot = True
		self.block.opacity = position
		self.tab.opacity = position
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def x0BlockChanged(self, value):
		self.project.refreshPlot = True
		self.block.x0 = value*self.project.grid.sizeScalefactor()
		self.tab.x0 = value
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def x1BlockChanged(self, value):
		self.project.refreshPlot = True
		self.block.x1 = value*self.project.grid.sizeScalefactor()
		self.tab.x1 = value
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def y0BlockChanged(self, value):
		self.project.refreshPlot = True
		self.block.y0 = value*self.project.grid.sizeScalefactor()
		self.tab.y0 = value
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def y1BlockChanged(self, value):
		self.project.refreshPlot = True
		self.block.y1 = value*self.project.grid.sizeScalefactor()
		self.tab.y1 = value
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def z0BlockChanged(self, value):
		self.project.refreshPlot = True
		self.block.z0 = value*self.project.grid.sizeScalefactor()
		self.tab.z0 = value
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def z1BlockChanged(self, value):
		self.project.refreshPlot = True
		self.block.z1 = value*self.project.grid.sizeScalefactor()
		self.tab.z1 = value
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def nameBlockChanged(self, name):
		self.block.name = name 				# Change name of block
		self.currentNode.setText(name)		# And adjust the name in the treeview

	def materialBlockChanged(self, index):
		self.project.refreshPlot = True
		if not(self.hold):
			self.tab.material = self.project.materials[index]
			self.block.material = self.project.materials[index]
			self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def createWidget(self):
		self.name = self.block.name
		self.tab = QtGui.QWidget()
		x0 = QtGui.QDoubleSpinBox()
		x1 = QtGui.QDoubleSpinBox()
		y0 = QtGui.QDoubleSpinBox()
		y1 = QtGui.QDoubleSpinBox()
		z0 = QtGui.QDoubleSpinBox()
		z1 = QtGui.QDoubleSpinBox()

		maximum = 1e3
		minimum = -1e3
		x0.setRange(minimum, maximum)
		x1.setRange(minimum, maximum)
		y0.setRange(minimum, maximum)
		y1.setRange(minimum, maximum)
		z0.setRange(minimum, maximum)
		z1.setRange(minimum, maximum)
		x0.setDecimals(3)
		x1.setDecimals(3)
		y0.setDecimals(3)
		y1.setDecimals(3)
		z0.setDecimals(3)
		z1.setDecimals(3)

		grid = self.project.grid
		x0.setValue(self.block.x0/grid.sizeScalefactor())
		x1.setValue(self.block.x1/grid.sizeScalefactor())
		y0.setValue(self.block.y0/grid.sizeScalefactor())
		y1.setValue(self.block.y1/grid.sizeScalefactor())
		z0.setValue(self.block.z0/grid.sizeScalefactor())
		z1.setValue(self.block.z1/grid.sizeScalefactor())

		x0.valueChanged.connect(self.x0BlockChanged)
		x1.valueChanged.connect(self.x1BlockChanged)
		y0.valueChanged.connect(self.y0BlockChanged)
		y1.valueChanged.connect(self.y1BlockChanged)
		z0.valueChanged.connect(self.z0BlockChanged)
		z1.valueChanged.connect(self.z1BlockChanged)

		x0Label = QtGui.QLabel("x0 ["+self.project.grid.sizeUnit+"]:")
		x1Label = QtGui.QLabel("x1 ["+self.project.grid.sizeUnit+"]:")
		y0Label = QtGui.QLabel("y0 ["+self.project.grid.sizeUnit+"]:")
		y1Label = QtGui.QLabel("y1 ["+self.project.grid.sizeUnit+"]:")
		z0Label = QtGui.QLabel("z0 ["+self.project.grid.sizeUnit+"]:")
		z1Label = QtGui.QLabel("z1 ["+self.project.grid.sizeUnit+"]:")

		nameLabel = QtGui.QLabel("Name:")
		nameEntry = QtGui.QLineEdit(self.block.name)
		enableLabel = QtGui.QLabel("Enable:")
		enableCheckBox = QtGui.QCheckBox()

		if self.block.enable == True:
			enableCheckBox.setCheckState(QtCore.Qt.Checked)
		nameEntry.textChanged.connect(self.nameBlockChanged)

		materialLabel = QtGui.QLabel("Material:")
		opacityLabel = QtGui.QLabel("Opacity:")
		materialComboBox = QtGui.QComboBox()
		opacitySlider = QtGui.QSlider(QtCore.Qt.Horizontal)
		opacitySlider.setValue(self.block.opacity)
		materialComboBox.currentIndexChanged.connect(self.materialBlockChanged)

		enableCheckBox.stateChanged.connect(self.addBlockEnable)
		opacitySlider.sliderMoved.connect(self.opacityBlockChanged)

		self.hold = True
		index = -1
		for item in self.project.materials:
			materialComboBox.addItem(item.name)
			if self.block.material != None and self.block.material.name == item.name:
				index = materialComboBox.count()-1
		self.hold = False
		materialComboBox.setCurrentIndex(index)
		# x0Label.setFixedWidth(x0Label.minimumSizeHint().width())

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.setColumnStretch(3, 1)
		layout.addWidget(nameLabel, 0, 0, QtCore.Qt.AlignRight)
		layout.addWidget(nameEntry, 0, 1, 1, 3)

		layout.addWidget(enableLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(enableCheckBox, 1, 1)

		layout.addWidget(x0Label, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(x0, 2, 1)
		layout.addWidget(x1Label, 2, 2, QtCore.Qt.AlignRight)
		layout.addWidget(x1, 2, 3)

		layout.addWidget(y0Label, 3, 0, QtCore.Qt.AlignRight)
		layout.addWidget(y0, 3, 1)
		layout.addWidget(y1Label, 3, 2, QtCore.Qt.AlignRight)
		layout.addWidget(y1, 3, 3)

		layout.addWidget(z0Label, 4, 0, QtCore.Qt.AlignRight)
		layout.addWidget(z0, 4, 1)
		layout.addWidget(z1Label, 4, 2, QtCore.Qt.AlignRight)
		layout.addWidget(z1, 4, 3)

		layout.addWidget(materialLabel, 5, 0, QtCore.Qt.AlignRight)
		layout.addWidget(materialComboBox, 5, 1)

		layout.addWidget(opacityLabel, 6, 0, QtCore.Qt.AlignRight)
		layout.addWidget(opacitySlider, 6, 1)

		layout.addItem(spacer, 7, 2)
		self.tab.setLayout(layout)

		# self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def addBlock(self):
		name = 'Block'+str(len(self.project.geometries)+1)
		self.block = Block((0,0,0), (0,0,0), name)						# Create new block
		self.project.addGeometry(self.block)							# Add it to the geometry list
		projectNode = self.project.item
		projectNode.child(2).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(2).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(2).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(2).child(projectNode.child(2).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class Cylinder():
	def __init__(self, P1, P2, r, name='Cylinder', geometryType='Cylinder'):
		self.geometryType = geometryType
		self.name = name
		self.enable = True

		self.x0 = min(P1[0], P2[0])		
		self.x1 = max(P1[0], P2[0])

		self.y0 = min(P1[1], P2[1])
		self.y1 = max(P1[1], P2[1])

		self.z0 = min(P1[2], P2[2])
		self.z1 = max(P1[2], P2[2])

		self.r = abs(r)

	def pointInGeometry(P):
		return True

class Sphere():
	def __init__(self, P, r, name='Sphere', geometryType='Sphere'):
		self.geometryType = geometryType
		self.name = name
		self.enable = True

		self.x = P[0]
		self.y = P[1]
		self.z = P[2]
		self.r = abs(r)

	def pointInGeometry(P):
		return True

class Cone():
	def __init__(self, P1, P2, r, name='Cone', geometryType='Cone'):
		self.geometryType = geometryType
		self.name = name
		self.enable = True

		self.x0 = min(P1[0], P2[0])		
		self.x1 = max(P1[0], P2[0])

		self.y0 = min(P1[1], P2[1])
		self.y1 = max(P1[1], P2[1])

		self.z0 = min(P1[2], P2[2])
		self.z1 = max(P1[2], P2[2])

		self.r = abs(r)

	def pointInGeometry(P):
		return True

class Workplane():
	def __init__(self, name, P, angles, plotter2D, grid):
		self.type = 'Workplane'
		self.name = name
		self.P = P
		self.alpha, self.beta, self.gamma = angles
		self.enable = True
		self.plotter2D = plotter2D
		self.grid = grid
		self.geometries = []
		self.recompute = True	# (Re)compute rotation matrices

	def addPolygon2D(self, polygon2D):
		self.geometries.append(polygon2D)

	def addCircle(self, circle):
		self.geometries.append(circle)

	def addExtrusion(self, extrusion):
		self.geometries.append(extrusion)

	def setAlpha(self, value):
		self.alpha = value
		self.recompute = True

	def setBeta(self, value):
		self.beta = value
		self.recompute = True

	def setGamma(self, value):
		self.gamma = value
		self.recompute = True

	def computeRotationMatrices(self):
		alpha = self.alpha/180.0*pi
		beta = self.beta/180.0*pi
		gamma = self.gamma/180.0*pi

		Ralpha = [[1, 0, 0], [0, cos(alpha), -sin(alpha)], [0, sin(alpha), cos(alpha)]]
		Rbeta = [[cos(beta), 0, sin(beta)], [0, 1, 0], [-sin(beta), 0, cos(beta)]]
		Rgamma = [[cos(gamma), -sin(gamma), 0], [sin(gamma), cos(gamma), 0], [0, 0, 1]]
		self.R = dot(Ralpha, dot(Rbeta, Rgamma))

		Ralphainv = [[1, 0, 0], [0, cos(alpha), sin(alpha)], [0, -sin(alpha), cos(alpha)]]
		Rbetainv = [[cos(beta), 0, -sin(beta)], [0, 1, 0], [sin(beta), 0, cos(beta)]]
		Rgammainv = [[cos(gamma), sin(gamma), 0], [-sin(gamma), cos(gamma), 0], [0, 0, 1]]
		self.Rinv = dot(Rgammainv, dot(Rbetainv, Ralphainv))

	def uvwToxyz(self, points):		# Accepts an array of points
		if self.recompute == True:
			self.computeRotationMatrices()
			self.recompute = False

		for k in range(len(points)):
			points[k] = dot(self.R, points[k]) + self.P

		return points

	def xyzTouvw(self, point): 	# Accepts only a single point
		if self.recompute == True:
			self.computeRotationMatrices()
			self.recompute = False

		return dot(self.Rinv, array(point) - array(self.P))

	def updatePlot(self):
		self.plotter2D.updatePlot(self)

class WorkplaneWidget():
	def __init__(self, project, mayavi, plotter2D, view, model, workplane=None):
		self.project = project
		self.view = view
		self.workplane = workplane
		self.mayavi = mayavi
		self.plotter2D = plotter2D
		self.model = model
		self.grid = project.grid
		if workplane==None:
			self.addWorkplane()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def PxChanged(self, value):
		self.workplane.P[0] = value*self.grid.sizeScalefactor()
		self.project.refreshPlot = True
		self.mayavi.visualization.update_plot(self.project)

	def PyChanged(self, value):
		self.workplane.P[1] = value*self.grid.sizeScalefactor()
		self.project.refreshPlot = True
		self.mayavi.visualization.update_plot(self.project)

	def PzChanged(self, value):
		self.workplane.P[2] = value*self.grid.sizeScalefactor()
		self.project.refreshPlot = True
		self.mayavi.visualization.update_plot(self.project)

	def alphaChanged(self, value):
		self.workplane.setAlpha(value)
		self.project.refreshPlot = True
		self.mayavi.visualization.update_plot(self.project)

	def betaChanged(self, value):
		self.workplane.setBeta(value)
		self.project.refreshPlot = True
		self.mayavi.visualization.update_plot(self.project)

	def gammaChanged(self, value):
		self.workplane.setGamma(value)
		self.project.refreshPlot = True
		self.mayavi.visualization.update_plot(self.project)

	def nameWorkplaneChanged(self, name):
		self.workplane.name = name
		self.currentNode.setText(name)

	def createWidget(self):
		self.tab = QtGui.QWidget()
		Px = QtGui.QDoubleSpinBox()
		Py = QtGui.QDoubleSpinBox()
		Pz = QtGui.QDoubleSpinBox()
		alpha = QtGui.QDoubleSpinBox()
		beta = QtGui.QDoubleSpinBox()
		gamma = QtGui.QDoubleSpinBox()

		Px.valueChanged.connect(self.PxChanged)
		Py.valueChanged.connect(self.PyChanged)
		Pz.valueChanged.connect(self.PzChanged)
		alpha.valueChanged.connect(self.alphaChanged)
		beta.valueChanged.connect(self.betaChanged)
		gamma.valueChanged.connect(self.gammaChanged)

		maximum = 1e5
		minimum = -1e5
		Px.setRange(minimum, maximum)
		Py.setRange(minimum, maximum)
		Pz.setRange(minimum, maximum)
		alpha.setRange(-180, 180)
		beta.setRange(-180, 180)
		gamma.setRange(-180, 180)
		Px.setDecimals(3)
		Py.setDecimals(3)
		Pz.setDecimals(3)
		alpha.setDecimals(3)
		beta.setDecimals(3)
		gamma.setDecimals(3)

		Px.setValue(self.workplane.P[0]/self.grid.sizeScalefactor())
		Py.setValue(self.workplane.P[1]/self.grid.sizeScalefactor())
		Pz.setValue(self.workplane.P[2]/self.grid.sizeScalefactor())
		alpha.setValue(self.workplane.alpha)
		beta.setValue(self.workplane.beta)
		gamma.setValue(self.workplane.gamma)

		PxLabel = QtGui.QLabel("Px ["+self.grid.sizeUnit+"]:")
		PyLabel = QtGui.QLabel("Py ["+self.grid.sizeUnit+"]:")
		PzLabel = QtGui.QLabel("Pz ["+self.grid.sizeUnit+"]:")
		alphaLabel = QtGui.QLabel(u"\u03B1 [\u00B0]")
		betaLabel = QtGui.QLabel(u"\u03B2 [\u00B0]")
		gammaLabel = QtGui.QLabel(u"\u0263 [\u00B0]")

		nameLabel = QtGui.QLabel("Name:")
		nameEntry = QtGui.QLineEdit(self.workplane.name)
		nameEntry.textChanged.connect(self.nameWorkplaneChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.addWidget(nameLabel, 0, 0, QtCore.Qt.AlignRight)
		layout.addWidget(nameEntry, 0, 1, 1, 3)

		layout.addWidget(PxLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(Px, 1, 1)
		layout.addWidget(alphaLabel, 1, 2, QtCore.Qt.AlignRight)
		layout.addWidget(alpha, 1, 3)

		layout.addWidget(PyLabel, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(Py, 2, 1)
		layout.addWidget(betaLabel, 2, 2, QtCore.Qt.AlignRight)
		layout.addWidget(beta, 2, 3)

		layout.addWidget(PzLabel, 3, 0, QtCore.Qt.AlignRight)
		layout.addWidget(Pz, 3, 1)
		layout.addWidget(gammaLabel, 3, 2, QtCore.Qt.AlignRight)
		layout.addWidget(gamma, 3, 3)

		layout.addItem(spacer, 4, 1)
		self.tab.setLayout(layout)

	def addWorkplane(self):
		name = 'Workplane'+str(len(self.project.geometries)+1)
		self.workplane = Workplane(name, P=[0,0,0], angles=(0,0,0), plotter2D=self.plotter2D, grid=self.project.grid)
		self.project.addGeometry(self.workplane)				# Add it to the geometry list
		projectNode = self.project.item
		projectNode.child(2).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(2).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(2).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(2).child(projectNode.child(2).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class Polygon2D():
	def __init__(self, name, workplane):
		self.type = 'Polygon 2D'
		self.name = name
		self.workplane = workplane
		self.P = [[0,0]]
		self.material = None
		self.opacity = 99
		self.enable = True

	def addPoint(self, point):
		self.P.append(point)
		self.workplane.updatePlot()

	def addPointAt(self, point, index):
		self.P.insert(index, point)
		self.workplane.updatePlot()

	def deletePointAt(self, index):
		del self.P[index]

	def pointInGeometry(self, P):
		inside = False
		P = self.workplane.xyzTouvw(P)

		extrusions = (x for x in self.workplane.geometries if x.type == 'Extrusion')
		for extrusion in extrusions:
			if extrusion.height > P[2]:		# First check the height
				inside = True

		if inside==True:		# If the height is ok, then proceed in 2D
			i = 0
			j = len(self.P)-1
			c = False
			while i < len(self.P):
				if ((self.P[i][1] > P[1]) != (self.P[j][1] > P[1])) and \
					(P[0] < (self.P[j][0] - self.P[i][0])*(P[1] - self.P[i][1])/(self.P[j][1] - self.P[i][1]) + self.P[i][0]):
					c = not(c)
				j = i
				i += 1
			return c
		return False

	def returnIntersections(self, x, y): 	# scan line algorithm, returns intersection along z-line with given x and y
		minHeight = 0
		maxHeight = 0

class Polygon2DWidget():
	def __init__(self, project, workplaneIndex, view, model, polygon2D=None):
		self.project = project
		self.workplaneIndex = workplaneIndex
		self.workplane = self.project.geometries[self.workplaneIndex]
		self.view = view
		self.polygon2D = polygon2D
		self.model = model
		self.grid = self.project.grid
		self.hold = False

		if polygon2D == None:
			self.addPolygon2D()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def namePolygon2DChanged(self, name):
		self.polygon2D.name = name
		self.currentNode.setText(name)

	def addButtonClicked(self):
		selection = self.tableView.selectionModel()
		if selection.hasSelection() and len(selection.selectedIndexes()) > 0:
			pos = selection.selectedIndexes()
			self.polygon2D.addPointAt([0,0], pos[0].row()+1)
		else:
			self.polygon2D.addPoint([0,0])

		self.tableModel.updateView()

	def deleteButtonClicked(self):
		selection = self.tableView.selectionModel()
		if selection.hasSelection() and len(selection.selectedIndexes()) > 0:
			indexes = selection.selectedIndexes()
			rows = []
			for k in range(len(indexes)):
				rows.append(indexes[k].row())	# Can contain duplicates
			uniqueRows = sorted(set(rows), reverse=True)	# Sort reverse, because last row needs to be deleted first
			for row in uniqueRows:
				self.polygon2D.deletePointAt(row)

		self.tableModel.updateView()

	def materialPolygonChanged(self, index):
		if self.hold == False:
			self.polygon2D.material = self.project.materials[index]
			self.workplane.updatePlot()

	def opacityPolygonChanged(self, value):
		self.polygon2D.opacity = value
		self.workplane.updatePlot()

	def enablePolygon(self, state):
		self.polygon2D.enable = state

	def createWidget(self):
		self.tab = QtGui.QWidget()
		self.tableView = QtGui.QTableView()
		self.tableModel = MyTableModel(self.tab, self.polygon2D.P, \
			["u ["+self.grid.sizeUnit+"]", "v ["+self.grid.sizeUnit+"]"], self.grid)
		self.tableView.setModel(self.tableModel)
		self.tableView.resizeColumnsToContents()
		addButton = QtGui.QPushButton("+")
		deleteButton = QtGui.QPushButton("-")

		# self.tab.setMinimumWidth(0)
		nameLabel = QtGui.QLabel("Name:")
		enableLabel = QtGui.QLabel("Enable:")
		opacityLabel = QtGui.QLabel("Opacity:")
		materialLabel = QtGui.QLabel("Material:")
		nameEntry = QtGui.QLineEdit(self.polygon2D.name)
		nameEntry.resize(nameEntry.minimumSizeHint())

		materialComboBox = QtGui.QComboBox()
		opacitySlider = QtGui.QSlider(QtCore.Qt.Horizontal)
		opacitySlider.setValue(self.polygon2D.opacity)
		enableCheckBox = QtGui.QCheckBox()

		enableCheckBox.setChecked(self.polygon2D.enable)

		self.hold = True
		index = -1
		for item in self.project.materials:
			materialComboBox.addItem(item.name)
			if self.polygon2D.material != None and self.polygon2D.material.name == item.name:
				index = materialComboBox.count()-1
		self.hold = False
		materialComboBox.setCurrentIndex(index)

		nameEntry.textChanged.connect(self.namePolygon2DChanged)
		self.tableModel.dataChanged.connect(self.workplane.updatePlot)
		materialComboBox.currentIndexChanged.connect(self.materialPolygonChanged)
		enableCheckBox.stateChanged.connect(self.enablePolygon)
		opacitySlider.sliderMoved.connect(self.opacityPolygonChanged)
		addButton.clicked.connect(self.addButtonClicked)
		deleteButton.clicked.connect(self.deleteButtonClicked)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 2)
		layout.addWidget(nameLabel, 0, 0, QtCore.Qt.AlignRight)
		layout.addWidget(nameEntry, 0, 1, 1, 4)

		layout.addWidget(enableLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(enableCheckBox, 1, 1)

		layout.addWidget(self.tableView, 2, 1, 4, 4)
		layout.addWidget(addButton, 2, 5)
		layout.addWidget(deleteButton, 3, 5)

		layout.addWidget(opacityLabel, 6, 0, QtCore.Qt.AlignRight)
		layout.addWidget(opacitySlider, 6, 1)

		layout.addWidget(materialLabel, 7, 0, QtCore.Qt.AlignRight)
		layout.addWidget(materialComboBox, 7, 1)

		layout.addItem(spacer, 8, 1)
		self.tab.setLayout(layout)

	def addPolygon2D(self):
		name = 'Polygon'+str(len(self.project.geometries[self.workplaneIndex].geometries)+1)
		self.polygon2D = Polygon2D(name, self.workplane)
		self.workplane.addPolygon2D(self.polygon2D)				# Add it to the geometry list
		projectNode = self.project.item
		projectNode.child(2).child(self.workplaneIndex).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(2).child(self.workplaneIndex).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(2).child(self.workplaneIndex).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(2).child(self.workplaneIndex).child(projectNode.child(2).child(self.workplaneIndex).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class Circle():
	def __init__(self, name, workplane):
		self.type = 'Circle'
		self.name = name
		self.workplane = workplane
		self.center = [0,0]
		self.radius = 0
		self.material = None
		self.opacity = 99
		self.enable = True

	def pointInGeometry(self, P):
		inside = False
		P = self.workplane.xyzTouvw(P)

		extrusions = (x for x in self.workplane.geometries if x.type == 'Extrusion')
		for extrusion in extrusions:
			if extrusion.height > P[2]:		# First check the height
				inside = True

		if inside==True and (P[0]-self.center[0])**2 + (P[1]-self.center[1])**2 <= self.radius**2:
			return True
		return False

class CircleWidget():
	def __init__(self, project, workplaneIndex, view, model, circle=None):
		self.project = project
		self.workplaneIndex = workplaneIndex
		self.workplane = self.project.geometries[self.workplaneIndex]
		self.view = view
		self.circle = circle
		self.model = model
		self.grid = self.project.grid
		self.hold = False

		if circle == None:
			self.addCircle()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def nameCircleChanged(self, name):
		self.circle.name = name
		self.currentNode.setText(name)

	def centerUChanged(self, value):
		self.circle.center[0] = value*self.grid.sizeScalefactor()
		self.workplane.updatePlot()

	def centerVChanged(self, value):
		self.circle.center[1] = value*self.grid.sizeScalefactor()
		self.workplane.updatePlot()

	def radiusChanged(self, value):
		self.circle.radius = value*self.grid.sizeScalefactor()
		self.workplane.updatePlot()

	def materialCircleChanged(self, index):
		if self.hold == False:
			self.circle.material = self.project.materials[index]
			self.workplane.updatePlot()

	def opacityCircleChanged(self, value):
		self.circle.opacity = value
		self.workplane.updatePlot()

	def enableCircle(self, state):
		if state:
			self.circle.enable = True
		else:
			self.circle.enable = False

	def createWidget(self):
		self.tab = QtGui.QWidget()

		centerU = QtGui.QDoubleSpinBox()
		centerV = QtGui.QDoubleSpinBox()
		radius = QtGui.QDoubleSpinBox()

		centerU.setRange(-1E4,1E4)
		centerV.setRange(-1E4,1E4)
		radius.setRange(0,1E4)
		centerU.setDecimals(3)
		centerV.setDecimals(3)
		radius.setDecimals(3)

		centerU.setValue(self.circle.center[0]/self.grid.sizeScalefactor())
		centerV.setValue(self.circle.center[1]/self.grid.sizeScalefactor())
		radius.setValue(self.circle.radius/self.grid.sizeScalefactor())

		nameLabel = QtGui.QLabel("Name:")
		enableLabel = QtGui.QLabel("Enable:")
		opacityLabel = QtGui.QLabel("Opacity:")
		materialLabel = QtGui.QLabel("Material:")
		centerULabel = QtGui.QLabel("Center u ["+self.grid.sizeUnit+"]:")
		centerVLabel = QtGui.QLabel("Center v ["+self.grid.sizeUnit+"]:")
		radiusLabel = QtGui.QLabel("Radius ["+self.grid.sizeUnit+"]:")
		nameEntry = QtGui.QLineEdit(self.circle.name)
		nameEntry.resize(nameEntry.minimumSizeHint())

		materialComboBox = QtGui.QComboBox()
		opacitySlider = QtGui.QSlider(QtCore.Qt.Horizontal)
		opacitySlider.setValue(self.circle.opacity)
		enableCheckBox = QtGui.QCheckBox()

		enableCheckBox.setChecked(self.circle.enable)

		self.hold = True
		index = -1
		for item in self.project.materials:
			materialComboBox.addItem(item.name)
			if self.circle.material != None and self.circle.material.name == item.name:
				index = materialComboBox.count()-1
		self.hold = False
		materialComboBox.setCurrentIndex(index)

		nameEntry.textChanged.connect(self.nameCircleChanged)
		centerU.valueChanged.connect(self.centerUChanged)
		centerV.valueChanged.connect(self.centerVChanged)
		radius.valueChanged.connect(self.radiusChanged)
		materialComboBox.currentIndexChanged.connect(self.materialCircleChanged)
		enableCheckBox.stateChanged.connect(self.enableCircle)
		opacitySlider.sliderMoved.connect(self.opacityCircleChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 2)
		layout.addWidget(nameLabel, 0, 0, QtCore.Qt.AlignRight)
		layout.addWidget(nameEntry, 0, 1, 1, 4)

		layout.addWidget(enableLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(enableCheckBox, 1, 1)

		layout.addWidget(centerULabel, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(centerU, 2, 1)

		layout.addWidget(centerVLabel, 3, 0, QtCore.Qt.AlignRight)
		layout.addWidget(centerV, 3, 1)

		layout.addWidget(radiusLabel, 4, 0, QtCore.Qt.AlignRight)
		layout.addWidget(radius, 4, 1)

		layout.addWidget(opacityLabel, 5, 0, QtCore.Qt.AlignRight)
		layout.addWidget(opacitySlider, 5, 1)

		layout.addWidget(materialLabel, 6, 0, QtCore.Qt.AlignRight)
		layout.addWidget(materialComboBox, 6, 1)

		layout.addItem(spacer, 7, 1)
		self.tab.setLayout(layout)

	def addCircle(self):
		name = 'Circle'+str(len(self.project.geometries[self.workplaneIndex].geometries)+1)
		self.circle = Circle(name, self.workplane)
		self.workplane.addCircle(self.circle)				# Add it to the geometry list
		projectNode = self.project.item
		projectNode.child(2).child(self.workplaneIndex).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(2).child(self.workplaneIndex).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(2).child(self.workplaneIndex).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(2).child(self.workplaneIndex).child(projectNode.child(2).child(self.workplaneIndex).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class Extrusion():
	def __init__(self, name, workplane):
		self.type = 'Extrusion'
		self.name = name
		self.workplane = workplane
		self.height = 0
		self.enable = True

	def adjustHeight(value):
		self.height = value

class ExtrusionWidget():
	def __init__(self, project, workplaneIndex, view, model, extrusion=None):
		self.project = project
		self.workplaneIndex = workplaneIndex
		self.workplane = self.project.geometries[self.workplaneIndex]
		self.view = view
		self.extrusion = extrusion
		self.model = model
		self.grid = self.project.grid

		if extrusion == None:
			self.addExtrusion()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def nameExtrusionChanged(self, name):
		self.extrusion.name = name
		self.currentNode.setText(name)

	def heightChanged(self, value):
		self.extrusion.height = value*self.grid.sizeScalefactor()

	def enableChanged(self, state):
		if state:
			self.extrusion.enable = True
		else:
			self.extrusion.enable = False

	def createWidget(self):
		self.tab = QtGui.QWidget()
		
		nameLabel = QtGui.QLabel("Name:")
		enableLabel = QtGui.QLabel("Enable:")
		heightLabel = QtGui.QLabel("Height ["+self.grid.sizeUnit+"]:")
		nameEntry = QtGui.QLineEdit(self.extrusion.name)
		height = QtGui.QDoubleSpinBox()
		enable = QtGui.QCheckBox()
		enable.setChecked(self.extrusion.enable)
		nameEntry.resize(nameEntry.minimumSizeHint())
		height.setRange(-1e5, 1e5)
		height.setDecimals(3)

		height.setValue(self.extrusion.height/self.grid.sizeScalefactor())

		nameEntry.textChanged.connect(self.nameExtrusionChanged)
		height.valueChanged.connect(self.heightChanged)
		enable.stateChanged.connect(self.enableChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 2)
		layout.addWidget(nameLabel, 0, 0, QtCore.Qt.AlignRight)
		layout.addWidget(nameEntry, 0, 1, 1, 2)

		layout.addWidget(heightLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(height, 1, 1)

		layout.addWidget(enableLabel, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(enable, 2, 1)

		layout.addItem(spacer, 3, 1)
		self.tab.setLayout(layout)

	def addExtrusion(self):
		name = 'Extrusion'+str(len(self.project.geometries[self.workplaneIndex].geometries)+1)
		self.extrusion = Extrusion(name, self.workplane)
		self.workplane.addExtrusion(self.extrusion)				# Add it to the geometry list
		projectNode = self.project.item
		projectNode.child(2).child(self.workplaneIndex).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(2).child(self.workplaneIndex).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(2).child(self.workplaneIndex).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(2).child(self.workplaneIndex).child(projectNode.child(2).child(self.workplaneIndex).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class Material():
	def __init__(self, name='Copper', mur=1, epsr=1, sigma=5.8e7, color=(50,50,250)):
		self.name = name
		self.mur = mur
		self.epsr = epsr
		self.sigma = sigma
		self.color = color

class MaterialWidget():
	def __init__(self, project, mayavi, view, model, material=None):
		self.project = project
		self.view = view
		self.material = material
		self.mayavi = mayavi
		self.model = model
		if material==None:
			self.addMaterial()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def epsrMaterialChanged(self, value):
		self.material.epsr = value

	def murMaterialChanged(self, value):
		self.material.mur = value 

	def sigmaMaterialChanged(self, value):
		self.material.sigma = value

	def nameMaterialChanged(self, name):
		self.material.name = name 		# Change name of material
		self.currentNode.setText(name)		# And adjust the name in the treeview

	def createWidget(self):
		self.tab = QtGui.QWidget()
		epsr = QtGui.QDoubleSpinBox()
		mur = QtGui.QDoubleSpinBox()
		sigma = QtGui.QDoubleSpinBox()

		maximum = 1e10
		minimum = 0
		epsr.setRange(minimum, maximum)
		mur.setRange(minimum, maximum)
		sigma.setRange(minimum, maximum)
		epsr.setDecimals(3)
		mur.setDecimals(3)
		sigma.setDecimals(3)

		epsr.setValue(self.material.epsr)
		mur.setValue(self.material.mur)
		sigma.setValue(self.material.sigma)

		epsr.valueChanged.connect(self.epsrMaterialChanged)
		mur.valueChanged.connect(self.murMaterialChanged)
		sigma.valueChanged.connect(self.sigmaMaterialChanged)

		epsrLabel = QtGui.QLabel(u"\u03B5\u1D63")
		murLabel = QtGui.QLabel(u"\u03BC\u1D63")
		sigmaLabel = QtGui.QLabel(u"\u03C3")

		nameLabel = QtGui.QLabel("Name:")
		nameEntry = QtGui.QLineEdit(self.material.name)
		nameEntry.textChanged.connect(self.nameMaterialChanged)

		colorLabel = QtGui.QLabel("Color:")
		self.colorButton = QtGui.QPushButton("")
		self.colorButton.clicked.connect(self.pickColor)
		self.colorButton.setMaximumWidth(30)
		self.colorButton.setMaximumHeight(20)
		color = QtGui.QColor(self.material.color[0], self.material.color[1], self.material.color[2])
		self.colorButton.setStyleSheet('background-color: %s;' %color.name())

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.addWidget(nameLabel, 0, 0, QtCore.Qt.AlignRight)
		layout.addWidget(nameEntry, 0, 1)

		layout.addWidget(epsrLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(epsr, 1, 1)

		layout.addWidget(murLabel, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(mur, 2, 1)

		layout.addWidget(sigmaLabel, 3, 0, QtCore.Qt.AlignRight)
		layout.addWidget(sigma, 3, 1)

		layout.addWidget(colorLabel, 4, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.colorButton, 4, 1)

		layout.addItem(spacer, 5, 1)
		self.tab.setLayout(layout)

	def pickColor(self):
		colorPicker = QtGui.QColorDialog()
		if colorPicker.exec_():
			self.colorButton.setStyleSheet('background-color: %s;' %colorPicker.currentColor().name())
			self.tab.color = colorPicker.currentColor().getRgb()[0:3]
			self.material.color = self.tab.color
			self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def addMaterial(self):
		name = 'Material'+str(len(self.project.materials)+1)
		self.material = Material(name, mur=1, epsr=1, sigma=0, color=(50, 50, 250))
		self.project.addMaterial(self.material)				# Add it to the material list
		projectNode = self.project.item
		projectNode.child(3).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(3).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(3).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(3).child(projectNode.child(3).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class PortPlane():
	def __init__(self, name, project, pos, plane, plotter2D):
		self.type = 'Port Plane'
		self.name = name
		self.pos = pos
		self.plane = plane 		# xy, yz or xz
		self.plotter2D = plotter2D
		self.project = project
		self.grid = project.grid
		self.sources = []

	def addSource(self, source):
		self.sources.append(source)

	def updatePlot(self):
		self.plotter2D.updatePlot(self)

class PortPlaneWidget():
	def __init__(self, project, mayavi, plotter2D, view, model, portPlane=None):
		self.project = project
		self.view = view
		self.portPlane = portPlane
		self.mayavi = mayavi
		self.plotter2D = plotter2D
		self.model = model
		self.grid = project.grid
		if portPlane==None:
			self.addPortPlane()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def planeChanged(self, index):
		if index == 0:
			self.portPlane.plane = 'xy'
		elif index == 1:
			self.portPlane.plane = 'yz'
		elif index == 2:
			self.portPlane.plane = 'zx'
		self.project.refreshPlot = True
		self.mayavi.visualization.update_plot(self.project)

	def posChanged(self, value):
		self.portPlane.pos = value*self.grid.sizeScalefactor()
		self.project.refreshPlot = True
		self.mayavi.visualization.update_plot(self.project)

	def namePortPlaneChanged(self, name):
		self.portPlane.name = name
		self.currentNode.setText(name)

	def createWidget(self):
		self.tab = QtGui.QWidget()
		name = QtGui.QLineEdit(self.project.name)
		pos = QtGui.QDoubleSpinBox()
		plane = QtGui.QComboBox()

		plane.addItems(('xy-plane', 'yz-plane', 'zx-plane'))
		nameLabel = QtGui.QLabel("Name:")
		posLabel = QtGui.QLabel("Position ["+self.grid.sizeUnit+"]:")
		planeLabel = QtGui.QLabel("Plane:")

		pos.setRange(-1e4, 1e4)
		pos.setDecimals(3)

		name.setText(self.portPlane.name)
		pos.setValue(self.portPlane.pos/self.grid.sizeScalefactor())
		plane.setCurrentIndex(ord(self.portPlane.plane[0])-ord('x'))

		name.textChanged.connect(self.namePortPlaneChanged)
		pos.valueChanged.connect(self.posChanged)
		plane.currentIndexChanged.connect(self.planeChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.addWidget(nameLabel, 0, 0, QtCore.Qt.AlignRight)
		layout.addWidget(name, 0, 1, 1, 3)

		layout.addWidget(posLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(pos, 1, 1)

		layout.addWidget(planeLabel, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(plane, 2, 1)

		layout.addItem(spacer, 3, 1)
		self.tab.setLayout(layout)
		self.mayavi.visualization.update_plot(self.project)

	def addPortPlane(self):
		name = 'Port plane'+str(len(self.project.excitations)+1)
		self.portPlane = PortPlane(name, self.project, 0, 'xy', self.plotter2D)
		self.project.addExcitation(self.portPlane)						# Add it to the excitation list
		projectNode = self.project.item
		projectNode.child(4).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(4).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(4).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(4).child(projectNode.child(4).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class WaveguidePort():
	def __init__(self, name, pos, portPlane):
		self.name = name
		self.type = 'Waveguide Port'
		self.pos = pos 	# ((x0, x1), (y0, y1))
		self.portPlane = portPlane
		self.field = 'Ex'
		self.frequency = 1E9 
		self.modulation = 'cos(2*pi*f*t)'
		self.modes = []		# Contains the mode
		self.nmodes = 4		# Number of modes to compute
		self.mode = 0		# mode picker
		self.grid = self.portPlane.project.grid
		self.w = []		# Eigenvalues	
		self.v = [] 	# Eigenvectors

		self.Ex = []
		self.Ey = []
		self.Hz = []

	def assignMaterials(self):
		workplanes = (x for x in self.portPlane.project.geometries if x.type == 'Workplane')
		for workplane in workplanes:
			extrusions = (x for x in workplane.geometries if x.type == 'Extrusion')
			for extrusion in extrusions:
				geometries = (x for x in workplane.geometries if (x.type != 'Extrusion') and (x.material != None))
				for geometry in geometries:
					if geometry.type == 'Circle':
						resolution = 100
						center = geometry.center
						radius = geometry.radius
						dalpha = 2*pi/resolution
						points = [[radius*cos(k*dalpha)+center[0], radius*sin(k*dalpha)+center[1]] for k in range(resolution)]

					elif geometry.type == 'Polygon 2D':
						points = geometry.P

					length = len(points)
					bottom = workplane.uvwToxyz([[points[k][0], points[k][1], 0] for k in range(length)])
					top = workplane.uvwToxyz([[points[k][0], points[k][1], extrusion.height] for k in range(length)])

					cutplane = Cutplane(orientation=self.portPlane.plane)
					intersection = cutplane.computeIntersection(self.portPlane.pos, bottom, top)

					if self.portPlane.plane == 'xy':
						for j in range(self.jmin, self.jmax):
							crossings = cutplane.rayCasting(intersection, self.grid.distanceY[j])

							for k in range(len(crossings)-1):
								i1 = (int)(crossings[k]/amax(self.grid.dx))
								i2 = (int)(crossings[k+1]/amax(self.grid.dx))
								for i in range(max(self.imin, i1), min(self.imax, i2)):
									index = (i-self.imin)*(self.sizeY-1)+j-self.jmin
									self.epsX[index, index] = geometry.material.epsr - 1j*geometry.material.sigma/(2*pi*self.frequency*self.grid.eps0)

									index = (i-self.imin)*self.sizeY+j-self.jmin
									self.epsY[index, index] = geometry.material.epsr - 1j*geometry.material.sigma/(2*pi*self.frequency*self.grid.eps0)

									index = (i-self.imin)*(self.sizeY-1)+j-self.jmin
									self.epsZinv[index, index] = 1.0/(geometry.material.epsr - 1j*geometry.material.sigma/(2*pi*self.frequency*self.grid.eps0))

	def calculateModes(self):
		self.k0 = 2*pi*self.frequency/self.grid.c

		if self.portPlane.plane == 'xy':
			a = (int)(self.grid.cellsX*self.pos[0][0]/self.grid.sizeX)
			b = (int)(self.grid.cellsX*self.pos[0][1]/self.grid.sizeX)
			self.imin = min(a, b)
			self.imax = max(a, b)
			a = (int)(self.grid.cellsY*self.pos[1][0]/self.grid.sizeY)
			b = (int)(self.grid.cellsY*self.pos[1][1]/self.grid.sizeY)
			self.jmin = min(a, b)
			self.jmax = max(a, b)
			self.dx = self.grid.dx[self.imin:self.imax]
			self.dy = self.grid.dy[self.jmin:self.jmax]
			self.distanceX = self.grid.distanceX
			self.distanceY = self.grid.distanceY

		self.sizeX = self.imax - self.imin
		self.sizeY = self.jmax - self.jmin


		self.epsX = sparse.eye((self.sizeX)*(self.sizeY-1))
		self.epsY = sparse.eye((self.sizeX-1)*self.sizeY)
		self.epsZinv = sparse.eye((self.sizeX-1)*(self.sizeY-1))

		self.epsX = sparse.csr_matrix(self.epsX, dtype=complex)
		self.epsY = sparse.csr_matrix(self.epsY, dtype=complex)
		self.epsZinv = sparse.csr_matrix(self.epsZinv, dtype=complex)

		self.assignMaterials()

		diag1 = 1.0/self.dx
		diag2 = -1.0/self.dx
		DAx = sparse.spdiags([diag1, diag2], array([0, -1]), self.sizeX, self.sizeX-1)

		diag1 = 1.0/self.dy
		diag2 = -1.0/self.dy
		DAy = sparse.spdiags([diag1, diag2], array([0, -1]), self.sizeY, self.sizeY-1)

		diag1 = 1.0/self.dx
		diag2 = -1.0/self.dx
		DBx = sparse.spdiags([diag1, diag2], array([0, -1]), self.sizeX, self.sizeX-1)

		diag1 = 1.0/self.dy
		diag2 = -1.0/self.dy
		DBy = sparse.spdiags([diag1, diag2], array([0, -1]), self.sizeY, self.sizeY-1)

		diag1 = 1.0/self.dx
		diag2 = -1.0/self.dx
		DCx = sparse.spdiags([diag1, diag2], array([1, 0]), self.sizeX-1, self.sizeX)

		diag1 = 1.0/self.dy
		diag2 = -1.0/self.dy
		DCy = sparse.spdiags([diag1, diag2], array([1, 0]), self.sizeY-1, self.sizeY)

		diag1 = 1.0/self.dx
		diag2 = -1.0/self.dx
		DDx = sparse.spdiags([diag1, diag2], array([1, 0]), self.sizeX-1, self.sizeX)

		diag1 = 1.0/self.dy
		diag2 = -1.0/self.dy
		DDy = sparse.spdiags([diag1, diag2], array([1, 0]), self.sizeY-1, self.sizeY)

		Ay = sparse.kron(sparse.eye(self.sizeX-1), DAy)
		Ax = sparse.kron(DAx, sparse.eye(self.sizeY-1))
		By = sparse.kron(sparse.eye(self.sizeX), DBy)
		Bx = sparse.kron(DBx, sparse.eye(self.sizeY))

		Cy = sparse.kron(sparse.eye(self.sizeX), DCy)
		Cx = sparse.kron(DCx, sparse.eye(self.sizeY))
		Dy = sparse.kron(sparse.eye(self.sizeX-1), DDy)
		Dx = sparse.kron(DDx, sparse.eye(self.sizeY-1))

		Qxx = self.k0**2*self.epsX + Cy*By + self.k0**(-2)*Ax*self.epsZinv*(Dx*Cy*By - Dy*Cx*By + self.k0**2*Dx*self.epsX)
		Qyy = self.k0**2*self.epsY + Cx*Bx + self.k0**(-2)*Ay*self.epsZinv*(-Dx*Cy*Bx + Dy*Cx*Bx + self.k0**2*Dy*self.epsY)
		Qxy = -Cy*Bx + self.k0**(-2)*Ax*self.epsZinv*(-Dx*Cy*Bx + Dy*Cx*Bx + self.k0**2*Dy*self.epsY)
		Qyx = -Cx*By + self.k0**(-2)*Ay*self.epsZinv*(Dx*Cy*By - Dy*Cx*By + self.k0**2*Dx*self.epsX)

		Q = sparse.bmat([[Qxx, Qxy], [Qyx, Qyy]])
		Q = sparse.csc_matrix(Q)

		self.w, self.v = eigs(Q, sigma=self.k0**2, which='LM', k=self.nmodes)

		# T = Q*self.v[:,0]
		# print self.v[10,0], T[10]

		r = real(self.w > 0)*2-1		# Real part needs to be > 0
		i = imag(self.w > 0)*2-1		# And imag part < 0
		self.w = real(self.w*r) - 1j*imag(self.w*i)

		print self.w

		ex = real(self.v[0:(self.sizeX)*(self.sizeY-1), :])
		ey = real(self.v[(self.sizeX)*(self.sizeY-1):self.sizeX*(self.sizeY-1) + (self.sizeX-1)*self.sizeY, :])
		ez = -1j*self.epsZinv*((-Dx*Cy*By + Dy*Cx*By - self.k0**2*Dx*self.epsX)*ex + (Dx*Cy*Bx - Dy*Cx*Bx - self.k0**2*Dy*self.epsY)*ey)/(self.k0**2*sqrt(self.w))

		hx = 1j*(-1j*sqrt(self.w)*ey + Ay*ez)/(2*pi*self.frequency*self.grid.mu0)
		hy = 1j*(1j*sqrt(self.w)*ex - Ax*ez)/(2*pi*self.frequency*self.grid.mu0)
		hz = 1j*(-By*ex + Bx*ey)/(2*pi*self.frequency*self.grid.mu0)

		io.savemat('data', {'Ay':Ay, 'Ax':Ax, 'By':By, 'Bx':Bx, 'Cy':Cy, 'Cx':Cx, 'Dy':Dy, 'Dx':Dx, 'Q':Q})

		# ex2 = (-1j*sqrt(self.w)*hy + Cy*hz)/(1j*2*pi*self.frequency*self.grid.eps0)
		# ey2 = (1j*sqrt(self.w)*hx - Cx*hz)/(1j*2*pi*self.frequency*self.grid.eps0)
		# err = ex-ex2
		# err2 = ey-ey2
		# for k in range(4):
		# 	print max(err[:, k]/ex2[:,k]), max(err2[:, k]/ey2[:,k])

		# for k in range(4):
		# 	print amax(abs(ey[:,k]))/amax(abs(hx[:,k])), amax(abs(ex[:,k]))/amax(abs(hy[:,k]))

		# print angle(ex[100][0])*180/pi, angle(ex[200][0])*180/pi
		# print angle(ey[100][0])*180/pi, angle(ey[200][0])*180/pi
		# print angle(ez[100][0])*180/pi, angle(ez[200][0])*180/pi
		# print angle(hx[100][0])*180/pi, angle(hx[200][0])*180/pi
		# print angle(hy[100][0])*180/pi, angle(hy[200][0])*180/pi
		# print angle(hz[100][0])*180/pi, angle(hz[200][0])*180/pi

		self.Ex = []
		self.Ey = []
		self.Ez = []
		self.Hx = []
		self.Hy = []
		self.Hz = []
		for n in range(self.nmodes):
			self.Ex.append(reshape(real(ex[:, n]), (self.sizeX, self.sizeY-1)).astype('f4'))
			self.Ey.append(reshape(real(ey[:, n]), (self.sizeX-1, self.sizeY)).astype('f4'))
			self.Ez.append(reshape(imag(ez[:, n]), (self.sizeX-1, self.sizeY-1)).astype('f4'))
			self.Hx.append(reshape(real(hx[:, n]), (self.sizeX-1, self.sizeY)).astype('f4'))
			self.Hy.append(reshape(real(hy[:, n]), (self.sizeX, self.sizeY-1)).astype('f4'))
			self.Hz.append(reshape(imag(hz[:, n]), (self.sizeX, self.sizeY)).astype('f4'))

		self.portPlane.updatePlot()

class WaveguideWidget():
	def __init__(self, project, view, model, excitationIndex, waveguidePort=None):
		self.project = project
		self.view = view
		self.model = model
		self.excitationIndex = excitationIndex
		self.portPlane = project.excitations[excitationIndex]
		self.waveguidePort = waveguidePort
		self.grid = self.project.grid

		if waveguidePort == None:
			self.addWaveguidePort()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def nameChanged(self, name):
		self.waveguidePort.name = name
		self.currentNode.setText(name)

	def modulationChanged(self):
		self.waveguidePort.modulation = self.modulation.toPlainText()

	def showFieldChanged(self, index):
		if index == 0:
			self.waveguidePort.field = 'Ex'
		elif index == 1:
			self.waveguidePort.field = 'Ey'
		elif index == 2:
			self.waveguidePort.field = 'Ez'
		self.portPlane.updatePlot()

	def frequencyChanged(self, value):
		self.waveguidePort.frequency = value*self.grid.freqScalefactor()

	def u0Changed(self, value):
		self.waveguidePort.pos[0][0] = value*self.grid.sizeScalefactor()
		self.portPlane.updatePlot()

	def u1Changed(self, value):
		self.waveguidePort.pos[0][1] = value*self.grid.sizeScalefactor()
		self.portPlane.updatePlot()

	def v0Changed(self, value):
		self.waveguidePort.pos[1][0] = value*self.grid.sizeScalefactor()
		self.portPlane.updatePlot()

	def v1Changed(self, value):
		self.waveguidePort.pos[1][1] = value*self.grid.sizeScalefactor()
		self.portPlane.updatePlot()

	def nmodesChanged(self, value):
		self.waveguidePort.nmodes = (int)(value)
		self.portPlane.updatePlot()

	def modeSelectorChanged(self, value):
		self.waveguidePort.mode = (int)(value)
		if len(self.waveguidePort.w) != 0:
			self.neffLabel.setText("neff: "+str(real(sqrt(self.waveguidePort.w[value]))/self.waveguidePort.k0)+ \
				" + j"+str(imag(sqrt(self.waveguidePort.w[value]))/self.waveguidePort.k0))
		self.portPlane.updatePlot()

	def calculateModesClicked(self):
		self.modeSelector.setRange(0, self.waveguidePort.nmodes-1)
		self.waveguidePort.calculateModes()
		self.portPlane.updatePlot()

	def createWidget(self):
		self.tab = QtGui.QWidget()
		name = QtGui.QLineEdit(self.waveguidePort.name)
		self.frequency = QtGui.QDoubleSpinBox()
		self.modulation = QtGui.QPlainTextEdit()
		self.showField = QtGui.QComboBox()
		self.u0 = QtGui.QDoubleSpinBox()
		self.u1 = QtGui.QDoubleSpinBox()
		self.v0 = QtGui.QDoubleSpinBox()
		self.v1 = QtGui.QDoubleSpinBox()
		self.nmodes = QtGui.QDoubleSpinBox()
		self.modeSelector = QtGui.QSlider(QtCore.Qt.Horizontal)
		self.calculateModes = QtGui.QPushButton("Calculate modes")

		nameLabel = QtGui.QLabel("Name:")
		frequencyLabel = QtGui.QLabel("Frequency ["+self.grid.freqUnit+"]:")
		modulationLabel = QtGui.QLabel("Modulation:")
		showFieldLabel = QtGui.QLabel("Show field:")
		u0Label = QtGui.QLabel(u"u\u2080 ["+self.grid.sizeUnit+"]:")
		u1Label = QtGui.QLabel(u"u\u2081 ["+self.grid.sizeUnit+"]:")
		v0Label = QtGui.QLabel(u"v\u2080 ["+self.grid.sizeUnit+"]:")
		v1Label = QtGui.QLabel(u"v\u2081 ["+self.grid.sizeUnit+"]:")
		nmodeLabel = QtGui.QLabel(u"n\u00B0 of modes:")
		modeSelectorLabel = QtGui.QLabel("mode:")
		self.neffLabel = QtGui.QLabel("neff: ")
		self.showField.addItems(('Ex', 'Ey', 'Ez'))

		self.frequency.setRange(0, 1E15)
		self.nmodes.setRange(0, 100)
		self.u0.setRange(-1E4, 1E4)
		self.u1.setRange(-1E4, 1E4)
		self.v0.setRange(-1E4, 1E4)
		self.v1.setRange(-1E4, 1E4)
		self.u0.setDecimals(3)
		self.u1.setDecimals(3)
		self.v0.setDecimals(3)
		self.v1.setDecimals(3)
		self.frequency.setDecimals(3)
		self.nmodes.setDecimals(0)
		self.modeSelector.setRange(0, self.waveguidePort.nmodes-1)

		font = self.modulation.document().defaultFont()
		fontMetrics = QtGui.QFontMetrics(font)
		textSize = fontMetrics.size(0, "Hz")
		textHeight = textSize.height()+11		# Constant might need some tweaking
		self.modulation.setMaximumHeight(textHeight)
		self.modulation.setMinimumWidth(100)

		self.frequency.setValue(self.waveguidePort.frequency/self.grid.freqScalefactor())
		self.modulation.setPlainText(self.waveguidePort.modulation)
		self.showField.setCurrentIndex(ord(self.waveguidePort.field[1])-ord('x'))
		self.u0.setValue(self.waveguidePort.pos[0][0]/self.grid.sizeScalefactor())
		self.u1.setValue(self.waveguidePort.pos[0][1]/self.grid.sizeScalefactor())
		self.v0.setValue(self.waveguidePort.pos[1][0]/self.grid.sizeScalefactor())
		self.v1.setValue(self.waveguidePort.pos[1][1]/self.grid.sizeScalefactor())
		self.nmodes.setValue(self.waveguidePort.nmodes)
		self.modeSelector.setValue(self.waveguidePort.mode)

		self.u0.valueChanged.connect(self.u0Changed)
		self.u1.valueChanged.connect(self.u1Changed)
		self.v0.valueChanged.connect(self.v0Changed)
		self.v1.valueChanged.connect(self.v1Changed)
		self.frequency.valueChanged.connect(self.frequencyChanged)
		self.nmodes.valueChanged.connect(self.nmodesChanged)
		self.modeSelector.valueChanged.connect(self.modeSelectorChanged)
		self.calculateModes.clicked.connect(self.calculateModesClicked)
		self.modulation.textChanged.connect(self.modulationChanged)
		self.showField.currentIndexChanged.connect(self.showFieldChanged)
		name.textChanged.connect(self.nameChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.addWidget(nameLabel, 0, 0, QtCore.Qt.AlignRight)
		layout.addWidget(name, 0, 1, 1, 3)

		layout.addWidget(frequencyLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.frequency, 1, 1, 1, 3)

		layout.addWidget(modulationLabel, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.modulation, 2, 1, 1, 3)

		layout.addWidget(showFieldLabel, 3, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.showField, 3, 1, 1, 3)

		layout.addWidget(u0Label, 4, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.u0, 4, 1)

		layout.addWidget(u1Label, 5, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.u1, 5, 1)

		layout.addWidget(v0Label, 4, 2, QtCore.Qt.AlignRight)
		layout.addWidget(self.v0, 4, 3)

		layout.addWidget(v1Label, 5, 2, QtCore.Qt.AlignRight)
		layout.addWidget(self.v1, 5, 3)

		layout.addWidget(nmodeLabel, 6, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.nmodes, 6, 1)

		layout.addWidget(modeSelectorLabel, 6, 2, QtCore.Qt.AlignRight)
		layout.addWidget(self.modeSelector, 6, 3)

		layout.addWidget(self.neffLabel, 7, 0, 1, 3, QtCore.Qt.AlignRight)

		layout.addWidget(self.calculateModes, 8, 1, 1, 2)

		layout.addItem(spacer, 9, 1)
		self.tab.setLayout(layout)

	def addWaveguidePort(self):
		name = 'Waveguide Port'+str(len(self.portPlane.sources)+1)
		self.waveguidePort = WaveguidePort(name, [[0, 0], [0, 0]], self.portPlane)
		self.portPlane.addSource(self.waveguidePort)				# Add it to the geometry list
		projectNode = self.project.item
		projectNode.child(4).child(self.excitationIndex).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(4).child(self.excitationIndex).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(4).child(self.excitationIndex).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(4).child(self.excitationIndex).child(projectNode.child(4).child(self.excitationIndex).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class MyTableModel(QtCore.QAbstractTableModel):
	dataChanged = QtCore.Signal()

	def __init__(self, parent, mylist, header, grid, *args):
		QtCore.QAbstractTableModel.__init__(self, parent, *args)
		self.mylist = mylist
		self.header = header
		self.grid = grid

	def rowCount(self, parent):
		return len(self.mylist)

	def columnCount(self, parent):
		if len(self.mylist) > 0:
			return len(self.mylist[0])
		else:
			return 0

	def data(self, index, role):
		if not index.isValid():
			return None
		elif role != QtCore.Qt.DisplayRole:
			return None
		return self.mylist[index.row()][index.column()]

	def headerData(self, col, orientation, role):
		if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole and col < len(self.header):
			return self.header[col]
		return None

	def flags(self, index):
		return QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsEnabled

	def data(self, index, role):
		if index.isValid():
			if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
				return str(self.mylist[index.row()][index.column()]/self.grid.sizeScalefactor())
		return None

	def setData(self, index, value, role):
		if index.isValid():
			if role == QtCore.Qt.EditRole:
				self.mylist[index.row()][index.column()] = (float)(value)*self.grid.sizeScalefactor()
				self.dataChanged.emit()
				return True
		return False

	def updateView(self):
		self.layoutChanged.emit()

class Project():
	def __init__(self, name, plotter):
		self.name = name
		self.item = QtGui.QStandardItem(name)
		self.item.appendRow(QtGui.QStandardItem('Grid'))
		self.item.appendRow(QtGui.QStandardItem('Sensor'))
		self.item.appendRow(QtGui.QStandardItem('Geometry'))
		self.item.appendRow(QtGui.QStandardItem('Materials'))
		self.item.appendRow(QtGui.QStandardItem('Excitation'))
		self.item.appendRow(QtGui.QStandardItem('Solver'))
		self.item.appendRow(QtGui.QStandardItem('Results'))

		self.mayavi = plotter
		self.grid = Grid()
		self.solver = Solver(self)
		self.result = Result(self, self.mayavi)
		self.sensors = []
		self.geometries = []
		self.materials = []
		self.excitations = []
		self.plots = []

		self.refreshPlot = False		# Rerender the plot, rather than updating it

	def addSensor(self, sensor):
		self.sensors.append(sensor)

	def addGeometry(self, geometry):
		self.geometries.append(geometry)

	def addMaterial(self, material):
		self.materials.append(material)

	def addExcitation(self, excitation):
		self.excitations.append(excitation)

	def addPlot(self, plot):
		self.plots.append(plot)

class Solver():
	def __init__(self, project):
		self.sampleDistance = 2
		self.grid = project.grid
		self.project = project
		self.threadCount = cpu_count()
		self.progressBar = QtGui.QProgressBar()

		self.Ex = None
		self.Ey = None
		self.Ez = None
		self.Hx = None
		self.Hy = None
		self.Hz = None

	def allocateMemory(self):
		dtype = 'f4'
		self.WBEx = zeros((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)	# Work buffers
		self.WBEy = zeros((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.WBEz = zeros((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.WBHx = zeros((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.WBHy = zeros((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.WBHz = zeros((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)

		self.Ex = zeros((self.grid.timesteps/self.sampleDistance, self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)	# Output buffers
		self.Ey = zeros((self.grid.timesteps/self.sampleDistance, self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.Ez = zeros((self.grid.timesteps/self.sampleDistance, self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.Hx = zeros((self.grid.timesteps/self.sampleDistance, self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.Hy = zeros((self.grid.timesteps/self.sampleDistance, self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.Hz = zeros((self.grid.timesteps/self.sampleDistance, self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)

		self.CEx = ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CEy = ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CEz = ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)

		self.CExdy =  self.grid.dt/(self.grid.eps0*amax(self.grid.dy))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CExdz = -self.grid.dt/(self.grid.eps0*amax(self.grid.dz))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CEydz =  self.grid.dt/(self.grid.eps0*amax(self.grid.dz))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CEydx = -self.grid.dt/(self.grid.eps0*amax(self.grid.dx))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CEzdx =  self.grid.dt/(self.grid.eps0*amax(self.grid.dx))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CEzdy = -self.grid.dt/(self.grid.eps0*amax(self.grid.dy))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)

		self.CHxdz =  self.grid.dt/(self.grid.mu0*amax(self.grid.dz))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CHxdy = -self.grid.dt/(self.grid.mu0*amax(self.grid.dy))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CHydx =  self.grid.dt/(self.grid.mu0*amax(self.grid.dx))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CHydz = -self.grid.dt/(self.grid.mu0*amax(self.grid.dz))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CHzdy =  self.grid.dt/(self.grid.mu0*amax(self.grid.dy))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)
		self.CHzdx = -self.grid.dt/(self.grid.mu0*amax(self.grid.dx))*ones((self.grid.cellsX, self.grid.cellsY, self.grid.cellsZ), dtype=dtype)

	def defineMaterials(self):
		workplanes = (x for x in self.project.geometries if x.type == 'Workplane')
		for workplane in workplanes:
			extrusions = (x for x in workplane.geometries if x.type == 'Extrusion')
			for extrusion in extrusions:
				geometries = (x for x in workplane.geometries if (x.type != 'Extrusion') and (x.material != None))
				for geometry in geometries:
					if geometry.type == 'Circle':
						resolution = 100
						center = geometry.center
						radius = geometry.radius
						dalpha = 2*pi/resolution
						points = [[radius*cos(k*dalpha)+center[0], radius*sin(k*dalpha)+center[1]] for k in range(resolution)]

					elif geometry.type == 'Polygon 2D':
						points = geometry.P

					length = len(points)
					bottom = workplane.uvwToxyz([[points[k][0], points[k][1], 0] for k in range(length)])
					top = workplane.uvwToxyz([[points[k][0], points[k][1], extrusion.height] for k in range(length)])

					for k in range(len(self.grid.distanceZ)):
						cutplane = Cutplane(orientation='xy')
						intersection = cutplane.computeIntersection(self.grid.distanceZ[k], bottom, top)

						for j in range(len(self.grid.distanceY)):
							crossings = cutplane.rayCasting(intersection, self.grid.distanceY[j])

							for m in range(len(crossings)-1):
								i1 = (int)(crossings[m]/amax(self.grid.dx))
								i2 = (int)(crossings[m+1]/amax(self.grid.dx))
								for i in range(max(0, i1), min(len(self.grid.distanceX), i2)):
									self.CEx[i][j][k] = (2*self.grid.eps0*geometry.material.epsr - geometry.material.sigma*self.grid.dt)/(2*self.grid.eps0*geometry.material.epsr + geometry.material.sigma*self.grid.dt)
									self.CExdy[i][j][k] =  2*self.grid.dt/(2*self.grid.eps0*geometry.material.epsr + geometry.material.sigma*self.grid.dt)/self.grid.dy[j]
									self.CExdz[i][j][k] = -2*self.grid.dt/(2*self.grid.eps0*geometry.material.epsr + geometry.material.sigma*self.grid.dt)/self.grid.dz[k]

									self.CEy[i][j][k] = (2*self.grid.eps0*geometry.material.epsr - geometry.material.sigma*self.grid.dt)/(2*self.grid.eps0*geometry.material.epsr + geometry.material.sigma*self.grid.dt)
									self.CEydz[i][j][k] =  2*self.grid.dt/(2*self.grid.eps0*geometry.material.epsr + geometry.material.sigma*self.grid.dt)/self.grid.dz[k]
									self.CEydx[i][j][k] = -2*self.grid.dt/(2*self.grid.eps0*geometry.material.epsr + geometry.material.sigma*self.grid.dt)/self.grid.dx[i]

									self.CEz[i][j][k] = (2*self.grid.eps0*geometry.material.epsr - geometry.material.sigma*self.grid.dt)/(2*self.grid.eps0*geometry.material.epsr + geometry.material.sigma*self.grid.dt)
									self.CEzdx[i][j][k] =  2*self.grid.dt/(2*self.grid.eps0*geometry.material.epsr + geometry.material.sigma*self.grid.dt)/self.grid.dx[i]
									self.CEzdy[i][j][k] = -2*self.grid.dt/(2*self.grid.eps0*geometry.material.epsr + geometry.material.sigma*self.grid.dt)/self.grid.dy[j]

									self.CHxdz[i][j][k] =  self.grid.dt/(self.grid.mu0*geometry.material.mur*self.grid.dz[k])
									self.CHxdy[i][j][k] = -self.grid.dt/(self.grid.mu0*geometry.material.mur*self.grid.dy[j])

									self.CHydx[i][j][k] =  self.grid.dt/(self.grid.mu0*geometry.material.mur*self.grid.dx[i])
									self.CHydz[i][j][k] = -self.grid.dt/(self.grid.mu0*geometry.material.mur*self.grid.dz[k])

									self.CHzdy[i][j][k] =  self.grid.dt/(self.grid.mu0*geometry.material.mur*self.grid.dy[j])
									self.CHzdx[i][j][k] = -self.grid.dt/(self.grid.mu0*geometry.material.mur*self.grid.dx[i])

	def run(self):
		self.progressBar.setRange(0, self.grid.timesteps/self.sampleDistance)

		start = time()
		self.allocateMemory()
		self.defineMaterials()
		print "Preprocessing took %f seconds"%(time() - start)

		endX = self.grid.cellsX
		endY = self.grid.cellsY
		endZ = self.grid.cellsZ

		start = time()
		f = Field()
		f.setOutputBuffer_cpp(self.Ex, self.Ey, self.Ez, self.Hx, self.Hy, self.Hz)
		portPlanes = (x for x in self.project.excitations if x.type == 'Port Plane')
		for plane in portPlanes:
			ports = (x for x in plane.sources if x.type == 'Waveguide Port')
			for port in ports:
				modE = sin(2*pi*linspace(0, self.grid.timesteps-1, self.grid.timesteps)*self.grid.dt*port.frequency, dtype='f4')
				modH = sin(2*pi*linspace(0, self.grid.timesteps-1, self.grid.timesteps)*self.grid.dt*port.frequency - sqrt(real(port.w[port.mode]))*amax(self.grid.dz)/2.0, dtype='f4')
				f.addWaveguideSourceXY_cpp(port.imin, port.imax, port.jmin, port.jmax, (int)(port.portPlane.pos/amax(self.grid.dz)), \
					modE, modH, port.Ex[port.mode], port.Ey[port.mode], port.Hx[port.mode], port.Hy[port.mode])

		f.update_cpp(self.sampleDistance, self.grid.timesteps, self.grid.dt,  \
			self.WBEx, self.WBEy, self.WBEz, self.WBHx, self.WBHy, self.WBHz, self.CExdy, self.CExdz, \
			self.CEydx, self.CEydz, self.CEzdx, self.CEzdy, self.CHxdy, self.CHxdz, self.CHydx, self.CHydz, self.CHzdx, self.CHzdy, \
			self.CEx, self.CEy, self.CEz)

		elapsed_time = time() - start
		print "%f seconds taken"%(elapsed_time)

class SolverWidget():
	def __init__(self, project):
		self.project = project
		self.createWidget()

	def createWidget(self):
		self.tab = QtGui.QWidget()

		self.startButton = QtGui.QPushButton()
		self.threadCount = QtGui.QDoubleSpinBox()
		self.sampleDistance = QtGui.QDoubleSpinBox()

		self.threadCount.setRange(0, 100)
		self.sampleDistance.setRange(0, 1E5)

		self.threadCount.setDecimals(0)
		self.sampleDistance.setDecimals(0)

		self.threadCount.setValue(self.project.solver.threadCount)
		self.sampleDistance.setValue(self.project.solver.sampleDistance)

		self.startButton.setText("Run")
		threadCountLabel = QtGui.QLabel(u"n\u00B0 threads:")
		sampleDistanceLabel = QtGui.QLabel("Sample Distance:")

		self.startButton.clicked.connect(self.startClicked)
		self.threadCount.valueChanged.connect(self.threadCountChanged)
		self.sampleDistance.valueChanged.connect(self.sampleDistanceChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.addWidget(sampleDistanceLabel, 0, 0)
		layout.addWidget(self.sampleDistance, 0, 1)

		layout.addWidget(threadCountLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.threadCount, 1, 1)

		layout.addWidget(self.startButton, 2, 0)
		layout.addWidget(self.project.solver.progressBar, 2, 1)

		layout.addItem(spacer, 3, 0)
		self.tab.setLayout(layout)

	def startClicked(self):
		self.project.solver.run()

	def threadCountChanged(self, value):
		self.project.solver.threadCount = (int)(value)

	def sampleDistanceChanged(self, value):
		self.project.solver.sampleDistance = (int)(value)

class Result():
	def __init__(self, project, mayavi):
		self.grid = project.grid
		self.showFields = False
		self.project = project
		self.mayavi = mayavi
		self.timestep = 0
		self.createTab()

	def timeSliderChanged(self, value):
		self.timestep = value
		self.timeLabel.setText("Time: "+str(value*self.grid.dt/self.grid.timeScaleFactor()*self.project.solver.sampleDistance)[0:4]+" "+self.grid.freqTimeConversion[self.grid.freqUnit])
		if self.showFields:
			self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def showFieldsChanged(self, state):
		if state:
			self.showFields = True
			self.mayavi.visualization.update_plot(self.project)		# Update the 3D view
		else:
			self.showFields = False

	def createTab(self):
		self.tab = QtGui.QWidget()

		self.time = QtGui.QSlider(QtCore.Qt.Horizontal)
		# print (self.grid.timesteps-1)/self.project.solver.sampleDistance
		self.time.setRange(0, (self.grid.timesteps-1)/self.project.solver.sampleDistance)
		self.time.valueChanged.connect(self.timeSliderChanged)
		self.time.setMinimumWidth(100)


		showFieldsLabel = QtGui.QLabel("Show Fields: ")
		showFieldsCheckBox = QtGui.QCheckBox()
		if self.showFields == True:
			showFieldsCheckBox.setChecked(True)
		showFieldsCheckBox.stateChanged.connect(self.showFieldsChanged)

		self.time.setValue(self.timestep)
		# self.timeLabel = QtGui.QLabel("Time: 0 "+self.grid.freqTimeConversion[self.grid.freqUnit])
		self.timeLabel = QtGui.QLabel("Time: "+str(self.timestep*self.grid.dt/self.grid.timeScaleFactor()*self.project.solver.sampleDistance)[0:4]+" "+self.grid.freqTimeConversion[self.grid.freqUnit])

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.addWidget(self.timeLabel, 0, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.time, 0, 1)

		layout.addWidget(showFieldsLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(showFieldsCheckBox, 1, 1)

		layout.addItem(spacer, 3, 0)
		self.tab.setLayout(layout)

class ScalarCutplane():
	def __init__(self, name):
		self.name = name
		self.type = 'Scalar Cutplane'
		self.orientation = 'x_axes'
		self.minimum = -1e-5
		self.maximum = 1e-5
		self.data = 'Hz'
		self.position = 0		# Position along orientation
		self.autoScale = False
		self.plotted = False	# The cutplane has not yet been plotted
		self.plot = None

class ScalarCutplaneWidget():
	def __init__(self, project, mayavi, view, model, scalarCutplane=None):
		self.project = project
		self.mayavi = mayavi
		self.grid = self.project.grid
		self.view = view
		self.model = model
		self.scalarCutplane = scalarCutplane
		if scalarCutplane == None:
			self.addScalarCutplane()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def orientationChanged(self, index):
		self.project.refreshPlot = True
		if index == 0:
			self.scalarCutplane.orientation = 'x_axes'
		elif index == 1:
			self.scalarCutplane.orientation = 'y_axes'
		elif index == 2:
			self.scalarCutplane.orientation = 'z_axes'

	def dataChanged(self):
		self.scalarCutplane.data = self.data.toPlainText()

	def minimumChanged(self, value):
		self.project.refreshPlot = True
		self.scalarCutplane.minimum = value

	def maximumChanged(self, value):
		self.project.refreshPlot = True
		self.scalarCutplane.maximum = value

	def positionChanged(self, value):
		self.project.refreshPlot = True
		self.scalarCutplane.position = value*self.grid.sizeScalefactor()

	def nameChanged(self, name):
		self.scalarCutplane.name = name
		self.currentNode.setText(name)

	def autoScaleChanged(self, state):
		self.project.refreshPlot = True
		if state:
			self.scalarCutplane = True
			self.minimum.setEnabled(False)
			self.maximum.setEnabled(False)
		else:
			self.scalarCutplane = False
			self.minimum.setEnabled(True)
			self.maximum.setEnabled(True)

	def updateClicked(self):
		self.project.result.showFields = True
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def createWidget(self):
		self.tab = QtGui.QWidget()

		orientation = QtGui.QComboBox()
		self.data = QtGui.QPlainTextEdit()
		self.minimum = QtGui.QDoubleSpinBox()
		self.maximum = QtGui.QDoubleSpinBox()
		self.position = QtGui.QDoubleSpinBox()
		autoScale = QtGui.QCheckBox()
		update = QtGui.QPushButton("Update")
		name = QtGui.QLineEdit(self.project.name)

		update.setMaximumWidth(80)

		self.maximum.setRange(-1e4, 1e4)
		self.minimum.setRange(-1e4, 1e4)
		self.position.setRange(-1e4, 1e4)
		self.maximum.setDecimals(6)
		self.minimum.setDecimals(6)
		self.position.setDecimals(3)
		name.setText(self.scalarCutplane.name)

		self.minimum.setValue(self.scalarCutplane.minimum)
		self.maximum.setValue(self.scalarCutplane.maximum)
		self.position.setValue(self.scalarCutplane.position/self.grid.sizeScalefactor())
		self.data.setPlainText(self.scalarCutplane.data)
		autoScale.setChecked(self.scalarCutplane.autoScale)
		
		font = self.data.document().defaultFont()
		fontMetrics = QtGui.QFontMetrics(font)
		textSize = fontMetrics.size(0, "Hz")
		textHeight = textSize.height()+11		# Constant might need some tweaking
		self.data.setMaximumHeight(textHeight)
		self.data.setMinimumWidth(100)

		orientationLabel = QtGui.QLabel("Orientation:")
		dataLabel = QtGui.QLabel("Data:")
		minimumLabel = QtGui.QLabel("Min:")
		maximumLabel = QtGui.QLabel("Max:")
		positionLabel = QtGui.QLabel("Position ["+self.grid.sizeUnit+"]:")
		autoScaleLabel = QtGui.QLabel("Auto scale:")
		orientation.addItems(('x-axis', 'y-axis', 'z-axis'))
		nameLabel = QtGui.QLabel("Name:")

		orientation.setCurrentIndex(ord(self.scalarCutplane.orientation[0])-ord('x'))

		orientation.currentIndexChanged.connect(self.orientationChanged)
		self.data.textChanged.connect(self.dataChanged)
		self.minimum.valueChanged.connect(self.minimumChanged)
		self.maximum.valueChanged.connect(self.maximumChanged)
		autoScale.stateChanged.connect(self.autoScaleChanged)
		self.position.valueChanged.connect(self.positionChanged)
		update.clicked.connect(self.updateClicked)
		name.textChanged.connect(self.nameChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.addWidget(nameLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(name, 1, 1)

		layout.addWidget(dataLabel, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.data, 2, 1)

		layout.addWidget(orientationLabel, 3, 0, QtCore.Qt.AlignRight)
		layout.addWidget(orientation, 3, 1)

		layout.addWidget(positionLabel, 4, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.position, 4, 1)

		layout.addWidget(minimumLabel, 5, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.minimum, 5, 1)

		layout.addWidget(maximumLabel, 6, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.maximum, 6, 1)

		layout.addWidget(autoScaleLabel, 7, 0, QtCore.Qt.AlignRight)
		layout.addWidget(autoScale, 7, 1)

		layout.addWidget(update, 8, 1)

		layout.addItem(spacer, 9, 1)
		self.tab.setLayout(layout)

	def addScalarCutplane(self):
		name = 'Plot'+str(len(self.project.plots)+1)
		self.scalarCutplane = ScalarCutplane(name)
		self.project.addPlot(self.scalarCutplane)				# Add it to the material list
		projectNode = self.project.item
		projectNode.child(6).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(6).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(6).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(6).child(projectNode.child(6).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class ScalarVolume():
	def __init__(self, name):
		self.name = name
		self.type = 'Scalar Volume'
		self.minimum = 0
		self.maximum = 1e-5
		self.data = 'sqrt(Hx**2+Hy**2+Hz**2)'
		self.autoScale = False
		self.plotted = False
		self.plot = None

class ScalarVolumeWidget():
	def __init__(self, project, mayavi, view, model, scalarVolume=None):
		self.project = project
		self.mayavi = mayavi
		self.grid = self.project.grid
		self.view = view
		self.model = model
		self.scalarVolume = scalarVolume
		if scalarVolume == None:
			self.addScalarVolume()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def dataChanged(self):
		self.scalarVolume.data = self.data.toPlainText()

	def minimumChanged(self, value):
		self.project.refreshPlot = True
		self.scalarVolume.minimum = value

	def maximumChanged(self, value):
		self.project.refreshPlot = True
		self.scalarVolume.maximum = value

	def nameChanged(self, name):
		self.scalarVolume.name = name
		self.currentNode.setText(name)

	def autoScaleChanged(self, state):
		self.project.refreshPlot = True
		if state:
			self.scalarVolume.autoScale = True
			self.minimum.setEnabled(False)
			self.maximum.setEnabled(False)
		else:
			self.scalarVolume.autoScale = False
			self.minimum.setEnabled(True)
			self.maximum.setEnabled(True)

	def updateClicked(self):
		self.project.result.showFields = True
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def createWidget(self):
		self.tab = QtGui.QWidget()

		self.data = QtGui.QPlainTextEdit()
		self.minimum = QtGui.QDoubleSpinBox()
		self.maximum = QtGui.QDoubleSpinBox()
		autoScale = QtGui.QCheckBox()
		update = QtGui.QPushButton("Update")
		name = QtGui.QLineEdit(self.project.name)

		update.setMaximumWidth(80)

		self.maximum.setRange(0, 1e4)
		self.minimum.setRange(0, 1e4)
		self.maximum.setDecimals(6)
		self.minimum.setDecimals(6)
		name.setText(self.scalarVolume.name)

		self.minimum.setValue(self.scalarVolume.minimum)
		self.maximum.setValue(self.scalarVolume.maximum)
		self.data.setPlainText(self.scalarVolume.data)
		autoScale.setChecked(self.scalarVolume.autoScale)
		
		font = self.data.document().defaultFont()
		fontMetrics = QtGui.QFontMetrics(font)
		textSize = fontMetrics.size(0, "Hz")
		textHeight = textSize.height()+11		# Constant might need some tweaking
		self.data.setMaximumHeight(textHeight)
		self.data.setMinimumWidth(100)

		dataLabel = QtGui.QLabel("Data:")
		minimumLabel = QtGui.QLabel("Min:")
		maximumLabel = QtGui.QLabel("Max:")
		autoScaleLabel = QtGui.QLabel("Auto scale:")
		nameLabel = QtGui.QLabel("Name:")

		self.data.textChanged.connect(self.dataChanged)
		self.minimum.valueChanged.connect(self.minimumChanged)
		self.maximum.valueChanged.connect(self.maximumChanged)
		autoScale.stateChanged.connect(self.autoScaleChanged)
		update.clicked.connect(self.updateClicked)
		name.textChanged.connect(self.nameChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.addWidget(nameLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(name, 1, 1)

		layout.addWidget(dataLabel, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.data, 2, 1)

		layout.addWidget(minimumLabel, 3, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.minimum, 3, 1)

		layout.addWidget(maximumLabel, 4, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.maximum, 4, 1)

		layout.addWidget(autoScaleLabel, 5, 0, QtCore.Qt.AlignRight)
		layout.addWidget(autoScale, 5, 1)

		layout.addWidget(update, 6, 1)

		layout.addItem(spacer, 7, 1)
		self.tab.setLayout(layout)

	def addScalarVolume(self):
		name = 'Plot'+str(len(self.project.plots)+1)
		self.scalarVolume = ScalarVolume(name)
		self.project.addPlot(self.scalarVolume)				# Add it to the material list
		projectNode = self.project.item
		projectNode.child(6).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(6).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(6).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(6).child(projectNode.child(6).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class VectorVolume():
	def __init__(self, name):
		self.name = name
		self.type = 'Vector Volume'
		self.minimum = 0
		self.maximum = 1e-5
		self.xComponent = 'Hx'
		self.yComponent = 'Hy'
		self.zComponent = 'Hz'
		self.autoScale = False
		self.plotted = False	# The cutplane has not yet been plotted
		self.maskPoints = 10
		self.plot = None		# Contains the plot

class VectorVolumeWidget():
	def __init__(self, project, mayavi, view, model, vectorVolume=None):
		self.project = project
		self.mayavi = mayavi
		self.grid = self.project.grid
		self.view = view
		self.model = model
		self.vectorVolume = vectorVolume
		if vectorVolume == None:
			self.addVectorVolume()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def dataChanged(self):
		self.vectorVolume.data = self.data.toPlainText()

	def minimumChanged(self, value):
		self.project.refreshPlot = True
		self.vectorVolume.minimum = value

	def maximumChanged(self, value):
		self.project.refreshPlot = True
		self.vectorVolume.maximum = value

	def nameChanged(self, name):
		self.vectorVolume.name = name
		self.currentNode.setText(name)

	def autoScaleChanged(self, state):
		self.project.refreshPlot = True
		if state:
			self.vectorVolume.autoScale = True
			self.minimum.setEnabled(False)
			self.maximum.setEnabled(False)
		else:
			self.vectorVolume.autoScale = False
			self.minimum.setEnabled(True)
			self.maximum.setEnabled(True)

	def xComponentChanged(self):
		self.vectorVolume.xComponent = self.xComponent.toPlainText()

	def yComponentChanged(self):
		self.vectorVolume.yComponent = self.yComponent.toPlainText()

	def zComponentChanged(self):
		self.vectorVolume.zComponent = self.zComponent.toPlainText()

	def updateClicked(self):
		self.project.result.showFields = True
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def maskPointsChanged(self, value):
		self.project.refreshPlot = True
		self.vectorVolume.maskPoints = (int)(value)

	def createWidget(self):
		self.tab = QtGui.QWidget()

		self.xComponent = QtGui.QPlainTextEdit()
		self.yComponent = QtGui.QPlainTextEdit()
		self.zComponent = QtGui.QPlainTextEdit()
		self.minimum = QtGui.QDoubleSpinBox()
		self.maximum = QtGui.QDoubleSpinBox()
		self.maskPoints = QtGui.QDoubleSpinBox()
		autoScale = QtGui.QCheckBox()
		update = QtGui.QPushButton("Update")
		name = QtGui.QLineEdit(self.project.name)

		update.setMaximumWidth(80)

		self.maximum.setRange(0, 1e4)
		self.minimum.setRange(0, 1e4)
		self.maskPoints.setRange(0, 1E4)
		self.maximum.setDecimals(6)
		self.minimum.setDecimals(6)
		self.maskPoints.setDecimals(0)
		name.setText(self.vectorVolume.name)

		self.minimum.setValue(self.vectorVolume.minimum)
		self.maximum.setValue(self.vectorVolume.maximum)
		self.maskPoints.setValue(self.vectorVolume.maskPoints)
		self.xComponent.setPlainText(self.vectorVolume.xComponent)
		self.yComponent.setPlainText(self.vectorVolume.yComponent)
		self.zComponent.setPlainText(self.vectorVolume.zComponent)
		autoScale.setChecked(self.vectorVolume.autoScale)
		
		font = self.xComponent.document().defaultFont()
		fontMetrics = QtGui.QFontMetrics(font)
		textSize = fontMetrics.size(0, "Hz")
		textHeight = textSize.height()+11		# Constant might need some tweaking
		self.xComponent.setMaximumHeight(textHeight)
		self.xComponent.setMinimumWidth(100)
		self.yComponent.setMaximumHeight(textHeight)
		self.yComponent.setMinimumWidth(100)
		self.zComponent.setMaximumHeight(textHeight)
		self.zComponent.setMinimumWidth(100)

		xLabel = QtGui.QLabel("x comp.:")
		yLabel = QtGui.QLabel("y comp.:")
		zLabel = QtGui.QLabel("z comp.:")
		minimumLabel = QtGui.QLabel("Min:")
		maximumLabel = QtGui.QLabel("Max:")
		autoScaleLabel = QtGui.QLabel("Auto scale:")
		maskPointsLabel = QtGui.QLabel("Mask points:")
		nameLabel = QtGui.QLabel("Name:")

		self.xComponent.textChanged.connect(self.xComponentChanged)
		self.yComponent.textChanged.connect(self.yComponentChanged)
		self.zComponent.textChanged.connect(self.zComponentChanged)
		self.minimum.valueChanged.connect(self.minimumChanged)
		self.maximum.valueChanged.connect(self.maximumChanged)
		autoScale.stateChanged.connect(self.autoScaleChanged)
		update.clicked.connect(self.updateClicked)
		name.textChanged.connect(self.nameChanged)
		self.maskPoints.valueChanged.connect(self.maskPointsChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.addWidget(nameLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(name, 1, 1)

		layout.addWidget(xLabel, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.xComponent, 2, 1)
		layout.addWidget(yLabel, 3, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.yComponent, 3, 1)
		layout.addWidget(zLabel, 4, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.zComponent, 4, 1)

		layout.addWidget(minimumLabel, 5, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.minimum, 5, 1)

		layout.addWidget(maximumLabel, 6, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.maximum, 6, 1)

		layout.addWidget(maskPointsLabel, 7, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.maskPoints, 7, 1)

		layout.addWidget(autoScaleLabel, 8, 0, QtCore.Qt.AlignRight)
		layout.addWidget(autoScale, 8, 1)

		layout.addWidget(update, 9, 1)

		layout.addItem(spacer, 10, 1)
		self.tab.setLayout(layout)

	def addVectorVolume(self):
		name = 'Plot'+str(len(self.project.plots)+1)
		self.vectorVolume = VectorVolume(name)
		self.project.addPlot(self.vectorVolume)				# Add it to the material list
		projectNode = self.project.item
		projectNode.child(6).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(6).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(6).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(6).child(projectNode.child(6).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class VectorCutplane():
	def __init__(self, name):
		self.name = name
		self.type = 'Vector Cutplane'
		self.minimum = 0
		self.maximum = 1e-5
		self.xComponent = 'Hx'
		self.yComponent = 'Hy'
		self.zComponent = 'Hz'
		self.autoScale = False
		self.plotted = False	# The cutplane has not yet been plotted
		self.maskPoints = 1
		self.orientation = 'x_axes'
		self.plot = None		# Contains the plot	

class VectorCutplaneWidget():
	def __init__(self, project, mayavi, view, model, vectorCutplane=None):
		self.project = project
		self.mayavi = mayavi
		self.grid = self.project.grid
		self.view = view
		self.model = model
		self.vectorCutplane = vectorCutplane
		if vectorCutplane == None:
			self.addVectorCutplane()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def dataChanged(self):
		self.vectorCutplane.data = self.data.toPlainText()

	def minimumChanged(self, value):
		self.project.refreshPlot = True
		self.vectorCutplane.minimum = value

	def maximumChanged(self, value):
		self.project.refreshPlot = True
		self.vectorCutplane.maximum = value

	def nameChanged(self, name):
		self.vectorCutplane.name = name
		self.currentNode.setText(name)

	def autoScaleChanged(self, state):
		self.project.refreshPlot = True
		if state:
			self.vectorCutplane.autoScale = True
			self.minimum.setEnabled(False)
			self.maximum.setEnabled(False)
		else:
			self.vectorCutplane.autoScale = False
			self.minimum.setEnabled(True)
			self.maximum.setEnabled(True)

	def xComponentChanged(self):
		self.vectorCutplane.xComponent = self.xComponent.toPlainText()

	def yComponentChanged(self):
		self.vectorCutplane.yComponent = self.yComponent.toPlainText()

	def zComponentChanged(self):
		self.vectorCutplane.zComponent = self.zComponent.toPlainText()

	def orientationChanged(self, index):
		self.project.refreshPlot = True
		if index == 0:
			self.vectorCutplane.orientation = 'x_axes'
		elif index == 1:
			self.vectorCutplane.orientation = 'y_axes'
		elif index == 2:
			self.vectorCutplane.orientation = 'z_axes'

	def updateClicked(self):
		self.project.result.showFields = True
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def maskPointsChanged(self, value):
		self.project.refreshPlot = True
		self.vectorCutplane.maskPoints = (int)(value)

	def createWidget(self):
		self.tab = QtGui.QWidget()

		self.xComponent = QtGui.QPlainTextEdit()
		self.yComponent = QtGui.QPlainTextEdit()
		self.zComponent = QtGui.QPlainTextEdit()
		self.minimum = QtGui.QDoubleSpinBox()
		self.maximum = QtGui.QDoubleSpinBox()
		self.maskPoints = QtGui.QDoubleSpinBox()
		orientation = QtGui.QComboBox()
		autoScale = QtGui.QCheckBox()
		update = QtGui.QPushButton("Update")
		name = QtGui.QLineEdit(self.project.name)

		update.setMaximumWidth(80)

		self.maximum.setRange(0, 1e4)
		self.minimum.setRange(0, 1e4)
		self.maskPoints.setRange(0, 1E4)
		self.maximum.setDecimals(6)
		self.minimum.setDecimals(6)
		self.maskPoints.setDecimals(0)
		name.setText(self.vectorCutplane.name)

		self.minimum.setValue(self.vectorCutplane.minimum)
		self.maximum.setValue(self.vectorCutplane.maximum)
		self.maskPoints.setValue(self.vectorCutplane.maskPoints)
		orientation.setCurrentIndex(ord(self.vectorCutplane.orientation[0])-ord('x'))
		self.xComponent.setPlainText(self.vectorCutplane.xComponent)
		self.yComponent.setPlainText(self.vectorCutplane.yComponent)
		self.zComponent.setPlainText(self.vectorCutplane.zComponent)
		autoScale.setChecked(self.vectorCutplane.autoScale)
		
		font = self.xComponent.document().defaultFont()
		fontMetrics = QtGui.QFontMetrics(font)
		textSize = fontMetrics.size(0, "Hz")
		textHeight = textSize.height()+11		# Constant might need some tweaking
		self.xComponent.setMaximumHeight(textHeight)
		self.xComponent.setMinimumWidth(100)
		self.yComponent.setMaximumHeight(textHeight)
		self.yComponent.setMinimumWidth(100)
		self.zComponent.setMaximumHeight(textHeight)
		self.zComponent.setMinimumWidth(100)

		xLabel = QtGui.QLabel("x comp.:")
		yLabel = QtGui.QLabel("y comp.:")
		zLabel = QtGui.QLabel("z comp.:")
		minimumLabel = QtGui.QLabel("Min:")
		maximumLabel = QtGui.QLabel("Max:")
		autoScaleLabel = QtGui.QLabel("Auto scale:")
		maskPointsLabel = QtGui.QLabel("Mask points:")
		nameLabel = QtGui.QLabel("Name:")
		orientationLabel = QtGui.QLabel("Orientation:")
		orientation.addItems(('x-axis', 'y-axis', 'z-axis'))

		self.xComponent.textChanged.connect(self.xComponentChanged)
		self.yComponent.textChanged.connect(self.yComponentChanged)
		self.zComponent.textChanged.connect(self.zComponentChanged)
		self.minimum.valueChanged.connect(self.minimumChanged)
		self.maximum.valueChanged.connect(self.maximumChanged)
		autoScale.stateChanged.connect(self.autoScaleChanged)
		update.clicked.connect(self.updateClicked)
		name.textChanged.connect(self.nameChanged)
		self.maskPoints.valueChanged.connect(self.maskPointsChanged)
		orientation.currentIndexChanged.connect(self.orientationChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.addWidget(nameLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(name, 1, 1)

		layout.addWidget(xLabel, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.xComponent, 2, 1)
		layout.addWidget(yLabel, 3, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.yComponent, 3, 1)
		layout.addWidget(zLabel, 4, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.zComponent, 4, 1)

		layout.addWidget(minimumLabel, 5, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.minimum, 5, 1)

		layout.addWidget(maximumLabel, 6, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.maximum, 6, 1)

		layout.addWidget(orientationLabel, 7, 0, QtCore.Qt.AlignRight)
		layout.addWidget(orientation, 7, 1)

		layout.addWidget(maskPointsLabel, 8, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.maskPoints, 8, 1)

		layout.addWidget(autoScaleLabel, 9, 0, QtCore.Qt.AlignRight)
		layout.addWidget(autoScale, 9, 1)

		layout.addWidget(update, 10, 1)

		layout.addItem(spacer, 11, 1)
		self.tab.setLayout(layout)

	def addVectorCutplane(self):
		name = 'Plot'+str(len(self.project.plots)+1)
		self.vectorCutplane = VectorCutplane(name)
		self.project.addPlot(self.vectorCutplane)				# Add it to the material list
		projectNode = self.project.item
		projectNode.child(6).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(6).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(6).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(6).child(projectNode.child(6).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class Contour():
	def __init__(self, name):
		self.name = name
		self.type = 'Contour'
		self.minimum = -1e-5
		self.maximum = 1e-5
		self.data = 'Hz'
		self.contours = 4
		self.opacity = 0.5
		self.autoScale = False
		self.plotted = False	# The cutplane has not yet been plotted
		self.plot = None		# Contains the plot

class ContourWidget():
	def __init__(self, project, mayavi, view, model, contour=None):
		self.project = project
		self.mayavi = mayavi
		self.grid = self.project.grid
		self.view = view
		self.model = model
		self.contour = contour
		if contour == None:
			self.addContour()
		self.currentNode = self.model.itemFromIndex(self.view.currentIndex())
		self.createWidget()

	def dataChanged(self):
		self.contour.data = self.data.toPlainText()

	def minimumChanged(self, value):
		self.project.refreshPlot = True
		self.contour.minimum = value

	def maximumChanged(self, value):
		self.project.refreshPlot = True
		self.contour.maximum = value

	def positionChanged(self, value):
		self.project.refreshPlot = True
		self.contour.position = value*self.grid.sizeScalefactor()

	def nameChanged(self, name):
		self.contour.name = name
		self.currentNode.setText(name)

	def autoScaleChanged(self, state):
		self.project.refreshPlot = True
		if state:
			self.contour = True
			self.minimum.setEnabled(False)
			self.maximum.setEnabled(False)
		else:
			self.contour = False
			self.minimum.setEnabled(True)
			self.maximum.setEnabled(True)

	def opacityChanged(self, value):
		self.project.refreshPlot = True
		self.contour.opacity = value/99.0

	def contoursChanged(self, value):
		self.project.refreshPlot = True
		self.contour.contours = (int)(value)

	def updateClicked(self):
		self.project.result.showFields = True
		self.mayavi.visualization.update_plot(self.project)		# Update the 3D view

	def createWidget(self):
		self.tab = QtGui.QWidget()

		self.data = QtGui.QPlainTextEdit()
		self.minimum = QtGui.QDoubleSpinBox()
		self.maximum = QtGui.QDoubleSpinBox()
		self.contours = QtGui.QDoubleSpinBox()
		autoScale = QtGui.QCheckBox()
		update = QtGui.QPushButton("Update")
		name = QtGui.QLineEdit(self.project.name)
		opacity = QtGui.QSlider(QtCore.Qt.Horizontal)

		update.setMaximumWidth(80)

		self.maximum.setRange(-1e4, 1e4)
		self.minimum.setRange(-1e4, 1e4)
		self.contours.setRange(0, 1e4)
		self.maximum.setDecimals(6)
		self.minimum.setDecimals(6)
		self.contours.setDecimals(0)
		opacity.setRange(0, 100)

		self.minimum.setValue(self.contour.minimum)
		self.maximum.setValue(self.contour.maximum)
		self.contours.setValue(self.contour.contours)
		self.data.setPlainText(self.contour.data)
		autoScale.setChecked(self.contour.autoScale)
		name.setText(self.contour.name)
		opacity.setValue(self.contour.opacity*100)
		self.contours.setValue(self.contour.contours)
		
		font = self.data.document().defaultFont()
		fontMetrics = QtGui.QFontMetrics(font)
		textSize = fontMetrics.size(0, "Hz")
		textHeight = textSize.height()+11		# Constant might need some tweaking
		self.data.setMaximumHeight(textHeight)
		self.data.setMinimumWidth(100)

		dataLabel = QtGui.QLabel("Data:")
		minimumLabel = QtGui.QLabel("Min:")
		maximumLabel = QtGui.QLabel("Max:")
		autoScaleLabel = QtGui.QLabel("Auto scale:")
		nameLabel = QtGui.QLabel("Name:")
		opacityLabel = QtGui.QLabel("Opacity:")
		contoursLabel = QtGui.QLabel(u"n\u00B0 contours:")

		self.data.textChanged.connect(self.dataChanged)
		self.minimum.valueChanged.connect(self.minimumChanged)
		self.maximum.valueChanged.connect(self.maximumChanged)
		autoScale.stateChanged.connect(self.autoScaleChanged)
		update.clicked.connect(self.updateClicked)
		name.textChanged.connect(self.nameChanged)
		opacity.valueChanged.connect(self.opacityChanged)
		self.contours.valueChanged.connect(self.contoursChanged)

		spacer = QtGui.QSpacerItem(10, 50, hPolicy=QtGui.QSizePolicy.Minimum, vPolicy=QtGui.QSizePolicy.Expanding)

		layout = QtGui.QGridLayout()
		layout.setColumnStretch(1, 1)
		layout.addWidget(nameLabel, 0, 0, QtCore.Qt.AlignRight)
		layout.addWidget(name, 0, 1)

		layout.addWidget(dataLabel, 1, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.data, 1, 1)

		layout.addWidget(opacityLabel, 2, 0, QtCore.Qt.AlignRight)
		layout.addWidget(opacity, 2, 1)

		layout.addWidget(contoursLabel, 3, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.contours, 3, 1)

		layout.addWidget(minimumLabel, 4, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.minimum, 4, 1)

		layout.addWidget(maximumLabel, 5, 0, QtCore.Qt.AlignRight)
		layout.addWidget(self.maximum, 5, 1)

		layout.addWidget(autoScaleLabel, 6, 0, QtCore.Qt.AlignRight)
		layout.addWidget(autoScale, 6, 1)

		layout.addWidget(update, 7, 1)

		layout.addItem(spacer, 8, 1)
		self.tab.setLayout(layout)

	def addContour(self):
		name = 'Plot'+str(len(self.project.plots)+1)
		self.contour = Contour(name)
		self.project.addPlot(self.contour)								# Add it to the material list
		projectNode = self.project.item
		projectNode.child(6).appendRow(QtGui.QStandardItem(name))		# Add it to the treeview
		self.view.expand(projectNode.child(6).index())					# Expand the treeview at the new node

		self.view.selectionModel().setCurrentIndex(projectNode.child(6).index(), QtGui.QItemSelectionModel.Deselect)
		self.view.selectionModel().setCurrentIndex(projectNode.child(6).child(projectNode.child(6).rowCount()-1).index(), QtGui.QItemSelectionModel.Select)

class FDTD(QtGui.QWidget):
	def __init__(self):
		super(FDTD, self).__init__()
		self.initUI()

	def initUI(self):
		layout = QtGui.QGridLayout()
		self.setLayout(layout)

		# Create treeview and assign signal <-> slots
		self.view = QtGui.QTreeView()
		self.model = QtGui.QStandardItemModel()
		self.view.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
		self.view.setModel(self.model)
		self.view.setHeaderHidden(True)
		self.view.setUniformRowHeights(True)
		self.view.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
		self.view.selectionModel().currentChanged.connect(self.newSelection)
		self.hold = False	# Variable used for tracking changes in the selected row of view

		self.view.setFirstColumnSpanned(0, self.view.rootIndex(), True)
		self.projects = []

		self.view.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
		self.view.customContextMenuRequested.connect(self.openMenu)

		self.plotter2D = Plotter2D()

		plotter3D = QtGui.QWidget()
		self.mayavi = MayaviQWidget(plotter3D)
		self.newProject("Untitled")

		# Setup the layout
		self.hSplitter = QtGui.QSplitter(QtCore.Qt.Horizontal, self)
		self.hSplitter.addWidget(self.view)
		self.hSplitter.addWidget(self.mayavi)
		self.hSplitter.addWidget(self.plotter2D.widget)

		layout.addWidget(self.hSplitter, 0, 0)

		self.menubar = QtGui.QMenuBar()
		fileMenu = self.menubar.addMenu('&File')
		exit = QtGui.QAction( 'SomethingElse', self )
		fileMenu.addAction(exit)

		self.resize(1200, 700)
		self.center()
		self.setWindowTitle('Untitled')		# Later on replaced by the name of the project
		# plotter.show()
		self.show()

		self.plotter2D.widget.hide()			# After show, because show adds info about size apparently, which is needed when unhiding

	def center(self):
		qr = self.frameGeometry()
		cp = QtGui.QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

	def newProject(self, name="Untitled"):
		self.projects.append(Project(name, self.mayavi))
		self.model.appendRow(self.projects[-1].item)

	def deleteProject(self, index):
		del self.projects[index]
		self.model.removeRow(index)		# Probably incorrect, because the childs aren't removed, only the parent

	def openMenu(self, position):
		if len(self.projects) > 0:
			hierarchyLevel = 0
			node = self.view.currentIndex()
			selectedProject = 0
			while not(node.parent().data() is None):
				node = node.parent()
				selectedProject = node.row()
				hierarchyLevel += 1

			if hierarchyLevel == 0:
				menu = QtGui.QMenu()
				deleteAction = menu.addAction("Delete")
				action = menu.exec_(self.view.mapToGlobal(position))
				if action == deleteAction:
					self.deleteProject(selectedProject)

			elif hierarchyLevel == 1:
				if self.view.currentIndex().data() == "Grid":
					print "right click on Grid"
				elif self.view.currentIndex().data() == "Sensor":
					print "sensor", selectedProject
				elif self.view.currentIndex().data() == "Geometry":
					menu = QtGui.QMenu()
					block = menu.addAction("Block")
					workplane = menu.addAction("Workplane")
					action = menu.exec_(self.view.mapToGlobal(position))
					if action == block:
						self.createBlockWidget(selectedProject)
					elif action == workplane:
						self.createWorkplaneWidget(selectedProject)
				elif self.view.currentIndex().data() == "Materials":
					menu = QtGui.QMenu()
					newMaterial = menu.addAction("New Material")
					action = menu.exec_(self.view.mapToGlobal(position))
					if action == newMaterial:
						self.createMaterialWidget(selectedProject)
				elif self.view.currentIndex().data() == "Excitation":
					menu = QtGui.QMenu()
					addPortPlane = menu.addAction("Add Port Plane")
					action = menu.exec_(self.view.mapToGlobal(position))
					if action == addPortPlane:
						self.createPortPlaneWidget(selectedProject)
				# elif self.view.currentIndex().data() == "Solver":
				elif self.view.currentIndex().data() == "Results":
					menu = QtGui.QMenu()
					scalarVolume = menu.addAction("Scalar volume plot")
					scalarCutplane = menu.addAction("Scalar cutplane")
					vectorVolume = menu.addAction("Vectorial volume plot")
					vectorCutplane = menu.addAction("Vectorial cutplane")
					contour = menu.addAction("Contour plot")

					action = menu.exec_(self.view.mapToGlobal(position))
					if action == scalarVolume:
						self.createScalarVolumeWidget(selectedProject)
					elif action == scalarCutplane:
						self.createScalarCutplaneWidget(selectedProject)
					elif action == vectorVolume:
						self.createVectorVolumeWidget(selectedProject)
					elif action == vectorCutplane:
						self.createVectorCutplane(selectedProject)
					elif action == contour:
						self.createContourWidget(selectedProject)

			elif hierarchyLevel == 2:
				if self.view.currentIndex().parent().row() == 2:
					menu = QtGui.QMenu()
					if self.projects[selectedProject].geometries[self.view.currentIndex().row()].type == 'Workplane':
						addPolygon = menu.addAction("Add Polygon")
						addCircle = menu.addAction("Add Circle")
						addExtrusion = menu.addAction("Add Extrusion")
					deleteGeometry = menu.addAction("Delete Geometry")
					action = menu.exec_(self.view.mapToGlobal(position))
					if action == deleteGeometry:
						del self.projects[selectedProject].geometries[self.view.currentIndex().row()]
						self.projects[selectedProject].item.child(2).removeRow(self.view.currentIndex().row())
						self.mayavi.visualization.update_plot(self.projects[selectedProject])
					if action == addPolygon:
						self.createPolygon2DWidget(selectedProject, self.view.currentIndex().row())
					if action == addCircle:
						self.createCircleWidget(selectedProject, self.view.currentIndex().row())
					if action == addExtrusion:
						self.createExtrusionWidget(selectedProject, self.view.currentIndex().row())
				elif self.view.currentIndex().parent().row() == 3:
					menu = QtGui.QMenu()
					deleteMaterial = menu.addAction("Delete Material")
					action = menu.exec_(self.view.mapToGlobal(position))
					if action == deleteMaterial:
						self.projects[selectedProject].refreshPlot = True
						del self.projects[selectedProject].materials[self.view.currentIndex().row()]
						self.projects[selectedProject].item.child(3).removeRow(self.view.currentIndex().row())
				elif self.view.currentIndex().parent().row() == 4:
					menu = QtGui.QMenu()
					if self.projects[selectedProject].excitations[self.view.currentIndex().row()].type == 'Port Plane':
						addWaveguidePort = menu.addAction("Add Waveguide Port")
					deleteGeometry = menu.addAction("Delete Excitation")
					action = menu.exec_(self.view.mapToGlobal(position))
					if action == deleteGeometry:
						del self.projects[selectedProject].excitations[self.view.currentIndex().row()]
						self.projects[selectedProject].item.child(4).removeRow(self.view.currentIndex().row())
						self.mayavi.visualization.update_plot(self.projects[selectedProject])
					elif action == addWaveguidePort:
						self.createWaveguidePort(selectedProject, self.view.currentIndex().row())
				elif self.view.currentIndex().parent().row() == 6:
					menu = QtGui.QMenu()
					deletePlot = menu.addAction("Delete Plot")
					action = menu.exec_(self.view.mapToGlobal(position))
					if action == deletePlot:
						self.projects[selectedProject].refreshPlot = True
						del self.projects[selectedProject].plots[self.view.currentIndex().row()]
						self.projects[selectedProject].item.child(6).removeRow(self.view.currentIndex().row())

			elif hierarchyLevel == 3:
				menu = QtGui.QMenu()
				deleteGeometry = menu.addAction("Delete Geometry")
				action = menu.exec_(self.view.mapToGlobal(position))
				if action == deleteGeometry:
					self.projects[selectedProject].refreshPlot = True
					del self.projects[selectedProject].geometries[self.view.currentIndex().parent().row()].geometries[self.view.currentIndex().row()]
					self.projects[selectedProject].item.child(2).child(self.view.currentIndex().parent().row()).removeRow(self.view.currentIndex().row())


	def newSelection(self, selected, deselected):		# If selection changed from the qtreeview
		plot = '3D'

		hierarchyLevel = 0
		node = self.view.currentIndex()
		while not(node.parent().data() is None):
			node = node.parent()
			hierarchyLevel += 1

		selectedProject = 0			# Must be outside self.hold block
		selectedMenu = 0
		selectedSubMenu = 0
		# if hierarchyLevel == 2:
		# 	selectedProject = self.view.currentIndex().parent().parent().row()
		# 	selectedGeometry = self.view.currentIndex().row()
		# 	if self.projects[selectedProject].geometries[selectedGeometry].type == 'Workplane':
		# 		plot = '2D'
		if hierarchyLevel == 3:
			selectedProject = self.view.currentIndex().parent().parent().parent().row()
			selectedMenu = self.view.currentIndex().parent().parent().row()
			selectedSubMenu = self.view.currentIndex().parent().row()
			if selectedMenu == 2 and self.projects[selectedProject].geometries[selectedSubMenu].type == 'Workplane':
				plot = '2D'
			if selectedMenu == 4 and self.projects[selectedProject].excitations[selectedSubMenu].type == 'Port Plane':
				plot = '2D'

		if plot == '3D':
			if self.mayavi.isVisible() == False:
				self.projects[selectedProject].refreshPlot = True
				self.mayavi.visualization.update_plot(self.projects[selectedProject])
			self.mayavi.show()
			self.plotter2D.widget.hide()
		else:
			if self.plotter2D.widget.isVisible() == False:
				if selectedMenu == 2:
					self.plotter2D.updatePlot(self.projects[selectedProject].geometries[selectedSubMenu])
				elif selectedMenu == 4:
					self.plotter2D.updatePlot(self.projects[selectedProject].excitations[selectedSubMenu])
			self.mayavi.hide()
			self.plotter2D.widget.show()

		if self.hold == False:							# When creating a new geometry the selection is also changed, don't interrupt this process
			hierarchyLevel = 0
			node = self.view.currentIndex()
			while not(node.parent().data() is None):
				node = node.parent()
				hierarchyLevel += 1

			if hierarchyLevel == 1:
				selectedProject = self.view.currentIndex().parent().row()
				if self.view.currentIndex().row() == 0:		# Grid selected
					self.createGridWidget(selectedProject)
				if self.view.currentIndex().row() == 5:		# Solver selected
					self.createSolverWidget(selectedProject)
				if self.view.currentIndex().row() == 6:		# Result selected
					self.createResultWidget(selectedProject)

			if hierarchyLevel == 2:
				selectedProject = self.view.currentIndex().parent().parent().row()
				if self.view.currentIndex().parent().row() == 2:		# if Geometry selected
					if self.hold == False:
						selectedGeometry = self.view.currentIndex().row()
						selection = self.projects[selectedProject].geometries[selectedGeometry]
						if selection.type == 'Block':
							self.createBlockWidget(selectedProject, selection)
						if selection.type == 'Workplane':
							self.createWorkplaneWidget(selectedProject, selection)
				elif self.view.currentIndex().parent().row() == 3:
					if self.hold == False:
						selectedMaterial = self.view.currentIndex().row()
						self.createMaterialWidget(selectedProject, self.projects[selectedProject].materials[selectedMaterial])
				elif self.view.currentIndex().parent().row() == 4:
					if self.hold == False:
						selectedExcitation = self.view.currentIndex().row()
						self.createPortPlaneWidget(selectedProject, self.projects[selectedProject].excitations[selectedExcitation])
				elif self.view.currentIndex().parent().row() == 6:
					selectedPlot = self.view.currentIndex().row()
					selection = self.projects[selectedProject].plots[selectedPlot]
					if self.hold == False:
						if selection.type == 'Scalar Cutplane':
							self.createScalarCutplaneWidget(selectedProject, selection)
						elif selection.type == 'Scalar Volume':
							self.createScalarVolumeWidget(selectedProject, selection)
						elif selection.type == 'Vector Volume':
							self.createVectorVolumeWidget(selectedProject, selection)
						elif selection.type == 'Vector Cutplane':
							self.createVectorCutplane(selectedProject, selection)
						elif selection.type == 'Contour':
							self.createContourWidget(selectedProject, selection)

			if hierarchyLevel == 3:
				selectedProject = self.view.currentIndex().parent().parent().parent().row()
				selectedMenu = self.view.currentIndex().parent().parent().row()
				selectedSubMenu = self.view.currentIndex().parent().row()
				if self.hold == False:
					if selectedMenu == 2:
						selection = self.projects[selectedProject].geometries[selectedSubMenu].geometries[self.view.currentIndex().row()]
						if selection.type == 'Polygon 2D':
							self.createPolygon2DWidget(selectedProject, selectedSubMenu, selection)
						if selection.type == 'Circle':
							self.createCircleWidget(selectedProject, selectedSubMenu, selection)
						if selection.type == 'Extrusion':
							self.createExtrusionWidget(selectedProject, selectedSubMenu, selection)
					if selectedMenu == 4:
						selection = self.projects[selectedProject].excitations[selectedSubMenu].sources[self.view.currentIndex().row()]
						if selection.type == 'Waveguide Port':
							self.createWaveguidePort(selectedProject, selectedSubMenu, selection)

	def addTab(self, tab, name):
		if self.hSplitter.count() == 3:
			self.geometrySettings = QtGui.QTabWidget()
			self.geometrySettings.setTabsClosable(True)
			self.geometrySettings.tabCloseRequested.connect(self.closeTab)

			self.hSplitter.insertWidget(1, self.geometrySettings)

		self.geometrySettings.removeTab(0)
		self.geometrySettings.addTab(tab, name)

	def closeTab(self, currentIndex):
		self.geometrySettings.widget(currentIndex).deleteLater()
		self.geometrySettings.removeTab(currentIndex)
		self.hSplitter.widget(1).deleteLater()

	def createScalarVolumeWidget(self, selectedProject, scalarVolume=None):
		self.hold = True
		self.scalarVolumeWidget = ScalarVolumeWidget(self.projects[selectedProject], self.mayavi, self.view, self.model, scalarVolume)
		self.addTab(self.scalarVolumeWidget.tab, self.projects[selectedProject].name+"/Results/"+self.scalarVolumeWidget.scalarVolume.name)
		self.hold = False

	def createScalarCutplaneWidget(self, selectedProject, scalarCutplane=None):
		self.hold = True
		self.scalarCutplaneWidget = ScalarCutplaneWidget(self.projects[selectedProject], self.mayavi, self.view, self.model, scalarCutplane)
		self.addTab(self.scalarCutplaneWidget.tab, self.projects[selectedProject].name+"/Results/"+self.scalarCutplaneWidget.scalarCutplane.name)
		self.hold = False

	def createVectorVolumeWidget(self, selectedProject, vectorVolume=None):
		self.hold = True
		self.vectorVolumeWidget = VectorVolumeWidget(self.projects[selectedProject], self.mayavi, self.view, self.model, vectorVolume)
		self.addTab(self.vectorVolumeWidget.tab, self.projects[selectedProject].name+"/Results/"+self.vectorVolumeWidget.vectorVolume.name)
		self.hold = False

	def createVectorCutplane(self, selectedProject, vectorCutplane=None):
		self.hold = True
		self.vectorCutplaneWidget = VectorCutplaneWidget(self.projects[selectedProject], self.mayavi, self.view, self.model, vectorCutplane)
		self.addTab(self.vectorCutplaneWidget.tab, self.projects[selectedProject].name+"/Results/"+self.vectorCutplaneWidget.vectorCutplane.name)
		self.hold = False

	def createContourWidget(self, selectedProject, contour=None):
		self.hold = True
		self.contourWidget = ContourWidget(self.projects[selectedProject], self.mayavi, self.view, self.model, contour)
		self.addTab(self.contourWidget.tab, self.projects[selectedProject].name+"/Results/"+self.contourWidget.contour.name)
		self.hold = False

	def createBlockWidget(self, selectedProject, block=None):
		self.hold = True
		self.blockWidget = BlockWidget(self.projects[selectedProject], self.mayavi, self.view, self.model, block)
		self.addTab(self.blockWidget.tab, self.projects[selectedProject].name+"/Geometry/"+self.blockWidget.block.name)
		self.hold = False

	def createWorkplaneWidget(self, selectedProject, workplane=None):
		self.hold=True
		self.workplaneWidget = WorkplaneWidget(self.projects[selectedProject], self.mayavi, self.plotter2D, self.view, self.model, workplane)
		self.addTab(self.workplaneWidget.tab, self.projects[selectedProject].name+"/Geometry/"+self.workplaneWidget.workplane.name)
		self.hold = False

	def createPolygon2DWidget(self, selectedProject, workplaneIndex, polygon2D=None):
		self.hold=True
		self.polygon2DWidget = Polygon2DWidget(self.projects[selectedProject], workplaneIndex, self.view, self.model, polygon2D)
		self.addTab(self.polygon2DWidget.tab, self.projects[selectedProject].name+"/Geometry/"+\
			self.projects[selectedProject].geometries[workplaneIndex].name+"/"+self.polygon2DWidget.polygon2D.name)
		self.hold=False

	def createCircleWidget(self, selectedProject, workplaneIndex, circle=None):
		self.hold = True
		self.circleWidget = CircleWidget(self.projects[selectedProject], workplaneIndex, self.view, self.model, circle)
		self.addTab(self.circleWidget.tab, self.projects[selectedProject].name+"/Geometry/"+\
			self.projects[selectedProject].geometries[workplaneIndex].name+"/"+self.circleWidget.circle.name)
		self.hold = False

	def createExtrusionWidget(self, selectedProject, workplaneIndex, extrusion=None):
		self.hold = True
		self.extrusion2DWidget = ExtrusionWidget(self.projects[selectedProject], workplaneIndex, self.view, self.model, extrusion)
		self.addTab(self.extrusion2DWidget.tab, self.projects[selectedProject].name+"/Geometry/"+\
			self.projects[selectedProject].geometries[workplaneIndex].name+"/"+self.extrusion2DWidget.extrusion.name)
		self.hold = False

	def createMaterialWidget(self, projectIndex, material=None):
		self.hold = True
		self.materialWidget = MaterialWidget(self.projects[projectIndex], self.mayavi, self.view, self.model, material)
		self.addTab(self.materialWidget.tab, self.projects[projectIndex].name+"/Material/"+self.materialWidget.material.name)
		self.hold = False

	def createPortPlaneWidget(self, projectIndex, portPlane=None):
		self.hold = True
		self.portPlaneWidget = PortPlaneWidget(self.projects[projectIndex], self.mayavi, self.plotter2D, self.view, self.model, portPlane)
		self.addTab(self.portPlaneWidget.tab, self.projects[projectIndex].name+"/Excitation/"+self.portPlaneWidget.portPlane.name)
		self.hold = False

	def createWaveguidePort(self, projectIndex, excitationIndex, waveguidePort=None):
		self.hold = True
		self.waveguideWidget = WaveguideWidget(self.projects[projectIndex], self.view, self.model, excitationIndex, waveguidePort)
		self.addTab(self.waveguideWidget.tab, self.waveguideWidget.waveguidePort.name)
		self.hold = False

	def createResultWidget(self, projectIndex):
		self.projects[projectIndex].result.createTab()
		self.addTab(self.projects[projectIndex].result.tab, self.projects[projectIndex].name+"/Results")

	def createGridWidget(self, projectIndex):
		self.gridWidget = GridWidget(self.projects[projectIndex])
		self.addTab(self.gridWidget.tab, self.projects[projectIndex].name+"/Grid")

	def createSolverWidget(self, projectIndex):
		self.solverWidget = SolverWidget(self.projects[projectIndex])
		self.addTab(self.solverWidget.tab, self.projects[projectIndex].name+"/Solver")

if __name__ == '__main__':
	app = QtGui.QApplication.instance()
	window = FDTD()
	app.exec_()