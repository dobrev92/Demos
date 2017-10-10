#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
#from Gaussian import *

class NearestInterp:
	data = 0
	dim = 0

	def __init__(self, dataSet):
		self.data = dataSet
		self.dim = np.shape(self.data)[1]

	def interpolate(self, fInput):
		ret = 0
		prev_index = 0
		next_index = 0

		#seek previous and next indices
		for i in range(self.dim):
			if (self.data[0, i] > fInput):
				break
			prev_index = i
			next_index = i + 1

		if ((fInput - self.data[0, prev_index]) < (self.data[0, next_index] - fInput)):
			ret = self.data[1, prev_index]
		else:
			ret = self.data[1, next_index]

		return ret

#Linear interpolation(for equally spaced data)
class LinInterp:
	data = 0
	dim = 0

	def __init__(self, dataSet):
		self.data = dataSet
		self.dim = np.shape(self.data)[1]		

	def interpolate(self, fInput):
		ret = 0
		prev_index = 0
		next_index = 0

		#seek previous and next indices
		for i in range(self.dim):
			if (self.data[0, i] > fInput):
				break
			prev_index = i
			next_index = i + 1

		ret = self.data[1, prev_index] + (fInput - self.data[0, prev_index]) * (self.data[1, next_index] - self.data[1, prev_index]) / (self.data[0, next_index] - self.data[0, prev_index])
		return ret

#Lagrangian polynomial interpolation
class PolyInterp:
	data = 0

	def __init__(self, dataSet):
		self.data = dataSet

	def interpolate(self, fInput):
		sumation = 0
		product = 1
		dataSetSize = np.shape(self.data)[1]
		for j in range(dataSetSize):
			product = 1
			for i in range(dataSetSize):
				if (j != i):
					product *= (fInput - self.data[0][i]) / (self.data[0][j]- self.data[0][i])
			sumation += self.data[1][j] * product
		return sumation
	
#Cubic spline interpolation(for equally spaced data)
class SplineInterp:
	data = 0
	#difference between two points
	delta = 0
	#dimension of input data
	dim = 0
	#second order derivatives
	deriv = 0

	def __init__(self, dataSet):
		dim = 0
		mat = 0
		self.data = dataSet
		self.delta = self.data[0, 1] - self.data[0, 0]

		self.dim = np.shape(self.data)[1]
		mat = np.zeros((self.dim - 2, self.dim - 2))

		#free run-out condition
		mat[0, 0] = 2 / 3
		mat[0, 1] = 1 / 6
		mat[self.dim - 2 - 1, self.dim - 2 - 1] = 2 / 3
		mat[self.dim - 2 - 1, self.dim - 2 - 2] = 1 / 6

		#fill the matrix
		for i in range(self.dim - 4):
			mat[i + 1, i + 0] = 1 / 6
			mat[i + 1, i + 1] = 2 / 3
			mat[i + 1, i + 2] = 1 / 6

		s = np.zeros(self.dim - 2)
		for i in range(self.dim - 2):
			s[i] = (self.data[1, i + 2] - 2 * self.data[1, i + 1] + self.data[1, i]) / self.delta ** 2

		#solve the linear system
		mat = np.asmatrix(mat)
		inv = np.linalg.inv(mat)
		x = np.dot(inv, s)

		#fill the derivatives
		self.deriv = np.zeros(self.dim)
		for i in range(self.dim - 2):
			self.deriv[i + 1] = x[0, i]

	def interpolate(self, fInput):
		ret = 0
		term1 = 0
		term2 = 0
		term3 = 0
		prev_index = 0
		next_index = 0

		#seek previous and next indices
		for i in range(self.dim):
			if (self.data[0, i] > fInput):
				break
			prev_index = i
			next_index = i + 1

		term1 = (self.deriv[prev_index] / 6) * ((((self.data[0, next_index] - fInput) ** 3) / self.delta) - self.delta * (self.data[0, next_index] - fInput))
		term2 = (self.deriv[next_index] / 6) * ((((fInput - self.data[0, prev_index]) ** 3) / self.delta) - self.delta * (fInput - self.data[0, prev_index]))
		term3 = self.data[1, prev_index] * ((self.data[0, next_index] - fInput) / self.delta) + self.data[1, next_index] * ((fInput - self.data[0, prev_index]) / self.delta)
		ret = term1 + term2 + term3

		return ret

#main
def main():
	np.set_printoptions(precision=2, suppress=True, linewidth=120)
	dataSetInput = np.arange(0, 14, 1)
	dataSetOutput = np.sin(dataSetInput) #normalDistrib(5, 0.5, dataSetInput)
	dataSet = np.array([dataSetInput, dataSetOutput])
	test_set = np.arange(0, 12, 0.001)
	actual = np.sin(test_set) #normalDistrib(5, 0.5, test_set)
	test_set_size = np.shape(test_set)[0]
	out = np.zeros(test_set_size)
	interp = SplineInterp(dataSet)
	for x in range(test_set_size):
		single_out = interp.interpolate(test_set[x])
		out[x] = single_out
	
	#plot the result
	plt.plot(test_set, out)
	plt.plot(dataSetInput, dataSetOutput)
	#plt.plot(test_set, actual)
	plt.grid(True)
	plt.show()

main()
	
