import math
import numpy as np

class halton():
	
	def __init__(self, num = None, dim = None, indices = None, bases = None):
		"""
		Halton quasi-random number generator.
		
		Example usage:
		
		.. code-block:: python
		
		   sequence = halton(10, 3)
		   out = sequence.evaluate()
		   
		Here out is the coordinates of the first 10 points in the 
		three-dimensional Halton sequence with bases 2, 3 and 4.
		   
		Parameters
		----------
		num : int
			Number of samples for the Halton sequence, starting at index 0.
		dim : int
			Dimensions of the Halton sequence, starting at 1D using base 2.
		indices : numpy array or list
			Specific indices of the Halton sequence.
		bases : numpy array or list
			Specific bases of the Halton sequence.
			
		Methods:
		--------
		evaluate :
			Returns a numpy array of values of the Halton sequences in bases
			and at indices provided in the class initialization.
		evaluateone(i, b) :
			Returns the value of the Halton sequence of dimensional basis ``b`` 
			at index ``i``.
			
		Author: Jonas Svenstrup Hansen, November 2017 jonas.svenstrup@gmail.com
		"""
		if ((num != None) or (dim != None)) and ((indices != None) or (bases != None)):
			raise ValueError('Parameters num and dim or indices and bases cannot be supplied simultaneously.')
		elif (num == None) and (dim == None) and (indices == None) and (bases == None):
			raise ValueError('Parameters num and dim or indices and bases must be supplied.')
		elif ((indices != None) and (bases != None)):
			self.indices = np.asarray(indices)
			self.bases   = np.asarray(bases)
		elif (num != None) and (dim != None):
			self.indices = np.asarray(range(num))
			self.bases   = np.asarray(range(dim)) + 2
		else:
			raise ValueError('Unknown input error. Parameters num and dim or indices and bases must be supplied.')
		
	def evaluate(self):
		"""
		Returns
		-------
		result : numpy array
			Values of the Halton sequence for indices (rows) and bases (cols).
		"""
		result = np.zeros([self.indices.size, self.bases.size])
		for (i, idx) in enumerate(self.indices):
			for (j, base) in enumerate(self.bases):
				result[i,j] = self.evaluateone(idx, base)
		return result
	
	def evaluateone(self, i, b):
		"""
		Evaluate single value Halton sequence.
		
		Parameters
		----------
		i : int
			Index value.
		b : int
			Base value. Must be greater than 0.
		
		Notes
		-----
		Implemented from the pseudocode in:
		
			https://en.wikipedia.org/wiki/Halton_sequence
		"""
		if b == 0:
			raise ValueError('Parameter b is 0. b must be greater than 0.')
		
		f = 1
		r = 0
		
		while i > 0:
			f = f/b
			r = r + f * (i % b)
			i = math.floor(i/b)
			
		return r


if __name__ == "__main__":
	import random
	import matplotlib.pyplot as plt
	plt.close('all')
	
	maxdim = 1e3
	bases = [2, 3]
	index = np.arange(0, maxdim, 1)
	
	# Halton sequence of maxdim points:
	sequence = halton(index, bases)
	coords_quasi = sequence.values
	
	# Random (Mersenne Twister) sequence of maxdim points:
	coords_pseudo = np.asarray([[random.random(),random.random()] for i in range(1000)])
	
	coordstypes = [coords_quasi, coords_pseudo]
	names = ['Quasirandom (Halton)', 'Pseudorandom (Mersenne Twister)']
	
	fig, axes = plt.subplots(ncols=2, sharex=True, sharey=True)
	for (ax, coords, name) in zip(axes, coordstypes, names):
		ax.plot(coords[:,0],coords[:,1], linestyle='None', marker='.',
			markerfacecolor='black', markeredgecolor='None')
		ax.set(adjustable='box-forced', aspect='equal')
		ax.set_title(name)
		ax.set_xlim([0, 1])
		ax.set_ylim([0, 1])
	plt.show()