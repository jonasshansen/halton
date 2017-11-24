import math
import numpy as np

class halton(object):
	
	def __init__(self, indices, bases):
		"""
		Halton quasi-random number generator.
		
		Parameters
		----------
		indices : numpy array
			Index values for the Halton sequence.
		bases : numpy array
			Dimensional bases for the Halton sequence. Default is 2.
			
		Methods:
		--------
		evaluate :
			Returns
		evaluateone(index, base) :
			Returns coordinate
			
		Author: Jonas Svenstrup Hansen, November 2017, jonas.svenstrup@gmail.com
		"""
		self.indices = np.asarray(indices)
		self.bases   = np.asarray(bases)
		self.values = self.evaluate()
		
	def evaluate(self):
		"""
		Returns
		-------
		result : numpy array
			Values of the Halton sequence for indices (rows) and bases (cols).
		"""
		result = np.zeros([len(self.indices), len(self.bases)])
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
	
