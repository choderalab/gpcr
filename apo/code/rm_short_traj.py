import mdtraj as md
import os
import glob

filenames = glob.glob('trj*')

lengths = []

for filename in filenames:
	f =  md.open(filename)
	lengths.append(len(f))
	f.close()

for k in enumerate(lengths):
	if k[1] > 20:
		print k
		os.symlink('../inverse_agonist_b2ar_trajectories' % filenames[k[0]], '../inverse_agonist_b2ar_trajectories/%s' % filenames[k[0]])


