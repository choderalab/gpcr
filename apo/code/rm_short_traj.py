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
		os.symlink('../agonist_b2ar/%s' % filenames[k[0]], '../agonist_b2ar_processed/%s' % filenames[k[0]])


