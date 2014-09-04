import glob
import mixtape.featurizer, mixtape.tica, mixtape.cluster, mixtape.markovstatemodel
import mdtraj as md 
import sklearn.pipeline, sklearn.externals.joblib
import mixtape.utils


n_iter = 1000

n_choose = 100
stride = 1
lag_time = 1

PDB =  md.load_pdb('../../GPCR_NatureChemistry/reference-structures/apo_snapshot.pdb')

filenames = glob.glob("../../dcd_trajectories/apo_b2ar_processed/trj*")

train = [md.load(filename, top=PDB) for filename in filenames[::2]]
for i in range(10):
	
	featurizer = sklearn.externals.joblib.load("./featurizer%d-%d.job" % (i,n_choose))

	tica_optimizer = mixtape.selector.TICAOptimizer(featurizer, train, lag_time=lag_time)
	tica_optimizer.optimize(n_iter, train)


	sklearn.externals.joblib.dump(tica_optimizer.featurizer, "./featurizer%d-%d.job" % (i+1, n_choose), compress=True)


