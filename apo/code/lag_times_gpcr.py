import glob
import mixtape.featurizer, mixtape.tica, mixtape.cluster, mixtape.markovstatemodel
import mdtraj as md 
import sklearn.pipeline, sklearn.externals.joblib
import mixtape.utils
from parameters import build_full_featurizer
from mixtape import ghmm, selector, subset_featurizer, selector

n_iter = 1000

n_choose = 100
lag_times = [2,3,4,5]

PDB =  md.load_pdb('../../../GPCRexacycle/GPCR_NatureChemistry/reference-structures/apo_snapshot.pdb')

featurizer = build_full_featurizer(PDB, n_choose)

for lag_time in lag_times:
	
	filenames = glob.glob("../../../GPCRexacycle/dcd_trajectories/apo_b2ar_%d_stride/trj*" % lag_time)
	print ' loading apo_b2ar_%d_stride trajectories' % lag_time

	train = [md.load(filename, top=PDB) for filename in filenames[::2]]
	

	tica_optimizer = mixtape.selector.TICAOptimizer(featurizer, train, lag_time=lag_time)
	tica_optimizer.optimize(n_iter, train)


	sklearn.externals.joblib.dump(tica_optimizer.featurizer, "../joblib_dump/featurizer_lt_%d-%d.job" % (lag_time, n_choose), compress=True)


