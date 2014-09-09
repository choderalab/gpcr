import glob
import mixtape.featurizer, mixtape.tica, mixtape.cluster, mixtape.markovstatemodel
import mdtraj as md 
import sklearn.pipeline, sklearn.externals.joblib
import mixtape.utils
from mixtape import ghmm, subset_featurizer, selector

n_iter = 1000

n_choose = 100
stride = 1
lag_time = 1

PDB =  md.load_pdb('../../../GPCRexacycle/GPCR_NatureChemistry/reference-structures/apo_snapshot.pdb')

filenames = glob.glob("../../../GPCRexacycle/dcd_trajectories/apo_b2ar_processed/trj*")

print 'loading trajectorie'
train = [md.load(filename, top=PDB) for filename in filenames[::2]]
print 'starting featurizer'
	
featurizer = sklearn.externals.joblib.load("../joblib_dump/featurizer9-%d.job" % n_choose)
print 'starting optimization'
tica_optimizer = mixtape.selector.TICAOptimizer(featurizer, train, lag_time=lag_time)
tica_optimizer.optimize(n_iter, train)


sklearn.externals.joblib.dump(tica_optimizer.featurizer, "../joblib_dump/featurizer10-%d.job" % n_choose, compress=True)


