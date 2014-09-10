import glob
import mixtape.featurizer, mixtape.tica, mixtape.feature_selection 
import mdtraj as md 
import sklearn.pipeline, sklearn.externals.joblib
import mixtape.utils
from mixtape import subset_featurizer

n_iter = 1000
n_choose = 10
stride = 1
lag_time = 1
n_components = 3

PDB =  md.load_pdb('../../../GPCRexacycle/GPCR_NatureChemistry/reference-structures/active_crystal_reference.pdb')

filenames = glob.glob("../../../GPCRexacycle/dcd_trajectories/agonist_b2ar_processed/trj*")

print 'loading trajectories'
train = [md.load(filename, top=PDB) for filename in filenames[::2]]
print 'starting featurizer'
	
featurizer = mixtape.subset_featurizer(trajectories[0][0], n_choose)
print 'starting optimization'
model = mixtape.tica.tICA(lage_time = lag_time, n_components=n_components)
tica_optimizer = mixtape.feature_selection.Optimizer(featurizer, model, n_iter)

featurizer = tica_optimizer.optimize(train)


sklearn.externals.joblib.dump(featurizer, "../joblib_dump/featurizer%d-%d.job" % (n_components,n_choose),  compress=True)


