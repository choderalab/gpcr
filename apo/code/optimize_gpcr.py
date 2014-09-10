import glob
import mdtraj as md 
import mixtape.subset_featurizer, mixtape.tica, mixtape.feature_selection
import sklearn.pipeline, sklearn.externals.joblib

n_iter = 1000
n_choose = 2
n_components = 2
stride = 1
lag_time = 1

PDB =  md.load_pdb('../../../GPCRexacycle/GPCR_NatureChemistry/reference-structures/apo_snapshot.pdb')

filenames = glob.glob("../../../GPCRexacycle/dcd_trajectories/apo_b2ar_processed/trj*")

print 'loading trajectorie'
train = [md.load(filename, top=PDB) for filename in filenames[::2]]
print 'starting featurizer'

featurizer = mixtape.subset_featurizer.guess_featurizers(train[0][0], n_choose)
model = mixtape.tica.tICA(lag_time=lag_time, n_components=n_components)
tica_optimizer = mixtape.feature_selection.Optimizer(featurizer, model, n_iter)
	
print 'starting optimization'

featurizer = tica_optimizer.optimize(train)


sklearn.externals.joblib.dump(tica_optimizer.featurizer, "../joblib_dump/featurizer%d-%d.job" % (n_components, n_choose), compress=True)


