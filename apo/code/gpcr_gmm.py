import mixtape.featurizer, mixtape.tica, mixtape.cluster, mixtape.markovstatemodel, mixtape.ghmm
import numpy as np
import mdtraj as md
import sklearn.pipeline, sklearn.externals.joblib
import mixtape.utils
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import glob

n_choose = 10
stride = 1
lag_time = 1
iteration = 0

PDB =  md.load_pdb('../../../GPCRexacycle/GPCR_NatureChemistry/reference-structures/apo_snapshot.pdb')
filenames = glob.glob('../../../GPCRexacycle/dcd_trajectories/apo_b2ar_processed/trj*')

trajectories = [md.load(filename, top=PDB) for filename in filenames]
train = trajectories[0::2]
test = trajectories[1::2]

featurizer = sklearn.externals.joblib.load("../joblib_dump/featurizer%d-%d.job" % (iteration, n_choose))

n_components = 3
n_states = 3
tica = mixtape.tica.tICA(n_components=n_components, lag_time=lag_time)
subsampler = mixtape.utils.Subsampler(lag_time=lag_time)
msm = mixtape.markovstatemodel.MarkovStateModel(n_timescales=n_components)
cluster = mixtape.cluster.GMM(n_components=n_states, covariance_type='full')
feature_pipeline = sklearn.pipeline.Pipeline([("features", featurizer), ('tica', tica)])
cluster_pipeline = sklearn.pipeline.Pipeline([("features", featurizer), ('tica', tica), ("cluster", cluster)])
pipeline = sklearn.pipeline.Pipeline([("features", featurizer), ('tica', tica), ("subsampler", subsampler), ("cluster", cluster), ("msm", msm)])

pipeline.fit(train)
pipeline.score(train), pipeline.score(test)


X_all = feature_pipeline.transform(trajectories)
q = np.concatenate(X_all)

S_all = cluster_pipeline.transform(trajectories)
covars_ = cluster.covars_.diagonal(axis1=1, axis2=2)
for n in range(3):
	for i, j in [(0,1)]:
	    	fig = plt.figure()
    		plt.hexbin(q[:,i],q[:,j], bins = 'log')
    		plt.errorbar(cluster.means_[:, i], cluster.means_[:,j], xerr=covars_[:,i]**0.5, yerr=covars_[:,j]**0.5,fmt='kx', linewidth=4)
    		for k,data in enumerate(X_all):
        		if S_all[k][0]==n:
            			plt.plot(data[:,i], data[:,j], 'k')
    		fig.savefig('../figures/gpcr_tics_%d_state%d.pdf'%(n_choose, n))


ind = msm.draw_samples(S_all, 3)
samples = mixtape.utils.map_drawn_samples(ind, trajectories)

for i in range(n_states):
    for k, t in enumerate(samples[i]):
        t.save("../pdbs/state%d-%d-%d.pdb" % (i, k, n_choose))

