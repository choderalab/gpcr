
import mixtape.featurizer, mixtape.tica, mixtape.cluster, mixtape.markovstatemodel, mixtape.ghmm
import numpy as np
import mdtraj as md
import sklearn.pipeline, sklearn.externals.joblib
import mixtape.utils
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import glob

n_components = 3
n_states = 3
lag_time = 1

ACTIVE = md.load_pdb('../../../GPCRexacycle/GPCR_NatureChemistry/reference-structure/active-union.pdb')
INACTIVE = md.load_pdb('../../../GPCRexacycle/GPCR_NatureChemistry/reference-structure/inactive-union.pdb')
APO = md.load_pdb('../../../GPCRexacycle/GPCR_NatureChemistry/reference-structure/apo_snapshot.pdb')

trajectories = [md.load(filename, top=APO) for filename in filenames]
train = trajectories[0::2]
test = trajectories[1::2]

#atom selection
I121_F282 = table[(table['resSeq'] == 121) | (table['resSeq'] == 282)].index  
NPxxY = {'resSeq': [322,323,324,325,326,327]}
NPxxY_region = table[table.isin(NPxxY).any(1)].index
R131_L272 = table[(table['resSeq'] == 131) & (table['name'] == 'CA') | (table['resSeq'] == 272) & (table['name'] == 'CA')].index
binding_pocket = {'resSeq': [109, 110, 113, 114, 117, 118, 191, 192, 193, 195, 199, 200, 203, 204, 207, 286, 289, 290, 293, 305, 308, 309, 312, 316]}
binding_pocket_RMSD = table[table.isin(binding_pocket).any(1)].index

# setting features for featurizer 
distance = np.array([R131_L272])
Hel3_Hel6_dist = mixtape.featurizer.AtomPairsFeaturizer(distance)
NPxxY_RMSD_active = mixtape.featurizer.RMSDFeaturizer(ACTIVE_union, NPxxY_region)
NPxxY_RMSD_inactive = mixtape.featurizer.RMSDFeaturizer(INACTIVE_union, NPxxY_region)
Connector_RMSD_active = mixtape.featurizer.RMSDFeaturizer(ACTIVE_union, I121_F282)
Connector_RMSD_inactive = mixtape.featurizer.RMSDFeaturizer(INACTIVE_union, I121_F282)
bind_pock_RMSD_active = mixtape.featurizer.RMSDFeaturizer(ACTIVE_union, binding_pocket_RMSD)
bind_pock_RMSD_inactive = mixtape.featurizer.RMSDFeaturizer(INACTIVE_union, binding_pocket_RMSD)

# setting transformers
union = mixtape.featurizer.TrajFeatureUnion([("Helical_distance", Hel3_Hel6_dist), ("NPxxYrmsd_active", NPxxY_RMSD_inactive),("NPxxYrmsd_inactive", NPxxYrmsd_inactive), ("connecector_active", Connector_RMSD_active),("connector_inactive", connector_RMSD_inactive), ("binding pocket_active", bind_pock_RMSD_active), ("bind_pock_in", bind_pock_RMSD_inactive)])


tica = mixtape.tica.tICA(n_components=n_components, lag_time=lag_time)
subsampler = mixtape.utils.Subsampler(lag_time=lag_time)
cluster = mixtape.cluster.GMM(n_components=n_states, covariance_type='diag')
feature_pipeline = sklearn.pipeline.Pipeline([('features', union), ('tica', tica)])
cluster_pipeline = sklearn.pipeline.Pipeline([('features', union), ('tica', tica), ('cluster', cluster)])
pipeline = sklearn.pipeline.Pipeline([("features", union), ('tica', tica), ("subsampler", subsampler), ("cluster", cluster), ("msm", msm)])

pipeline.fit(train)
print pipeline.score(train), pipeline.score(test)

X_all = feature_pipeline.transform(trajectories)
q = np.concatenate(X_all)

S_all = cluster_pipeline.transform(trajectories)
covars_ = cluster.covars_.diagonal(axis1=0, axis2=1)
for n in range(3):
	for i, j in [(0,1)]:
	    	fig = plt.figure()
    		plt.hexbin(q[:,i],q[:,j], bins = 'log')
    		plt.errorbar(cluster.means_[:, i], cluster.means_[j], xerr=covars_[i]**0.5, yerr=covars_[:,j]**0.5,fmt='kx', linewidth=4)
    		for k,data in enumerate(X_all):
        		if S_all[k][0]==n:
            			plt.plot(data[:,i], data[:,j], 'k')
    		fig.savefig('../figures/gpcr_tics_%d_state%d.pdf'%(n_choose, n))


ind = msm.draw_samples(S_all, 3)
samples = mixtape.utils.map_drawn_samples(ind, trajectories)

for i in range(n_states):
	for k, t in enumerate(samples[i]):
		t.save("../pdbs/state%d-%d.pdb" % (i, k)

#Plot 2D hexbins of metrics
targets = {'Hel_3-6_dist':[], 'NPxxY_RMSD_act': [], 'NPxxy_RMSD_inact':[], 'Conn_RMSD_act':[], 'Conn_RMSD_inact':[], 'bind_RMSD_act':[], 'bind_RMSD_inact':[]}

targetsi_all = ['Hel_3-6_dist_all', 'NPxxY_RMSD_act_all', 'NPxxy_RMSD_inact_all', 'Conn_RMSD_act_all', 'Conn_RMSD_inact_all', 'bind_RMSD_act_all', 'bind_RMSD_inact_all']

for i, target in enumerate(zip(targets.iteritems(), targets_all)):
	target[0][-1].append([union.transformer_list[i][-1].partial_transform(t) for t in trajectories])
	target[1][-1] = np.concatenate(target[0][-1]:)
        fig = plt.figure()
	plt.hist(target[	
fig = plt.figure()
plt.hist(Hel_dist_all, bin=100)
plt.title('Helix 3-6 distance')
plt.xlabel('distance')
plt.ylabel('counts')
fig.savefig('../figures/hist_%s.pdb)
