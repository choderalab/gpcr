import mdtraj as md
import msmbuilder as msmb
import glob
import itertools
import numpy as np
import mixtape.subset_featurizer

figure_path = "/home/kyleb/src/kyleabeauchamp/tICABenchmark/figures/"

def load_trajectories(stride=1):
        trj0 = md.load("./system.subset.pdb")
        filenames = sorted(glob.glob("./Trajectories/*.h5"), key=msmb.utils.keynat)
        if len(filenames) > 0:
            trajectories = [md.load(filename, stride=stride) for filename in filenames]
            filenames = np.array(filenames)
        else:
            filenames = ["./Trajectories/1am7_%d.dcd" % i for i in range(15)]
            trajectories = [md.load(filename, top=trj0) for filename in filenames]
        return trj0, trajectories, filenames

n_tics = 5

n_em_iter = 10
n_init = 10


def build_full_featurizer(trj0, n_choose, rmsd_trajectory=None):
    atom_indices, pair_indices = mixtape.subset_featurizer.get_atompair_indices(trj0, keep_atoms=None, exclude_atoms=np.array([]))
    
    atom_featurizer1m = mixtape.subset_featurizer.SubsetAtomPairs(pair_indices, trj0, exponent=-1.0)
    atom_featurizer2m = mixtape.subset_featurizer.SubsetAtomPairs(pair_indices, trj0, exponent=-2.0)
    atom_featurizer1p = mixtape.subset_featurizer.SubsetAtomPairs(pair_indices, trj0, exponent= 1.0)
    
    prod_featurizer = mixtape.subset_featurizer.SubsetProductFeaturizer(mixtape.subset_featurizer.SubsetAtomPairs(pair_indices, trj0, exponent=-1.0), mixtape.subset_featurizer.SubsetAtomPairs(pair_indices, trj0, exponent=-1.0))
    
    cosphi = mixtape.subset_featurizer.SubsetCosPhiFeaturizer(trj0)
    sinphi = mixtape.subset_featurizer.SubsetSinPhiFeaturizer(trj0)
    cospsi = mixtape.subset_featurizer.SubsetCosPsiFeaturizer(trj0)
    sinpsi = mixtape.subset_featurizer.SubsetSinPsiFeaturizer(trj0)
    if rmsd_trajectory is None:
        featurizer = mixtape.subset_featurizer.SubsetFeatureUnion([("pairs1m", atom_featurizer1m), ("pairs2m", atom_featurizer2m), ("pairs1p", atom_featurizer1p), ("prod1m", prod_featurizer), ("cosphi", cosphi), ("sinphi", sinphi), ("cospsi", cospsi), ("sinpsi", sinpsi)])
    else:
        rmsd_featurizer = mixtape.subset_featurizer.SubsetRMSDFeaturizer(rmsd_trajectory)
        featurizer = mixtape.subset_featurizer.SubsetFeatureUnion([("pairs1m", atom_featurizer1m), ("pairs2m", atom_featurizer2m), ("pairs1p", atom_featurizer1p), ("prod1m", prod_featurizer), ("cosphi", cosphi), ("sinphi", sinphi), ("cospsi", cospsi), ("sinpsi", sinpsi), ("rmsd", rmsd_featurizer)])
    
    subsets = [[] for i in range(featurizer.n_featurizers)]
    subsets[0] = np.random.randint(0, atom_featurizer1m.n_max - 1, n_choose)
    featurizer.subsets = subsets    

    return featurizer
