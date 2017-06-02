"""
Routines to support machine learning for genotyping.
"""


import numpy as np

from sklearn.model_selection import StratifiedKFold

# Labels of features used in training
GT_FEATURES = (
    'SVTYPE', 'SVLEN', 'REF_COUNT', 'ALT_COUNT',
    'REF_CLIP', 'ALT_CLIP', 'PRIMARY_CLIP',
    'N_INSERT', 'INSERT_LOWER', 'INSERT_UPPER'
)


class GtModel:
    def __init__(self, model, scaler):
        self.model = model
        self.scaler = scaler

    def genotype(self, features):
        pass

def features_table_to_unscaled_array(features_table):
    pass

def stratify_folds(labels, k):
    """
    Return an n-x-2 matrix with one row for each feature (n). The first column is the row index of the feature array
    (X) for a feature, and the second column is a fold-set it belongs to (0 through k - 1).

    `k` Folds are generated over the labels preserving the proportion of each label in each fold. `labels` should
    include the learned value (genotype), but it can have other information concatenated to it. For example,
    'HET-DEL-InTRF' will preserve the proportion of each genotype, SV type, and TRF members in each fold.

    :param labels: An array with one element for each variant (row of X).
    :param k: Number of folds.

    :return: An n-x-2 matrix of feature indices and fold membership.
    """

    if k < 1:
        raise ValueError('Number of folds, k, must be positive: {}'.format(k))

    # Split indices into k folds
    skf_test_model = StratifiedKFold(shuffle=True, n_splits=k)

    skf_test = skf_test_model.split(np.zeros(len(labels)), labels)
    skf_test = [skf_set[1] for skf_set in skf_test]

    for cv_set in range(k):
        skf_test[cv_set] = np.expand_dims(skf_test[cv_set], axis=1)
        skf_test[cv_set] = np.append(skf_test[cv_set], np.repeat(cv_set, skf_test[cv_set].shape[0]).T[:, np.newaxis],
                                     axis=1)

    # Create one n x 2 array
    #  * where n is the number of samples
    #  * column 1 is the variant index
    #  * column 2 is the set the variant belongs to

    fold_array = np.concatenate(skf_test, axis=0)

    # Sort by variant index
    fold_array = fold_array[fold_array[:, 0].argsort(), :]

    # Return array
    return fold_array


def cv_set_iter(fold_array):
    """
    Return an iterator over fold sets. Each element is a two-element tuple for each fold of the cross-validation set.
    The first tuple element is a list of indices for the training set, and the last tuple is a list of indices for the
    test set.

    :param fold_array: A 2-dim array with variants as rows. The first column is a variant index (row index of the
        feature array, X) and the second column is a k-fold set it belongs to (typically 0 to k - 1).

    :return: In iterator of tuples containing two one-dimensional numpy arrays. Each tuple has a set of indices to
        train on (first element) and a set of indices to test on (second element). The indices are indexes to rows
        of the feature array (X).
    """
    for cv_set in sorted(set(fold_array[:, 1])):
        yield (
            fold_array[fold_array[:, 1] != cv_set, 0],
            fold_array[fold_array[:, 1] == cv_set, 0]
        )
