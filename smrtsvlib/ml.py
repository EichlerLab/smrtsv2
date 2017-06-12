"""
Routines to support machine learning for genotyping.
"""


import numpy as np
import pandas as pd

from sklearn.externals import joblib
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from sklearn.svm import SVC


# Labels of features used in training
GT_FEATURES = (
    'SVTYPE', 'SVLEN', 'REF_COUNT', 'ALT_COUNT',
    'REF_CLIP', 'ALT_CLIP', 'PRIMARY_CLIP',
    'N_INSERT', 'INSERT_LOWER', 'INSERT_UPPER'
)

# Genotype labels to integer map
GT_SVTYPE_TO_NUMERIC = {
    'INS': 0.0,
    'DEL': 1.0
}

# Genotype labels
GT_LABELS = ['HOM_REF', 'HET', 'HOM_ALT']


class GtModel:
    """
    Makes predictions on data with a trained model.
    """

    def __init__(self, predictor, scaler):
        """
        Create an object for genotype predictions.

        :param predictor: Trained predictor. If string, then it must be a path to a `joblib` file containing the
            predictor object.
        :param scaler: Scaler fit to the data used for training. If string, then it must be a path to a `joblib` file
            containing the scaler object.
        """

        # Check arguments
        if predictor is None:
            raise ValueError('Cannot load genotyper predictor `None`')

        if scaler is None:
            raise ValueError('Cannot load feature scaler `None`')

        if isinstance(predictor, str):
            predictor = joblib.load(predictor)

        if isinstance(scaler, str):
            scaler = joblib.load(scaler)

        if not isinstance(predictor, SVC):
            raise ValueError('Predictor must be class sklearn.svm.SVC: Found "{}"'.format(type(predictor)))

        if not isinstance(scaler, StandardScaler):
            raise ValueError(
                'Scaler must be class sklearn.preprocessing.StandardScaler: Found "{}"'.format(type(scaler))
            )

        # Set fields
        self.predictor = predictor
        self.scaler = scaler

    def genotype(self, features_table):
        """
        Get genotype labels for each variant.

        :param features_table: A name to the features table, a table of loaded features from the features table, or
            genotype features extracted and scaled.

        :return: An array with a genotype calls for each variant.
        """

        # Get features
        if type(features_table) != np.ndarray:
            X = features_to_array(features_table, self.scaler)
        else:
            X = features_table

        # Predict
        return self.predictor.predict(X)

    def density(self, features_table):
        """
        Get density estimation for each variant.

        :param features_table: A name to the features table, a table of loaded features from the features table, or
            genotype features extracted and scaled.

        :return: A DataFrame with one row for each variant and a column with a density estimation for each genotype
            call. Columns are labeled by genotype call.
        """

        # Get features
        if type(features_table) != np.ndarray:
            X = features_to_array(features_table, self.scaler)
        else:
            X = features_table

        # Predict
        df_density = pd.DataFrame(self.predictor.predict_proba(X), columns=self.predictor.classes_)
        df_density = df_density.loc[:, GT_LABELS]

        return df_density

    def genotype_and_density(self, features_table):
        """
        Get most likely genotype calls and density esitmation for each variant.

        :param features_table: A name to the features table, a table of loaded features from the features table, or
            genotype features extracted and scaled.

        :return: A tuple of genotypes (first element) and densities (second element).
        """

        # Get features
        if type(features_table) != np.ndarray:
            X = features_to_array(features_table, self.scaler)
        else:
            X = features_table

        # Return genotype and density
        return self.genotype(X), self.density(X)


def get_cv_score_table(clf):
    """
    Get a table (DataFrame) of CV parameters and scores for each combination.

    :param clf: Cross-validation object (GridSearchCV)
    :return:
    """

    # Create data frame
    df = pd.DataFrame(list(clf.cv_results_['params']))

    # Add test scores
    df['rank'] = clf.cv_results_['rank_test_score']
    df['test_mean'] = clf.cv_results_['mean_test_score']
    df['test_sd'] = clf.cv_results_['std_test_score']

    # Add scores over training data
    df['train_mean'] = clf.cv_results_['mean_train_score']
    df['train_sd'] = clf.cv_results_['std_train_score']

    # Add time metrics (s)
    df['fit_time_mean'] = clf.cv_results_['mean_fit_time']
    df['fit_time_sd'] = clf.cv_results_['std_fit_time']

    df['score_time_mean'] = clf.cv_results_['mean_score_time']
    df['score_time_sd'] = clf.cv_results_['std_score_time']

    return df


def features_to_array(features_table, scaler):
    """
    Get a scaled feature array and a scaler object as a two element tuple (in that order).

    The input DataFrame must contain columns in list `GT_FEATURES`. All other columns are ignored.

    `scaler` must be `None` when features are being prepared for training. The scaler returned by this function when
    setting up for training must be given to the function to transform data for analysis by the model. Therefore,
    this scaler should be serialized along with the model to predict labels for the genotyper.

    :param features_table: Pandas DataFrame of a feature table or a string path of where the feature table can be
        loaded from.
    :param scaler: Scaler object if one already exists. If `None`, a scaler is created.

    :return: A tuple of the scaled feature array ready for model training or evaluation (first element) and the scaler
        used to scale the data (second element). If the `scaler` argument is not `None`, then the secord element is this
        scaler, and if it is `None`, then the second argument is the scaler that was fit to this set of features.
    """

    # Check arguments
    X = features_to_unscaled_matrix(features_table)
    return scaler.transform(X)


def features_to_unscaled_matrix(features_table):
    """
    Read features and transform to a matrix of features with the correct columns and in the correct order. This
    function does not scale the features, which must be done before training or predicting.

    :param features_table: Pandas DataFrame of a feature table or a string path of where the feature table can be
        loaded from.

    :return: Numeric feature matrix, unscaled.
    """

    # Check arguments
    if features_table is None:
        raise ValueError('Cannot convert features table: None')

    if isinstance(features_table, str):
        features_table = pd.read_table(features_table, header=0)

    if not isinstance(features_table, pd.DataFrame):
        raise ValueError(
            'Argument "features_table" must be a Pandas DataFrame or a string path to a features file that can be '
            'loaded into a DataFrame: Found type "{}"'.format(type(features_table)))

    # Load
    X = features_table[list(GT_FEATURES)].copy()

    # Cast all features to float64
    X['SVTYPE'] = X['SVTYPE'].apply(lambda label: GT_SVTYPE_TO_NUMERIC[label])  # SVTYPE label numeric representation
    X = X.astype(np.float64)

    # Return feature matrix
    return X


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
