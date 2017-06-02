
rule gt_learn_test_cv:
    input:
        tab='test/cv/tune/set_{cv_set}.tab'
    output:
        tab='test/cv/test/cv_scores_{cv_set}.tab'
    run:
        pass

# gt_learn_test_split_tune_k
#
# Split a test set into k tuning sets. K-fold cross-validation is performed on these sets to optimize hyper-parameters
# and then used to train a model for the test set.
rule gt_learn_test_split_tune_k:
    input:
        tab='test/cv/sets.tab',
        features='features/features.tab'
    output:
        tab='test/cv/tune/set_{cv_set}.tab'
    params:
        k=4
    run:

        # Read table for this set
        set_index = pd.read_table(input.tab, header=0)

        set_index = set_index.loc[set_index['test_set'] == int(wildcards.cv_set), :]
        set_index = set(set_index['variant'])

        # Get features
        # Read data
        features = pd.read_table(input.features, header=0)

        # Get feature labels (augmented with TRF annotation for stratification)
        labels = features.apply(lambda row: '{}-{}'.format(row['CALL'], 'TRF' if row['TRF'] else 'NoTRF'), axis=1)
        labels = labels[set_index]

        # Split indices into k folds
        skf_tune_model = StratifiedKFold(shuffle=True, n_splits=params.k)

        skf_tune_model = skf_tune_model.split(np.zeros(labels.shape[0]), labels)
        skf_tune = [skf_set[1] for skf_set in skf_tune_model]

        for cv_set in range(params.k):
            skf_tune[cv_set] = np.expand_dims(skf_tune[cv_set], axis=1)
            skf_tune[cv_set] = np.append(skf_tune[cv_set], np.repeat(cv_set, skf_tune[cv_set].shape[0]).T[:, np.newaxis], axis=1)

        # Create one n x 2 array
        #  * where n is the number of samples
        #  * column 1 is the variant index
        #  * column 2 is the set the variant belongs to

        fold_array = np.concatenate(skf_tune, axis=0)

        # Sort by variant index
        fold_array = fold_array[fold_array[:, 0].argsort(), :]

        # Write
        pd.DataFrame(fold_array, columns=('variant', 'tune_set')).to_csv(output.tab, sep='\t', index=False)


### TRF Annotations ###

# gt_learn_model_anno_trf
#
# Get TRF annotations from UCSC.
rule gt_learn_model_anno_trf:
    output:
        bed=temp('model/anno/trf/trf.bed')
    params:
        trf_url=CONFIG['trf_url']
    shell:
        """wget -qO- {params.trf_url} | """
        """gunzip | """
        """awk -vOFS="\\t" '{{print $2, $3, $4}}' | """
        """bedtools merge """
        """>{output.bed}"""


# Part of gt_learn_model_annotate:
        # Add TRF annotations
        trf_set = set(shell(
            """bedtools intersect -a {input.features} -b {input.trf_bed} -wa -u -f {params.trf_overlap} | """
            """awk '{{print $5}}'""",
            iterable=True
        ))

        trf_set = {int(trf) for trf in trf_set}

        df['TRF'] = df['INDEX'].apply(lambda index: index in trf_set)


### ML ###


class KFoldGridSearchModel:

    def __init__(self, model, param_grid, cv, gt_callable, refit=True):

        # Set argument fields
        self.model = model
        self.param_grid = param_grid
        self.cv = list(cv)  # Unroll so it may be iterated over for each combination of parameters
        self.gt_callable = gt_callable
        self.refit = refit

        # Initialize fit parameters
        self.grid_table = None

    def fit(self, X, y):
        """
        Perform a grid-search cross-validation on features `X` and labels `y`.

        :param X: Features matrix (variants x features)
        :param y: Labels array (1-dim list of labels, one for each row in `X`).
        """

        # Reset CV parameters
        self.reset()

        # Get grid-search parameters
        self.grid_table = pd.DataFrame(
            list(itertools.product(*self.param_grid.values())),
            columns=self.param_grid.keys()
        )

        # Iterate over each set of parameters in the grid
        for row in self.grid_table:
            grid_index = row[0]
            model_params = row[1]

            # Iterate over each fold (as the test set, rest as training set)
            cv_index = 0

            for index_train, index_test in self.cv:
                self.model.set_params(**model_params).fit(X[index_train, :], y[index_train])

                model_predict = self.model.predict(X[index_test, :])

                f1_score(y[index_test], model_predict, average='weighted')
                precision_score(y[index_test], model_predict, average='weighted')
                recall_score(y[index_test], model_predict, average='weighted')
                accuracy_score(y[index_test], model_predict)
                confusion_matrix(y[index_test], model_predict)

                cv_index += 1

        # Return
        return

    def reset(self):
        """
        Reset CV stats.
        """

        self.grid_table = None

        pass
