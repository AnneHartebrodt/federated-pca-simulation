def prepare_data_iris():
    original = pd.read_csv(
        filepath_or_buffer='https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data',
        header=None,
        sep=',')
    r.seed(a=12)
    original = original.sample(frac=1, random_state=12).reset_index(
        drop=True)  # randomize dataframe to get rid of the order
    original.columns = ['sepal_len', 'sepal_wid', 'petal_len', 'petal_wid', 'class']
    original.dropna(how="all", inplace=True)  # drops the empty line at file-end
    labels = original.iloc[:, 4].values
    orig = sc.array(original)[:, 0:4]
    orig = orig.astype(float)
    return orig, label