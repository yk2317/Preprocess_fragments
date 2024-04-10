import musical
import pandas as pd
import pickle
import numpy as np
import sys

def run_nmf(X, min_n, max_n):
    return musical.DenovoSig(X, min_n_components=min_n, max_n_components=max_n, init='nndsvdar',
                             method='nmf', n_replicates=20, ncpu=10, conv_test_freq=1000,
                             conv_test_baseline='min-iter', max_iter=100000, bootstrap=True,
                             tol=1e-8, verbose=1, normalize_X=True)

def run_mvnmf(X, min_n, max_n):
    return musical.DenovoSig(X, min_n_components=min_n, max_n_components=max_n, init='nndsvdar',
                             method='mvnmf', n_replicates=20, ncpu=10,
                             mvnmf_lambda_tilde_grid=np.array([1e-10, 1e-9, 5e-9, 1e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1]),
                             mvnmf_hyperparameter_method='single', conv_test_freq=1000,
                             conv_test_baseline='min-iter', max_iter=100000, bootstrap=True,
                             tol=1e-8, verbose=1, normalize_X=True)

def run_analysis(input_file, output_file, method, min_n, max_n):
    X = pd.read_csv(input_file, index_col=0)

    if method == 'NMF':
        model = run_nmf(X, min_n, max_n)
    elif method == 'MVNMF':
        model = run_mvnmf(X, min_n, max_n)
    else:
        raise ValueError("Invalid method specified. Choose 'NMF' or 'MVNMF'.")

    model.fit()

    with open(output_file, 'wb') as f:
        pickle.dump(model, f, pickle.HIGHEST_PROTOCOL)

def main():
    if len(sys.argv) < 6:
        print("Usage: python script.py <input_file> <output_file> <method> <min_n> <max_n>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    method = sys.argv[3]
    min_n = int(sys.argv[4])
    max_n = int(sys.argv[5])

    run_analysis(input_file, output_file, method, min_n, max_n)

if __name__ == "__main__":
    main()

