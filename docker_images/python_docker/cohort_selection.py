import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_ann')
    parser.add_argument('-o', '--output', dest='output_ann')
    parser.add_argument('selections',
                        nargs='+',
                        help='Enter as "a=b".')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    d = {}
    for s in args.selections:
        k, v = [x.strip() for x in s.split('=')]
        d[k] = v

    df = pd.read_table(args.input_ann, index_col=0)
    for k, v in d.items():
        df = df.loc[df[k] == v]

    df.to_csv(args.output_ann, sep='\t')
