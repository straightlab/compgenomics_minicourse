import sys
import pandas as pd

def load_data(filename):
    df = pd.read_csv(filename, sep='\t', header=None, names=['species','N'])
    return df

def species_summary(df):
    out = df.groupby('species').sum().sort_values('N', ascending=False)
    return out


if __name__ == '__main__':
    args = sys.argv[1:]
    filename = sys.argv[1]
    df = load_data(filename)
    out = species_summary(df)
    print(out)