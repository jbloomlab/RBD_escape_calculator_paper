"""``snakemake`` file that runs analysis."""


import pandas as pd


rule all:
    input:
        'results/LucasIwasaki_data.csv'

rule get_LucasIwasaki_data:
    """Get data from https://www.nature.com/articles/s41586-021-04085-y."""
    output: csv='results/LucasIwasaki_data.csv'
    params:
        neut_data='https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-04085-y/MediaObjects/41586_2021_4085_MOESM1_ESM.xlsx',
        # RBD mutations from Extended Data Table 2
        # of https://www.nature.com/articles/s41586-021-04085-y 
        variant_muts = {'Ancestral (A)': '',
                        'B.1.526a': 'S477N',
                        'B.1.526b': 'S477N',
                        'B.1.1.7a': 'N501Y',
                        'B.1.517': 'N501T',
                        'B.1.526c': 'E484K',
                        'B.1.617.2': 'L452R T478K',
                        'R.1': 'E484K',
                        'B.1.427': 'L452R',
                        'B.1.526d': 'L452R',
                        'B.1.429': 'L452R',
                        'B.1.525': 'E484K',
                        'B.1.617.1': 'L452R E484Q',
                        'P.1': 'K417T E484K N501Y',
                        'B.1': 'E484K N501T',
                        'B.1.1.7b': 'E484K N501Y',
                        'B.1.351a': 'K417N E484K N501Y',
                        'B.1.351b': 'K417N E484K N501Y',
                        }
    run:
        df = pd.read_excel(params.neut_data)
        variant_cols = {'IC50 ' + variant: muts for variant, muts
                        in params.variant_muts.items()}
        assert set(df.columns).issuperset(variant_cols)
        df = (df
              .query('`Previous Exposure` in ["(+)", "(-)"]')
              .assign(prior_exposure=lambda x: x['Previous Exposure'].map({'(+)': True,
                                                                           '(-)': False}),
                      has_all_neuts=lambda x: x[list(variant_cols)].apply(lambda s: pd.notnull(s).all(),
                                                                          axis=1),
                      )
              .query('has_all_neuts == True')
              .melt(id_vars=['Sample_ID', 'prior_exposure'],
                    value_vars=variant_cols,
                    var_name='variant',
                    value_name='log10_IC50',
                    )
              .assign(variant=lambda x: x['variant'].str.replace('IC50 ', ''),
                      mutations=lambda x: x['variant'].map(params.variant_muts),
                      )
              )
        df.to_csv(output.csv, index=False)
