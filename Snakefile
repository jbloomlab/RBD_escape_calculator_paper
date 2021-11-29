"""``snakemake`` file that runs analysis."""


import io

import pandas as pd

import pdfplumber


rule all:
    input:
        'results/neut_studies/LucasIwasaki_data.csv',
        'results/neut_studies/UriuSato_data.csv'

rule get_UriuSato_data:
    """Get data from https://www.nejm.org/doi/full/10.1056/NEJMc2114706."""
    input: pdf='downloads/UriuSato_supp.pdf'
    output: csv='results/neut_studies/UriuSato_data.csv'
    params:
        # RBD mutations from Table S5 of https://www.nejm.org/doi/full/10.1056/NEJMc2114706
        variant_muts = {
                        'Parental': '',
                        'Alpha': 'N501Y',
                        'Beta': 'K417N E484K N501Y',
                        'Gamma': 'K417T E484K N501Y',
                        'Delta': 'L452R T478K',
                        'Epsilon': 'L452R',
                        'Lambda': 'L452Q F490S',
                        'Mu': 'R346K E484K N501Y',
                        }
    run:
        dfs = []
        with pdfplumber.open(input.pdf) as pdf:
            for group, page in [('convalescent', 10), ('vaccinated', 11)]:
                text = pdf.pages[page].extract_text()
                text = text[text.index('Donor'): text.index('Geometric')]
                text = text.replace(' a ', ' ').replace(' b ', ' ').replace('Donor ID', 'sample')
                dfs.append(pd.read_csv(io.StringIO(text), sep=' ').assign(group=group))
        df = (pd.concat(dfs, ignore_index=True)
              .melt(id_vars=['sample', 'group'],
                    value_vars=params.variant_muts,
                    value_name='IC50',
                    var_name='variant')
              .assign(RBD_mutations=lambda x: x['variant'].map(params.variant_muts),
                      IC50=lambda x: x['IC50'].str.replace(',', '').str.replace('<', '').astype(float),
                      variant=lambda x: pd.Categorical(x['variant'], params.variant_muts, ordered=True),
                      )
              .sort_values('variant')
              .assign(fold_change=lambda x: x['IC50'] / x.groupby(['sample', 'group'])['IC50'].transform('first'))
              )
        df.to_csv(output.csv, index=False)

rule get_LucasIwasaki_data:
    """Get data from https://www.nature.com/articles/s41586-021-04085-y."""
    output: csv='results/neut_studies/LucasIwasaki_data.csv'
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
              .assign(group=lambda x: x['Previous Exposure'].map({'(+)': 'vaccinated_then_infected',
                                                                  '(-)': 'vaccinated'}),
                      has_all_neuts=lambda x: x[list(variant_cols)].apply(lambda s: pd.notnull(s).all(),
                                                                          axis=1),
                      )
              .query('has_all_neuts == True')
              .rename(columns={'Sample_ID': 'sample'})
              .melt(id_vars=['sample', 'group'],
                    value_vars=variant_cols,
                    var_name='variant',
                    value_name='log10_IC50',
                    )
              .assign(variant=lambda x: pd.Categorical(x['variant'].str.replace('IC50 ', ''),
                                                       params.variant_muts, ordered=True),
                      IC50=lambda x: x['log10_IC50'].map(lambda y: 10**y),
                      RBD_mutations=lambda x: x['variant'].map(params.variant_muts),
                      )
              .sort_values('variant')
              .assign(fold_change=lambda x: x['IC50'] / x.groupby(['sample', 'group'])['IC50'].transform('first'))
              )
        df.to_csv(output.csv, index=False)
