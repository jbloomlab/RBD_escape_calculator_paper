"""``snakemake`` file that runs analysis."""


import io
import urllib.request
import shutil

import pandas as pd

import pdfplumber


rule all:
    input:
        'results/neut_studies/neut_studies.html',
        'results/neut_studies/neut_studies.png',
        'results/neut_studies/neut_studies_mean.png',
        'results/variants/variants.html',
        'results/variants/variants.png',
        'results/variants/variants.tex',
        'docs/neut_studies.html',
        'docs/variants.html',
        'docs/mini_example_escape_calc.html',
        'docs/escape_calc_chart.html',

rule get_html:
    params:
        'https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/bioRxiv_v1/docs/_includes/mini_example_escape_calc.html',
        'https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/bioRxiv_v1/docs/_includes/escape_calc_chart.html'
    output:
        ['docs/mini_example_escape_calc.html',
         'docs/escape_calc_chart.html']
    run:
        for url_in, f_out in zip(params, output):
            print(f"Getting {url_in=} to {f_out=}")
            urllib.request.urlretrieve(url_in, f_out)

rule cp_html:
    input:
        ['results/neut_studies/neut_studies.html',
         'results/variants/variants.html']
    output:
        ['docs/neut_studies.html',
         'docs/variants.html']
    run:
        for f_in, f_out in zip(input, output):
            shutil.copy(f_in, f_out)

rule plot_variants:
    input: csv='data/variant_RBD_muts.csv'
    output:
        html='results/variants/variants.html',
        png='results/variants/variants.png',
        table='results/variants/variants.tex',
    notebook: 'plot_variants.py.ipynb'

rule plot_neut_studies:
    input:
        expand("results/neut_studies/{study}_data.csv",
               study=['LucasIwasaki', 'UriuSato', 'WangHo'])
    output:
        html='results/neut_studies/neut_studies.html',
        png='results/neut_studies/neut_studies.png',
    params: metric="sum of mutations at site"
    notebook: 'plot_neut_studies.py.ipynb'

rule plot_neut_studies_mean:
    input:
        expand("results/neut_studies/{study}_data.csv",
               study=['LucasIwasaki', 'UriuSato', 'WangHo'])
    output:
        html='results/neut_studies/neut_studies_mean.html',
        png='results/neut_studies/neut_studies_mean.png',
    params: metric="mean of mutations at site"
    notebook: 'plot_neut_studies.py.ipynb'

rule get_WangHo_data:
    """Get data from https://www.nature.com/articles/s41586-021-03398-2."""
    input:
        convalescent='data/WangHo_Fig3b.txt',
        vaccinated='data/WangHo_Fig4b.txt',
    output: csv='results/neut_studies/WangHo_data.csv'
    run:
        dfs = []
        for group, path in [('convalescent', input.convalescent), ('vaccinated', input.vaccinated)]:
            dfs.append(
                pd.read_csv(path, sep=' ')
                .melt(id_vars='RBD_mutations',
                      var_name='sample',
                      value_name='fold_change',
                      )
                .assign(fold_change=lambda x: x['fold_change'].where(x['fold_change'] > 0,
                                                                     1 / x['fold_change'].abs()),
                        group=group,
                        )
                )
        pd.concat(dfs).to_csv(output.csv, index=False)

rule get_UriuSato_data:
    """Get data from https://www.nejm.org/doi/full/10.1056/NEJMc2114706."""
    input: pdf='data/UriuSato_supp.pdf'
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
                    value_name='IC50_str',
                    var_name='variant')
              .assign(RBD_mutations=lambda x: x['variant'].map(params.variant_muts),
                      IC50_str=lambda x: x['IC50_str'].astype(str),
                      IC50=lambda x: x['IC50_str'].str.replace(',', '').str.replace('<', '').astype(float),
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
              .assign(group=lambda x: x['Previous Exposure'].map({'(+)': 'infected and vaccinated',
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
