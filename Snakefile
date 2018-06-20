# Snakefile

# Perform a phyloHiC analysis

# Author: Sylvain PULICANI <pulicani@lirmm.fr>
# Created on: June 19, 2018
# Last modified on: June 19, 2018

from itertools import combinations
from os.path import dirname
from string import Template

shell.prefix("source /beegfs/home/pulicani/miniconda3beegfs/bin/activate py36;")

THRESHOLDS = ['00', '10', '25', '50', '75', '90']
ADJACENCIES = ['all', 'none', 'and', 'xor']

bindir = config['bin']
orthologs = config['orthologs']
name = config['experiment']
method = config['method']
resolutions = config['resolutions']
n = config['nb_replicates']


rule all:
    input:
        expand('{rootdir}/{res}/{name}/threshold_{th}_adj_{adj}/trees_{method}_{rep}/all_replicates.phylip',
               rootdir=config['rootdir'],
               res=resolutions,
               name=name,
               th=THRESHOLDS,
               adj=ADJACENCIES,
               method=method,
               rep=n)


rule cat_replicates:
    input:
        expand('{{rootdir}}/{{res}}/{name}/threshold_{{th}}_adj_{{adj}}/trees_{{type}}_{{rep}}/replicates_{n}.phylip',
               n=range(n), name=name)
    output:
        '{rootdir}/{res}/{name}/threshold_{th}_adj_{adj}/trees_{type}_{rep}/all_replicates.phylip'
    shell:
        'cat {input} > {output}'


rule make_bootstrap:
    input:
        [Template('{rootdir}/{res}/{name}/threshold_{th}_adj_{adj}/${d1}_${d2}_values.tsv.gz').substitute(d1=d1, d2=d2)
         for d1, d2 in combinations(config['datasets'], 2)]
    output:
        expand('{{rootdir}}/{{res}}/{{name}}/threshold_{{th}}_adj_{{adj}}/trees_{{method}}_{{rep}}/replicates_{n}.phylip',
               n=range(n))
    params:
        union='-u' if config['method'] == 'union' else '',
        outdir=lambda wildcards, output: dirname(output[0])
    shell:
        '{bindir}/bootstrap.py {params.union} {orthologs} {params.outdir} {n} {input}'


rule get_stats:
    input:
        expand('{{rootdir}}/{{res}}/{{name}}/pairs/{dataset}.tsv.gz',
               dataset=config['datasets'])
    output:
        '{rootdir}/{res}/{name}/pairs/stats.tsv'
    shell:
        "{bindir}/pairsStats {input} > {output}"



def get_threshold(threshold, stat):
    with open(stat, 'r') as f:
        lines = [l for l in f.read().split('\n') if l]
    d = dict(zip(lines[0].replace('%', '').split('\t'),
                 lines[-1].split('\t')))
    if threshold == '50':
        return '-t ' + d['Median']
    elif threshold == '00':
        return '-t ' + '0.0'
    else:
        return '-t ' + d.get(threshold, '0.0')

rule join_pairs:
    input:
        s = rules.get_stats.output,
        d1 = '{rootdir}/{res}/{name}/pairs/{dataset1}.tsv.gz',
        d2 = '{rootdir}/{res}/{name}/pairs/{dataset2}.tsv.gz'
    output:
        '{rootdir}/{res}/{name}/threshold_{th}_adj_{adj}/{dataset1}_{dataset2}_values.tsv.gz'
    # params:
    #     threshold=get_threshold
    run:
        s = '{}/{}/{}/pairs/stats.tsv'.format(
            wildcards.rootdir, wildcards.res, wildcards.name)
        # print(s)
        # t = get_threshold(wildcards.th, input.s)
        t = get_threshold(wildcards.th, s)
        shell("{bindir}/join_pairs.py {t} -a {wildcards.adj} {orthologs} {input.d1} {input.d2} {output}")




def make_pairs_input(wildcards):
    return (config['datasets'][wildcards.dataset]['genes'],
            config['datasets'][wildcards.dataset]['hic'][wildcards.res])

rule make_pairs:
    input:
        make_pairs_input
    output:
        '{rootdir}/{res}/{name}/pairs/{dataset}.tsv.gz'
    params:
        scramble='-s' if config['scramble'] else ''
    shell:
        "{bindir}/make_pairs.py -N {params.scramble} {input} {output}"
