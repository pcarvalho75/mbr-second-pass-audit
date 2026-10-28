"""First-pass support within quantity terciles, with Wilson intervals (carried from analysis/revision5_analysis.py)."""
import numpy as np
import pandas as pd
from mbraudit.singlecell.intervals import fraction


def wilson(k,n):
    from scipy.stats import norm
    if n==0:return np.nan,np.nan
    z=float(norm.ppf(.975));p=k/n;z2=z*z;den=1+z2/n
    center=(p+z2/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z2/(4*n*n))/den
    return max(0.,center-half),min(1.,center+half)


def original_bin(count,denominator):
    if count==0:return 'zero'
    if count*5<=denominator:return 'positive_0_20'
    if count*5<=2*denominator:return 'positive_20_40'
    if count*5<=3*denominator:return 'positive_40_60'
    return 'positive_60_100'


def support_intervals(reports,c2,output):
    from mbraudit.singlecell.design import H6_RUNS,M6,H16,M16
    from mbraudit.singlecell.accounting import load_quantity,filtered_pair,quantity_terciles,support_bin,species_record
    cohorts=[('homogeneous_6',list(H6_RUNS)),('mixed_6',list(M6)),('homogeneous_16',list(H16)),('mixed_16',list(M16)),('all_cells_64',['all_cells_64'])]
    schemes={'current_three_positive':(support_bin,['zero','low_positive','middle','high']),
        'original_four_positive':(original_bin,['zero','positive_0_20','positive_20_40','positive_40_60','positive_60_100'])}
    records=[];thresholds=[]
    for cohort,aliases in cohorts:
        for alias in aliases:
            data=c2 if alias=='all_cells_64' else load_quantity(reports,alias)
            a,b=filtered_pair(data,'F1');selected=b.loc[~b.direct].copy()
            selected['quantity_tercile'],cuts=quantity_terciles(a,selected)
            for row in cuts:thresholds.append({'cohort':cohort,'pool':alias,**row})
            denominator=len(data['cells'])-1
            for scheme,(classifier,labels) in schemes.items():
                selected['bin']=[classifier(int(s),denominator) for s in selected.support]
                for tercile in ['all','low','middle','high','quantity_unavailable']:
                    part=selected if tercile=='all' else selected[selected.quantity_tercile.eq(tercile)]
                    for label in labels:
                        frame=part[part['bin'].eq(label)]
                        records.append({'cohort':cohort,'pool':alias,'filter':'F1','bin_scheme':scheme,
                            'quantity_tercile':tercile,'support_bin':label,**species_record(frame.species)})
            print('Support cells and original-bin sensitivity',alias,'complete.',flush=True)
    perpool=pd.DataFrame(records)
    keys=['cohort','filter','bin_scheme','quantity_tercile','support_bin']
    fields=['observations','human','yeast','species_exclusive_denominator','shared','shared_IL','unknown']
    pooled=perpool.groupby(keys,sort=False,as_index=False)[fields].sum()
    pooled['observed_yeast_fraction']=fraction(pooled.yeast.to_numpy(),pooled.species_exclusive_denominator.to_numpy())
    bounds=[wilson(int(k),int(n)) for k,n in zip(pooled.yeast,pooled.species_exclusive_denominator)]
    pooled['wilson_low95']=[x[0] for x in bounds];pooled['wilson_high95']=[x[1] for x in bounds]
    for name in ['observed_yeast_fraction','wilson_low95','wilson_high95']:pooled[name+'_pct']=100*pooled[name]
    pooled['n']=pooled.species_exclusive_denominator
    pooled['point_label']=[f'{int(k)}/{int(n)}' for k,n in zip(pooled.yeast,pooled.n)]
    pooled['zero_yeast_label']=[f'no yeast among {int(n)}' if k==0 and n>0 else 'not estimable: n=0' if n==0 else '' for k,n in zip(pooled.yeast,pooled.n)]
    pooled['interval_interpretation']='pointwise Wilson count interval; does not account for recipient/precursor dependence'
    perpool.to_csv(output/'support_per_pool.tsv',sep='\t',index=False,na_rep='NA')
    pooled.to_csv(output/'support_wilson.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(thresholds).to_csv(output/'quantity_cutoffs.tsv',sep='\t',index=False,na_rep='NA')
    contrasts=[]
    for context,part in pooled.groupby(keys[:-1],sort=False):
        indexed=part.set_index('support_bin')
        reference='high' if context[2]=='current_three_positive' else 'positive_60_100'
        ref=indexed.loc[reference]
        for label,row in indexed.iterrows():
            if label==reference:continue
            if ref.n==0 or row.n==0:value=np.nan;text='not estimable: n=0'
            elif ref.yeast==0:value=np.nan;text=f'no yeast among {int(ref.n)} high-support observations'
            else:
                value=row.observed_yeast_fraction/ref.observed_yeast_fraction
                text=f'{value:.3f}'
            contrasts.append(dict(zip(keys[:-1],context),support_bin=label,reference_bin=reference,
                numerator_yeast=int(row.yeast),numerator_n=int(row.n),reference_yeast=int(ref.yeast),reference_n=int(ref.n),
                fraction_ratio=value,ratio_display=text))
    pd.DataFrame(contrasts).to_csv(output/'support_contrasts.tsv',sep='\t',index=False,na_rep='NA')
    primary=pooled[pooled.cohort.eq('all_cells_64')&pooled.quantity_tercile.isin(['low','middle','high'])].copy()
    primary['plot_point']=primary.n.gt(0)&primary.support_bin.ne('zero')
    primary[primary.bin_scheme.eq('current_three_positive')].to_csv(output/'figure3_current.tsv',sep='\t',index=False,na_rep='NA')
    primary[primary.bin_scheme.eq('original_four_positive')].to_csv(output/'figure3_original_bins.tsv',sep='\t',index=False,na_rep='NA')
