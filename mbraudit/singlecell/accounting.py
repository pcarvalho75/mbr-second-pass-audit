"""Classes retained, added and lost; represented entities, completeness and species composition.

Carried from analysis/revision3_analysis.py (coverage and loss helpers), analysis/revision2_analysis.py
(workflow_definitions) and analysis/revision4_analysis.py. Bodies are unchanged; only in-body imports point to
the new modules.
"""
from collections import defaultdict
import json

import numpy as np
import pandas as pd


FILTERS=('F3','F4','F0')


LEVELS=('precursor','peptide','protein_group')


BINS=('zero','low_positive','middle','high')


TERCILES=('low','middle','high')


SEED=20261018


DRAWS=2000


def event_keys(frame):return set(zip(frame.Run,frame.key))


def maps_from_union(first,final):
    from mbraudit.singlecell.reports import with_prec
    from mbraudit.singlecell.reports import annotations,group_units
    keys=sorted(set(first.key)|set(final.key))
    annotation=annotations(with_prec(first),with_prec(final),keys)
    units=group_units(keys,annotation,'protein_group')
    groups={key:units[key] if annotation[key]['groups'] else None for key in keys}
    peptides={key:annotation[key]['stripped'] for key in keys}
    audit={'precursor_keys':len(keys),'keys_without_group':sum(not row['groups'] for row in annotation.values()),
           'keys_with_multiple_group_strings':sum(len(row['groups'])>1 for row in annotation.values()),
           'exact_group_strings':len({g for row in annotation.values() for g in row['groups']}),
           'fixed_group_components':len({g for g in groups.values() if g is not None})}
    return peptides,groups,audit


def entity_frame(frame,level):
    column={'precursor':'key','peptide':'peptide','protein_group':'protein_group'}[level]
    result=frame[['Run',column,'Precursor.Quantity']].rename(columns={column:'entity'})
    return result[result.entity.notna()].copy()


def coverage(first,final,cells,filter_name):
    output=[]
    for level in ['precursor','peptide','protein_group']:
        a,b=entity_frame(first,level),entity_frame(final,level)
        ar,br=set(a.entity),set(b.entity)
        common=ar&br
        ap=a[np.isfinite(a['Precursor.Quantity'])&a['Precursor.Quantity'].gt(0)].drop_duplicates(['Run','entity'])
        bp=b[np.isfinite(b['Precursor.Quantity'])&b['Precursor.Quantity'].gt(0)].drop_duplicates(['Run','entity'])
        row={'filter':filter_name,'level':level,'cells':len(cells),'first_represented_entities':len(ar),
             'final_represented_entities':len(br),'common_entities':len(common),
             'first_entity_cell_ids':len(a.drop_duplicates(['Run','entity'])),
             'final_entity_cell_ids':len(b.drop_duplicates(['Run','entity'])),
             'first_positive_quantity_entity_cells':len(ap),'final_positive_quantity_entity_cells':len(bp),
             'first_all_row_denominator':len(ar)*len(cells),'final_all_row_denominator':len(br)*len(cells),
             'common_row_denominator':len(common)*len(cells),
             'first_positive_common_entity_cells':int(ap.entity.isin(common).sum()),
             'final_positive_common_entity_cells':int(bp.entity.isin(common).sum()),
             'first_unknown_entity_events':len(first)-len(a),'final_unknown_entity_events':len(final)-len(b)}
        for label,num,den in [('first_all_row_completeness',len(ap),row['first_all_row_denominator']),
                              ('final_all_row_completeness',len(bp),row['final_all_row_denominator']),
                              ('first_common_row_completeness',row['first_positive_common_entity_cells'],row['common_row_denominator']),
                              ('final_common_row_completeness',row['final_positive_common_entity_cells'],row['common_row_denominator'])]:
            row[label]=num/den if den else np.nan
        output.append(row)
    return output


def representation(lost,final,filter_name):
    peptide_cell=defaultdict(set);peptide_cohort=defaultdict(set)
    group_cell=defaultdict(set);group_cohort=defaultdict(set)
    group_peptides_cell=defaultdict(set);group_peptides_cohort=defaultdict(set)
    for row in final.itertuples():
        peptide_cell[row.Run,row.peptide].add(row.key);peptide_cohort[row.peptide].add(row.key)
        if pd.notna(row.protein_group):
            group_cell[row.Run,row.protein_group].add(row.key);group_cohort[row.protein_group].add(row.key)
            group_peptides_cell[row.Run,row.protein_group].add(row.peptide);group_peptides_cohort[row.protein_group].add(row.peptide)
    records=[]
    for row in lost.itertuples():
        cell_peptide=peptide_cell[row.Run,row.peptide];all_peptide=peptide_cohort[row.peptide]
        known=pd.notna(row.protein_group)
        cell_group=group_cell[row.Run,row.protein_group] if known else set()
        all_group=group_cohort[row.protein_group] if known else set()
        records.append({'filter':filter_name,'cell':row.Run,'modified_sequence':row.key[0],'charge':row.key[1],
            'peptide':row.peptide,'protein_group':row.protein_group,'species':row.species,
            'loss_category':row.loss_category,'first_run_q':np.nan,
            'first_quantity':np.nan,
            'original_quantity_percentile':row.original_quantity_percentile,
            'original_first_pass_recurrence':row.original_first_pass_recurrence,
            'filtered_first_pass_recurrence':row.filtered_first_pass_recurrence,
            'peptide_in_final_cell':bool(cell_peptide),'peptide_in_final_cell_other_precursor':bool(cell_peptide-{row.key}),
            'peptide_in_final_cohort':bool(all_peptide),'peptide_in_final_cohort_other_precursor':bool(all_peptide-{row.key}),
            'protein_group_annotation_known':known,
            'protein_group_in_final_cell':bool(cell_group) if known else pd.NA,
            'protein_group_in_final_cell_other_precursor':bool(cell_group-{row.key}) if known else pd.NA,
            'protein_group_in_final_cell_other_peptide':bool(group_peptides_cell[row.Run,row.protein_group]-{row.peptide}) if known else pd.NA,
            'protein_group_in_final_cohort':bool(all_group) if known else pd.NA,
            'protein_group_in_final_cohort_other_precursor':bool(all_group-{row.key}) if known else pd.NA,
            'protein_group_in_final_cohort_other_peptide':bool(group_peptides_cohort[row.protein_group]-{row.peptide}) if known else pd.NA})
    result=pd.DataFrame(records)
    if len(result):
        # Explicit column access avoids namedtuple renaming of dotted report fields.
        result['first_run_q']=lost['Q.Value'].to_numpy()
        result['first_quantity']=lost['Precursor.Quantity'].to_numpy()
    return result


def loss_entities(first,final,lost,filter_name):
    result=[]
    for level in ['precursor','peptide','protein_group']:
        a,b,l=entity_frame(first,level),entity_frame(final,level),entity_frame(lost,level)
        aset,bset,lset=set(a.entity),set(b.entity),set(l.entity)
        ac,bc,lc=[set(zip(f.Run,f.entity)) for f in [a,b,l]]
        result.append({'filter':filter_name,'level':level,'first_entities':len(aset),'final_entities':len(bset),
            'entities_touched_by_lost_observations':len(lset),'touched_entities_still_in_final_cohort':len(lset&bset),
            'first_entities_wholly_absent_final_cohort':len(aset-bset),
            'cell_entities_touched_by_lost_observations':len(lc),'touched_cell_entities_still_represented':len(lc&bc),
            'first_entity_cells_absent_final':len(ac-bc),'first_entity_cells':len(ac),'final_entity_cells':len(bc),
            'lost_observations_without_entity_mapping':len(lost)-len(l)})
    return result


def workflow_definitions():
    from mbraudit.singlecell.design import H16,M16
    return [('first_pass','F0',['all_cells_64'],True),
            ('all_cells_F0','F0',['all_cells_64'],False),
            ('all_cells_F3','F3',['all_cells_64'],False),
            ('all_cells_F4','F4',['all_cells_64'],False),
            ('homogeneous_16_F0','F0',list(H16),False),
            ('homogeneous_16_F4','F4',list(H16),False),
            ('mixed_16_F0','F0',list(M16),False),
            ('mixed_16_F4','F4',list(M16),False)]


def species_record(labels):
    values=pd.Series(list(labels),dtype='object')
    h=int(values.eq('target').sum());e=int(values.eq('entrapment').sum());n=h+e
    return {'observations':len(values),'human':h,'yeast':e,'species_exclusive_denominator':n,
            'shared':int(values.eq('shared').sum()),'shared_IL':int(values.eq('shared_IL').sum()),
            'unknown':int(values.eq('unknown').sum()),'observed_yeast_fraction':e/n if n else np.nan}


def fraction_ratio(numerator,denominator):
    if not np.isfinite(numerator) or not np.isfinite(denominator):return np.nan
    if denominator==0:return np.inf if numerator>0 else np.nan
    return numerator/denominator


def collapse_labels(labels):
    values=set(labels)
    if len(values)==1:return next(iter(values))
    if 'unknown' in values:return 'unknown'
    if 'shared_IL' in values:return 'shared_IL'
    return 'shared'


def fixed_species(first,final,registry):
    both=pd.concat([first,final],ignore_index=True)
    maps={}
    for level,col in [('precursor','key'),('peptide','peptide')]:
        maps[level]=both.groupby(col,sort=False).species.agg(collapse_labels).to_dict()
    from mbraudit.species import classify
    groups={};audit=[]
    for group,rows in both.groupby('protein_group',sort=False):
        assert group.startswith('group:')
        annotation=classify(group[6:],registry.target,registry.entrap)
        labels=set(rows.species)
        if annotation=='entrapment' and 'entrapment' not in labels:
            label='shared_IL' if 'shared_IL' in labels else 'unknown'
        elif annotation=='target' and 'entrapment' in labels:label='shared'
        else:label=annotation
        groups[group]=label
        audit.append({'protein_group':group,'fixed_species':label,'annotation_species':annotation,
            'any_exclusive_yeast_constituent':bool('entrapment' in labels),
            'any_ambiguous_constituent':bool(labels-{'target','entrapment'}),
            'all_constituents_exclusive_human':labels=={'target'}})
    maps['protein_group']=groups
    return maps,pd.DataFrame(audit)


def load_quantity(reports,alias):
    from mbraudit.species import accept_rows
    from mbraudit.singlecell.reports import KEY
    from mbraudit.singlecell.reports import report_dir
    loaded=reports.load(alias);out={}
    for stage,file in [('first','report-first-pass.parquet'),('final','report.parquet')]:
        frame=loaded[stage].copy()
        raw=accept_rows(pd.read_parquet(report_dir(alias)/file,
            columns=['Run',*KEY,'Decoy','Q.Value','Precursor.Quantity']),targets_only=True)
        assert not raw.duplicated(['Run','key']).any()
        frame=frame.merge(raw[['Run','key','Precursor.Quantity']],on=['Run','key'],how='left',validate='one_to_one',sort=False)
        assert len(frame)==len(raw)
        out[stage]=frame
    out['cells']=sorted(set(out['first'].Run)|set(out['final'].Run))
    return out


def filtered_pair(loaded,filt):
    from mbraudit.filters import filter_mask
    from mbraudit.singlecell.reports import QCOLS
    a=loaded['first'].loc[filter_mask(loaded['first'][QCOLS].fillna(2.),filt)].copy()
    b=loaded['final'].loc[filter_mask(loaded['final'][QCOLS].fillna(2.),filt)].copy()
    own=set(zip(a.Run,a.key));counts=a.groupby('key').Run.nunique().to_dict()
    b['direct']=[(cell,key) in own for cell,key in zip(b.Run,b.key)]
    b['support']=[counts.get(key,0)-int(direct) for key,direct in zip(b.key,b.direct)]
    assert b.support.between(0,len(loaded['cells'])-1).all()
    return a,b


def group_representation(frame,final):
    known=frame.protein_group.notna()
    pairs=set(zip(final.Run,final.protein_group))
    triples=set(zip(final.Run,final.protein_group,final.peptide))
    count=final.groupby(['Run','protein_group']).peptide.nunique().to_dict()
    any_group=np.array([(c,g) in pairs for c,g in zip(frame.Run,frame.protein_group)])
    other=np.array([count.get((c,g),0)-int((c,g,p) in triples)>0 for c,g,p in zip(frame.Run,frame.protein_group,frame.peptide)])
    n=int(known.sum())
    return {'group_annotation_evaluable_observations':n,'same_cell_group_in_final':int((known&any_group).sum()),
            'same_cell_group_in_final_fraction':float(any_group[known].mean()) if n else np.nan,
            'same_cell_group_other_peptide_in_final':int((known&other).sum()),
            'same_cell_group_other_peptide_in_final_fraction':float(other[known].mean()) if n else np.nan}


def positive_pairs(frame,col):
    rows=frame.loc[np.isfinite(frame['Precursor.Quantity'])&frame['Precursor.Quantity'].gt(0)&frame[col].notna()]
    return set(zip(rows.Run,rows[col]))


def entity_composition(first,final,maps,filt):
    rows=[];ledger=[];human=[]
    for level,col in [('precursor','key'),('peptide','peptide'),('protein_group','protein_group')]:
        aa=set(first[col].dropna());bb=set(final[col].dropna())
        absent=aa-bb;present=aa&bb
        records={}
        for state,keys in [('wholly_absent',absent),('still_represented',present)]:
            record=species_record(maps[level][key] for key in sorted(keys))
            rows.append({'filter':filt,'level':level,'state':state,**record});records[state]=record
        ratio=fraction_ratio(records['wholly_absent']['observed_yeast_fraction'],records['still_represented']['observed_yeast_fraction'])
        for row in rows[-2:]:row['absent_vs_still_yeast_fraction_ratio']=ratio
        for key in sorted(aa):
            ledger.append({'filter':filt,'level':level,'entity':str(key),'species':maps[level][key],
                'wholly_absent_final':key in absent})
        if level=='protein_group':
            ac=set(zip(first.Run,first[col]));bc=set(zip(final.Run,final[col]));lost=ac-bc
            ap=positive_pairs(first,col);bp=positive_pairs(final,col)
            human.append({'filter':filt,'human_groups_wholly_absent':sum(maps[level][g]=='target' for g in absent),
                'human_group_cell_identifications_absent':sum(maps[level][g]=='target' for c,g in lost),
                'human_groups_with_cell_identifications_absent':len({g for c,g in lost if maps[level][g]=='target'}),
                'human_positive_quantity_group_cells_absent':sum(maps[level][g]=='target' for c,g in ap-bp),
                'all_group_cell_identifications_absent':len(lost),'all_positive_quantity_group_cells_absent':len(ap-bp),
                'first_human_groups':sum(maps[level][g]=='target' for g in aa)})
    return rows,ledger,human


def global_tail(first,final):
    allkeys=set(first.key);finalkeys=set(final.key)
    good=set(first.loc[first['Global.Q.Value'].le(.01),'key'])
    failed=allkeys-good;absent=allkeys-finalkeys
    return {'first_F0_keys':len(allkeys),'final_F0_keys':len(finalkeys),'wholly_absent_F0_keys':len(absent),
        'first_keys_without_any_Global_Q_at_or_below_0_01':len(failed),
        'globally_failing_keys_wholly_absent_final':len(failed&absent),
        'globally_failing_keys_present_final':len(failed&finalkeys),
        'wholly_absent_keys_with_a_global_pass':len(absent&good),
        'share_absent_keys_without_global_pass':len(failed&absent)/len(absent) if absent else np.nan,
        'net_represented_key_decrease':len(allkeys)-len(finalkeys),
        'share_net_decrease_without_global_pass':len(failed&absent)/(len(allkeys)-len(finalkeys)) if len(allkeys)!=len(finalkeys) else np.nan}


def all_cell_accounting(reports,loaded,output):
    from mbraudit.singlecell.accounting import maps_from_union,event_keys,coverage,loss_entities,representation
    from mbraudit.singlecell.reports import KEY
    first,final=loaded['first'],loaded['final']
    peptides,groups,mapping=maps_from_union(first,final)
    assert mapping['keys_without_group']==mapping['keys_with_multiple_group_strings']==0
    for frame in [first,final]:
        frame['peptide']=frame.key.map(peptides);frame['protein_group']=frame.key.map(groups)
    maps,group_audit=fixed_species(first,final,reports.registry)
    group_audit.to_csv(output/'group_species_audit.tsv',sep='\t',index=False)
    pd.DataFrame([mapping]).to_csv(output/'mapping_audit.tsv',sep='\t',index=False)
    pd.DataFrame([global_tail(first,final)]).to_csv(output/'global_q_tail.tsv',sep='\t',index=False)
    library=reports.library();library=library[library.Decoy.eq(0)]
    librarykeys=set(zip(library[KEY[0]],library[KEY[1]].astype(int)))
    positive=first['Precursor.Quantity'].where(np.isfinite(first['Precursor.Quantity'])&first['Precursor.Quantity'].gt(0))
    first['original_quantity_percentile']=positive.groupby(first.Run).rank(method='average',pct=True)
    first['original_first_pass_recurrence']=first.key.map(first.groupby('key').Run.nunique())
    accounts=[];cover=[];classes=[];compositions=[];entity_ledger=[];human=[];losses=[];lost_ledgers=[]
    fp0=event_keys(first);fin0=event_keys(final)
    for filt in FILTERS:
        a,b=filtered_pair(loaded,filt)
        ak,bk=event_keys(a),event_keys(b)
        retained=a.loc[[(c,k) in bk for c,k in zip(a.Run,a.key)]].copy()
        lost=a.loc[[(c,k) not in bk for c,k in zip(a.Run,a.key)]].copy()
        only=b.loc[~b.direct].copy()
        lost['loss_category']=np.where(lost.key.isin(librarykeys),'library_present_absent_final','library_excluded')
        accounts.append({'filter':filt,'primary':filt=='F3','first_observations':len(a),'final_observations':len(b),
            'retained_observations':len(retained),'second_pass_only_observations':len(only),'lost_observations':len(lost),
            'lost_library_present':int(lost.key.isin(librarykeys).sum()),'lost_library_excluded':int((~lost.key.isin(librarykeys)).sum()),
            'second_pass_only_with_F0_firstpass':sum((c,k) in fp0 for c,k in zip(only.Run,only.key)),
            'lost_with_F0_final':sum((c,k) in fin0 for c,k in zip(lost.Run,lost.key)),
            'net_observation_change_fraction':len(b)/len(a)-1})
        cover.extend(coverage(a,b,loaded['cells'],filt));losses.extend(loss_entities(a,b,lost,filt))
        a['filtered_first_pass_recurrence']=a.key.map(a.groupby('key').Run.nunique())
        lost['filtered_first_pass_recurrence']=lost.key.map(a.groupby('key').Run.nunique())
        lost_ledgers.append(representation(lost,b,filt))
        sets={'first_pass':a,'final':b,'retained':retained,'second_pass_only':only,'lost':lost,
            'lost_library_present':lost[lost.key.isin(librarykeys)],'lost_library_excluded':lost[~lost.key.isin(librarykeys)]}
        for name,frame in sets.items():classes.append({'filter':filt,'class':name,**species_record(frame.species),**group_representation(frame,b)})
        co,le,hu=entity_composition(a,b,maps,filt)
        compositions.extend(co);entity_ledger.extend(le);human.extend(hu)
    tables={'accounting.tsv':accounts,'represented_coverage.tsv':cover,'class_species.tsv':classes,
        'entity_species.tsv':compositions,'entity_ledger.tsv':entity_ledger,'human_group_losses.tsv':human,'loss_entities.tsv':losses}
    for name,rows in tables.items():pd.DataFrame(rows).to_csv(output/name,sep='\t',index=False,na_rep='NA')
    pd.concat(lost_ledgers,ignore_index=True).to_csv(output/'lost_observation_representation.tsv',sep='\t',index=False,na_rep='NA')
    print('All-cell accounting and entity composition complete.',flush=True)


def support_bin(count,denominator):
    if count==0:return 'zero'
    if count*5<=denominator:return 'low_positive'
    if count*2<=denominator:return 'middle'
    return 'high'


def quantity_terciles(first,selected):
    finite=first.loc[np.isfinite(first['Precursor.Quantity'])&first['Precursor.Quantity'].gt(0)]
    cut=finite.groupby('Run')['Precursor.Quantity'].quantile([1/3,2/3]).unstack()
    n=finite.groupby('Run').size()
    records=[{'cell':cell,'first_positive_quantity_observations':int(n[cell]),
        'lower_tercile_quantity':row.iloc[0],'upper_tercile_quantity':row.iloc[1]} for cell,row in cut.iterrows()]
    # Second-pass-only observations have no same-filter first-pass quantity.
    values=selected['Precursor.Quantity'].to_numpy()
    lower=selected.Run.map(cut.iloc[:,0]).to_numpy();upper=selected.Run.map(cut.iloc[:,1]).to_numpy()
    valid=np.isfinite(values)&(values>0)&np.isfinite(lower)&np.isfinite(upper)
    terciles=np.full(len(selected),'quantity_unavailable',dtype=object)
    terciles[valid]=np.where(values[valid]<=lower[valid],'low',np.where(values[valid]<=upper[valid],'middle','high'))
    return terciles,records


def support_strata(reports,loaded_cache,output):
    from mbraudit.singlecell.design import H6_RUNS,M6,H16,M16
    cohorts=[('homogeneous_6',list(H6_RUNS)),('mixed_6',list(M6)),
        ('homogeneous_16',list(H16)),('mixed_16',list(M16)),('all_cells_64',['all_cells_64'])]
    parts=[];thresholds=[];totals=[]
    for cohort,aliases in cohorts:
        for alias in aliases:
            loaded=loaded_cache.get(alias)
            if loaded is None:loaded=load_quantity(reports,alias);loaded_cache[alias]=loaded
            for filt in FILTERS:
                a,b=filtered_pair(loaded,filt);selected=b.loc[~b.direct].copy()
                denom=len(loaded['cells'])-1
                selected['support_bin']=[support_bin(int(k),denom) for k in selected.support]
                selected['quantity_tercile'],cuts=quantity_terciles(a,selected)
                for row in cuts:thresholds.append({'cohort':cohort,'pool':alias,'filter':filt,**row})
                original_own=set(zip(loaded['first'].Run,loaded['first'].key))
                totals.append({'cohort':cohort,'pool':alias,'filter':filt,'second_pass_only_observations':len(selected),
                    'quantity_unavailable':int(selected.quantity_tercile.eq('quantity_unavailable').sum()),
                    'same_filter_support_denominator':denom,
                    'second_pass_only_with_F0_firstpass':sum((c,k) in original_own for c,k in zip(selected.Run,selected.key))})
                for tercile in ['all',*TERCILES,'quantity_unavailable']:
                    subset=selected if tercile=='all' else selected[selected.quantity_tercile.eq(tercile)]
                    for label in BINS:
                        frame=subset[subset.support_bin.eq(label)]
                        parts.append({'cohort':cohort,'pool':alias,'filter':filt,'quantity_tercile':tercile,
                            'support_bin':label,**species_record(frame.species)})
            print('Support strata',cohort,alias,flush=True)
    keys=['cohort','filter','quantity_tercile','support_bin']
    counts=['observations','human','yeast','species_exclusive_denominator','shared','shared_IL','unknown']
    perpool=pd.DataFrame(parts)
    pooled=perpool.groupby(keys,sort=False,as_index=False)[counts].sum()
    pooled['observed_yeast_fraction']=pooled.yeast.div(pooled.species_exclusive_denominator.where(pooled.species_exclusive_denominator.gt(0)))
    pooled['quantity_source']='final quantity compared with same-filter recipient first-pass quantity cutoffs'
    perpool.to_csv(output/'support_quantity_per_pool.tsv',sep='\t',index=False,na_rep='NA')
    pooled.to_csv(output/'support_quantity_bins.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(thresholds).to_csv(output/'quantity_tercile_thresholds.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(totals).to_csv(output/'support_pool_totals.tsv',sep='\t',index=False,na_rep='NA')
    ratios=[];gradients=[]
    for context,rows in pooled.groupby(['cohort','filter','quantity_tercile'],sort=False):
        indexed=rows.set_index('support_bin')
        for label in ['zero','low_positive','middle']:
            numerator=indexed.loc[label];denominator=indexed.loc['high']
            ratios.append(dict(zip(['cohort','filter','quantity_tercile'],context),comparison=label+'_vs_high',
                numerator_yeast=int(numerator.yeast),numerator_denominator=int(numerator.species_exclusive_denominator),
                denominator_yeast=int(denominator.yeast),denominator_denominator=int(denominator.species_exclusive_denominator),
                observed_yeast_fraction_ratio=fraction_ratio(numerator.observed_yeast_fraction,denominator.observed_yeast_fraction)))
        values=indexed.loc[['low_positive','middle','high'],'observed_yeast_fraction'].to_numpy()
        defined=bool(np.isfinite(values).all())
        gradients.append(dict(zip(['cohort','filter','quantity_tercile'],context),
            all_positive_bins_evaluable=defined,positive_bins_nonincreasing=bool(np.all(np.diff(values)<=0)) if defined else pd.NA,
            low_positive_above_high=bool(values[0]>values[2]) if np.isfinite(values[[0,2]]).all() else pd.NA,
            low_positive_fraction=values[0],middle_fraction=values[1],high_fraction=values[2],
            interpretation='descriptive point gradient; no causal or abundance-independent filter validation'))
    pd.DataFrame(ratios).to_csv(output/'support_quantity_ratios.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(gradients).to_csv(output/'support_quantity_gradient.tsv',sep='\t',index=False,na_rep='NA')


def workflows(reports,loaded_cache,output):
    from mbraudit.singlecell.bootstrap import bootstrap,interval,FACTOR
    from mbraudit.singlecell.accounting import workflow_definitions
    from mbraudit.singlecell.reports import natural
    ref={f:filtered_pair(loaded_cache['all_cells_64'],f)[1] for f in FILTERS}
    rows=[];long=[];perpool=[]
    cells=sorted(loaded_cache['all_cells_64']['cells'],key=natural)
    for name,filt,aliases,isfirst in workflow_definitions():
        frames=[];alternatives=[]
        for alias in aliases:
            loaded=loaded_cache.get(alias)
            if loaded is None:loaded=load_quantity(reports,alias);loaded_cache[alias]=loaded
            a,b=filtered_pair(loaded,filt)
            frame=a if isfirst else b
            if isfirst:frame=frame.copy();frame['direct']=True
            frames.append(frame)
            if not isfirst:
                original_own=set(zip(loaded['first'].Run,loaded['first'].key))
                alternative=frame.loc[[(c,k) not in original_own for c,k in zip(frame.Run,frame.key)]]
                alternatives.append(alternative)
            perpool.append({'workflow':name,'filter':filt,'pool':alias,'accepted_observations':len(frame),
                **species_record(frame.loc[~frame.direct,'species'])})
        frame=pd.concat(frames,ignore_index=True);selected=frame.loc[~frame.direct]
        assert not frame.duplicated(['Run','key']).any()
        previous=pd.concat(alternatives,ignore_index=True) if alternatives else selected
        record=species_record(selected.species)
        row={'workflow':name,'filter':filt,'accepted_observations':len(frame),'second_pass_only_observations':len(selected),
            'second_pass_only_relative_to_F0_firstpass':len(previous),'reclassified_from_F0_firstpass':len(selected)-len(previous),
            'represented_precursors':frame.key.nunique(),'represented_protein_groups':frame['Protein.Group'].dropna().replace('',np.nan).nunique(),
            'accepted_observation_retention_same_filter':len(frame)/len(ref[filt]),
            'represented_precursor_retention_same_filter':frame.key.nunique()/ref[filt].key.nunique(),
            **record}
        for scheme in ['recipient','precursor','two_way']:
            values=bootstrap([selected],cells,'fdp',scheme,draws=DRAWS,seed=SEED,batch=100)
            fractions=values/FACTOR;stats=interval(fractions)
            if selected.empty:stats['interval_status']='not applicable: no second-pass-only observations'
            label=name+'_'+scheme
            pd.DataFrame({'draw':np.arange(1,DRAWS+1),'observed_yeast_fraction':fractions}).to_csv(output/('draws_workflow_'+label+'.tsv'),sep='\t',index=False,na_rep='NA')
            long.append({'workflow':name,'filter':filt,'scheme':scheme,'point':record['observed_yeast_fraction'],
                **stats,'seed':SEED,'recipients':len(cells),'precursor_blocks':selected.key.nunique()})
            for key,value in stats.items():row[scheme+'_'+key]=value
        rows.append(row)
        print('Workflow',name,'second-pass-only',len(selected),flush=True)
    pd.DataFrame(rows).to_csv(output/'workflow_tradeoffs.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(long).to_csv(output/'workflow_intervals.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(perpool).to_csv(output/'workflow_per_pool.tsv',sep='\t',index=False,na_rep='NA')


def class_ratios_and_scales(output,candidate_path):
    candidate=pd.read_csv(candidate_path,sep='\t')
    empirical=candidate.loc[candidate.library.eq('empirical')&candidate.unit.eq('precursor')&
        candidate.charge.astype(str).eq('all')&candidate.classification.eq('I_as_L_exclusive')].iloc[0]
    empirical_ratio=float(empirical.r)
    ratios=[];scales=[]
    classes=pd.read_csv(output/'class_species.tsv',sep='\t')
    for filt,rows in classes.groupby('filter',sort=False):
        ix=rows.set_index('class');den=ix.loc['retained']
        for name in ['second_pass_only','lost','lost_library_present','lost_library_excluded']:
            row=ix.loc[name]
            ratios.append({'filter':filt,'numerator_class':name,'denominator_class':'retained',
                'numerator_yeast':row.yeast,'numerator_denominator':row.species_exclusive_denominator,
                'denominator_yeast':den.yeast,'denominator_denominator':den.species_exclusive_denominator,
                'observed_yeast_fraction_ratio':fraction_ratio(row.observed_yeast_fraction,den.observed_yeast_fraction)})
    pd.DataFrame(ratios).to_csv(output/'class_yeast_ratios.tsv',sep='\t',index=False,na_rep='NA')
    # Candidate-ratio scales apply to precursor observations, not distinct protein-group units.
    tables=[('class_species.tsv','class'),('workflow_tradeoffs.tsv','workflow'),('support_quantity_bins.tsv','support')]
    for name,domain in tables:
        frame=pd.read_csv(output/name,sep='\t')
        for source in frame.to_dict('records'):
            record={key:source[key] for key in ['filter','class','level','state','workflow','cohort','quantity_tercile','support_bin'] if key in source}
            record.update({'domain':domain,'human':source['human'],'yeast':source['yeast'],
                'species_exclusive_denominator':source['species_exclusive_denominator'],
                'observed_yeast_fraction':source['observed_yeast_fraction'],
                'historical_r':.276123,'empirical_library_r':empirical_ratio,
                'conditional_combined_historical':source['observed_yeast_fraction']*(1+1/.276123),
                'conditional_combined_empirical':source['observed_yeast_fraction']*(1+1/empirical_ratio),
                'interpretation':'two candidate-opportunity assumptions; neither is calibrated precursor FDP'})
            scales.append(record)
    pd.DataFrame(scales).to_csv(output/'conditional_scales.tsv',sep='\t',index=False,na_rep='NA')


def figures(output):
    account=pd.read_csv(output/'accounting.tsv',sep='\t')
    classes=pd.read_csv(output/'class_species.tsv',sep='\t')
    flows=classes.loc[classes['class'].isin(['first_pass','retained','second_pass_only','lost_library_present','lost_library_excluded'])].copy()
    flows['observed_yeast_pct']=100*flows.observed_yeast_fraction
    flows['same_cell_group_retained_pct']=100*flows.same_cell_group_in_final_fraction
    flows['same_cell_group_other_peptide_pct']=100*flows.same_cell_group_other_peptide_in_final_fraction
    flows.to_csv(output/'figure_flow.tsv',sep='\t',index=False,na_rep='NA')
    cover=pd.read_csv(output/'represented_coverage.tsv',sep='\t')
    for col in ['first_all_row_completeness','final_all_row_completeness','first_common_row_completeness','final_common_row_completeness']:cover[col+'_pct']=100*cover[col]
    cover.to_csv(output/'figure_completeness.tsv',sep='\t',index=False,na_rep='NA')
    entity=pd.read_csv(output/'entity_species.tsv',sep='\t');entity['observed_yeast_pct']=100*entity.observed_yeast_fraction
    entity.to_csv(output/'figure_entity_species.tsv',sep='\t',index=False,na_rep='NA')
    support=pd.read_csv(output/'support_quantity_bins.tsv',sep='\t');support['observed_yeast_pct']=100*support.observed_yeast_fraction
    support.to_csv(output/'figure_support_quantity.tsv',sep='\t',index=False,na_rep='NA')
    overview=account.merge(cover[cover.level.eq('precursor')],on='filter',validate='one_to_one')
    overview.to_csv(output/'central_accounting.tsv',sep='\t',index=False,na_rep='NA')
