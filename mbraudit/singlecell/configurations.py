"""Side-by-side tables of two search configurations (carried from analysis/revision6_sensitivity.py)."""
import numpy as np
import pandas as pd


def read(path):return pd.read_csv(path,sep='\t',float_precision='round_trip')


def paired_tables(original,norefine,keys,columns):
    assert not original.duplicated(keys).any() and not norefine.duplicated(keys).any()
    assert set(map(tuple,original[keys].to_numpy()))==set(map(tuple,norefine[keys].to_numpy()))
    a=original[keys+columns].rename(columns={c:'original_'+c for c in columns})
    b=norefine[keys+columns].rename(columns={c:'norefine_'+c for c in columns})
    return a.merge(b,on=keys,validate='one_to_one')


def empirical_opportunities(reports,original_snapshot,output):
    from mbraudit.singlecell.reports import report_dir
    from mbraudit.singlecell.candidates import label
    previous=read(original_snapshot)
    old=previous[previous.library.eq('empirical')&previous.unit.eq('precursor')&previous.charge.astype(str).eq('all')&previous.classification.eq('I_as_L_exclusive')].iloc[0]
    path=report_dir('all_cells_64')/'report-lib.parquet'
    raw=pd.read_parquet(path,columns=['Modified.Sequence','Precursor.Charge','Stripped.Sequence','Protein.Ids','Decoy'])
    raw=raw[raw.Decoy.eq(0)].drop(columns='Decoy').drop_duplicates()
    assert not raw.duplicated(['Modified.Sequence','Precursor.Charge']).any()
    labeled=label(raw,reports.registry,'returned empirical library')
    counts=labeled.species.value_counts()
    h,e=int(counts.get('target',0)),int(counts.get('entrapment',0))
    oldrow={'configuration':'original_default','total_candidates':int(old.total_candidates),'human_only':int(old.human_only),
        'yeast_only':int(old.yeast_only),'shared':int(old.shared_accessions),'shared_IL':int(old.shared_IL),
        'unknown':int(old.unknown),'r':float(old.r)}
    newrow={'configuration':'disabled_refinement','total_candidates':len(labeled),'human_only':h,'yeast_only':e,
        'shared':int(counts.get('shared',0)),'shared_IL':int(counts.get('shared_IL',0)),'unknown':int(counts.get('unknown',0)),
        'r':e/h if h else np.nan}
    assert sum(newrow[x] for x in ['human_only','yeast_only','shared','shared_IL','unknown'])==len(labeled)
    pd.DataFrame([oldrow,newrow]).to_csv(output/'empirical_library_opportunities.tsv',sep='\t',index=False)
    labeled.to_csv(output/'returned_empirical_candidate_species.tsv',sep='\t',index=False)
    return oldrow['r'],newrow['r']


def comparison_tables(old4,old5,new4,new5,output,oldr,newr):
    from mbraudit.species import combined_factor
    from mbraudit.singlecell.accounting import fraction_ratio
    a=read(old4/'accounting.tsv');b=read(new4/'accounting.tsv')
    cols=['first_observations','final_observations','retained_observations','second_pass_only_observations',
        'lost_observations','lost_library_present','lost_library_excluded','net_observation_change_fraction',
        'second_pass_only_with_F3_firstpass','lost_with_F3_final']
    accounting=paired_tables(a,b,['filter'],cols)
    for col in cols:accounting['change_'+col]=accounting['norefine_'+col]-accounting['original_'+col]
    accounting.to_csv(output/'accounting_comparison.tsv',sep='\t',index=False)
    covcols=['first_represented_entities','final_represented_entities','common_entities',
        'first_all_row_completeness','final_all_row_completeness','first_common_row_completeness','final_common_row_completeness',
        'first_positive_common_entity_cells','final_positive_common_entity_cells','common_row_denominator']
    coverage=paired_tables(read(old4/'represented_coverage.tsv'),read(new4/'represented_coverage.tsv'),['filter','level'],covcols)
    coverage.to_csv(output/'coverage_comparison.tsv',sep='\t',index=False)
    entitycols=['observations','human','yeast','species_exclusive_denominator','shared','shared_IL','unknown',
        'observed_yeast_fraction','absent_vs_still_yeast_fraction_ratio']
    paired_tables(read(old4/'entity_species.tsv'),read(new4/'entity_species.tsv'),['filter','level','state'],entitycols).to_csv(output/'entity_species_comparison.tsv',sep='\t',index=False,na_rep='NA')
    hu=['human_groups_wholly_absent','human_group_cell_identifications_absent','human_positive_quantity_group_cells_absent',
        'human_groups_with_cell_identifications_absent']
    paired_tables(read(old4/'human_group_losses.tsv'),read(new4/'human_group_losses.tsv'),['filter'],hu).to_csv(output/'human_group_loss_comparison.tsv',sep='\t',index=False)
    classes=[]
    for directory,r in [(old4,oldr),(new4,newr)]:
        frame=read(directory/'class_species.tsv')
        retained=frame[frame['class'].eq('retained')].set_index('filter').observed_yeast_fraction.to_dict()
        frame['r_free_ratio_vs_retained']=[fraction_ratio(value,retained[filt]) for value,filt in zip(frame.observed_yeast_fraction,frame['filter'])]
        frame['historical_r']=.276123;frame['configuration_empirical_r']=r
        frame['conditional_historical']=frame.observed_yeast_fraction*combined_factor(.276123)
        frame['conditional_configuration_empirical']=frame.observed_yeast_fraction*combined_factor(r)
        classes.append(frame)
    fields=['observations','human','yeast','species_exclusive_denominator','observed_yeast_fraction','r_free_ratio_vs_retained',
        'historical_r','configuration_empirical_r','conditional_historical','conditional_configuration_empirical',
        'same_cell_group_in_final_fraction','same_cell_group_other_peptide_in_final_fraction']
    cls=paired_tables(*classes,['filter','class'],fields)
    cls['change_yeast_fraction_percentage_points']=100*(cls.norefine_observed_yeast_fraction-cls.original_observed_yeast_fraction)
    cls['yeast_fraction_fold_change']=[fraction_ratio(new,old) for new,old in zip(cls.norefine_observed_yeast_fraction,cls.original_observed_yeast_fraction)]
    cls.to_csv(output/'class_species_scales_comparison.tsv',sep='\t',index=False,na_rep='NA')
    fixedcols=['human','yeast','species_exclusive_denominator','observed_yeast_fraction','conditional_combined_historical',
        'conditional_combined_empirical','historical_r','empirical_library_r']
    fixed_old=read(old4/'conditional_scales.tsv');fixed_new=read(new4/'conditional_scales.tsv')
    paired_tables(fixed_old[fixed_old.domain.eq('class')],fixed_new[fixed_new.domain.eq('class')],['filter','class'],fixedcols).to_csv(output/'frozen_fixed_scale_comparison.tsv',sep='\t',index=False,na_rep='NA')
    f_old=read(old5/'support_wilson.tsv');f_new=read(new5/'support_wilson.tsv')
    f_old=f_old[f_old.cohort.eq('all_cells_64')];f_new=f_new[f_new.cohort.eq('all_cells_64')]
    supportcols=['observations','human','yeast','n','observed_yeast_fraction','wilson_low95','wilson_high95','point_label','zero_yeast_label']
    support=paired_tables(f_old,f_new,['filter','bin_scheme','quantity_tercile','support_bin'],supportcols)
    support['change_yeast_fraction_percentage_points']=100*(support.norefine_observed_yeast_fraction-support.original_observed_yeast_fraction)
    support.to_csv(output/'support_wilson_comparison.tsv',sep='\t',index=False,na_rep='NA')
    ci=['point','low95','high95','draws','undefined_draws','infinite_draws','interval_status','yeast','species_exclusive_denominator','precursor_blocks']
    paired_tables(read(old5/'class_intervals.tsv'),read(new5/'class_intervals.tsv'),['class','scheme'],ci).to_csv(output/'class_interval_comparison.tsv',sep='\t',index=False,na_rep='NA')
    ri=['point','low95','high95','draws','undefined_draws','infinite_draws','interval_status','numerator_yeast','numerator_denominator','reference_yeast','reference_denominator','precursor_blocks']
    paired_tables(read(old5/'class_ratio_intervals.tsv'),read(new5/'class_ratio_intervals.tsv'),['contrast','scheme'],ri).to_csv(output/'class_ratio_interval_comparison.tsv',sep='\t',index=False,na_rep='NA')
    rc_old=read(old5/'support_contrasts.tsv');rc_new=read(new5/'support_contrasts.tsv')
    rc_old=rc_old[rc_old.cohort.eq('all_cells_64')];rc_new=rc_new[rc_new.cohort.eq('all_cells_64')]
    rc=['numerator_yeast','numerator_n','reference_yeast','reference_n','fraction_ratio','ratio_display']
    paired_tables(rc_old,rc_new,['filter','bin_scheme','quantity_tercile','support_bin','reference_bin'],rc).to_csv(output/'support_contrast_comparison.tsv',sep='\t',index=False,na_rep='NA')
    assessments=[]
    for configuration,ac,cv,cl,sp in [('original_default',a,read(old4/'represented_coverage.tsv'),classes[0],f_old),
        ('disabled_refinement',b,read(new4/'represented_coverage.tsv'),classes[1],f_new)]:
        for filt in ['F1','F2','F3']:
            account=ac[ac['filter'].eq(filt)].iloc[0];coverage_row=cv[cv['filter'].eq(filt)&cv.level.eq('precursor')].iloc[0]
            classrows=cl[cl['filter'].eq(filt)].set_index('class')
            assessments.append({'configuration':configuration,'filter':filt,'accepted_observations_gain':account.final_observations>account.first_observations,
                'represented_precursors_decrease':coverage_row.final_represented_entities<coverage_row.first_represented_entities,
                'common_row_completeness_gain':coverage_row.final_common_row_completeness>coverage_row.first_common_row_completeness,
                'second_pass_only_vs_retained_yeast_fraction_ratio':classrows.loc['second_pass_only','r_free_ratio_vs_retained'],
                'lost_vs_retained_yeast_fraction_ratio':classrows.loc['lost','r_free_ratio_vs_retained']})
        for bins,labels in [('current_three_positive',['low_positive','middle','high']),
            ('original_four_positive',['positive_0_20','positive_20_40','positive_40_60','positive_60_100'])]:
            rows=sp[sp.bin_scheme.eq(bins)&sp.quantity_tercile.eq('high')].set_index('support_bin').loc[labels]
            values=rows.observed_yeast_fraction.to_numpy()
            assessments.append({'configuration':configuration,'filter':'F1','support_bin_scheme':bins,
                'high_tercile_all_positive_bins_evaluable':bool(np.isfinite(values).all()),
                'high_tercile_positive_bins_nonincreasing':bool(np.all(np.diff(values)<=0)) if np.isfinite(values).all() else pd.NA,
                'high_tercile_low_and_high_both_have_yeast':bool(rows.yeast.iloc[0]>0 and rows.yeast.iloc[-1]>0),
                'high_tercile_low_vs_high_ratio':fraction_ratio(values[0],values[-1]),
                'high_tercile_low_yeast':int(rows.yeast.iloc[0]),'high_tercile_low_n':int(rows.n.iloc[0]),
                'high_tercile_high_yeast':int(rows.yeast.iloc[-1]),'high_tercile_high_n':int(rows.n.iloc[-1])})
    pd.DataFrame(assessments).to_csv(output/'pattern_assessment.tsv',sep='\t',index=False,na_rep='NA')
