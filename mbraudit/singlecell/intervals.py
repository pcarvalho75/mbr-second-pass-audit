"""Joint two-way bootstrap intervals for class fractions and class ratios (carried from analysis/revision5_analysis.py)."""
import numpy as np
import pandas as pd


SEED=20261018


DRAWS=2000


BASE_CLASSES=('retained','second_pass_only','lost_library_present','lost_library_excluded')


SCHEMES=('recipient','precursor','two_way')


def fraction(e,n):
    e,n=np.broadcast_arrays(np.asarray(e,dtype=float),np.asarray(n,dtype=float))
    return np.divide(e,n,out=np.full(e.shape,np.nan),where=n>0)


def derive_counts(base):
    result=dict(base)
    result['lost']=(base['lost_library_present'][0]+base['lost_library_excluded'][0],
                    base['lost_library_present'][1]+base['lost_library_excluded'][1])
    result['first_pass']=(base['retained'][0]+result['lost'][0],base['retained'][1]+result['lost'][1])
    result['final']=(base['retained'][0]+base['second_pass_only'][0],base['retained'][1]+base['second_pass_only'][1])
    return result


def joint_draws(frames,cells,scheme,draws=DRAWS,seed=SEED,batch=100):
    """Reuse the fixed block/multiplicity code; share every draw across all classes."""
    from mbraudit.singlecell.bootstrap import matrices,multiplicities
    arrays,keys=matrices(frames,cells)
    rng=np.random.default_rng(seed)
    cw=multiplicities(rng,len(cells),draws) if scheme in ('recipient','two_way') else None
    output={name:[[],[]] for name in BASE_CLASSES}
    for start in range(0,draws,batch):
        stop=min(start+batch,draws)
        pw=multiplicities(rng,len(keys),stop-start) if scheme in ('precursor','two_way') else None
        for name,channels in zip(BASE_CLASSES,arrays):
            for index,channel in enumerate(channels[1:]):
                if scheme=='recipient':values=cw[start:stop]@channel.sum(axis=1)
                elif scheme=='precursor':values=pw@channel.sum(axis=0)
                elif scheme=='two_way':values=(cw[start:stop]@channel*pw).sum(axis=1)
                else:raise ValueError('Unknown scheme')
                output[name][index].extend(values.tolist())
    base={name:tuple(np.asarray(x) for x in values) for name,values in output.items()}
    return derive_counts(base),len(keys)


def class_intervals(reports,loaded,output):
    from mbraudit.singlecell.accounting import filtered_pair
    from mbraudit.singlecell.reports import KEY
    from mbraudit.singlecell.bootstrap import interval,ratio
    a,b=filtered_pair(loaded,'F1')
    final=set(zip(b.Run,b.key))
    retained=a.loc[[(c,k) in final for c,k in zip(a.Run,a.key)]]
    lost=a.loc[[(c,k) not in final for c,k in zip(a.Run,a.key)]]
    only=b.loc[~b.direct]
    lib=reports.library();lib=lib[lib.Decoy.eq(0)]
    library_keys=set(zip(lib[KEY[0]],lib[KEY[1]].astype(int)))
    frames=[retained,only,lost[lost.key.isin(library_keys)],lost[~lost.key.isin(library_keys)]]
    points=derive_counts({name:(int(f.species.eq('entrapment').sum()),int(f.species.isin(['target','entrapment']).sum())) for name,f in zip(BASE_CLASSES,frames)})
    points_rows=[{'class':name,'observations':len(f),'human':int(f.species.eq('target').sum()),
        'yeast':points[name][0],'species_exclusive_denominator':points[name][1],
        'observed_yeast_fraction':float(fraction(*points[name]))} for name,f in zip(BASE_CLASSES,frames)]
    pd.DataFrame(points_rows).to_csv(output/'class_points.tsv',sep='\t',index=False)
    intervals=[];contrasts=[]
    for scheme in SCHEMES:
        counts,nkeys=joint_draws(frames,loaded['cells'],scheme)
        ledger={'draw':np.arange(1,DRAWS+1)}
        for name,(e,n) in counts.items():
            values=fraction(e,n);stats=interval(values)
            ledger[name+'_yeast']=e;ledger[name+'_denominator']=n;ledger[name+'_fraction']=values
            intervals.append({'class':name,'scheme':scheme,'yeast':points[name][0],
                'species_exclusive_denominator':points[name][1],'point':float(fraction(*points[name])),
                **stats,'seed':SEED,'recipients':len(loaded['cells']),'precursor_blocks':nkeys,
                'bootstrap_universe':'joint union of four disjoint F1 observation classes'})
        for numerator in ['second_pass_only','lost']:
            values=ratio(counts[numerator][0],counts[numerator][1],counts['retained'][0],counts['retained'][1])
            point=float(ratio(*points[numerator],*points['retained']))
            label=numerator+'_vs_retained';ledger[label]=values
            contrasts.append({'contrast':label,'scheme':scheme,'point':point,
                'numerator_yeast':points[numerator][0],'numerator_denominator':points[numerator][1],
                'reference_yeast':points['retained'][0],'reference_denominator':points['retained'][1],
                **interval(values),'seed':SEED,'recipients':len(loaded['cells']),'precursor_blocks':nkeys,
                'pairing':'same recipient and precursor multiplicities shared across numerator and reference'})
        pd.DataFrame(ledger).to_csv(output/('joint_draws_'+scheme+'.tsv'),sep='\t',index=False,na_rep='NA')
        print('Joint class/ratio bootstrap',scheme,'complete;',nkeys,'precursor blocks.',flush=True)
    pd.DataFrame(intervals).to_csv(output/'class_intervals.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(contrasts).to_csv(output/'class_ratio_intervals.tsv',sep='\t',index=False,na_rep='NA')
