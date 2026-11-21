"""Independent read-only evidence audit; writes review JSON and diagnostic plot.
Does not edit transport inputs, inventories or their hashes.
"""
from pathlib import Path
import csv,json,re,hashlib,math
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
W=ROOT/'non_upload/operating_coefficients'

def main():
 rows=list(csv.DictReader((ROOT/'data/operating_states.csv').open()))
 coef=list(csv.DictReader((ROOT/'data/operating_coefficients.csv').open()))
 evidence=[];fig,axes=plt.subplots(2,1,figsize=(10,7))
 for r in rows:
  d=W/'openmc_cases'/r['case_id']
  with h5py.File(d/'statepoint.200.h5') as f:
   assert int(f['n_particles'][()])==40000 and int(f['n_inactive'][()])==50 and int(f['current_batch'][()])==200
   assert int(f['seed'][()])==int(r['seed'])
   np.testing.assert_allclose(f['k_combined'][()],[float(r['k_openmc']),float(r['k_openmc_sd'])],rtol=1e-14)
   x=np.array(f['k_generation']);e=np.array(f['entropy']) if 'entropy' in f else np.array([])
   a=x[50:];acf=float(np.corrcoef(a[:-1],a[1:])[0,1])
   axes[0].plot(np.arange(1,201),(x-a.mean())/a.std(ddof=1),alpha=.18,lw=.5)
   if len(e):axes[1].plot(np.arange(1,201),(e-e[50:].mean())/e[50:].std(ddof=1),alpha=.18,lw=.5)
   version=str(f.attrs.get('openmc_version','not stored'));git=str(f.attrs.get('git_sha1','not stored'))
  prov=json.loads((d/'run_provenance.json').read_text())
  text=(d/'openmc_stdout.log').read_text(errors='replace')
  warnings=[line.strip() for line in text.splitlines() if re.search(r'WARNING|ERROR',line)]
  assert all('Ignoring number of threads specified on command line' in line for line in warnings),(r['case_id'],warnings)
  evidence.append(dict(case=r['case_id'],lag1=acf,entropy=bool(len(e)),schedule=prov.get('schedule','not recorded in per-run record'),version=version,git=git,tallies=(d/'tallies.xml').exists(),warnings=warnings))
 for ax in axes:
  ax.axvline(50,color='black',ls='--');ax.set_xlim(1,200);ax.set_ylim(-6,6);ax.set_xlabel('Batch');ax.grid(alpha=.15)
 axes[0].set_ylabel('k-generation standardized to active batches')
 axes[1].set_ylabel('Entropy standardized to active batches')
 fig.tight_layout();fig.savefig(W/'source_diagnostic_overlay.png',dpi=150);plt.close(fig)
 # Local preserved copies and installed accepted listings must match exactly.
 checks=[]
 for prefix,n in [('OP_CONV_J2',21),('OP_M40_J2',11)]:
  for b in ('00p0','07p5','20p0'):
   name=f'{prefix}_BU{b}.result';p=W/'dragon_results'/name
   t=p.read_text();assert 'CONVERGENCE NOT REACHED' not in t and '>|VRF_COMPLETED' in t
   assert len(re.findall(r'^>\|OP_RESULT ',t,re.M))==n
   original=Path('/Users/alokjha/DRAGON/5.1/Dragon/Darwin_arm64')/name
   assert p.read_bytes()==original.read_bytes()
   checks.append(dict(file=name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),states=n))
 # Use Decimal arithmetic as an independent reciprocal-k evaluation.
 from decimal import Decimal,localcontext
 with localcontext() as ctx:
  ctx.prec=40
  for c in coef:
   rr=[r for r in rows if all(float(r[k])==float(c[k]) for k in ('burnup_GWd_tU','step_K')) and r['variable']==c['variable'] and r['density_mode']==c['density_mode']]
   lo=next(r for r in rr if float(r['sign'])<0);hi=next(r for r in rr if float(r['sign'])>0)
   for code in ('dragon','openmc'):
    val=Decimal(100000)*(1/Decimal(lo['k_'+code])-1/Decimal(hi['k_'+code]))/(2*Decimal(c['step_K']))
    assert abs(float(val)-float(c[code+'_pcm_K']))<1e-10
 steps=list(csv.DictReader((ROOT/'data/operating_step_convergence.csv').open()))
 density=list(csv.DictReader((ROOT/'data/operating_density_contributions.csv').open()))
 assert len(steps)==15 and len(density)==12
 for d in density:
  a=[c for c in coef if c['variable']==d['variable'] and all(float(c[k])==float(d[k]) for k in ('burnup_GWd_tU','step_K'))]
  fixed=next(c for c in a if c['density_mode']=='fixed');path=next(c for c in a if c['density_mode']=='pressure')
  for code in ('dragon','openmc'):
   assert abs(float(d[code+'_density_contribution_pcm_K'])-(float(path[code+'_pcm_K'])-float(fixed[code+'_pcm_K'])))<1e-10
  assert math.isclose(float(d['openmc_sd_pcm_K']),math.hypot(float(fixed['openmc_sd_pcm_K']),float(path['openmc_sd_pcm_K'])),rel_tol=1e-12)
 diagnostics=list(csv.DictReader((ROOT/'data/operating_source_convergence.csv').open()))
 zsummary={}
 for key in ('z_active_halves','z_inactive_tail','entropy_z_inactive_tail'):
  values=np.array([float(r[key]) for r in diagnostics]);values=values[np.isfinite(values)]
  zsummary[key]=dict(cases=len(values),mean=float(values.mean()),sd=float(values.std(ddof=1)),max_abs=float(max(abs(values))))
 out=dict(status='passed_scoped_numerical_review',source_diagnostic_summary=zsummary,states=63,decimal_reciprocal_checks=60,
  max_all_step_intercode_z=max(abs(float(c['difference_z_independent'])) for c in coef),
  max_step_difference_z=max(abs(float(c['openmc_step_change_pcm_K'])/float(c['openmc_step_change_sd_pcm_K'])) for c in steps),
  mean_lag1=float(np.mean([r['lag1'] for r in evidence])),
  schedule_missing_records=sum(r['schedule'].startswith('not recorded') for r in evidence),
  schedule_gap='Runtime schedule absent from original per-run records. Current runner records it prospectively; do not backfill an unobserved historical environment.',
  solver_default_source='Installed Dragon/src/FLUGPI.f:101 ROUMIN=1.0E-5; EPSOUT=ROUMIN at162; EPSUNK defaults EPSINR at431. Online manual differs.',
  accepted_listing_checks=checks,source_records=evidence,
  limits='Finite stationarity/lag checks do not prove source convergence or physical correctness. No material diagnostic triggers observed; no runs discarded or rerun selectively.')
 (ROOT/'qa/operating_independent_review.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k not in ('accepted_listing_checks','source_records')},indent=2))
if __name__=='__main__':main()
