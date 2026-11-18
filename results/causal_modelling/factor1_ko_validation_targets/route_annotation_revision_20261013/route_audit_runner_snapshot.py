import csv,json,collections,hashlib,copy
from pathlib import Path
ROOT=Path('/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored')
b=ROOT/'results/causal_modelling/factor1_ko_validation_targets/alternating_prioritization';review=json.load(open(ROOT/'support/cosmos_mouse/factor1_ko_validation_review.json'));pa=json.load(open(b/'path_audit.json'))
base=ROOT/'results/causal_modelling/factor1_tf_to_lipid_expanded_species';mapping={r['node_id']:r['feature_id']for r in csv.DictReader((base/'inputs/input_mapping.csv').open())};weights={r['feature_id']:float(r['weight'])for r in csv.DictReader((base/'inputs/input_mapping.csv').open())};labels={r['node_id']:r for r in csv.DictReader((base/'inputs/nodes.csv').open())};tf={r['source']:r for r in csv.DictReader((ROOT/review['frozen_tf_activity_file']).open())}
seed_files={layer:base/('currency_v2_expanded_layers'+('4'if layer=='layers4'else'6')+'_lower_thresholds')/'upstream_input.csv'for layer in ['layers4','layers6']};seeds={layer:{r['source'].removeprefix('Gene__'):float(r['score'])for r in csv.DictReader(fn.open())}for layer,fn in seed_files.items()};assert seeds['layers4']==seeds['layers6'];assert len(seeds['layers4'])==49
edges={};attrs={}
for layer,path in review['networks'].items():
 d=ROOT/path;edges[layer]={(r['source'],r['target']):r for r in csv.DictReader((d/'network_labeled_SIF.csv').open())};attrs[layer]={r['node_id']:r for r in csv.DictReader((d/'network_labeled_ATT.csv').open())}
records={}
for fn in [ROOT/'results/causal_modelling/cosmos_mouse_pkn/lipid_rules_reactions.jsonl',ROOT/'results/causal_modelling/cosmos_mouse_pkn/rhea_reactions.jsonl',ROOT/'results/causal_modelling/cosmos_mouse_pkn_expanded_species/expanded_species_reactions.jsonl']:
 for l in fn.open():
  r=json.loads(l);records[r['reaction_id']]=r
warnings=[{'edge':'Nr2f1 -> Rxra','pmid':'10900149','url':'https://pubmed.ncbi.nlm.nih.gov/10900149/','issue':'COUP-TFI/RXR heterodimers bind an ALDH2 promoter element and suppress RXR-associated reporter induction; abstract does not establish context-free positive NR2F1-to-RXRA causal activation','status':'Raw positive fitted edge retained for topology audit, but its sign/context is questioned; do not present as verified activating relay'}, {'edge':'Nr4a2 -> Rxra -> Pld1','pmid':'17433303','url':'https://pubmed.ncbi.nlm.nih.gov/17433303/','issue':'NR4A2/RXRA association is a different receptor-complex context from vitaminD/VDR/RXRA induction of PLD1 in human HaCaT cells; concatenation needs unmodeled partner/ligand transfer','status':'Model-contiguous route, context-specific mechanistic hypothesis'}, {'edge':'Hoxa10 -> Ednra','pmid':'20371740','url':'https://pubmed.ncbi.nlm.nih.gov/20371740/','issue':'Source transcription evidence is human endometrial stromal cells, not hepatocytes; Hoxa10 expression is weak/absent in historical hepatic culture controls','status':'Six-layer fitted topology with hepatic transfer uncertain'}]
chemistry={}
def reaction_info(rid):
 if rid in chemistry:return
 r=records[rid]
 def parts(side):return [{k:p.get(k)for k in ['id','name','stoichiometry','compartment','measurement_feature_id','template_chemical_id','conditional_acyl_chain','unmeasured']}for p in r[side]]
 chemistry[rid]={'reaction_id':rid,'resource':r['resource'],'supporting_reaction_id':r.get('supporting_reaction_id'),'master_reaction_id':r.get('master_reaction_id'),'substrates':parts('substrates'),'products':parts('products'),'equation':r.get('equation'),'source_template_equation':r.get('source_template_equation'),'direction':r['direction'],'rule_family':r.get('rule_family'),'matched_feature_ids':r.get('matched_feature_ids'),'assumptions':r.get('assumptions',[]),'match_details':r.get('match_details'),'graph_semantics':'Catalyst capacity -> reaction instance -> product; complete joint substrates and stoichiometry are not enforced in MOON','annotated_catalysts':[{'gene':p.get('gene_symbol'),'uniprot':p.get('uniprot'),'source_url':p.get('source_url'),'rhea':p.get('annotated_rhea_id')}for p in r.get('provenance',[])if p.get('gene_symbol')]}
def route_record(layer,gene,feature,p,kind,target):
 genes=[n[6:]for n in p['nodes']if n.startswith('Gene__')];rids=[];ge=[]
 for n in p['nodes']:
  if n.startswith('Enzyme__'):
   rid=labels[n]['reaction_id'];reaction_info(rid);rids.append(rid)
 for x,y,s in zip(p['nodes'],p['nodes'][1:],p['signs']):
  if x.startswith('Gene__')and y.startswith('Gene__'):
   e=edges[layer].get((x,y));ge.append({'source':x[6:],'target':y[6:],'sign':s,'evidence_kind':'source_corrected_unscored_prefix'if e is None else'frozen_graph_edge','mechanism':e.get('mechanisms')if e else'proposed transcription from corrected primary source','references':json.loads(e['references_json'])if e else p['source_corrected_prefix']['source_ids']})
 activity=float(tf[gene]['oriented_activity'])if gene in tf else None
 warnings_here=[w for w in warnings if (w['edge']=='Nr2f1 -> Rxra'and'Nr2f1'in genes)or(w['edge']=='Nr4a2 -> Rxra -> Pld1'and'Nr4a2'in genes and'Pld1'in genes)or(w['edge']=='Hoxa10 -> Ednra'and'Hoxa10'in genes and'Ednra'in genes)]
 return {'network':layer,'route_kind':kind,'source_estimated_TF':gene,'source_TF_activity':activity,'actually_supplied_upstream_TF':gene in seeds[layer],'KO_target':target,'KO_target_gene_index':genes.index(target),'genes':genes,'nodes':p['nodes'],'edge_signs':p['signs'],'source_levels':p['levels'],'gene_edges':ge,'reaction_ids':rids,'readout':feature,'readout_factor1_weight':weights[feature],'route_sign':p['path_sign'],'factor1_source_sign_agrees':activity*p['path_sign']*weights[feature]>0 if activity is not None else None,'conditional_source_TF_KO_direction':'down'if p['path_sign']==1else'up','conditional_local_target_KO_direction':'down','source_hypothesis_id':p.get('source_hypothesis_id'),'minimum_full_path_inferred_absolute_score':p['minimum_full_path_inferred_absolute_score'],'minimum_frozen_portion_inferred_absolute_score':p['minimum_frozen_portion_inferred_absolute_score'],'mechanistically_fully_verified':False,'source_warnings':warnings_here}
all_sets=[]
for item in review['sets']:
 wanted={f for gr in item['primary_readout_groups']for f in review['readout_groups'][gr]};branches=[]
 for target in item['targets']:
  frozen={};corrected={};partial={}
  for layer in ['layers4','layers6']:
   rr=[]
   for seed in seeds[layer]:
    for feature,rows in pa['frozen_screened'][layer].get(seed,{}).items():
     if feature not in wanted:continue
     for p in rows:
      if'Gene__'+target in p['nodes']:rr.append(route_record(layer,seed,feature,p,'frozen_seed_to_lipid_path',target))
   frozen[layer]=rr
   rr=[]
   for source,ep in pa['source_corrected_extensions'][layer].items():
    for feature,rows in ep.items():
     if feature not in wanted:continue
     for p in rows:
      if'Gene__'+target in p['nodes']:rr.append(route_record(layer,source,feature,p,'source_corrected_unscored_estimated_TF_extension',target))
   corrected[layer]=rr
   rr=[]
   for feature,rows in pa['frozen_screened'][layer].get(target,{}).items():
    if feature not in wanted:continue
    for p in rows:
     for n in p['nodes']:
      if n.startswith('Enzyme__'):reaction_info(labels[n]['reaction_id'])
     rr.append({'nodes':p['nodes'],'edge_signs':p['signs'],'readout':feature,'reaction_ids':[labels[n]['reaction_id']for n in p['nodes']if n.startswith('Enzyme__')],'route_kind':'local_effector_to_lipid_frozen_subpath','no_invented_TF_prefix':True})
   partial[layer]=rr
  representative='layers4'if frozen['layers4']else('layers6'if frozen['layers6']else None)
  if target=='Pla2g4a':gap='No surviving fitted full TF prefix: Gata6->Mapk1 excluded. Mapk1->Pla2g4a local kinase/lipase relationship remains, but no alternate supplied-TF prefix is present in either selected graph.'
  elif target=='Plpp1':gap='No surviving fitted full TF prefix: Ppara->Plpp1 excluded. Source-corrected Pparg->Plpp1 is available only as an unscored extension from an estimated TF below the original upstream cutoff.'
  elif target=='Plcb1':gap='No surviving four-layer supplied-TF prefix after Gata6->Mapk1 exclusion; six-layer Hoxa10 receptor/G-protein routes remain as context-dependent topology.'
  elif target in {'Cebpb','Pparg'}:gap='Target is represented by a source-corrected unscored TF prefix, not by a fitted selected-network TF-to-lipid path.'
  else:gap=None
  branches.append({'target':target,'preferred_frozen_network':representative,'full_frozen_TF_route_exists':representative is not None,'full_frozen_seed_routes':frozen,'source_corrected_estimated_TF_extensions':corrected,'local_frozen_effector_subpaths':partial,'gap_or_qualification':gap})
 all_sets.append({'rank':item['rank'],'set_id':item['set_id'],'targets':item['targets'],'category':item['category'],'primary_readouts':sorted(wanted),'branches':branches,'compensation':item['compensation'],'conditional_gates':item['conditional_gates']})
seed_inventory=[]
for g,a in seeds['layers4'].items():seed_inventory.append({'gene':g,'activity':a,'actually_supplied':True,'present_in_selected_layers4':'Gene__'+g in attrs['layers4'],'present_in_selected_layers6':'Gene__'+g in attrs['layers6'],'surviving_screened_readouts_layers4':list(pa['frozen_screened']['layers4'].get(g,{})),'surviving_screened_readouts_layers6':list(pa['frozen_screened']['layers6'].get(g,{}))})
res={'analysis':'read_only_full_TF_to_lipid_route_audit','source_warning_ledger':warnings,'supplied_upstream_TF_count':49,'supplied_TFs':seed_inventory,'estimated_but_not_supplied_TFs':[{ 'gene':g,'activity':float(tf[g]['oriented_activity']),'in_original_mapped_PKN':tf[g]['in_pkn'],'reason':'Below abs(activity)>2 cutoff; no fitted upstream input constraint'}for g in ['Rxra','Pparg']],'sets':all_sets,'complete_reactions':chemistry,'global_limits':['Original fitted graphs, scores and KO ranks unchanged; no new MOON run','No repaired upstream path is silently added; source-corrected extensions explicitly unscored','Graph positive metabolic signs mean permissive formation/participation, not concentration control','Complete chemical substrates/co-products shown in records are not joint-availability/flux constraints in MOON','A model-contiguous TF prefix can connect receptor/heterodimer observations from different ligand/cell contexts; its direction is not thereby experimentally established']}
Path('/private/tmp/ko_full_route_audit.json').write_text(json.dumps(res,indent=2)+'\n')
print('Completed',len(all_sets),'sets;',len(chemistry),'complete reaction records')
for s in all_sets:
 print(s['rank'],s['targets'],[(x['target'],len(x['full_frozen_seed_routes']['layers4']),len(x['full_frozen_seed_routes']['layers6']),len(x['source_corrected_estimated_TF_extensions']['layers4']))for x in s['branches']])
