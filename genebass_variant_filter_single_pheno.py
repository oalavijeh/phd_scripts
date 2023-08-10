import hail as hl

def bhr_genebass_variant_to_gene_lof_syn(var_type, upper_bound, lower_bound, name):
    #Load variants and phenotype files
    genebass_variant = hl.read_matrix_table('gs://ukbb-exome-public/500k/results/variant_results.mt')
    genebass_variant = genebass_variant.rename({'AF.Cases': 'AF_Cases', 'AF.Controls': 'AF_Controls'})
    genebass_variant = genebass_variant.filter_cols(genebass_variant.phenocode == "132532")
    #Filter variants and count
    genebass_variant = genebass_variant.filter_rows((genebass_variant.call_stats.AF < upper_bound) & (genebass_variant.call_stats.AF >= lower_bound), keep=True)
    genebass_variant = genebass_variant.filter_rows(genebass_variant.annotation == var_type)
    nvar = genebass_variant.count()[0]
    
    #Add BHR-specific parameters
    genebass_variant = genebass_variant.annotate_entries(af_overall = ((genebass_variant.n_cases*genebass_variant.AF_Cases) + (genebass_variant.n_controls*genebass_variant.AF_Controls))/(genebass_variant.n_cases + genebass_variant.n_controls),
                    prevalence = genebass_variant.n_cases/(genebass_variant.n_cases + genebass_variant.n_controls))
    genebass_variant = genebass_variant.annotate_entries(beta_binary = ((2*genebass_variant.prevalence*(genebass_variant.AF_Cases - genebass_variant.af_overall))/hl.sqrt(2*genebass_variant.af_overall*(1-genebass_variant.af_overall)*genebass_variant.prevalence*(1-genebass_variant.prevalence))))
    
    genebass_variant = genebass_variant.annotate_entries(variant_variance = hl.if_else(genebass_variant.trait_type == "continuous",
                                                            ((1/(genebass_variant.SE*hl.sqrt(genebass_variant.n_eff)))**2),
                                                            2*genebass_variant.af_overall*(1-genebass_variant.af_overall)))
                                                                                                                                                    
    genebass_variant = genebass_variant.annotate_entries(beta_per_allele = hl.if_else(genebass_variant.trait_type == "continuous",
                                                          genebass_variant.BETA,
                                                          genebass_variant.beta_binary/(hl.sqrt(genebass_variant.variant_variance))))
                                                                                                                                                               
                 
    #Export gene summary statistics file
    genebass_variant.entries().export('~/dec_bhr_ms_variant_ss_400k_final_thin_withnullburden_'+str(var_type)+'_nvar'+str(nvar)+'_low'+str(lower_bound)+'_high'+str(upper_bound)+'_'+str(name)+'.txt.bgz')
     
     
#variant frequency and function filters
upper = [1e-5,0.0001,0.001]
lower = [0,1e-5,0.0001]
names = ['group1','group2','group3']
var_type = ['pLoF']

for variant_type in range(len(var_type)):
    for grouping in range(len(names)):
        bhr_genebass_variant_to_gene_lof_syn(var_type[variant_type], upper[grouping], lower[grouping], names[grouping])
