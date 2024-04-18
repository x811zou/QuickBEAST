import argparse
import multiprocessing
import os
import sys
import subprocess
import tempfile
import shutil
import pandas as pd
from scipy.stats import norm
import numpy as np
from scipy import stats
from datetime import datetime, timedelta
import time
import logging
from multiprocessing import Pool
from numpy.random import uniform
from scipy.stats import t
from scipy.optimize import minimize
# Configure logging
logging.basicConfig(level=logging.INFO)

"""
python /home/scarlett/github/QuickBEAST/calculate_p_value_from_qb_mode.py test_data/test_new_pis test_data/qb_test_new_pis
python /home/scarlett/github/QuickBEAST/calculate_p_value_from_qb_mode.py /data2/simulation/parametrized/ASE_0.05_error/g-1000_h-3_d-30_t-1.txt /data2/stan/quickBEAST/a8.789625_b8.789625/lambda0.04545/parametrized/ASE_0.05_error/withcache/g-1000_h-3_d-30_t-1.txt
python /home/scarlett/github/QuickBEAST/calculate_p_value_from_qb_mode.py /data2/NA12878/beastie/runModel_phased_even101/chr1-22_alignBiasp0.05_s0.7_a0.05_sinCov0_totCov1_W1000K1000/tmp/NA12878_real_alignBiasafterFilter.phasedByshapeit2.cleaned.modelinput.tsv /data2/NA12878/qb/NA12878_highestcoverage.tsv
""" 

def calculate_norm_2sided_pvalue(loc, scale, sample_X):
    # Get the cumulative probability of X (left tail)
    left_tail_p = stats.norm.cdf(sample_X, loc, scale)
    
    # Get the survival function value for X (right tail)
    right_tail_p = stats.norm.sf(sample_X, loc, scale)
    
    # Return the double-sided p-value
    return 2 * min(left_tail_p, right_tail_p)

def calculate_t_2sided_pvalue(df, loc, scale, sample_X):
    # Get the cumulative probability of X (left tail)
    left_tail_p = stats.t.cdf(sample_X, df, loc, scale)
    
    # Get the survival function value for X (right tail)
    right_tail_p = stats.t.sf(sample_X, df, loc, scale)
    
    # Return the double-sided p-value
    return 2 * min(left_tail_p, right_tail_p)

def calculate_st_2sided_pvalue(a, loc, scale, sample_X):
    # Get the cumulative probability of X (left tail)
    left_tail_p = stats.skewnorm.cdf(sample_X, a, loc, scale)
    
    # Get the survival function value for X (right tail)
    right_tail_p = stats.skewnorm.sf(sample_X, a, loc, scale)
    
    # Return the double-sided p-value
    return 2 * min(left_tail_p, right_tail_p)

def generate_fields(geneID,M, D, theta,switching_error = 0.05):
    # calculate probability for binomial distribution
    p = theta / (1.0 + theta)
    # calculate alternative and reference read counts for each het
    alt_counts = np.random.binomial(D, p, M)
    ref_counts = D - alt_counts
    # construct output fields
    switched = False
    fields = [geneID, M]
    
    for counts in zip(alt_counts, ref_counts):
        if switched:
            fields.append(counts[1])
            fields.append(counts[0])
        else:
            fields.append(counts[0])
            fields.append(counts[1])
        if np.random.uniform() <= switching_error:
            switched = not switched
    
    fields += [switching_error] * (M - 1)  # Phasing error (We assume all pairs to be -1)

    return fields

def run_qb_parallel(genes):
    parameter = 8.789625

    # prepare files
    temp_dir = tempfile.mkdtemp()
    input_file_path = os.path.join(temp_dir, f'input.txt')
    output_file_path = os.path.join(temp_dir, f"output.txt")

    # conversion.write_genes_to_quickbeast_input_file(genes, input_file_path)
    with open(input_file_path, 'w') as file:
        for gene in genes:
            file.write('\t'.join(map(str, gene)) + '\n')
    
    #run qb
    try:
        subprocess.run([f"./QuickBEAST --alpha {parameter} --beta {parameter} --mean --mode -f {input_file_path} --fixMaxHetSite > {output_file_path}"], check=True, shell=True)
        # logging.info(f"QB for {input_file_path} executed successfully")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running QB for {input_file_path}: {e}")

    #read and convert results
    qb = pd.read_csv(output_file_path, delimiter="\t", header=0)
    qb.columns = ['geneID', 'qb_mean', 'qb_var', 'qb_mode']
    # cleanup
    shutil.rmtree(temp_dir)

    return qb

def format_geneID(geneID):
    # Find the index of the last period in the string
    last_period_index = geneID.rfind(".")
    
    # If there is a period in the string, remove characters after it
    if last_period_index != -1:
        return geneID[:last_period_index]
    
    # If there is no period, return the original string
    return geneID

def calculate_time(start_t, end_t):
    elapsed_time_seconds = end_t - start_t
    elapsed_time_formatted = str(timedelta(seconds=elapsed_time_seconds))
    return elapsed_time_formatted

def simulate_null_genes(number_of_hets, average_read_depth_per_het, pi=0.05):
    mode_parameter_fit_st = "PASSED"
    mean_parameter_fit_st = "PASSED"
    # SIMULATE GENE INPUTS
    simulated_genes = []
    for i in range(1000):
        gene_id = f"gene_{i}"
        simulated_genes.append(generate_fields(gene_id, M=number_of_hets, D=average_read_depth_per_het, theta=1, switching_error = pi))

    # RUN qb simulated genes 
    simulated_gene_results = run_qb_parallel(simulated_genes)
    # compute summary statistics
    mode = simulated_gene_results["qb_mode"].tolist()
    mean = simulated_gene_results["qb_mean"].tolist()
    # df: degrees of freedom
    # loc: location
    # scale: scale of std
    # mode
    try:
        st_df_mode, st_loc_mode, st_scale_mode = stats.skewnorm.fit(mode)
    except Exception as e:
        mode_parameter_fit_st = "FAILED"
        st_df_mode, st_loc_mode, st_scale_mode = manual_fitting(number_of_hets,average_read_depth_per_het,e,mode)
    # mean
    try:
        st_df_mean, st_loc_mean, st_scale_mean = stats.skewnorm.fit(mean)
    except Exception as e:
        mean_parameter_fit_st = "FAILED"
        st_df_mean, st_loc_mean, st_scale_mean = manual_fitting(number_of_hets,average_read_depth_per_het,e,mean)
    
    return number_of_hets, average_read_depth_per_het*number_of_hets, st_df_mode, st_loc_mode, st_scale_mode, mode_parameter_fit_st, st_df_mean, st_loc_mean, st_scale_mean, mean_parameter_fit_st

def manual_fitting(number_of_hets,average_read_depth_per_het,e,values):
    print(f"An error occurred during skewnorm fitting: {e}")
    print(f"het-count: {number_of_hets}, read-depth per het: {average_read_depth_per_het}")
    # Define the negative log-likelihood function
    def neg_log_likelihood(params):
        a, loc, scale = params
        return -np.sum(np.log(stats.skewnorm.pdf(values, a, loc, scale)))
    
    # Initial parameter guesses
    initial_guess = [0, np.mean(values), np.std(values)]
    
    # Boundaries for the parameters
    param_bounds = [(-np.inf, np.inf), (np.min(values), np.max(values)), (1e-7, np.inf)]
    
    # Minimize the negative log-likelihood function
    result = minimize(neg_log_likelihood, initial_guess, bounds=param_bounds)
    
    if result.success:
        st_df, st_loc, st_scale = result.x
    else:
        raise ValueError(result.message)
    return st_df, st_loc, st_scale
        
def simulate_null_genes_helper(args):
    geneID, num_hets, total_count = args
    n_hets, total_readepth, st_df, st_loc, st_scale, mode_fitting_para, st_df_mean, st_loc_mean, st_scale_mean, mean_fitting_para = simulate_null_genes(num_hets, int(total_count / num_hets))
    return geneID, n_hets, total_readepth, st_df, st_loc, st_scale, mode_fitting_para, st_df_mean, st_loc_mean, st_scale_mean, mean_fitting_para

def run_null_simulations(input_genes, disable_cache):
    def cache_key(run):
        if disable_cache:
            return run[0]
        return f"h-{run[1]}_d-{run[2]}"
    
    null_simulation_cache = {}
    null_simulation_data = []
    runs_completed = 0
    simulation_start_t = time.time()
    with multiprocessing.Pool(12) as pool:
        all_runs = [(gene[0], gene[1], sum(gene[2:2+2*gene[1]])) for gene in input_genes]
        unique_runs = set([(cache_key(r), r[1], r[2]) for r in all_runs])
        logging.info(f"{datetime.now()} Starting {len(unique_runs)} unique gene simulation runs {len(all_runs)} total runs")
        for output in pool.imap_unordered(simulate_null_genes_helper, unique_runs):
            null_simulation_cache[output[0]] = output[1], output[2], output[3], output[4], output[5], output[6],output[7], output[8], output[9], output[10]

            runs_completed += 1
            if runs_completed % 100 == 0:
                genes_per_second = runs_completed / (time.time() - simulation_start_t)
                seconds_remaining = (len(unique_runs) - runs_completed) / genes_per_second
                logging.info(f"{datetime.now()} Completed {runs_completed} gene simulations ({round(genes_per_second, 2)} genes/s) [ETA: {round(seconds_remaining, 1)}s]")

    for r in all_runs:
        key = cache_key(r)
        n_hets, total_readepth, st_df, st_loc, st_scale, mode_fit_para, st_df_mean, st_loc_mean, st_scale_mean, mean_fit_para = null_simulation_cache[key]
        null_simulation_data.append((r[0], n_hets, total_readepth, st_df, st_loc, st_scale, mode_fit_para, st_df_mean, st_loc_mean, st_scale_mean, mean_fit_para))

    null_simulation_df = pd.DataFrame(null_simulation_data, columns=['geneID', 'n_hets', 'total_count', 'st_df' ,'st_loc' ,'st_scale', 'mode_st_parameter_fit', 'st_df_mean' ,'st_loc_mean' ,'st_scale_mean','mean_st_parameter_fit'])

    simulation_end_t = time.time()
    logging.info(f"Finished simulations in ${calculate_time(simulation_start_t, simulation_end_t)}")

    return null_simulation_df

def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("input_file_path", help="Input file path")
    argparser.add_argument("output_file_path", help="Output file path")
    argparser.add_argument("--disable-cache", help="Disable caching of null simulations", action="store_true", default=False)
    argparser.add_argument("--disable-null-simulation", help="Disable null simulations", action="store_true", default=False)
    argparser.add_argument("--fix-phasing-error", help="fix phasing error to 5%", action="store_true", default=False)
    args = argparser.parse_args()

    input_file_path = args.input_file_path
    output_file_path = args.output_file_path
    disable_cache = args.disable_cache
    disable_null_simulation = args.disable_null_simulation
    fix_phasing_error = args.fix_phasing_error
    start_t = time.time()

    # step1: run qB on input genes
    input_genes = []
    with open(input_file_path) as file:
        for line in file.readlines():
            fields = line.strip().split("\t")

            gene_id = format_geneID(fields[0])
            n_hets = int(fields[1])
            counts = map(int, fields[2:(n_hets*2+2)])

            unused = ''
            pis = []
            if n_hets > 1:
                #unused = fields[n_hets*2+2]
                #pis = map(float, fields[n_hets*2+2 + 1:])
                pis = map(float, fields[n_hets*2+2 + 0:])
                if len(fields) == n_hets*2+2:
                    # generate pis for gene input
                    pis = uniform(0.05, 0.05, n_hets-1)
            
            pis = map(lambda x: 0.05 if x == -1 else x, pis)

            if fix_phasing_error:
                pis = map(lambda x: 0.05, pis)

            gene = [
                gene_id,
                n_hets,
                *counts,
                *pis
            ]
            input_genes.append(gene)

    logging.info(f"Going to run quickbeast on {len(input_genes)} genes")
    gene_df = run_qb_parallel(input_genes)
    gene_df = gene_df.dropna(subset=['qb_mean'])
    qb_end_t = time.time()
    logging.info(f"Finished QB on input in ${calculate_time(start_t, qb_end_t)}")

    # step2: NULL simulation
    if not disable_null_simulation:
        null_simulation_df = run_null_simulations(input_genes, disable_cache)
        gene_df = pd.merge(gene_df, null_simulation_df, on="geneID")

        gene_df['mode_st_p_value'] = gene_df.apply(lambda row: calculate_st_2sided_pvalue(row['st_df'], row['st_loc'], row['st_scale'], row['qb_mode']), axis=1)
        gene_df['mean_st_p_value'] = gene_df.apply(lambda row: calculate_st_2sided_pvalue(row['st_df_mean'], row['st_loc_mean'], row['st_scale_mean'], row['qb_mean']), axis=1)

        columns_to_drop = ['st_df','st_loc','st_scale','st_df_mean','st_loc_mean','st_scale_mean']
        gene_df = gene_df.drop(columns=columns_to_drop)

    # Save the DataFrame as a TSV file
    gene_df.to_csv(output_file_path, sep='\t', index=False)

    complete_time_t = time.time()
    logging.info(f"Completed QB on input in ${calculate_time(start_t, complete_time_t)}")
    
# usage
if __name__ == "__main__":
    main()