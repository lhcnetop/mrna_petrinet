import os
import polars as pl
import re
import glob
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def extract_parameters_from_filename(filename):
    """Extract initial ribosomes and excess factor from filename."""
    # Pattern to match both formats:
    # 1. chains_marking_simulation_history_{ribosomes}_ribosomes_{excess}_excess_mrna.parquet (old format)
    # 2. chains_marking_simulation_history_{ribosomes}_ribosomes_{excess}_excess_mrna_{timestamp}.parquet (new format)
    pattern = r'chains_marking_simulation_history_(\d+)_ribosomes_([\d.]+)_excess_mrna(?:_\d{8}_\d{6})?\.parquet'
    match = re.search(pattern, filename)
    
    if match:
        initial_ribosomes = int(match.group(1))
        excess_factor = float(match.group(2))
        return initial_ribosomes, excess_factor
    else:
        return None, None

def analyze_simulation_results():
    """Analyze all simulation results in the chains marking experiment folder."""
    experiment_folder = "data/chains_marking_experiment"
    
    if not os.path.exists(experiment_folder):
        print(f"Experiment folder {experiment_folder} not found!")
        return
    
    # Get all parquet files in the experiment folder
    parquet_files = glob.glob(os.path.join(experiment_folder, "*.parquet"))
    
    if not parquet_files:
        print(f"No parquet files found in {experiment_folder}")
        return
    
    print(f"Found {len(parquet_files)} simulation files to analyze...")
    
    results = []
    
    for file_path in parquet_files:
        filename = os.path.basename(file_path)
        initial_ribosomes, excess_factor = extract_parameters_from_filename(filename)
        
        if initial_ribosomes is None or excess_factor is None:
            print(f"Could not extract parameters from filename: {filename}")
            continue
        
        try:
            # Read the parquet file
            df = pl.read_parquet(file_path)
            
            # Get the maximum value of p_preinsulin column
            if 'p_preinsulin' in df.columns:
                max_protein = df.select('p_preinsulin').max().item()
            else:
                max_protein = None
                print(f"Warning: 'p_preinsulin' column not found in {filename}")
            
            # Get simulation metadata
            total_steps = df.height
            final_step = df.tail(1) if df.height > 0 else None
            
            # Extract initial chains marking from the simulation data
            # Get the value of p_chainA_0 at step 0 (first row)
            initial_chains_marking = None
            if df.height > 0 and 'p_chainA_0' in df.columns:
                first_row = df.head(1)
                initial_chains_marking = first_row.select('p_chainA_0').item()
            
            results.append({
                'filename': filename,
                'initial_ribosomes': initial_ribosomes,
                'excess_factor': excess_factor,
                'initial_chains_marking': initial_chains_marking,
                'max_protein_production': max_protein,
                'total_simulation_steps': total_steps,
                'final_protein_count': final_step.select('p_preinsulin').item() if final_step is not None and 'p_preinsulin' in final_step.columns else None
            })
            
            print(f"Analyzed {filename}: {initial_ribosomes} initial ribosomes, {excess_factor:.2f} excess factor, max protein: {max_protein}")
            
        except Exception as e:
            print(f"Error analyzing {filename}: {e}")
    
    # Create results DataFrame
    results_df = pl.DataFrame(results)
    
    if results_df.height == 0:
        print("No valid results found!")
        return None
    
    # Group by parameters and calculate averages
    print(f"\nCalculating averages for {results_df.height} simulation runs...")
    
    # Group by initial_ribosomes and excess_factor, then calculate averages
    averaged_results = results_df.group_by(['initial_ribosomes', 'excess_factor']).agg([
        pl.col('max_protein_production').mean().alias('avg_max_protein_production'),
        pl.col('max_protein_production').median().alias('median_max_protein_production'),
        pl.col('max_protein_production').std().alias('std_max_protein_production'),
        pl.col('max_protein_production').min().alias('min_max_protein_production'),
        pl.col('max_protein_production').max().alias('max_max_protein_production'),
        pl.col('initial_chains_marking').mean().alias('avg_initial_chains_marking'),
        pl.col('total_simulation_steps').mean().alias('avg_total_simulation_steps'),
        pl.len().alias('number_of_runs')
    ])
    
    # Add percentage of max output protein column
    # Max possible protein = min(100, initial_chains_marking)
    averaged_results = averaged_results.with_columns([
        pl.min_horizontal(pl.col('avg_initial_chains_marking'), 100).alias('max_possible_protein'),
        (pl.col('avg_max_protein_production') / pl.min_horizontal(pl.col('avg_initial_chains_marking'), 100) * 100).alias('percentage_max_protein_production')
    ])
    
    # Sort by initial ribosomes and excess factor
    averaged_results = averaged_results.sort(['initial_ribosomes', 'excess_factor'])
    
    # Save both detailed and averaged results
    detailed_analysis_file = os.path.join("data", "chains_marking_detailed_analysis_results.parquet")
    averaged_analysis_file = os.path.join("data", "chains_marking_analysis_results.parquet")
    
    results_df.write_parquet(detailed_analysis_file)
    averaged_results.write_parquet(averaged_analysis_file)
    
    print(f"Detailed results saved to: {detailed_analysis_file}")
    print(f"Averaged results saved to: {averaged_analysis_file}")
    print(f"Found {averaged_results.height} unique parameter combinations")
    
    return averaged_results

if __name__ == "__main__":
    print("=== Starting Chains Marking Analysis ===")
    results = analyze_simulation_results()
    if results is not None:
        print(f"\n=== Analysis Complete ===")
        print(f"Found {results.height} unique parameter combinations")
        print(f"Results saved to: data/chains_marking_analysis_results.parquet")
        print(f"To create plots, run: python scripts/run_chains_marking_plots.py")
    else:
        print("\n=== Analysis Failed ===") 