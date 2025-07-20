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
        pl.col('max_protein_production').std().alias('std_max_protein_production'),
        pl.col('max_protein_production').min().alias('min_max_protein_production'),
        pl.col('max_protein_production').max().alias('max_max_protein_production'),
        pl.col('initial_chains_marking').mean().alias('avg_initial_chains_marking'),
        pl.col('total_simulation_steps').mean().alias('avg_total_simulation_steps'),
        pl.count().alias('number_of_runs')
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

def create_plots(results_df):
    """Create various plots from the analysis results."""
    if results_df.height == 0:
        print("No data to plot!")
        return
    
    print("\nCreating plots...")
    
    # 1. Heatmap of average max protein production vs initial ribosomes and excess factor
    pivot_df = results_df.pivot(
        values="avg_max_protein_production",
        index="initial_ribosomes", 
        on="excess_factor"
    )
    
    # Convert pivot to list format for heatmap
    heatmap_data = []
    for ribosomes in pivot_df.select("initial_ribosomes").to_series().to_list():
        row = []
        for col in pivot_df.columns[1:]:  # Skip initial_ribosomes column
            filtered_df = pivot_df.filter(pl.col("initial_ribosomes") == ribosomes).select(col)
            if filtered_df.height > 0:
                value = filtered_df.item()
            else:
                value = None  # or 0, depending on how you want to handle missing data
            row.append(value)
        heatmap_data.append(row)
    
    # Debug: Print some data to verify
    print(f"Heatmap data shape: {len(heatmap_data)} rows")
    if heatmap_data:
        print(f"First row length: {len(heatmap_data[0])}")
        print(f"Sample values: {heatmap_data[0][:3]}...")
    
    fig_heatmap = px.imshow(
        heatmap_data,
        title="Average Protein Production Heatmap - Ribosome Availability Experiment",
        labels=dict(x="Excess Factor", y="Initial Ribosomes", color="Avg Max Protein Production"),
        aspect="auto"
    )
    fig_heatmap.write_html("data/chains_marking_protein_production_heatmap.html")
    print("✓ Heatmap saved to data/chains_marking_protein_production_heatmap.html")
    
    # 2. 3D Surface plot - properly structured data
    # Create a pivot table for 3D surface
    pivot_3d = results_df.pivot(
        values="avg_max_protein_production",
        index="initial_ribosomes", 
        on="excess_factor"
    )
    
    # Extract data for 3D surface
    # Get unique values for x and y axes
    x_values = sorted(results_df.select("excess_factor").unique().to_series().to_list())
    y_values = sorted(results_df.select("initial_ribosomes").unique().to_series().to_list())
    z_values = []
    
    for ribosomes in y_values:
        row = []
        for excess in x_values:
            # Find the value for this ribosomes/excess combination
            filtered_df = results_df.filter(
                (pl.col("initial_ribosomes") == ribosomes) & 
                (pl.col("excess_factor") == excess)
            ).select("avg_max_protein_production")
            if filtered_df.height > 0:
                value = filtered_df.item()
            else:
                value = None  # or 0, depending on how you want to handle missing data
            row.append(value)
        z_values.append(row)
    
    fig_3d = go.Figure(data=[go.Surface(
        x=x_values,
        y=y_values,
        z=z_values,
        colorscale='Viridis'
    )])
    fig_3d.update_layout(
        title="3D Surface: Average Protein Production vs Ribosome Availability vs Excess Factor",
        scene=dict(
            xaxis_title="Excess Factor",
            yaxis_title="Initial Ribosomes", 
            zaxis_title="Avg Max Protein Production"
        )
    )
    fig_3d.write_html("data/chains_marking_protein_production_3d_surface.html")
    print("✓ 3D Surface plot saved to data/chains_marking_protein_production_3d_surface.html")
    
    # 3. Line plot: Average protein production vs excess factor for different ribosome counts
    fig_lines = px.line(
        results_df,
        x="excess_factor",
        y="avg_max_protein_production",
        color="initial_ribosomes",
        title="Average Protein Production vs Excess Factor by Ribosome Count",
        labels=dict(x="Excess Factor", y="Avg Max Protein Production", color="Initial Ribosomes")
    )
    fig_lines.write_html("data/chains_marking_protein_production_vs_excess_factor.html")
    print("✓ Line plot saved to data/chains_marking_protein_production_vs_excess_factor.html")
    
    # 4. Scatter plot with size and color based on percentage of max protein production
    fig_scatter = px.scatter(
        results_df,
        x="initial_ribosomes",
        y="excess_factor",
        size="percentage_max_protein_production",
        color="percentage_max_protein_production",
        title="Ribosome Availability Parameter Space Analysis - % of Max Protein Production",
        labels=dict(x="Initial Ribosomes", y="Excess Factor", size="% of Max Protein Production", color="% of Max Protein Production")
    )
    
    # Set log scale for x-axis (initial ribosomes)
    fig_scatter.update_xaxes(type="log")
    fig_scatter.write_html("data/chains_marking_parameter_space_analysis.html")
    print("✓ Scatter plot saved to data/chains_marking_parameter_space_analysis.html")
    
    # 5. Line plot: Average protein production vs ribosome availability for different ribosome counts
    fig_chains_lines = px.line(
        results_df,
        x="excess_factor",
        y="avg_max_protein_production",
        color="initial_ribosomes",
        title="Average Protein Production vs Excess Factor by Ribosome Count",
        labels=dict(x="Excess Factor", y="Avg Max Protein Production", color="Initial Ribosomes")
    )
    
    fig_chains_lines.write_html("data/chains_marking_protein_production_vs_excess_factor_by_ribosomes.html")
    print("✓ Protein production vs excess factor by ribosomes line plot saved to data/chains_marking_protein_production_vs_excess_factor_by_ribosomes.html")
    
    # 6. Line plot: Percentage of max protein production vs excess factor for different ribosome counts
    fig_percentage_lines = px.line(
        results_df,
        x="excess_factor",
        y="percentage_max_protein_production",
        color="initial_ribosomes",
        title="Percentage of Max Protein Production vs Excess Factor by Ribosome Count",
        labels=dict(x="Excess Factor", y="% of Max Protein Production", color="Initial Ribosomes")
    )
    
    fig_percentage_lines.write_html("data/chains_marking_percentage_protein_production_vs_excess_factor.html")
    print("✓ Percentage protein production vs excess factor line plot saved to data/chains_marking_percentage_protein_production_vs_excess_factor.html")
    
    print(f"\nAll plots saved to the data/ folder!")

if __name__ == "__main__":
    results = analyze_simulation_results()
    if results is not None:
        create_plots(results) 