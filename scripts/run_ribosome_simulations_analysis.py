import os
import polars as pl
import re
import glob
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def extract_parameters_from_filename(filename):
    """Extract ribosome count and excess factor from filename."""
    # Pattern to match: ribosome_simulation_history_{ribosomes}_ribosomes_{excess}_excess.parquet
    pattern = r'ribosome_simulation_history_(\d+)_ribosomes_([\d.]+)_excess\.parquet'
    match = re.search(pattern, filename)
    
    if match:
        ribosomes = int(match.group(1))
        excess_factor = float(match.group(2))
        return ribosomes, excess_factor
    else:
        return None, None

def analyze_simulation_results():
    """Analyze all simulation results in the experiment folder."""
    experiment_folder = "data/ribosome_experiment"
    
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
        ribosomes, excess_factor = extract_parameters_from_filename(filename)
        
        if ribosomes is None or excess_factor is None:
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
            
            results.append({
                'filename': filename,
                'ribosomes': ribosomes,
                'excess_factor': excess_factor,
                'max_protein_production': max_protein,
                'total_simulation_steps': total_steps,
                'final_protein_count': final_step.select('p_preinsulin').item() if final_step is not None and 'p_preinsulin' in final_step.columns else None
            })
            
            print(f"Analyzed {filename}: {ribosomes} ribosomes, {excess_factor:.2f} excess factor, max protein: {max_protein}")
            
        except Exception as e:
            print(f"Error analyzing {filename}: {e}")
    
    # Create results DataFrame
    results_df = pl.DataFrame(results)
    
    # Sort by ribosomes and excess factor
    results_df = results_df.sort(['ribosomes', 'excess_factor'])
    
    # Save analysis results to Parquet file in the root of the data folder
    analysis_file = os.path.join("data", "simulation_analysis_results.parquet")
    results_df.write_parquet(analysis_file)
    
    print(f"\nAnalysis complete! Results saved to: {analysis_file}")
        
    return results_df

def create_plots(results_df):
    """Create various plots from the analysis results."""
    if results_df.height == 0:
        print("No data to plot!")
        return
    
    print("\nCreating plots...")
    
    # 1. Heatmap of max protein production vs ribosomes and excess factor
    pivot_df = results_df.pivot(
        values="max_protein_production",
        index="ribosomes", 
        on="excess_factor"
    )
    
    # Convert pivot to list format for heatmap
    heatmap_data = []
    for ribosome in pivot_df.select("ribosomes").to_series().to_list():
        row = [ribosome]
        for col in pivot_df.columns[1:]:  # Skip ribosomes column
            value = pivot_df.filter(pl.col("ribosomes") == ribosome).select(col).item()
            row.append(value)
        heatmap_data.append(row)
    
    fig_heatmap = px.imshow(
        heatmap_data,
        title="Protein Production Heatmap",
        labels=dict(x="Excess Factor", y="Initial Ribosomes", color="Max Protein Production"),
        aspect="auto"
    )
    fig_heatmap.write_html("data/protein_production_heatmap.html")
    print("✓ Heatmap saved to data/protein_production_heatmap.html")
    
    # 2. 3D Surface plot - properly structured data
    # Create a pivot table for 3D surface
    pivot_3d = results_df.pivot(
        values="max_protein_production",
        index="ribosomes", 
        on="excess_factor"
    )
    
    # Extract data for 3D surface
    # Get unique values for x and y axes
    x_values = sorted(results_df.select("excess_factor").unique().to_series().to_list())
    y_values = sorted(results_df.select("ribosomes").unique().to_series().to_list())
    z_values = []
    
    for ribosome in y_values:
        row = []
        for excess in x_values:
            # Find the value for this ribosome/excess combination
            value = results_df.filter(
                (pl.col("ribosomes") == ribosome) & 
                (pl.col("excess_factor") == excess)
            ).select("max_protein_production").item()
            row.append(value)
        z_values.append(row)
    
    fig_3d = go.Figure(data=[go.Surface(
        x=x_values,
        y=y_values,
        z=z_values,
        colorscale='Viridis'
    )])
    fig_3d.update_layout(
        title="3D Surface: Protein Production vs Ribosomes vs Excess Factor",
        scene=dict(
            xaxis_title="Excess Factor",
            yaxis_title="Initial Ribosomes", 
            zaxis_title="Max Protein Production"
        )
    )
    fig_3d.write_html("data/protein_production_3d_surface.html")
    print("✓ 3D Surface plot saved to data/protein_production_3d_surface.html")
    
    # 3. Line plot: Protein production vs excess factor for different ribosome counts
    fig_lines = px.line(
        results_df,
        x="excess_factor",
        y="max_protein_production",
        color="ribosomes",
        title="Protein Production vs Excess Factor by Ribosome Count",
        labels=dict(x="Excess Factor", y="Max Protein Production", color="Initial Ribosomes")
    )
    fig_lines.write_html("data/protein_production_vs_excess_factor.html")
    print("✓ Line plot saved to data/protein_production_vs_excess_factor.html")
    
    # 4. Scatter plot with size based on protein production
    fig_scatter = px.scatter(
        results_df,
        x="ribosomes",
        y="excess_factor",
        size="max_protein_production",
        color="max_protein_production",
        title="Parameter Space Analysis",
        labels=dict(x="Initial Ribosomes", y="Excess Factor", size="Max Protein Production", color="Max Protein Production")
    )
    fig_scatter.write_html("data/parameter_space_analysis.html")
    print("✓ Scatter plot saved to data/parameter_space_analysis.html")
    
    print(f"\nAll plots saved to the data/ folder!")

if __name__ == "__main__":
    results = analyze_simulation_results()
    if results is not None:
        create_plots(results) 