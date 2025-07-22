#!/usr/bin/env python3
"""
Script to create plots from chains marking analysis results.
This script reads the analysis results and creates various visualizations.
"""

import os
import math
from re import template
import polars as pl
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def create_plots(results_df):
    """Create various plots from the analysis results."""
    if results_df.height == 0:
        print("No data to plot!")
        return
    
    # Filter data for excess factor <= 4.0, initial ribosomes <= 1000, and exclude specific ribosome values
    filtered_df = results_df.filter(
        (pl.col("excess_factor") <= 4.0) & 
        (pl.col("initial_ribosomes") <= 1000) &
        ~pl.col("excess_factor").is_in([0.2]) &
        ~pl.col("initial_ribosomes").is_in([60,80, 120, 140, 160, 170, 180, 190])
    )
    print(f"Filtered data: {filtered_df.height} parameter combinations (excess factor <= 4.0, ribosomes <= 1000, excluding 110,120,140,160,170)")
    
    if filtered_df.height == 0:
        print("No data to plot after filtering!")
        return
    
    # Debug: Check columns in filtered_df
    print(f"Filtered DataFrame columns: {filtered_df.columns}")
    print(f"Original DataFrame columns: {results_df.columns}")
    
    # Check if required columns exist
    required_columns = ["initial_ribosomes", "excess_factor", "avg_max_protein_production", "percentage_max_protein_production"]
    missing_columns = [col for col in required_columns if col not in filtered_df.columns]
    
    if missing_columns:
        print(f"Error: Missing columns in filtered DataFrame: {missing_columns}")
        print("Using original DataFrame instead...")
        filtered_df = results_df  # Fallback to original data
    
    # Ensure we have the right data types and columns for pivot
    try:
        # Test if pivot works with filtered data
        test_pivot = filtered_df.pivot(
            values="avg_max_protein_production",
            index="initial_ribosomes", 
            on="excess_factor"
        )
        print("✓ Pivot test successful with filtered data")
    except Exception as e:
        print(f"Pivot test failed with filtered data: {e}")
        print("Using original DataFrame for all plots...")
        filtered_df = results_df
    
    print("\nCreating plots...")
    
    # 1. Heatmap of median max protein production vs initial ribosomes and excess factor
    pivot_df = filtered_df.pivot(
        values="median_max_protein_production",
        index="initial_ribosomes", 
        on="excess_factor"
    )
    
    # Convert pivot to list format for heatmap
    heatmap_data = []
    for ribosomes in pivot_df.select("initial_ribosomes").to_series().to_list():
        row = []
        for col in pivot_df.columns[1:]:  # Skip initial_ribosomes column
            filtered_row = pivot_df.filter(pl.col("initial_ribosomes") == ribosomes).select(col)
            if filtered_row.height > 0:
                value = filtered_row.item()
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
        title="Median Protein Production Heatmap - Ribosome Availability Experiment",
        labels=dict(x="Excess Factor", y="Initial Ribosomes", color="Median Max Protein Production"),
        aspect="auto"
    )
    fig_heatmap.write_html("data/chains_marking_protein_production_heatmap.html")
    print("✓ Heatmap saved to data/chains_marking_protein_production_heatmap.html")
    
    # 2. 3D Surface plot - Improved version with better interpolation
    # Create pivot table for surface plot
    pivot_df = filtered_df.pivot(
        values="median_max_protein_production",
        index="initial_ribosomes", 
        on="excess_factor"
    )
    
    # Convert pivot to list format for surface plot with better handling of missing data
    surface_data = []
    for ribosomes in pivot_df.select("initial_ribosomes").to_series().to_list():
        row = []
        for col in pivot_df.columns[1:]:  # Skip initial_ribosomes column
            filtered_row = pivot_df.filter(pl.col("initial_ribosomes") == ribosomes).select(col)
            if filtered_row.height > 0:
                value = filtered_row.item()
                if value is None:
                    value = 0  # Replace None with 0
            else:
                value = 0  # Use 0 for missing values in surface plot
            row.append(value)
        surface_data.append(row)
    
    # Get x and y coordinates for surface plot
    x_coords = [float(col) for col in pivot_df.columns[1:]]
    y_coords = pivot_df.select("initial_ribosomes").to_series().to_list()
    
    # Debug: Print some info about the data
    print(f"Surface data shape: {len(surface_data)} x {len(surface_data[0]) if surface_data else 0}")
    print(f"X coordinates range: {min(x_coords):.2f} to {max(x_coords):.2f}")
    print(f"Y coordinates range: {min(y_coords)} to {max(y_coords)}")
    
    # Create improved surface plot with better settings
    fig_3d = go.Figure(data=[go.Surface(
        x=x_coords,
        y=y_coords,
        z=surface_data,
        colorscale='Viridis',
        opacity=0.9,
        connectgaps=True,  # Don't connect gaps to avoid artifacts
        showscale=True,
        hoverinfo='z'
    )])
    
    fig_3d.update_layout(
        title="3D Surface: Median Protein Production vs Ribosome Availability vs Excess Factor",
        scene=dict(
            xaxis_title="Excess Factor",
            yaxis_title="Initial Ribosomes", 
            zaxis_title="Median Max Protein Production",
            yaxis=dict(type="log"),
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5)  # Better viewing angle
            )
        ),
        width=800,
        height=600,
        template='none'
    )
    fig_3d.write_html("data/chains_marking_protein_production_3d_surface.html")
    print("✓ 3D Surface plot saved to data/chains_marking_protein_production_3d_surface.html")
    
    # 2b. Alternative: 3D Scatter plot for better visualization of sparse data
    fig_3d_scatter = go.Figure(data=[go.Scatter3d(
        x=filtered_df.select("excess_factor").to_series().to_list(),
        y=filtered_df.select("initial_ribosomes").to_series().to_list(),
        z=filtered_df.select("median_max_protein_production").to_series().to_list(),
        mode='markers',
        marker=dict(
            size=8,
            color=filtered_df.select("median_max_protein_production").to_series().to_list(),
            colorscale='Viridis',
            opacity=1,
            showscale=True
        ),
        text=[f"Ribosomes: {r}, Excess: {e:.2f}, Protein: {p:.1f}" 
              for r, e, p in zip(
                  filtered_df.select("initial_ribosomes").to_series().to_list(),
                  filtered_df.select("excess_factor").to_series().to_list(),
                  filtered_df.select("median_max_protein_production").to_series().to_list()
              )],
        hovertemplate='<b>%{text}</b><extra></extra>'
    )])
    
    fig_3d_scatter.update_layout(
        title="3D Scatter: Median Protein Production vs Ribosome Availability vs Excess Factor",
        scene=dict(
            xaxis_title="Excess Factor",
            yaxis_title="Initial Ribosomes", 
            zaxis_title="Median Max Protein Production",
            yaxis=dict(type="log"),
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5)
            )
        ),
        width=800,
        height=600
    )
    fig_3d_scatter.write_html("data/chains_marking_protein_production_3d_scatter.html")
    print("✓ 3D Scatter plot saved to data/chains_marking_protein_production_3d_scatter.html")
    
    
    
    # 3. Line plot: Median protein production vs excess factor for different ribosome counts
    fig_lines = px.line(
        filtered_df,
        x="excess_factor",
        y="median_max_protein_production",
        color="initial_ribosomes",
        title="Median Protein Production vs Excess Factor by Ribosome Count",
        labels=dict(x="Excess Factor", y="Median Max Protein Production", color="Initial Ribosomes")
    )
    # Set log scale for color axis (initial ribosomes)
    fig_lines.update_layout(coloraxis=dict(colorscale='Viridis'))
    fig_lines.write_html("data/chains_marking_protein_production_vs_excess_factor.html")
    print("✓ Line plot saved to data/chains_marking_protein_production_vs_excess_factor.html")
    
    # 4. Scatter plot with size and color based on percentage of max protein production
    fig_scatter = px.scatter(
        filtered_df,
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
    
    # 5. Line plot: Median protein production vs ribosome availability for different ribosome counts
    fig_chains_lines = px.line(
        filtered_df,
        x="excess_factor",
        y="median_max_protein_production",
        color="initial_ribosomes",
        title="Median Protein Production vs Excess Factor by Ribosome Count",
        labels=dict(x="Excess Factor", y="Median Max Protein Production", color="Initial Ribosomes")
    )
    
    # Set log scale for color axis (initial ribosomes)
    fig_chains_lines.update_layout(coloraxis=dict(colorscale='Viridis'))
    fig_chains_lines.write_html("data/chains_marking_protein_production_vs_excess_factor_by_ribosomes.html")
    print("✓ Protein production vs excess factor by ribosomes line plot saved to data/chains_marking_protein_production_vs_excess_factor_by_ribosomes.html")
    
    # 6. Line plot: Percentage of max protein production vs excess factor for different ribosome counts
    # Create gradient colors manually
    ribosome_values = sorted(filtered_df.select("initial_ribosomes").unique().to_series().to_list())
    min_ribosomes = min(ribosome_values)
    max_ribosomes = max(ribosome_values)
    
    fig_percentage_lines = go.Figure()
    
    for ribosomes in ribosome_values:
        # Filter data for this ribosome count
        ribosome_data = filtered_df.filter(pl.col("initial_ribosomes") == ribosomes)
        
        # Calculate color based on logarithmic position in range
        if min_ribosomes > 0:
            log_ratio = (math.log(ribosomes) - math.log(min_ribosomes)) / (math.log(max_ribosomes) - math.log(min_ribosomes))
        else:
            log_ratio = (ribosomes - min_ribosomes) / (max_ribosomes - min_ribosomes)
        
        # Create trace for this ribosome count
        fig_percentage_lines.add_trace(go.Scatter(
            x=ribosome_data.select("excess_factor").to_series().to_list(),
            y=ribosome_data.select("percentage_max_protein_production").to_series().to_list(),
            mode='lines+markers',
            name=f'{ribosomes} ribosomes',
            line=dict(color=f'rgb({int(255 * log_ratio)}, {int(100 + 155 * log_ratio)}, {int(255 - 255 * log_ratio)})'),
            marker=dict(size=6)
        ))
    
    fig_percentage_lines.update_layout(
        title="Percentage of Max Protein Production vs Excess Factor by Ribosome Count",
        xaxis_title="Excess Factor",
        yaxis_title="% of Max Protein Production",
        showlegend=True
    )
    fig_percentage_lines.write_html("data/chains_marking_percentage_protein_production_vs_excess_factor.html")
    print("✓ Percentage protein production vs excess factor line plot saved to data/chains_marking_percentage_protein_production_vs_excess_factor.html")
    
    print(f"\nAll plots saved to the data/ folder!")

def load_analysis_results():
    """Load the analysis results from the parquet file."""
    analysis_file = os.path.join("data", "chains_marking_analysis_results.parquet")
    
    if not os.path.exists(analysis_file):
        print(f"Error: Analysis results file not found: {analysis_file}")
        print("Please run the analysis script first: python scripts/run_chains_marking_analysis.py")
        return None
    
    try:
        results_df = pl.read_parquet(analysis_file)
        print(f"Loaded analysis results: {results_df.height} parameter combinations")
        return results_df
    except Exception as e:
        print(f"Error loading analysis results: {e}")
        return None

def main():
    """Main function to load results and create plots."""
    print("=== Chains Marking Plotting Script ===")
    
    # Load analysis results
    results_df = load_analysis_results()
    
    if results_df is not None:
        # Create plots
        create_plots(results_df)
        print("\n=== Plotting Complete ===")
    else:
        print("\n=== Plotting Failed ===")

if __name__ == "__main__":
    main() 