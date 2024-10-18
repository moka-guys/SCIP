import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from jinja2 import Template
import plotly.io as pio
import numpy as np

# Assume FL_SNPs, mean_pat, median_pat, IQR_Pat, mean_Fet, median_Fet, IQR_Fet, Total_count, Alt_count are calculated

# Create HTML header template
html_header_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{ report_name }}</title>
    <style>
        body { font-family: Arial, sans-serif; }
        h1 { color: #333; }
        .summary { margin-bottom: 20px; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        img { display: block; margin: 0 auto; }
    </style>
</head>
<body>
    <h1>{{ report_name }}</h1>
    <h2>Summary</h2>
    </body>
</html>
"""

# Create a basic HTML template using jinja2
html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{ report_name }}</title>
    <style>
        body { font-family: Arial, sans-serif; }
        h1 { color: #333; }
        .summary { margin-bottom: 20px; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        img { display: block; margin: 0 auto; }
    </style>
</head>
<body>
    <h2>{{ report_section_name }}</h2>

    <h3>Summary Statistics</h3>
    <div class="summary">
        <p><strong>Mean Paternal Fraction:</strong> {{ mean_pat }}</p>
        <p><strong>Median Paternal Fraction:</strong> {{ median_pat }}</p>
        <p><strong>IQR Paternal Fraction:</strong> {{ IQR_Pat }}</p>
        <p><strong>Mean Fetal Fraction:</strong> {{ mean_Fet }}</p>
        <p><strong>Median Fetal Fraction:</strong> {{ median_Fet }}</p>
        <p><strong>IQR Fetal Fraction:</strong> {{ IQR_Fet }}</p>
        <p><strong>Total Count:</strong> {{ Total_count }}</p>
        <p><strong>Alternative Count:</strong> {{ Alt_count }}</p>
        <p><strong>Alternative Ratio:</strong> {{ Alt_ratio }}</p>
    </div>

    <h3>FL_SNPs Table</h3>
    <div>{{ table_html }}</div>

    <h3>Graphs</h3>
    <div>{{ sprt_plot|safe }}</div>
    <div>{{ chr11_plot|safe }}</div>
</body>
</html>
"""

# Create a basic HTML template using jinja2 - for the summary of prediction results
pred_sum_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{ report_name }}</title>
    <style>
        body { font-family: Arial, sans-serif; }
        h1 { color: #333; }
        .summary { margin-bottom: 20px; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        img { display: block; margin: 0 auto; }
    </style>
</head>
<body>
    <div class="Summary">
        <p><strong>Predicted fetal genotype for the {{ allele }}:</strong> {{ prediction }}</p>
    </div>
    </body>
</html>
"""

# Generate Plotly graphs directly as HTML
def generate_sprt_plot(Total_count, Alt_count, d, g, d_wt, g_wt):
    
    # Define the x values
    x_vals = np.arange(50, 120001)

    # Calculate limits
    Upper_limit_graph_rmd = ((np.log(8) / x_vals) - np.log(d)) / np.log(g)
    Lower_limit_graph_rmd = ((np.log(1 / 8) / x_vals) - np.log(d)) / np.log(g)
    Upper_limit_graph_wt = ((np.log(8) / x_vals) - np.log(d_wt)) / np.log(g_wt)
    Lower_limit_graph_wt = ((np.log(1 / 8) / x_vals) - np.log(d_wt)) / np.log(g_wt)

    # Calculate the y-axis range
    g_range = (
        min(Upper_limit_graph_rmd.min(), Lower_limit_graph_rmd.min(),
            Upper_limit_graph_wt.min(), Lower_limit_graph_wt.min(),
            Alt_count / Total_count),
        max(Upper_limit_graph_rmd.max(), Lower_limit_graph_rmd.max(),
            Upper_limit_graph_wt.max(), Lower_limit_graph_wt.max(),
            Alt_count / Total_count)
    )

    # Create a Plotly figure
    fig = go.Figure()

    # Add traces for each line
    fig.add_trace(go.Scatter(x=x_vals, y=Upper_limit_graph_rmd, mode='lines', 
                             name='Upper Limit RMD', line=dict(color='red')))
    fig.add_trace(go.Scatter(x=x_vals, y=Lower_limit_graph_rmd, mode='lines', 
                             name='Lower Limit RMD', line=dict(color='pink')))
    fig.add_trace(go.Scatter(x=x_vals, y=Upper_limit_graph_wt, mode='lines', 
                             name='Upper Limit WT', line=dict(color='green')))
    fig.add_trace(go.Scatter(x=x_vals, y=Lower_limit_graph_wt, mode='lines', 
                             name='Lower Limit WT', line=dict(color='blue')))

    # Calculate the observed value
    observed_y = Alt_count / Total_count
    fig.add_trace(go.Scatter(x=[Total_count], y=[observed_y], mode='markers', 
                             name='Observed', marker=dict(color='blue', size=10)))

    # Update layout with title, axis labels, and range
    fig.update_layout(
        title="Modified SPRT",
        xaxis_title="Total number of counts",
        yaxis_title="Pr over-represented allele",
        xaxis=dict(range=[min(x_vals), max(x_vals)]),
        yaxis=dict(range=g_range),
        legend=dict(title="Legend")
    )
    # Return the HTML representation of the plot
    return pio.to_html(fig, full_html=False)



def generate_chr11_plot(median_pat, Alt_count,Total_count): #TODO how is it deciding what allele to show?
    # Define the x and y ranges
    xrange = np.arange(5225264, 5227272)
    yrange = [0, 1.2]

    # Define the HBB exon ranges
    HBB_exon_1 = np.arange(5227071, 5226929, -1)
    HBB_exon_2 = np.arange(5226799, 5226576, -1)
    HBB_exon_3 = np.arange(5225726, 5225463, -1)

    # Define the fetal variants
    Fetal_Het_alt = [((100 - median_pat) / 100)] * len(xrange)
    Fetal_Hom_alt = [(0.5 + (median_pat / 100))] * len(xrange)
    Fetal_Hom_ref = [(0.5 - (median_pat / 100))] * len(xrange)
    Fetal_Het_ref = [(median_pat / 100)] * len(xrange)
    Variant_of_interest = [(Alt_count / Total_count)]

    # Create a Plotly figure
    fig = go.Figure()

    # Add traces for each line
    fig.add_trace(go.Scatter(x=xrange, y=[0.5] * len(xrange), mode='lines', 
                             line=dict(color='black'), name='Chr11'))
    fig.add_trace(go.Scatter(x=HBB_exon_1, y=[0.5] * len(HBB_exon_1), mode='lines', 
                             line=dict(color='purple', width=10), name='HBB exon 1'))
    fig.add_trace(go.Scatter(x=HBB_exon_2, y=[0.5] * len(HBB_exon_2), mode='lines', 
                             line=dict(color='green', width=10), name='HBB exon 2'))
    fig.add_trace(go.Scatter(x=HBB_exon_3, y=[0.5] * len(HBB_exon_3), mode='lines', 
                             line=dict(color='orange', width=10), name='HBB exon 3'))
    fig.add_trace(go.Scatter(x=xrange, y=Fetal_Het_alt, mode='lines', 
                             line=dict(color='cadetblue', width=4), name='Fetal Het Alt'))
    fig.add_trace(go.Scatter(x=[5226925], y=Variant_of_interest, mode='markers', 
                             marker=dict(color='red', size=10), name='Variant of Interest'))

    # Update layout with title and axis labels
    fig.update_layout(
        title="Chromosome 11 Plot",
        xaxis_title="Position",
        yaxis_title="Fraction",
        yaxis=dict(range=yrange),
        legend=dict(title="Legend")
    )
    return pio.to_html(fig, full_html=False)

def generate_html_content(mean_pat, median_pat, IQR_Pat, mean_Fet, \
                          median_Fet, IQR_Fet, Total_count, Alt_count,FL_SNPs, allele,\
                            report_name, d, g, d_wt, g_wt):
    # Generate the graphs
    sprt_plot_html = generate_sprt_plot(Total_count,Alt_count,d, g, d_wt, g_wt)
    chr11_plot_html = generate_chr11_plot(median_pat, Alt_count, Total_count)

    # Save FL_SNPs DataFrame as HTML table
    table_html = FL_SNPs.to_html(index=False)

    # Fill in the HTML template
    html_content = Template(html_template).render(
        mean_pat=mean_pat,
        median_pat=median_pat,
        IQR_Pat=IQR_Pat,
        mean_Fet=mean_Fet,
        median_Fet=median_Fet,
        IQR_Fet=IQR_Fet,
        Total_count=Total_count,
        Alt_count=Alt_count,
        Alt_ratio=Alt_count / Total_count,
        table_html=table_html,
        sprt_plot=sprt_plot_html,
        chr11_plot=chr11_plot_html,
        report_section_name = "Report for the " + allele,
        report_name = report_name
    )

    return html_content

def generate_summary_html_content(report_name,allele, prediction):
    html_summary_content = Template(pred_sum_template).render(
        report_name = report_name,
        allele = allele,
        prediction = prediction
    )

    return html_summary_content

def generate_html_header(report_name):
    html_header = Template(html_header_template).render(
        report_name = report_name
    )
    return html_header


# # Write to an HTML file
# with open('FL_SNPs_report.html', 'w') as f:
#     f.write(html_content)

# print("Interactive report generated successfully: FL_SNPs_report.html")
