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
def generate_sprt_plot(Total_count, Alt_count):
    x_vals = list(range(50, 120001)) # TODO softcode these values? how calculated in fetal_get_pred.py?
    Upper_limit_graph_rmd = ((np.log(8) / x_vals) - np.log(0.5)) / np.log(0.8)
    Lower_limit_graph_rmd = ((np.log(1/8) / x_vals) - np.log(0.5)) / np.log(0.8)
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x_vals, y=Upper_limit_graph_rmd, mode='lines', name='Upper Limit RMD', line=dict(color='red')))
    fig.add_trace(go.Scatter(x=x_vals, y=Lower_limit_graph_rmd, mode='lines', name='Lower Limit RMD', line=dict(color='pink')))
    fig.add_trace(go.Scatter(x=[Total_count], y=[Alt_count / Total_count], mode='markers', name='Observed', marker=dict(color='blue')))
    
    fig.update_layout(title="Modified SPRT", xaxis_title="Total number of counts", yaxis_title="Pr over-represented allele")
    
    return pio.to_html(fig, full_html=False)

def generate_chr11_plot(): #TODO how is it deciding what allele to show?
    xrange = np.arange(5225264, 5227272)
    yrange = [0, 1.2]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=xrange, y=[0.5]*len(xrange), mode='lines', name='Chr11', line=dict(color='black')))
    
    fig.update_layout(title="Chromosome 11 Plot", xaxis_title="Position", yaxis_title="Fraction")
    
    return pio.to_html(fig, full_html=False)

def generate_html_content(mean_pat, median_pat, IQR_Pat, mean_Fet, \
                          median_Fet, IQR_Fet, Total_count, Alt_count,FL_SNPs, allele,\
                            report_name):
    # Generate the graphs
    sprt_plot_html = generate_sprt_plot(Total_count,Alt_count)
    chr11_plot_html = generate_chr11_plot()

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