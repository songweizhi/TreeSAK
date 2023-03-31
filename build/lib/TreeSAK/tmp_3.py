
import pandas as pd
import plotly.express as px


def plot_distribution(df_txt, output_plot):

    df = pd.read_table(df_txt, sep=',')
    fig = px.violin(df, y="Section", x="Value", color="Class", points=False, orientation="h")
    fig.update_traces(side="positive", fillcolor='rgba(0,0,0,0)', width=1.8)
    fig.update_traces(showlegend=True)
    fig.layout.template = "simple_white"
    # fig.layout.width = 700
    # fig.layout.height = 750
    # fig.update_xaxes(range=[40, 0])
    # fig.update_layout(margin_t=10, title_text='Demo', title_x=0.5)
    fig.write_image(output_plot)


df_txt      = '/Users/songweizhi/Desktop/demo.csv'
output_plot = '/Users/songweizhi/Desktop/demo.pdf'
plot_distribution(df_txt, output_plot)
