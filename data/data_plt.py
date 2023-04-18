import pickle
import numpy as np
import plotly.graph_objects as go
import pandas as pd
from datetime import datetime, timedelta
import os


start_date = datetime(2022, 10, 1)   # 起始日期
end_date = datetime(2023, 4, 3)   # 结束日期

date_list = []   # 存储日期的列表

# 循环遍历每一天，并将日期添加到列表中
while start_date < end_date:
    date_list.append(start_date.strftime('%Y%m%d'))
    start_date += timedelta(days=1)


fig = go.Figure()

GAP = 20
Flight = "CPA250"
for date in date_list:
    csvPath = "./flightaware_data/" + Flight + "/" + Flight + "_csv/" + Flight + "_" + date + ".csv"
    # print(csvPath)
    if os.path.exists(csvPath):
        df = pd.read_csv(csvPath)
        for i in np.arange(0, round((len(df) - GAP)), GAP):
            fig.add_trace(
                go.Scattergeo(
                    lon=[df.iloc[i, 1], df.iloc[i + GAP, 1]],
                    lat=[df.iloc[i, 0], df.iloc[i + GAP, 0]],
                    mode='lines',
                    line=dict(width=1, color='black', dash='solid')))

Flight = "CES7551"
GAP = 20
date_list_choose = ["20230107", "20230114", "20230115", "20230117",
                    "20230128", "20230204", "20230205", "20230304",
                    "20230307", "20230311", "20230318", "20230402"]
for date in date_list_choose:
    csvPath = "./flightaware_data/" + Flight + "/" + Flight + "_csv/" + Flight + "_" + date + ".csv"
    # print(csvPath)
    if os.path.exists(csvPath):
        df = pd.read_csv(csvPath)
        for i in np.arange(0, len(df) - GAP, GAP):
            fig.add_trace(
                go.Scattergeo(
                    lon=[df.iloc[i, 1], df.iloc[i + GAP, 1]],
                    lat=[df.iloc[i, 0], df.iloc[i + GAP, 0]],
                    mode='lines',
                    line=dict(width=1, color='red')))

fig.update_layout(
    title='Multiple Geo Subplots',
    geo=dict(
        scope='world',
        projection_type='orthographic',  # 'natural earth', 'orthographic'
        coastlinecolor='rgb(217, 217, 217)',
        showcountries=True,
        countrycolor='rgb(128, 128, 128)',
        showocean=True,
        oceancolor='rgb(12, 47, 94)',  # rgb(10, 88, 175)
        showland=True,
        landcolor='rgb(217, 217, 217)'
    )
)

fig.show()
