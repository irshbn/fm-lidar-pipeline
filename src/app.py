import pandas as pd
import numpy as np
import io
import base64
import datetime as dt
from datajoint.errors import DuplicateError
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, dash_table, html, dcc, Input, Output, State, callback
from src import my_schema as db

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

query = db.Experiment
df = pd.DataFrame(query.fetch())

### App layout ###
app = Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)

app.layout = html.Div(children=[
    html.H1('FM-LIDAR Pipeline Demo',
            style={'textAlign':'center'}),
    
    html.Div('Select experiment:'),

    dcc.Dropdown(df['exp_desc'],
                 id='exp-selector'),

    dcc.Tabs(
        id='upload-dashboard-tabs',
        children=[
            dcc.Tab(label='Dashboard', value='dashboard-tab'),
            dcc.Tab(label='Upload data', value='upload-tab')
        ],
        value='dashboard-tab'
    ),

    html.Div(id='tab-content')
])
              

### Waveform upload parser ###
def parse_wfmdata(contents: str, filename: str):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        blob = np.genfromtxt(io.StringIO(decoded.decode('utf-8')),delimiter=',',missing_values='')
        delta = blob[1,1]
        wfm = blob[3:,1]

        id_str, date, time = filename.removeprefix('Wfm').removesuffix('.csv').split(sep='_')

        id = int(id_str.removesuffix('.'))

        dd, mm, yyyy = date.split(sep='.')
        h, m, s = time.split(sep='-')
        timestamp=dt.datetime(int(yyyy), int(mm), int(dd),
                              int(h), int(m), int(s),
                              tzinfo=dt.timezone(dt.timedelta(hours=+3)))

    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])
    
    df = pd.DataFrame(dict(
        x=np.linspace(0,delta*len(wfm),len(wfm)),
        y=wfm
    ))

    meta_df = pd.DataFrame([dict(
        wfm_id=id,
        wfm_timestamp=timestamp,
        wfm_delta=delta
    )])

    return html.Div([
        dcc.Graph(
            figure=px.line(df, x='x', y='y',
                           title=f'Extracted waveform',
                           labels={"x": "Time (s)", "y": "Amplitude (a.u.)"},
                           render_mode='webgl'),
            id='wfm-graph'
        ),

        dash_table.DataTable(meta_df.to_dict('records'),
                             id='wfm-metadata'),

        dcc.ConfirmDialogProvider(
            children=html.Button('Push to database'),
            id='push-to-db',
            message='You are going to insert waveform data into database. Continue?'
        )
    ])


### Callbacks and controls ###
@callback(Output('output-data-upload', 'children', allow_duplicate=True),
          Input('upload-data', 'contents'),
          State('upload-data','filename'),
          prevent_initial_call=True)
def update_output(contents, filename):
    if contents is not None:
        return parse_wfmdata(contents, filename)


@callback(Output('upload-data','disabled'),
          Input('exp-selector','value'))
def disallow_upload(selected_value):
    if selected_value is None:
        return True
    else:
        return False


@callback(Output('output-data-upload', 'children'),
          Input('push-to-db','submit_n_clicks'),
          State('wfm-graph','figure'),
          State('wfm-metadata','data'),
          State('exp-selector','value'),
          State('output-data-upload', 'children'),
          prevent_initial_call=True)
def push_wfm(n_clicks, fig: go.Figure, metadata, exp_desc, status_quo):
    if n_clicks:
        id = (db.Experiment & f'exp_desc = "{exp_desc}"').fetch1('exp_id')
        data = fig['data'][0]['y']
        entry = dict(
            wfm_id=metadata[0]['wfm_id'],
            exp_id=id,
            wfm_timestamp=dt.datetime.strptime(metadata[0]['wfm_timestamp'],'%Y-%m-%dT%H:%M:%S%z'),
            wfm_delta=metadata[0]['wfm_delta'],
            wfm_data=data
        )
        try:
            db.Waveform.insert1(entry)
        except DuplicateError:
            return html.Div(html.H5('Entry with the same ID already exists. Try again.'))

        db.WaveformPreconditioned.populate()
        db.LeastSquaresFit.populate()

        return html.Div(html.H5('Success!'))

    return status_quo


@callback(Output('tab-content','children'),
          Input('upload-dashboard-tabs','value'),
          Input('exp-selector','value'))
def render_tab(selected_tab, exp_desc):
    if selected_tab == 'upload-tab':
        return [
            dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Waveform Files'),
                ' (WfmXX._DD.MM.YYYY_hh-mm-ss.csv)'
            ]),
            style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
            }),

            html.Div(id='output-data-upload')
        ]
    elif selected_tab == 'dashboard-tab':
        if exp_desc is None:
            return html.Div('Select an experiment to view results.')
        
        exp_id = (db.Experiment & f'exp_desc = "{exp_desc}"').fetch1('exp_id')
        
        query = db.LeastSquaresFit * db.Waveform & f'exp_id = {exp_id}'

        df = pd.DataFrame(query.fetch())
        df = df[df.columns.difference(['wfm_data','wfm_timestamp','wfm_delta','exp_id'])]
        df['lsq_error95p'] = df['lsq_sigma'] * 2

        fig = px.line(df, 
                      x='wfm_id', 
                      y='lsq_frequency',
                      error_y='lsq_error95p',
                      title='Frequency Graph',
                      labels=dict(
                          wfm_id='Index',
                          lsq_frequency='Frequency (Hz)'
                      ))

        return [
            dcc.Graph(figure=fig),

            dash_table.DataTable(df.to_dict('records'))
        ]


if __name__ == '__main__':
    app.run(debug=True)