#!/usr/bin/env python

#Import packages
import pandas as pd
import numpy as np
from os import listdir,path
import warnings
import math
import csv
import argparse
import gzip
from sys import exit
import plotly.graph_objects as go
import dash
from dash import html, dcc
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import plotly.express as px
from dash import Dash, dcc, html, Input, Output, State, callback, dash_table

#Function to process all files uploaded and selections selected
def process_logic(df_csv, df_xlsx, fcount_value, out_dir_value, chip_dim_value, timepoint_value, instrument_dropdown):
    # Set up 
    ifc = chip_dim_value
    global instrument_type
    instrument_type = instrument_dropdown # EP1 or BM (Biomark)
    fcount = fcount_value # experiment file that will be used
    tgap = 1 # time gap between mixing of reagents (end of chip loading) and t0 image in minutes

    if instrument_type == 'EP1':
        global count_tp
        count_tp = 1 # End point run
    else:
        count_tp = timepoint_value # number of timepoints, standard for 2h techdev is 25 atm

    # Define variables based on inputs
    global exp_name
    exp_name = fcount + '_' + str(ifc)
    #csv_file = fcount + '_' + ifc + '_' + instrument_type +'_rawdata.csv' # csv file with raw data from Fluidigm RT PCR software
    global csv_file
    csv_file= df_csv

    # Write a path to that excel sheet
    layout_file = df_xlsx # uses 192_assignment.xlsx for layout

    global out_folder
    out_folder = out_dir_value # actual directory is called output!!!!!!!!!!!!!!!!!!

    # Definition of functions
    # Create a dictonary for timepoints, creates a key/value pair (t1/4), (t2/9), (t3/14) etc
    time_assign = {}
    for cycle in range(1,38):
        tpoint = "t" + str(cycle)
        time_assign[tpoint] = tgap + 3 + (cycle-1) * 5
    # calculate the real timing of image
    # used for image and axis labeling 
    def gettime(tname):
        realt = time_assign[tname] # a key is used on the dictonary to get the value, which is stored in realt 
        return (realt)             # returns the value
    
    global probe_df, reference_df, bkgd_ref_df, bkgd_probe_df

#With iloc, it reads the files starting at 0 for the row numbering, so the numbers are 1 less than the original arg parse, because the argparse reads the csv as a path, starting from 1 for the row numbers
   #Create dataframes per results of interest from selected instrument type and ifc 
    if (instrument_type == 'BM'):        
        if (ifc == 92):
            probe_df = csv_file.iloc[18448:27666, : ] #header = 18449, nrows = 9216 # FAM #See how to extract the row and header with iloc and remove read_csv
            reference_df = csv_file.iloc[9230:18448, :] #header = 9231, nrows = 9216) # ROX
            bkgd_ref_df = csv_file.loc[27666:36884, :] #header = 27667, nrows = 9216)
            bkgd_probe_df = csv_file.iloc[36884:46101, :] #header = 36885, nrows = 9216)
        if (ifc == 192):
            #probe_df = pd.read_csv(csv_file,header = 9233, nrows = 4608)
            #reference_df = pd.read_csv(csv_file, header = 4623, nrows = 4608)
            #bkgd_ref_df = pd.read_csv(csv_file, header = 13843, nrows = 4608)
            #bkgd_probe_df = pd.read_csv(csv_file,header = 18453, nrows = 4608)
            probe_df = csv_file.iloc[9237:13846, :] #sep=",", header = 9238, nrows = 4608, skip_blank_lines=False)
            reference_df = csv_file.iloc[4626:9235, :] #header = 4627, nrows = 4608, skip_blank_lines=False)
            bkgd_ref_df = csv_file.iloc[13848:18457, :] #sep=",", header = 13849, nrows = 4608, skip_blank_lines=False)
            bkgd_probe_df = csv_file.iloc[18459:23068, :] #sep=",", header = 18460, nrows = 4608, skip_blank_lines=False)
    elif (instrument_type == 'EP1'):
        if (ifc == 192):
            probe_df = csv_file.iloc[9237:13847, :] #header = 9238, nrows = 4608)
            reference_df = csv_file.iloc[4627:9237, :] #header = 4628, nrows = 4608)
            bkgd_ref_df = csv_file.ioc[18457:23067, :] #header = 18458, nrows = 4608)
            bkgd_probe_df = csv_file.iloc[23067:27677, :] #header = 23068, nrows = 4608)

    # List of dataframes to process
    dataframes = [probe_df, reference_df, bkgd_ref_df, bkgd_probe_df]
    
    # Loop through each dataframe to reindex them with Chamber ID
    for i, df in enumerate(dataframes):
        if not df.empty:
            new_columns = df.iloc[0].tolist()  # Extract the header row and set it as column names
            df.columns = new_columns
            df = df.drop(df.index[0])  # Drop the first row (header row) b/c now its doubled 
            df.set_index('Chamber ID', inplace=True)  # Set 'Chamber ID' as the index
            dataframes[i] = df

    probe_df, reference_df, bkgd_ref_df, bkgd_probe_df = dataframes

    # rename column names
    probe_df.columns = ['t' + str(col) for col in probe_df.columns]
    reference_df.columns = ['t' + str(col) for col in reference_df.columns]
    bkgd_ref_df.columns = ['t' + str(col) for col in bkgd_ref_df.columns]
    bkgd_probe_df.columns = ['t' + str(col) for col in bkgd_probe_df.columns]
    
    
    # Substract the background from the probe and reference data
    probe_bkgd_substracted = probe_df.astype(int).subtract(bkgd_probe_df.astype(int)) #convert to an int to be able to perform subtraction
    ref_bkgd_substracted = reference_df.astype(int).subtract(bkgd_ref_df.astype(int)) #convert to an int to be able to perform subtraction

    # Normalize the probe signal with the reference dye signal
    signal_df = pd.DataFrame(probe_bkgd_substracted/ref_bkgd_substracted) # a table is created with the normalized values

    # reset index
    signal_df = signal_df.reset_index()
    
    # split Chamber ID into SampleID and AssayID
    splitassignment = signal_df['Chamber ID'].str.split("-",n=1,expand=True) 
    signal_df["sampleID"] = splitassignment[0]
    signal_df["assayID"] = splitassignment[1]

    #set index again to Chamber ID
    signal_df = signal_df.set_index('Chamber ID')
  
    # Remove decimal points from column names in signal_df
    new_columns_signal_df = {col: col.replace('.0', '') for col in signal_df.columns}
    signal_df.rename(columns=new_columns_signal_df, inplace=True)
    
    sampleID_list = signal_df.sampleID.unique()
    assayID_list = signal_df.assayID.unique()


    
    # Save csv
    signal_out_csv_1 = f"{out_folder}/{exp_name}_{instrument_type}_1_signal_bkgdsubtracted_norm_{str(count_tp)}.csv"
    signal_df.to_csv(signal_out_csv_1)
    # Create two dictionaries that align the IFC wells to the sample and assay names
    #samples_layout_wo_string = pd.read_excel(path.join('',layout_file),sheet_name='layout_samples')
    
    
    samples_layout_wo_string = layout_file['layout_samples']
    samples_layout = samples_layout_wo_string.map(str)
    
    assays_layout_wo_string = layout_file['layout_assays']
    assays_layout = assays_layout_wo_string.map(str)
    
    assays = layout_file['assays']
    assays_zip = zip(assays_layout.values.reshape(-1),assays.values.reshape(-1))
    assays_dict = dict(assays_zip)
    
    samples = layout_file['samples']

    #function for appending numbers to repeated NTC/blank samples
    from itertools import count
    
    counter = count(1)

    def append_num_to_ntc(x):
        if x in ('NTC', 'blank'):
            return f"{x}_{next(counter)}"
        return x

    #updating samples using the function
    samples = samples.map(lambda x: append_num_to_ntc(x))
    # Create a dictionary with sample numbers and their actual sample name
    samples_zip = zip(samples_layout.values.reshape(-1),samples.values.reshape(-1))
    samples_dict = dict(samples_zip)

    # Map assay and sample names
    signal_df['assay'] = signal_df['assayID'].map(assays_dict)
    signal_df['sample'] = signal_df['sampleID'].map(samples_dict)
 
    # Save csv
    signal_out_csv_2 = f"{out_folder}/ {exp_name}_{instrument_type}_2_signal_bkgdsubtracted_norm_named_{str(count_tp)}.csv"
    #signal_df.to_csv(path.join(out_folder, exp_name+'_' +instrument_type +'_2_signal_bkgdsubtracted_norm_named_' + str(count_tp) +'.csv'))
    signal_df.to_csv(signal_out_csv_2)

    # # Transform and summarize data for plotting

    #create list with timepoints
    # count_tp = 23 # in case you were wrong before
    t_names = []
    for i in range(1,count_tp+1):
        t_names.append(('t' + str(i)))

    # Create a list of all assays and samples
    # only indicate columns with unique assays. np.unique could be used on the list, but messes up our preferred the order
    
    global new_array
    if ifc == 96:
        new_array = np.stack(assays[['C1','C2','C3','C4','C5']].values,axis=-1)
    if ifc == 192:
        new_array = np.stack(assays[['C1','C2', 'C3']].values,axis=-1)
    
    global assay_list
    assay_list = np.concatenate(new_array).tolist()
    print('identified crRNAs: ',len(assay_list))

    # Do this for the samples
    if ifc == 96:
        new_array = np.stack(samples[['C1','C2','C3','C4','C5','C6','C7', 'C8','C9','C10','C11','C12']].values,axis=-1)
    if ifc == 192:
        new_array = np.stack(samples[['C1','C2','C3','C4','C5','C6', 'C7','C8','C9','C10','C11','C12', 'C13','C14','C15','C16','C17','C18', 'C19','C20','C21','C22','C23','C24']].values,axis=-1)
    #     new_array = np.stack(samples[['C1','C2','C3','C4','C5','C6',\
    #                                   'C7','C8','C9','C10','C11','C12']].values,axis=-1)
    
    global sample_list
    sample_list = np.concatenate(new_array).tolist()
    print('identified samples: ',len(sample_list))

    # Grouped medians 
    grouped_stuff = signal_df.groupby(['assay','sample']) #medians = signal_df.groupby(['assay', 'sample']).median(numeric_only= True).reset_index() do this instead when you rewrite the script
    medians = grouped_stuff.median(numeric_only = True)


    # Creating dataframes in a loop: 
    # https://stackoverflow.com/questions/30635145/create-multiple-dataframes-in-loop
    global med_frames
    med_frames = {}
    for name in t_names:
        #time_med = signal_df.groupby(['assay','sample']).median()[name].unstack()
        time_med = medians[name].unstack()
        time_med.index.names=['']
        time_med.columns.names=['']
        med_frames[name] = time_med

    # Write results

    # Write TSV, with one row per (timepoint, sample (target), assay (guide)) triplet
    # The value written is: {median across replicates[(probe signal - probe background) / (reference signal - reference background)]}
    # for different time points

    #with gzip.open(path.join(out_folder, exp_name+ '_'+ instrument_type + '_merged.tsv.gz'), 'wt') as fw:
    with gzip.open(f"{out_folder}/{exp_name}_{instrument_type}_merged.tsv.gz", 'wt') as file:    
        def write_row(row):
            file.write('\t'.join(str(x) for x in row) + '\n')
        header = ['timepoint', 'minute', 'guide', 'target', 'value']
        write_row(header)

        global tp
        for tp in med_frames.keys():
            rt = gettime(tp)
            for target in sample_list:
                for guide in assay_list:
                    v = med_frames[tp][target][guide]
                    row = [tp, rt, guide, target, v]
                    write_row(row)
    return med_frames[tp][sample_list].reindex(assay_list)
    
#CALCULATE THE THRESHOLD
def NTC_logic(df_csv, df_xlsx, fcount_value, out_dir_value, chip_dim_value, timepoint_value, instrument_dropdown):
        
    #Create a dataframe with assays and samples for specific time point
    threshold = process_logic(df_csv, df_xlsx, fcount_value, out_dir_value, chip_dim_value, timepoint_value, instrument_dropdown) #med_frames[tp][sample_list]#.reset_index()

    #Filter threshold dataframe sample names by samples and NTC
    global threshold_with_ntc 
    threshold_with_ntc = threshold.filter(like='NTC', axis = 1)
    threshold_without_ntc = threshold.filter(regex= '^(?!.*NTC).*$', axis = 1) #regex function that selects columns that do not contain NTC

    #Create a loop to divide every column in threshold_without_ntc by each column in threshold_with_ntc, creating 51 new dataframes
    # Create a dictionary to store the result data frames
    threshold_data_frames = {}

    # Iterate over columns in threshold_with_NTC
    for col_name, col_values in threshold_with_ntc.items():
        # Perform division for each column in threshold_without_NTC by the current column in threshold_with_NTC
        result_df = threshold_without_ntc.div(col_values, axis=0)  # Use axis=0 to divide row-wise
        
        # Store the result data frame in the dictionary with a meaningful key
        key = f"Sample_divided_by_{col_name}"
        threshold_data_frames[key] = result_df #To print results in dictionary, do threshold_data_frames["Sample_divided_by_"col_name"]
    #print(threshold_data_frames["Sample_divided_by_InfA NTC"])

    #Melt all the data frames in the dictionary
    melted_threshold = [] #empty list to store melted dataframe

    for key, df in threshold_data_frames.items():
        melted_dictionary = pd.melt(df, ignore_index=False, #this preserves the assay list as the index
                                    var_name = 'Samples',
                                    value_name=f'{key}')
        melted_threshold.append(melted_dictionary)

    # Concatenate all melted data frames into a single data frame
    global combined_melted_threshold
    combined_melted_threshold = pd.concat(melted_threshold, axis=1)

    #Drop duplicated Samples columns that occured after concatenating
    combined_melted_threshold = combined_melted_threshold.loc[:, ~combined_melted_threshold.columns.duplicated()]
    combined_melted_threshold.reset_index(inplace=True) #Makes the assay index a column

    #Rename the new Assay column
    combined_melted_threshold.rename(columns={"" : "Assay"}, inplace= True)

    #Change the name of the NTC columns to drop preceeding 'sample' name
    # Specify the prefix to be removed
    prefix_to_remove = 'Sample_divided_by_'

    # Remove the specified prefix from column names
    new_column_names = [col.replace(prefix_to_remove, '') for col in combined_melted_threshold.columns]

    # Create a dictionary to map old column names to new column names
    columns_mapping = {old_col: new_col for old_col, new_col in zip(combined_melted_threshold.columns, new_column_names)}

    # Rename the columns using the created mapping
    combined_melted_threshold.rename(columns=columns_mapping, inplace=True)

    #Save to a csv
    threshold_out_csv = f"{out_folder}/ {exp_name}_{instrument_type}_threshold_{str(count_tp)}.csv"
    combined_melted_threshold.to_csv(threshold_out_csv)
    return combined_melted_threshold

    
def threshold_logic(df_csv, df_xlsx, fcount_value, out_dir_value, chip_dim_value, timepoint_value, instrument_dropdown):
    #Extract the dilution factor from the sample name in the sample column and create two new columns
    combined_melted_threshold_heat = NTC_logic(df_csv, df_xlsx, fcount_value, out_dir_value, chip_dim_value, timepoint_value, instrument_dropdown)

    extract = combined_melted_threshold_heat['Samples'].str.split("-",n=1,expand=True) #Sample is split into two after '-'. n=1 means split into two indicies 
    combined_melted_threshold_heat['Sample ID'] = extract[0] #First indicies
    combined_melted_threshold_heat['Dilution Factor'] = extract[1] #Second indicies

    #filter out the ntc columns to use as the values for pivoting
    ntc_columns2 = combined_melted_threshold_heat.filter(like='NTC').columns
  
    #Pivot the sample ID column into a wide format, keeping the Dilution Factor and Assay columns intact
    pivot_threshold_heat = pd.pivot(combined_melted_threshold_heat,
                                    index= ['Assay', 'Dilution Factor'],
                                    columns = 'Sample ID',
                                    values = ntc_columns2

    )
    pivot_threshold = pivot_threshold_heat.reset_index()

    #Melt the dataframe for better graphing
    melted_df = pd.melt(combined_melted_threshold_heat, id_vars=['Assay', 'Samples'], 
                        value_vars=ntc_columns2, 
                        var_name='NTC', value_name='Value')
  
    #Turn the new_array (Holds the map for the samples per column from the excel file) array into a dataframe for better manipulatation 
    threshold_array = new_array
    threshold_array_df = pd.DataFrame(threshold_array) #creates a shape of [24,8]

    #Transpose dataframe to create a shape of [8,24] to match the excel map
    threshold_array_df = threshold_array_df.T
    threshold_array_df.columns = [f'Column{i}' for i in range(1, threshold_array_df.shape[1] + 1)] #name the columns 1-24
    #print(threshold_array_df)

    #Melt dataframe for easier merging, and remove all sample rows with NTC
    threshold_array_melt = pd.melt(threshold_array_df,
                                var_name='Column Number', value_name='Samples')
    #print(threshold_array_melt)
    threshold_array_without_ntc = threshold_array_melt[~threshold_array_melt['Samples'].str.contains('NTC')]
    #print(threshold_array_without_ntc)

    #Merge the melted_df and the threshold_array_without_ntc dataframes together for graphing
    #merged_threshold_df
    global merged_threshold_df
    merged_threshold_df = pd.merge(melted_df, threshold_array_without_ntc,
                                    how= 'left',
                                    on= 'Samples') 
    return merged_threshold_df

    
#Create Dash App
import base64
import datetime
import io

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css'] #These are external CSS files that can be applied to the Dash web application to enhance its styling. 

carmen = dash.Dash(__name__, external_stylesheets= external_stylesheets) #initilazies the Dash app and applies the style sheet

#Layout for File upload and criteria selection

carmen.layout = html.Div([
    html.Div(children=[
    html.H1(children='CARMEN', # Site Heading
            style={'textAlign': 'center'},), # Align to center

    html.Div(children='''Multiplexed CRISPR-based microfluidic platform for clinical testing of respiratory viruses and
                         identification of SARS-CoV-2 variants'''), # Site Sub-Heading
        ]),
    
    html.Hr(), #horizontal line
    
    html.Div(children='''
        Upload CSV file. Be sure '.csv' is in the file name.
    '''),
    
    html.Div([
        dcc.Upload(
        id='upload-csv-data',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Upload CSV File')
        ]),
        style={
            'width': '20%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'solid',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        # Allow multiple files to be uploaded
        multiple=False
    ),
    html.Div(id='output-data-csv-upload'), #container where the parsed data lives

    html.Div(children='''
        Upload excel file. Be sure '.xlsx' is in the file name.
    '''),
    
    dcc.Upload(
        id='upload-xlsx-data',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select xlsx Files')
        ]),
        style={
            'width': '20%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'solid',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        # Allow multiple files to be uploaded
        multiple=False
    ),
    html.Div(id='output-data-xlsx-upload'), #container where the parsed data lives
    
    # Add dcc.Store components to store data and file names
    dcc.Store(id='csv-data-store'),
    dcc.Store(id='xlsx-data-store'),
    dcc.Store(id='csv-filename-store'),
    
    ]),

    html.Hr(), #horizontal line
    
    html.Div(children='''
        Select ifc chip dimmension.
    '''),
    
    html.Div([
        dcc.Dropdown( #dropdown to selct from
            id='chip-dim-dropdown',
            options=[{'label': '96', 'value': 96},
                     {'label': '192', 'value': 192}],
            multi=False,
            placeholder="select ifc chip-dim",
            style={'marginBottom': '30px'}),

    html.Div(children='''
        Select desired timepoint.
    '''),
    
    dcc.Dropdown( #dropdown to select from
            id='timepoint-dropdown',
            options=[{'label': i, 'value': i} for i in range(1, 14)],
            multi=False,
            placeholder="select timepoint",
            style={'marginBottom': '30px'}),
    
    html.Div(children='''
        Select instrument type.
    '''),
    
    dcc.Dropdown( #dropdown to selct from
        id='instrument-dropdown',
        options=[{'label': 'Biomek', 'value': 'BM'},
                    {'label': 'EP1', 'value': 'EP1'}],
        multi=False,
        placeholder="select instrument type")

    ]),
    
    html.Hr(), #horizontal line

    html.Div([
    html.Div(children='''
        Enter the name of the new file.
    '''),
    
    html.Div([
        dcc.Input( # empty text box to write in
            id='fcount-textarea',
            value='Enter experiment file name',
            type='text'
        )
    ], style={'marginBottom': '30px'}),  # Adding some margin for spacing

    html.Div(children='''
        Enter where the outputs should be saved. EX: 'output'.
    '''),
    
    html.Div([
        dcc.Input( # empty text box to write in
            id='out_dir-textarea',
            value='Enter output directory',
            type='text'
        )
        ])
    ]),
    
    html.Hr(), #horizontal line
    
    html.Button(id='Submit-Data',  # Submit all options
                children='SUBMIT DATA',
                n_clicks=0,
                style={'fontSize': '20px',
                       'width': '20%',
                       'height': '60px'}
    ),

    # Add placeholders for the heatmaps
    html.Div([
        dcc.Graph(id='fam-heatmap'),
        dcc.Graph(id='ntc-heatmap')
    ]),
    
    
    # Add Dropdown placeholders for Threshold Heatmap   
    html.Div([
        html.Div(children='''
        Select Sample.
    '''),
        
        dcc.Dropdown(id='sample-dropdown',
                 placeholder="Select sample",
                 style={'marginBottom': '30px'}),
    
        html.Div(children='''
        Select NTC.
    '''),
        
        dcc.Dropdown(id='NTC-dropdown',
                 placeholder="Select NTC"),
    
        # Add Placeholder for Threshold Heatmap
        dcc.Graph(id='threshold-heatmap')]),

    # # Placeholder for filtered threshold values table
    # html.Div([
    #     html.H3(f'Table of Sample Threshold Passing or Failing Values per Guide RNA'),
    #     dash_table.DataTable(id='threshold-values-table',
    #     page_size=15, # Display only 15 rows at a time
    #     style_table={'border': '2px solid black'},  # Add border to the table
    #     style_cell={'border': '1px solid grey'}, # Add border to the cells
    #         )  
    #     ])
    
    #Placeholder for Guide RND (Assay) filtering
    html.Div([
        html.Div(children='''
        Select Guide RNA.
    '''),
        
        dcc.Dropdown(id='assay-dropdown',
                 placeholder="Select Guide RNA",
                 style={'marginBottom': '30px'}),
    
        html.Div(children='''
        Select NTC.
    '''),
        
        dcc.Dropdown(id='NTC2-dropdown',
                 placeholder="Select NTC"),
    ]),
    
    #Placeholder for Threshold Values Table    
    html.Div([
        html.H3(f'Table of Sample Threshold Values per Guide RNA'),
        dash_table.DataTable(
            id='threshold-table',
            columns=[],
            data=[],
            page_size=15,
            style_table={'border': '2px solid black'},  # Add border to the table
            style_cell={'border': '1px solid grey'}, # Add border to the cells
            tooltip_delay=0,
            tooltip_duration=None
            )
        ])
    ])

# Function to parse data from uploaded CSV file
def parse_csv(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
            return df.to_dict('records')
    except Exception as e:
        print(e)
        return "There was an error reading your file"

# Function to parse data from uploaded xlsx file
def parse_excel(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        if 'xls' in filename:
            xls = pd.ExcelFile(io.BytesIO(decoded))
            sheet_names = xls.sheet_names
            sheet_data = {}
            for sheet in sheet_names:
                df = pd.read_excel(xls, sheet_name=sheet)
                sheet_data[sheet] = df.to_dict('records')
            return sheet_data
    except Exception as e:
        print(e)
        return "There was an error reading your file"
    
# Callback for handling csv upload and processing
@carmen.callback(
    Output('csv-data-store', 'data'),
    Output('output-data-csv-upload', 'children'),
    Input('upload-csv-data', 'contents'),
    State('upload-csv-data', 'filename'))

def csv_process(csv_contents, csv_filename):
    #Check if csv has been uploaded
    global df_csv
    if csv_contents is not None:
        df_csv = parse_csv(csv_contents, csv_filename) #read file like this
        if df_csv:
            return df_csv, html.Div([
                html.H5(csv_filename)
            ])
    return None, html.Div()

# Callback for handling xlsx upload and processing
@carmen.callback(
    Output('xlsx-data-store', 'data'),
    Output('output-data-xlsx-upload', 'children'),
    Input('upload-xlsx-data', 'contents'),
    State('upload-xlsx-data', 'filename'))

def xlsx_process(contents_xlsx, filename_xlsx):
    global df_xlsx
    if contents_xlsx is not None:
        df_xlsx = parse_excel(contents_xlsx, filename_xlsx) #read file like this
        if df_xlsx:
            return df_xlsx, html.Div([
                html.H5(filename_xlsx),
                html.H6('Sheets: ' + ', '.join(df_xlsx.keys()))
            ])
    return None, html.Div()

# Callback for handling processing of all uploaded files and dropdown selections to produce FAM and NTC heatmap
@carmen.callback(
    Output('fam-heatmap', 'figure'),
    Output('ntc-heatmap', 'figure'),
    Output('sample-dropdown', 'options'),
    Output('sample-dropdown', 'value'),
    Output('NTC-dropdown', 'options'),
    Output('NTC-dropdown', 'value'),
    Output('assay-dropdown', 'options'),
    Output('assay-dropdown', 'value'),
    Output('NTC2-dropdown', 'options'),
    Output('NTC2-dropdown', 'value'),
    # Output('threshold-values-table', 'data'),
    # Output('threshold-values-table', 'columns'),
    Input('Submit-Data', 'n_clicks'),
    State('csv-data-store', 'data'),
    State('xlsx-data-store', 'data'),
    State('fcount-textarea', 'value'),
    State('out_dir-textarea', 'value'),
    State('chip-dim-dropdown', 'value'),
    State('timepoint-dropdown', 'value'),
    State('instrument-dropdown', 'value')
)

def button_click(n_clicks, csv_data, xlsx_data, fcount_value, out_dir_value, chip_dim_value, timepoint_value, instrument_dropdown):
    # Empty lists to hold figures
    fam_heatmap = {}
    ntc_heatmap = {}
    
    # Empty lists to hold dropdown values
    sample_options = []
    sample_value = ''
    ntc_options = []
    ntc_value = ''
    assay_options = []
    assay_value = ''
    ntc2_options = []
    ntc2_value = ''
    
    # Check if the submit button is clicked    
    if n_clicks != 0: 
    # Modify processing logic based on dropdown values selected
        if csv_data is not None and xlsx_data is not None:
            # Convert csv_data back to DataFrame
            df_csv = pd.DataFrame(csv_data)
            
            # Convert xlsx_data back to DataFrame dictionary
            df_xlsx = {sheet: pd.DataFrame(data) for sheet, data in xlsx_data.items()}
            
            #Process logic
            processed_data = process_logic(df_csv, df_xlsx, fcount_value, out_dir_value, chip_dim_value, timepoint_value, instrument_dropdown) #but apply the modifications based on what the user selects for these options
            processed_data2 = NTC_logic(df_csv, df_xlsx, fcount_value, out_dir_value, chip_dim_value, timepoint_value, instrument_dropdown)
            processed_data3 = threshold_logic(df_csv, df_xlsx, fcount_value, out_dir_value, chip_dim_value, timepoint_value, instrument_dropdown)
            #Check if processed_data is not None or empty
            if processed_data is not None and processed_data2 is not None and processed_data3 is not None:
                # Generate and update the FAM heatmap
                fam_heatmap = generate_fam_heatmap(processed_data)  # function created later in script
                # Generate and update the NTC heatmap
                ntc_heatmap = generate_ntc_heatmap(processed_data2)  # function created later in script
                
                # Update sample-dropdown and ntc-dropdown for threshold heatmap
                sample_options= [{'label': sample, 'value': sample} for sample in merged_threshold_df['Column Number'].unique()]
                sample_value = sample_options[0]['value'] if sample_options else ''

                ntc_options = [{'label': ntc, 'value': ntc} for ntc in merged_threshold_df['NTC'].unique()]
                ntc_value = ntc_options[0]['value'] if ntc_options else ''
                
                # Table for Threshold Values
                #Update dropdown values for threshold table
                assay_options= [{'label': assay, 'value': assay} for assay in merged_threshold_df['Assay'].unique()]
                assay_value = assay_options[0]['value'] if assay_options else ''

                ntc2_options = [{'label': ntc2, 'value': ntc2} for ntc2 in merged_threshold_df['NTC'].unique()]
                ntc2_value = ntc2_options[0]['value'] if ntc2_options else ''

                # threshold_table_data = threshold_filter.to_dict('records')
                # threshold_table_columns = [{'name': i, 'id': i} for i in threshold_filter.columns]
                return [fam_heatmap, ntc_heatmap, sample_options, sample_value, ntc_options, ntc_value, assay_options, assay_value, ntc2_options, ntc2_value]

#Vizualize the data

#Plotly Heatmaps of FAM normalized data, with another heatmap with NTC FAM Normalized data
def generate_fam_heatmap(norm_df):  
    
    #Dataframe to be used
    norm_df = med_frames[tp][sample_list].reindex(assay_list)    

    # Heatmap of FAM normalized data
    fig1 = go.Figure(data=go.Heatmap(
        x=norm_df.columns.tolist(),
        y=norm_df.index.tolist(),
        z=norm_df,
        colorscale='Reds',
    ))

    fig1.update_layout(
        xaxis_title='Samples',
        yaxis_title="Assay",
        title=f'Heatmap of FAM Normalized Data for {csv_file} at Time-Point {count_tp}',
        autosize=False,
        width=2800,
        height=800,
        margin=dict(l=50, r=50, b=50, t=175, pad=5)
    )
    return fig1

def generate_ntc_heatmap(NTC_fam):
    
    #Dataframe to be used
    NTC_fam = threshold_with_ntc # Renamed this dataframe so i don't alter changes to the main dataframe, also for ease of tracking different dataframes and what they do
    
    # Heatmap of NTC only FAM normalized data
    fig2 = go.Figure(data=go.Heatmap(
        x=NTC_fam.columns.tolist(),
        y=NTC_fam.index.tolist(),
        z=NTC_fam,
        colorscale='reds',
    ))

    fig2.update_layout(
        xaxis_title='NTC',
        yaxis_title="Assay",
        title=f'Heatmap of FAM Normalized NTC Data for {csv_file} at Time-Point {count_tp}',
        autosize=False,
        width=2800,
        height=800,
        margin=dict(l=50, r=50, b=50, t=175, pad=5)
    )
    return fig2

# Callback for generating threshold heatmap
@carmen.callback(
    Output('threshold-heatmap', 'figure'),
    Input('sample-dropdown', 'value'),
    Input('NTC-dropdown', 'value'),
)

def generate_threshold_heatmap(sample_value, ntc_value):
    filtered_df = merged_threshold_df[(merged_threshold_df['Column Number'] == sample_value) & (merged_threshold_df['NTC'] == ntc_value)]

    # Create Threshold Heatmap
    threshold_fig = go.Figure(data=go.Heatmap(
        x=filtered_df['Assay'],
        y=filtered_df["Samples"],
        z=filtered_df['Value'],
        colorscale='Reds',
        ))

    threshold_fig.update_layout(
        xaxis_title='Assay',
        yaxis_title="Samples",
        title=f'Heatmap for Threshold Values for {csv_file} at Time-Point {count_tp}',
        autosize=False,
        width=1000,
        height=800,
         margin=dict(l=50, r=50, b=50, t=175, pad=5)
    )

    return threshold_fig

# Callback for filtering and displaying the table
@carmen.callback(
    Output('threshold-table', 'data'),
    Output('threshold-table', 'columns'),
    Output('threshold-table', 'tooltip_data'),
    Input('assay-dropdown', 'value'),
    Input('NTC2-dropdown', 'value')
)
def update_table(assay_value, ntc2_value):
    #threshold_filter = merged_threshold_df[threshold_filter['Assay'] == assay_value] & threshold_filter[threshold_filter['NTC'] == ntc2_value]

    threshold_filter = merged_threshold_df
     # Empty lists to hold table values
    threshold_table_data = []
    threshold_table_columns = []
    
    if assay_value:
        threshold_filter = threshold_filter[threshold_filter['Assay'] == assay_value]
    if ntc2_value:
        threshold_filter= threshold_filter[threshold_filter['NTC'] == ntc2_value]
    
    # Ensure 'Samples' and 'Column Number' are of appropriate data types
    threshold_filter['Samples'] = threshold_filter['Samples'].astype(str)
    threshold_filter['Column Number'] = threshold_filter['Column Number'].astype(str)
    
    #Pivot dataframe to match excel layout
    pivot_data = pd.pivot(threshold_filter,
                          columns = 'Column Number',
                          values = 'Samples')
    threshold_table_data = pivot_data.reset_index().to_dict('records')
    threshold_table_columns = [{'name': str(i), 'id': str(i)} for i in pivot_data.columns]

    #Incorporate tooltip data(hover data)
    tooltip_data = [
        {
            column: {'value': str(value), 'type': 'markdown'}
            for column, value in row.items()
        } for row in threshold_filter[['Value']].to_dict('records')
    ]
    return threshold_table_data, threshold_table_columns, tooltip_data

carmen.run_server(debug=True)