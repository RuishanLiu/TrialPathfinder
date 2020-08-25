from sqlalchemy import create_engine
import pandas as pd
import numpy as np
import datetime
import os

MISSING = 'missing' # string to indicate missingness of a feature
ICD10_CLASS1 = ['A00-B99', 'C00-D49', 'D50-D89', 'E00-E89', 'F01-F99', 'G00-G99', 
                'H00-H59', 'H60-H95', 'I00-I99', 'J00-J99', 'K00-K95', 'L00-L99', 
                'M00-M99', 'N00-N99', 'O00-O9A', 'P00-P96', 'Q00-Q99', 'R00-R99',
                'S00-T88', 'V00-Y99', 'Z00-Z99']


redshift_con = create_engine('postgresql://mdhrw:mdhRW2019!@phc-curated-data-01-uat-redshift-cluster.cgqkjd9bwok8.eu-central-1.redshift.amazonaws.com:5439/phcdap')

def load_data_aws(name, dataset='edm', n_data=None, start_data=0):
    '''
    Load data as pandas array.
    name: table name
    dataset: 'edm' or 'cgdb'
    n_data: number of data loaded (top n_data). Load all if None.
    '''
#     str_n = '' if n_data is None else 'top %d' % n_data
#     df = pd.read_sql_query('select %s * from mdh_%s.%s;' % (str_n, dataset, name), redshift_con)
    
    str_n = '' if n_data is None else ' LIMIT %d OFFSET %d' % (n_data, start_data)
    df = pd.read_sql_query('select * from mdh_%s.%s%s;' % (dataset, name, str_n), redshift_con)
    
#     str_n = '' if n_data is None else 'top %d' % n_data
#     query = "select %s * from (select patientid from mdh_%s.enhanced_advancednsclc) a inner join(select * from mdh_%s.%s) b on a.patientid = b.patientid;"  % (str_n, dataset, dataset, name) 
#     df = pd.read_sql_query(query, redshift_con)
#     df = df.iloc[:, 1:]
    
    return df


def load_data(name, dir_load='../../data', dataset='edm'):
    '''
    Load data as pandas array.
    name: table name
    dataset: 'edm' or 'cgdb'
    n_data: number of data loaded (top n_data). Load all if None.
    '''
    if dataset == 'edm':
        file_name = os.path.join(dir_load, 'edm/FLT_GNE_%s_20190630.txt' % name)
        try:
            # df = pd.read_csv(file_name, sep='\x01')
            chunksize = 1e6
            tfr = pd.read_csv(file_name, sep='\x01', chunksize=chunksize, iterator=True)
            df = pd.concat(tfr, ignore_index=True)
        except:
            df = pd.read_csv(file_name+'.gz', compression='gzip', sep='\x01')
    else:
        file_name = os.path.join(dir_load, 'cgdb/%s.txt' % name)
        df = pd.read_csv(file_name, sep='|')
    return df


def process_table(df, enhanced_advancednsclc=None):
    '''Processing loaded tables'''
    # All the column names to be lower case
    df.columns = [c_name.lower() for c_name in df.columns]
    # Select patients which is advanced
    if enhanced_advancednsclc is not None:
        df = df.loc[df['patientid'].isin(enhanced_advancednsclc['patientid'])]
    else:
        print('Enhanced_AdvancedNSCLC Table not loaded yet')
    # Change date time datatype
    for name_feature in df.columns:
        if 'date' in name_feature:
            df[name_feature] = pd.to_datetime(df[name_feature], errors = 'coerce')
    return df