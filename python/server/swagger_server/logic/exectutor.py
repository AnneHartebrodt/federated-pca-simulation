from flask import current_app, jsonify
from PCA import simulation_runner
from PCA import master
from server.swagger_server.logic import scheduler, config
from server.swagger_server.logic import InputFileError
import json
from server.swagger_server import numpy_encoder
import traceback
import numpy as np

'''
In this file the steps defined in the workflow file need to be implemented
The function names need to correspond to the names defined in the workflow_description file.
'''

def  setup(ismaster, params):
    '''
    This method is called by the Controller whenever setup is called. It runs all
    the necessary steps to Configure the app and the Scheduler
    :param ismaster:
    :return:
    '''
    current_app.logger.info('Configuring workflow...')
    current_app.PCA = master.Distributed_DP_PCA(current_app.logger)
    current_app.logger.info('... adding configuration file')
    current_app.Config = config.Configuration(ismaster, params, current_app.logger)
    current_app.logger.info('...setting up scheduler')
    current_app.Scheduler = scheduler.Scheduler(current_app.Config.get_workflow_definition())
    current_app.logger.info('... configuration done')

def init():
    '''
    Read local data and perform local PCA
    :return:
    '''

    #def data_import(self, filename, seed=11, nr_samples=None, header=None, rownames=None, sep='\t', outfile=None):
    current_app.logger.info('Executing step 1')
    current_app.logger.info('Data mounted to /mounted_local_folder/')
    params = current_app.Config.get_params()
    current_app.logger.info('Parameters:')
    current_app.logger.info(params)
    try:
        current_app.logger.info('Reading data ...')
        current_app.logger.info(params['datafile'])
        current_app.logger.info(params['header'])
        current_app.logger.info(params['sample_ids'])
        current_app.logger.info(params['separator'])
        current_app.logger.info(params['scaled_datafile'])
        scaled_data, sample_ids, variable_names = current_app.PCA.data_import(filename=params['datafile'], seed = 11, nr_samples=None, header=params['header'], rownames=params['sample_ids'], sep =params['separator'], outfile=params['scaled_datafile'])
        current_app.logger.info('... done!')
    #except:
    except:
        current_app.logger.error('Error reading local data!')
        current_app.logger.error(traceback.print_exc())
    try:
        current_app.logger.info('Performing local PCA ...')
        PC = current_app.PCA.local_PCA(scaled_data, params['epsilon'], params['delta'], noise=True, ndims=10)
        current_app.logger.info('... done!')
    except:
        current_app.logger.error('Error performing local PCA!')
        current_app.logger.error(traceback.print_exc())
    current_app.logger.info('Finished local PCA')
    return PC.tolist()


def step2():
    current_app.logger.info('Master only ')
    svd_list = []
    try:
        current_app.logger.info(current_app.inbox)
        incoming = current_app.inbox
        # the object is still double json encoded.
        # first level generic json encoding of Flask app
        # second level numpy encoded data object
        incoming = parse_json(incoming)
        current_app.logger.info('incoming data:')
        current_app.logger.info(incoming)
        svd_list=incoming
    except:
        current_app.logger.error(traceback.print_exc())

    try:
        local= np.array(current_app.outbox)
        current_app.logger.info('local data:')
        current_app.logger.info(local)
        svd_list.append(local)
    except:
        current_app.logger.error(traceback.print_exc())
    current_app.logger.info(svd_list)
    u,v = current_app.PCA.aggregate_partial_SVDs(svd_list= svd_list, ndims=4)
    return u.tolist()

def parse_json(jsonob):
    print(jsonob)
    if isinstance(jsonob, list):
        deparsed=[]
        for elem in jsonob:
            deparsed.append(np.array(json.loads(elem['content']['payload'])['data']))
    else:
        deparsed = np.array(json.loads(jsonob['payload'][0]['content']['payload'])['data'])
    return deparsed


def finalize():
    params = current_app.Config.get_params()
    current_app.logger.info('Finalize')
    if current_app.Config.get_master():
        with open(params['output'], 'a') as handle:
            handle.write(json.dumps(current_app.outbox))
    else:
        data = parse_json(current_app.inbox)
        with open(params['output'], 'a') as handle:
            handle.write(json.dumps(parse_json(current_app.inbox).tolist()))
    return 'Finished, nothing to tell'

