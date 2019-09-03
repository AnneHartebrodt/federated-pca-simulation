import connexion
from flask import current_app

from server.swagger_server.models.action_result import ActionResult  # noqa: E501
from server.swagger_server.models.status_response import StatusResponse  # noqa: E501
from server.swagger_server.logic import exectutor, scheduler
import pandas as pd
import json as json
import time

def data_get():  # noqa: E501
    """
    Get data from container. This function should only be called if program is
    in sync, broadcast or finished mode

     # noqa: E501
    :rtype: List[ActionResult]
    """
    current_app.logger.info('Calling data_get() -- docker status '+current_app.status)
    if current_app.status not in ['sync', 'broadcast', 'finished']:
        pass
    else:
        current_app.logger.info('Status: '+current_app.status)
        data=current_app.outbox
        current_app.logger.info('Returning outbox')
        current_app.logger.info(data)
        response = ActionResult(data=data)
        #current_app.Scheduler.execute_current_step()
        return response


def data_put(data=None):  # noqa: E501
    """Pass data to the container

    https://connexion.readthedocs.io/en/latest/request.html
    Zalando connexion> header parameters not passed to handler functions as parameters
    Header parameters have to be accessed via connexion.request.headers['<parameter>']

    Used to pass data from slave into the container # noqa: E501
    This function should be used to pass data to the container without triggering the
    next action. It can be used for asynchronous computation once feature is
    implemented.

    :param dataComplete: 0 - if data is to come, 1 - if all slaves has sent their data.
    Flag 1 needs to be used as a trigger for next action in the container
    :type dataComplete: int
    :param data: Additional data to pass to action (sync data), JSON format
    :type data: str

    :rtype: None
    """

    return 'do some magic'


def interrupt_post():  # noqa: E501
    """Stop operations in container

    The purpose of this call is to interrupt container operations and
    it should put container in sync status. # noqa: E501

    :rtype: None
    """
    # this should reset the state
    return 'do some magic!'


def log_get(level):  # noqa: E501
    """Get current logs from container

     # noqa: E501

    :param level: Log level:  * &#x60;info&#x60; - returns logs of info level and higher, in this case, everything  * &#x60;warning&#x60; - returns logs of warning level and higher (warning, error, fatal)  * &#x60;error&#x60; - returns logs of error level and higher (error, fatal)  * &#x60;fatal&#x60; - returns logs of fatal level 
    :type level: str

    :rtype: List[LogResponse]
    """
    return 'do some magic!'


def setup_post():  # noqa: E501
    """setup Docker for proper operation.

    This is used send configuration settings for the container # noqa: E501

    :param data: Additional configuration data, JSON format
    :type data: str

    :rtype: None
    """

    current_app.logger.info('Calling setup_post() -- docker status ' + current_app.status)
    ismaster = connexion.request.json['isMaster']
    params = connexion.request.json['params']
    params= parameter_parser(params)
    exectutor.setup(ismaster, params)
    response = StatusResponse(status=current_app.status)
    current_app.status = 'waiting'
    return response

def parameter_parser(params):
    '''
    Parse setup parameters into dictionary.
    >>> parameter_parser([{'key': 'number_type', 'value': 5}, {'key': 'text_type', 'value': 'Text'}, {'key': 'select', 'value': 'second'}, {'key': 'boolean', 'value': True}, {'key': 'checkbox', 'value': ['first', 'third']}])
    {'number_type': 5, 'text_type': 'Text', 'select': 'second', 'boolean': True, 'checkbox': ['first', 'third']}
    :param params: Input array containing an array of dictionaries with the parameters to pass to the setup call.
    :return: A dictionary containing the values
    '''
    current_app.logger.info('Parsing parameters')
    new = {}
    for p in params:
        new[p['key']]=p['value']
    return new

def start_post(data=None):  # noqa: E501
    """Starts container operations,
    p
     # noqa: E501
    :param data: Additional data to pass to action, JSON format
    :type data: str

    :rtype: None
    """
    current_app.logger.info('Calling start_post() -- docker status ' + current_app.status)
    try:
        # get request data
        current_app.logger.info('Putting request data')
        current_app.logger.info(connexion.request.json)
        data = connexion.request.json
    except:
        data = None
    if data!=None:
        # Save request data to inbox
        current_app.logger.info('Data is not None')
        current_app.inbox = data
    # Execute the next step as described in the workflow file
    # this also sets the pointer to the next step
    current_app.Scheduler.execute_current_step()
    return StatusResponse(current_app.status)



def status_get():  # noqa: E501
    """get status of Docker container

    Returns the current state of the Docker container  # noqa: E501
    This call should also return the current sync id, which is the syncId of the
    step that has to be performed next.
    :rtype: List[StatusResponse]
    """

    current_app.logger.info('Calling status_get() -- docker status ' + current_app.status)
    try:
        status = current_app.status
        try:
            sync = current_app.Scheduler.get_current_step()
        except:
            # Sync Id should not be below 0.
            sync = -1
        response = StatusResponse(status=status, sync_id=sync)
    except:
        response = StatusResponse(status='error')
    return response
