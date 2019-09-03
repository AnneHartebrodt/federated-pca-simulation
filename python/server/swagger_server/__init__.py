#!/usr/bin/env python3

import connexion

from server.swagger_server import encoder
import os as os
import logging
import logging.handlers
from flask.logging import default_handler

status = 'init'

#/mounted_local_folder
app = connexion.App(__name__, specification_dir='./swagger/')
app.app.json_encoder = encoder.JSONEncoder
app.add_api('swagger.yaml', arguments={'title': 'Controller To Docker API'})
app.app.inbox = None
app.app.outbox = None
app.app.status = 'init'

# set up internal logging for app
handler1 = logging.handlers.RotatingFileHandler(
        '/mounted_local_folder/applog1.txt',
        maxBytes=1024 * 1024)
handler1.setLevel(logging.DEBUG)


app.app.logger.removeHandler(default_handler)
app.app.logger.setLevel(logging.DEBUG)
app.app.logger.addHandler(handler1)
app.app.logger.info('Starting up!')
# request logging
handler = logging.handlers.RotatingFileHandler(
        '/mounted_local_folder/applogw.txt',
        maxBytes=1024 * 1024)
handler.setLevel(logging.DEBUG)
logging.getLogger('werkzeug').setLevel(logging.DEBUG)
logging.getLogger('werkzeug').addHandler(handler)

import server.swagger_server.controllers.default_controller
