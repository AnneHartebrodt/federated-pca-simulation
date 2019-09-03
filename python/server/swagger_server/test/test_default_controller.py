# coding: utf-8

from __future__ import absolute_import

from flask import json
from six import BytesIO

from server.swagger_server.models.action_result import ActionResult  # noqa: E501
from server.swagger_server.models.log_response import LogResponse  # noqa: E501
from server.swagger_server.models.status_response import StatusResponse  # noqa: E501
from server.swagger_server.test import BaseTestCase


class TestDefaultController(BaseTestCase):
    """DefaultController integration test stubs"""

    def test_data_get(self):
        """Test case for data_get

        get data from container
        """
        response = self.client.open(
            '/balazsorban/ControllerToDocker/1.0.1/data',
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_data_put(self):
        """Test case for data_put

        Pass data to the container
        """
        headers = [('dataComplete', 56),
                   ('data', 'data_example')]
        response = self.client.open(
            '/balazsorban/ControllerToDocker/1.0.1/data',
            method='PUT',
            headers=headers)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_interrupt_post(self):
        """Test case for interrupt_post

        Stop operations in container
        """
        response = self.client.open(
            '/balazsorban/ControllerToDocker/1.0.1/interrupt',
            method='POST')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_log_get(self):
        """Test case for log_get

        Get current logs from container
        """
        response = self.client.open(
            '/balazsorban/ControllerToDocker/1.0.1/log'.format(level='level_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_setup_post(self):
        """Test case for setup_post

        setup Docker for proper operation.
        """
        headers = [('data', 'data_example')]
        response = self.client.open(
            '/balazsorban/ControllerToDocker/1.0.1/setup',
            method='POST',
            headers=headers)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_start_post(self):
        """Test case for start_post

        Starts container operations
        """
        headers = [('data', 'data_example')]
        response = self.client.open(
            '/balazsorban/ControllerToDocker/1.0.1/start',
            method='POST',
            headers=headers,
            content_type='application/json')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_status_get(self):
        """Test case for status_get

        get status of Docker container
        """
        response = self.client.open(
            '/balazsorban/ControllerToDocker/1.0.1/status',
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
