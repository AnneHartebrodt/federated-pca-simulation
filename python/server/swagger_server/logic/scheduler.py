from server.swagger_server.logic import exectutor as exe
from flask import current_app


class Scheduler():

    def __init__(self, schedule):
        if len(schedule)==0:
            raise ValueError('Your internal Container workflow specification is empty. Please provide a valid specification!')
        flattened = []
        flattened = list(flattened.extend(item) for item in schedule)
        if(flattened.__contains__('paused')&flattened.__contains__('error')):
            raise ValueError('Your internal Container workflow specification should not containe error or paused states')


        self.schedule = schedule
        self.currentStep = 0

        """
        These are possible container states
            "init"
            - "waiting"
            - "sync"
            - "broadcast"
            - "running"
            - "paused"
            - "error"
            - "finished"
            """

    # schedule is currently of the form
    """
    """


    def get_current_step_index(self) -> int:
        '''
        Get numeric index of the current step
        :return: index
        '''
        return  self.currentStep

    def get_current_step(self):
        '''
        Get name of the current step
        :return: current step name
        '''
        if (len(self.schedule)>self.currentStep):
            return self.schedule[self.currentStep]['step']
        else:
            return self.schedule[-1]['step']

    def get_current_status(self)-> str:
        '''
        Get status docker container has to be in after succesful execution of
        the current step.
        :return: status
        '''
        if (len(self.schedule)>self.currentStep):
            return self.schedule[self.currentStep]['container_status']
        else:
            return self.schedule[-1]['container_status']

    def get_current_command(self) -> str:
        '''
        Return function name corresponding to function definition in executor file
        :return: function name
        '''
        if (len(self.schedule)>self.currentStep):
            return self.schedule[self.currentStep]['command']
        else:
            return self.schedule[-1]['command']

    def next_step(self):
        '''
        Move pointer to next step
        :return: No return
        '''
        if (len(self.schedule) > self.currentStep):
            self.currentStep = self.currentStep+1
        else:
            self.currentStep = len(self.schedule)-1

        current_app.logger.info(self.get_current_status())

    def execute_current_step(self):
        '''
        Execute the current step, add data to the outbox, if applicable, move
        pointer to next step
        :return: No return
        '''
        current_app.logger.info('Executing step:' +str(self.get_current_command()))
        current_app.status = 'running'
        result = getattr(exe, self.get_current_command())() # execute step
        if(result is not None):
            current_app.logger.info('New data in outbox')
            current_app.outbox = result
        current_app.logger.info('...done')
        current_app.status = self.get_current_status()
        self.next_step()  # move pointer to next step


