class Configuration:
    def __init__(self, ismaster, params, logger):
        self._ismaster = ismaster
        self.logger = logger
        self.workflow_definition = self.parse_configuration_file(params['workflow_description'], ismaster)
        self._params = params
        self.parse_params()


    def parse_configuration_file(self, workflow_definition_file, ismaster):
        '''
        Parse workflow configuration file which contains the description of
        the internal docker workflow to be performed. Linewise entries are of
        the form, no header:
        syncId;stateAfterExecutionOfStep;masterOrBoth;nameOfFuntion;Description
        Example :
        1;sync;both;init;start workflow
        2;broadcast;master;step2;waiting for data from slaves then running global PCA
        3;finished;both;finalize;saving data
        :param workflow_definition_file: Name of the workflow definition file
        :param ismaster: True if Docker container runs as master
        :return: workflow definition as list of dictionaries
        '''
        workflow_definition = []
        with open(workflow_definition_file, 'r') as handle:
            for line in handle.readlines():
                line = line.strip()
                ll = line.split(';')
                if(ll[2]=='master'):
                    if(ismaster):
                        d = {}
                        d['step'] = ll[0]
                        d['container_status'] = ll[1]
                        d['command'] = ll[3]
                        d['description'] = ll[4]
                        workflow_definition.append(d)
                elif(ll[2]=='slave'):
                    if(not ismaster):
                        d = {}
                        d['step'] = ll[0]
                        d['container_status'] = ll[1]
                        d['command'] = ll[3]
                        d['description'] = ll[4]
                        workflow_definition.append(d)
                else:
                    d = {}
                    d['step'] = ll[0]
                    d['container_status'] = ll[1]
                    d['command'] = ll[3]
                    d['description'] = ll[4]
                    workflow_definition.append(d)
        self.logger.info(workflow_definition)
        return workflow_definition


    def get_workflow_definition(self):
        '''
        Getter for workflow definition
        :return: workflow definition
        '''
        return self.workflow_definition

    def set_master(self, ismaster):
        '''
        Set whether docker container runs as master or slave
        :param ismaster:
        :return:
        '''
        self._ismaster = ismaster

    def get_master(self):
        '''
        Get exection modus of program (master or slave)
        :return:
        '''
        return self._ismaster

    def get_params(self):
        '''
        Get parsed parameters
        '''
        return self._params

    def parse_params(self):
        '''
        App specific function to parse the parameter object passed with the
        setup call.
        :return:
        '''
        if int(self._params['sample_ids']) != -1:
            self._params['sample_ids'] = int(self._params['sample_ids'])
        else:
            self._params['sample_ids']= None
        if  int(self._params['header']) !=-1:
            self._params['header'] = int(self._params['header'])
        else:
            self._params['header'] = None
        self._params['epsilon'] = float(self._params['epsilon'])
        self._params['delta'] = float(self._params['delta'])