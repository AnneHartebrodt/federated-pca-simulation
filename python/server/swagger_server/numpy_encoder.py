import json
import numpy as np

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.ndarray, np.float64)):
            content = {
                'ndarray' :
                    [[float(y) for y in x] for x in obj]
            }
            print('np.array')
            return content
        return json.JSONEncoder.default(self, obj)

class NumpyDecoder(json.JSONDecoder):
    """

    """

    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, obj):
        print(obj)
        print(type(obj))
        if 'ndarray' in obj:
            return np.asarray(obj.get('ndarray'))
        return json.JSONDecoder.default(self,obj)