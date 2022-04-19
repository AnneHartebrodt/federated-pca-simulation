import pandas as pd
import scipy as sc
import numpy as np
import os.path as op
import os
import re




def create_dataframe(basepath, suffix):
    data_list = []
    for d in os.walk(basepath):
        print(d)
        if len(d[2])>0 and len(d[1])==0:
            levels = d[0].replace( basepath, "").split(('/'))
            fff = list(filter(lambda x: re.search(re.compile(suffix), x), d[2]))
            for f in fff:
                current_data =  pd.read_csv(op.join(d[0],f), sep='\t', header=None)
                counter = 0
                for l in levels:
                    if l != "":
                        current_data["C"+str(counter)] = l
                        counter = counter+1
                print(f)
                current_data['filename']= f
                data_list.append(current_data)

    data_list = pd.concat(data_list, axis=0)
    print(op.join(basepath, 'summary'+suffix+'.tsv'))
    data_list.to_csv(op.join(basepath ,'summary'+suffix+'.tsv'), sep='\t', header=False, index=False)
    return data_list

if __name__ == '__main__':
    basepath = "/home/anne/Documents/featurecloud/pca/approximative-vertical/results/matrix"
    # create_dataframe(basepath=basepath, suffix='.angles.u')
    # create_dataframe(basepath=basepath, suffix='.angles.v')
    # df = create_dataframe(basepath=basepath, suffix='.mev.u')
    # df = create_dataframe(basepath=basepath, suffix='.mev.v')
    # #
    # df = create_dataframe(basepath=basepath, suffix='.eo')
    # df = create_dataframe(basepath=basepath, suffix='.transmission')
    df = create_dataframe(basepath=basepath, suffix='.transmission')