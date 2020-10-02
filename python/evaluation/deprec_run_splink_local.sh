{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {
    "Collapsed": "false"
   },
   "source": [
    "## run sPlink PCA"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "metadata": {
    "Collapsed": "false"
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Process is interrupted.\n"
     ]
    }
   ],
   "source": [
    "%%bash\n",
    "conda activate fed-gwas\n",
    "cd /work/Documents/gwas/data/hapmap/raw\n",
    "export PYTHONPATH=$PYTHONPATH:/work/Documents/gwas/fed-gwas-client\n",
    "for pop in ASW CEU CHB CHD GIH JPT LWK MEX MKK TSI YRI\n",
    "do\n",
    "mkdir /work/Documents/gwas/results/hapmap/$pop\n",
    "cd $pop; \n",
    "python3 /work/Documents/gwas/fed-gwas-client/client/project/my_test_client.py --dataset hapmap3_r1_b36_fwd.${pop}.thin --pcs 10 --maf 0.01 --resultdir ../../../../results/hapmap/${pop} --fileprefix hapmap3_r1_b36_fwd.${pop}.thin\n",
    "cd ..\n",
    "done"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {
    "Collapsed": "false"
   },
   "source": [
    "## Run sPLINK PCA to obtain eigenvectors"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {
    "Collapsed": "false",
    "pycharm": {
     "name": "#%%\n"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Process is terminated.\n"
     ]
    }
   ],
   "source": [
    "%%bash\n",
    "export PYTHONPATH=$PYTHONPATH:/work/Documents/gwas/fed-gwas-client\n",
    "cd /work/Documents/gwas/data/1000g/raw\n",
    "for e in {1..22} ; do\n",
    "cd chr${e}\n",
    "python3 /work/Documents/gwas/fed-gwas-client/client/project/my_test_client.py --dataset chr${e}.thin --pcs 10 --maf 0.01 --resultdir ../../../../results/1000g/chr${e}/ --fileprefix chr${e}\n",
    "cd ..\n",
    "done\n",
    "\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.6"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}