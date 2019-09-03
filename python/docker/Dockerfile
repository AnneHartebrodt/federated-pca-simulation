# https://medium.com/google-cloud/deploying-scipy-numpy-with-docker-205e9afac3b0
# [START docker]
FROM python:3
# [END docker]# [START system dependencies]
RUN apt-get update
RUN apt-get install -y \  
    libpng-dev \
    freetype* \
    libblas-dev \
    liblapack-dev \
    libatlas-base-dev \
    gfortran
# [END system dependencies]# [START app environment]
WORKDIR /app  
ADD ../server/requirements.txt /app/requirements.txt
RUN pip install --upgrade pip  
# [END app environment]# [START install & configure wheel]
RUN pip install wheel  
ENV PIP_WHEEL_DIR=/wheelhouse  
ENV WHEELHOUSE=/wheelhouse  
ENV PIP_FIND_LINKS=/wheelhouse  
# [START install & configure wheel]
RUN while read p; do pip wheel $p; done < requirements.txt