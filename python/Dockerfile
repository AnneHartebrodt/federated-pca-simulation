FROM datascience-container:latest

RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app
COPY server/requirements.txt /usr/src/app/

RUN while read p; \
    do pip install --no-index -f /wheelhouse $p; \
    done < requirements.txt

ADD . /usr/src/app
COPY . /usr/src/app

COPY server/workflow_description.tsv /usr/src/app


EXPOSE 80 3306 6379

#ENTRYPOINT ["python3"]

#CMD ["-m", "swagger_server"]
CMD python -m server.swagger_server > /mounted_local_folder/error.log