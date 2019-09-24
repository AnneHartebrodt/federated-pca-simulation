FROM datascience-container:latest

RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app
COPY requirements.txt /usr/src/app/

RUN while read p; \
    do pip install --no-index -f /wheelhouse $p; \
    done < requirements.txt

ADD simulation /usr/src/app
COPY simulation /usr/src/app



EXPOSE 80 3306 6379

ENTRYPOINT ["python3"]

CMD ["-m", "real_data", "-f", "$1", "-o", "$2"]