FROM python:3.6
RUN pip install pandas
RUN pip install matplotlib
RUN pip install plotly
RUN pip install jinja2
RUN mkdir /code /resources
RUN mkdir -p /html_report
WORKDIR /code
COPY ./code/*.py /code/
RUN chmod +x /code/*py
ENTRYPOINT ["python","/code/scip.py"]