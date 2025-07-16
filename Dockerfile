FROM python:3.9-slim

WORKDIR /app


COPY requirements.txt /app/

RUN pip install --upgrade pip && \
    pip install --no-cache-dir --upgrade --force-reinstall -r requirements.txt

COPY main.py /app
COPY pipeline.py /app
COPY /src /app/src


CMD ["python", "main.py"]