FROM ubuntu:jammy

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
                git \
                gfortran \
                python3 \
                python3-pip \
                ipython3 \
                && rm -rf /var/lib/apt/lists/*

WORKDIR /home/user/
ADD ./ ./capytaine/
RUN pip install ./capytaine
