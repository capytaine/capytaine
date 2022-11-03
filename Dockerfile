# Build with:
#  docker build -f Dockerfile -t capytaine .
#
# Run with:
#  docker run -it -v $(pwd):/home/user capytaine
#
FROM ubuntu:jammy

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
                git \
                gfortran \
                python3 \
                python3-pip \
                python3-numpy \
                ipython3 \
                python-is-python3 \
                && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/
RUN git clone --recurse-submodules https://github.com/capytaine/capytaine
RUN pip install ./capytaine

WORKDIR /home/user/
CMD ["ipython3"]
