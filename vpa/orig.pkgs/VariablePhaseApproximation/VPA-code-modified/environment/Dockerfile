# hash:sha256:21e6729b7f38726d944af2d56f3518376ec2f5e9a4a85c429d2130890be1c1b7
FROM registry.codeocean.com/codeocean/ubuntu:18.04.3

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gfortran=4:7.4.0-1ubuntu2.3 \
    && rm -rf /var/lib/apt/lists/*
