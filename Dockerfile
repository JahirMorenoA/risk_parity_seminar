FROM rocker/tidyverse

#RUN apt-get update

RUN install2.r --error --deps TRUE devtools
