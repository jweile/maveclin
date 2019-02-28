# Copyright (C) 2018  Jochen Weile, Roth Lab
#
# This file is part of MaveClin.
#
# MaveClin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MaveClin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with MaveClin.  If not, see <https://www.gnu.org/licenses/>.

# This file defines a docker image encapsulating a CGI-enabled apache server
# with a RestFul service and web-app for MaveClin with all its dependencies.
#
# A makefile is provided to assist with building and running

FROM ubuntu:18.04
MAINTAINER Jochen Weile <jochenweile@gmail.com>
ARG UPDATE_RPKG=0

#Force apt-get to skip confirmation questions
ENV DEBIAN_FRONTEND noninteractive

#Install apache httpd and enable CGI
RUN apt update && \
    apt -y dist-upgrade && \
    apt --no-install-recommends install -y apache2 && \
    a2enmod cgi headers

# Apache ports for HTTP and HTTPS
EXPOSE 80
EXPOSE 443
# EXPOSE 8004

#directory for software setup
RUN mkdir /setup
WORKDIR /setup

#install dependencies for app
RUN apt update && \
    apt --no-install-recommends install -y r-base \
    g++ make libcurl4-openssl-dev libssl-dev zlib1g-dev\
    libjson-c-dev libxml2 libxml2-dev cron sqlite3 
    #mariadb-client libmariadb-client-lgpl-dev
#    wget dssp clustalo  

#install R packages
#The UPDATE_RPKG argument doens't actually do anything in the script
#it just allows forcing docker to override the cache for this instruction
#by running $docker build --build-arg UPDATE_RPKG=`date +%s` -t jweile/maveclin .
COPY installDependencies.R /setup/
RUN Rscript installDependencies.R ${UPDATE_RPKG}


# Create maveclin's directory structure and grant access to apache
RUN mkdir -p /var/www/html/maveclin/httpdocs/imgcache &&\
    chmod 777 /var/www/html/maveclin/httpdocs/imgcache &&\
    mkdir -p /var/www/maveclin/logs &&\
    chmod 777 /var/www/maveclin/logs &&\
    mkdir -p /var/www/maveclin/cache/ &&\
    chmod 777 /var/www/maveclin/cache/

#set environment variables so the app knows were to write/read data
ENV MAVECLIN_CACHE /var/www/maveclin/cache/
ENV MAVECLIN_LOGS /var/www/maveclin/logs/
#Note: the apache user and cron cannot see these variables, so they must also be set in the
# apache config file (maveclin.conf)!!!
RUN printf "export MAVECLIN_CACHE=/var/www/maveclin/cache/\n">>/etc/apache2/envvars &&\
    printf "export MAVECLIN_LOGS=/var/www/maveclin/logs/\n">>/etc/apache2/envvars &&\
    printf "export MAVECLIN_CACHE=/var/www/maveclin/cache/\n">>/etc/environment &&\
    printf "export MAVECLIN_LOGS=/var/www/maveclin/logs/\n">>/etc/environment

#Move CGI scripts to staging location and grant access to apache
ADD cgi/ /var/www/html/maveclin/httpdocs/
# COPY maveclin_launcher.R /var/www/html/maveclin/
# COPY searchIndex.csv /var/www/maveclin/cache/
RUN chmod a+rx /var/www/html/maveclin/httpdocs/*.R #&&\
    # chmod a+rx /var/www/html/maveclin/*.R &&\
    # mkdir /var/www/html/maveclin/httpdocs/results/ &&\
    # chmod 777 /var/www/html/maveclin/httpdocs/results/ &&\
    # chmod 666 /var/www/maveclin/cache/searchIndex.csv

#Introduce apache configuration files for app and enable it
COPY security.conf /etc/apache2/conf-available/
COPY maveclin.conf /etc/apache2/sites-available/
RUN a2dissite 000-default.conf && a2ensite maveclin

# Print logs to stdout
# RUN \
#   ln -sf /proc/self/fd/1 /var/log/apache2/access.log && \
#   ln -sf /proc/self/fd/1 /var/log/apache2/error.log && \
#   # ln -sf /proc/self/fd/1 /var/www/maveclin/logs/exec.log && \
#   ln -sf /proc/self/fd/1 /var/www/maveclin/logs/error.log && \
#   ln -sf /proc/self/fd/1 /var/www/maveclin/logs/access.log

#setup startup script and daemon
# COPY startup.sh daemon.R sync.R /setup/
COPY startup.sh sync.R calibrationCron.R dbInit.sql /setup/

#run DB init script and allow www user to write to DB
RUN sqlite3 ${MAVECLIN_CACHE}maveclin.db < dbInit.sql &&\
    chmod 666 ${MAVECLIN_CACHE}maveclin.db


#setup cronjob for synchronization script
ADD crontab /etc/cron.d/maveclin-cron
RUN chmod 0644 /etc/cron.d/maveclin-cron &&\
    crontab /etc/cron.d/maveclin-cron &&\
    chmod a+x /setup/sync.R /setup/calibrationCron.R

#startup script starts apache and daemon
CMD bash /setup/startup.sh
# CMD cron && rm -rf /run/httpd/* && apachectl -DFOREGROUND
