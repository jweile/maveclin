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

COMMIT := $(shell git rev-parse HEAD|cut -c1-8)

build:	
	docker build --build-arg UPDATE_RPKG=`date +%s` -t jweile/maveclin .
	# docker build --build-arg UPDATE_RPKG=`date +%s` -t jweile/maveclin:$(COMMIT) .
	# docker tag jweile/maveclin:$(COMMIT) jweile/maveclin:latest

clean:
	docker stop maveclin&&docker container rm maveclin

run:
	# nohup docker run -t -p 80:80 --name maveclin jweile/maveclin:$(COMMIT) &
	nohup docker run -t -p 80:80 --name maveclin jweile/maveclin &

runwmount:
# 	nohup docker run -t -p 80:80 --name maveclin \
# 		--mount type=bind,source=/home/jweile/projects/,target=/projects,readonly \
# 		jweile/maveclin&:$(COMMIT) &
	nohup docker run -t -p 80:80 --name maveclin \
		--mount type=bind,source=/home/jweile/projects/,target=/projects,readonly \
		jweile/maveclin &

push:
	docker push jweile/maveclin:$(COMMIT)
	docker push jweile/maveclin:latest
