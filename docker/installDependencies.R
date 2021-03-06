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

install.packages(c("hash","httr","RJSONIO","devtools","RSQLite","flock"),repos="http://cran.utstat.utoronto.ca/")

library(devtools)

install_github("jweile/yogitools")
install_github("jweile/yogilog")
install_github("jweile/cgir")
install_github("VariantEffect/hgvsParseR")
install_github("VariantEffect/rapimave")
install_github("jweile/maveclin")
