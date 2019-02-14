# Copyright (C) 2018  Jochen Weile, Roth Lab
#
# This file is part of MaveVis.
#
# MaveVis is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MaveVis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with MaveVis.  If not, see <https://www.gnu.org/licenses/>.

library("maveclin")

context("persistence")

setup <- function() {
	# sqlfile <- system.file("testdata","dbInit.sql",package="maveclin")
	sqlfile <- "docker/dbInit.sql"
	stopifnot(file.exists(sqlfile))
	dbfile <- tempfile(fileext=".db")
	system(paste("sqlite3",dbfile,"<",sqlfile))
	return(dbfile)
}
teardown <- function(dbfile) {
	stopifnot(file.remove(dbfile))
}

test_that("scoresets",{
	
	dbfile <- setup()
	ps <- new.persistence.connection(dbfile)

	#check that isKnown() returns false where no entry exists
	expect_false(ps$isKnown("urn:mavedb:00000001-a-1"))
	#and that functions complain when the argument is invalid
	expect_error(ps$isKnown("urn:mavedb:00000001-a"))

	#check that adding a new scorest works
	ps$newScoreset("urn:mavedb:00000001-a-1","UBE2I","ENSG00000103275")
	expect_true(ps$isKnown("urn:mavedb:00000001-a-1"))
	expect_equal(ps$getStatus("urn:mavedb:00000001-a-1"),"new")
	#check that invalid scoreset data is rejected
	expect_error(ps$newScoreset("urn:mavedb:00000001-a-1","\"; DROP TABLE scoresets;","ENSG00000103275"))

	mafCutoff <- 0.2
	flip <- TRUE
	ps$setParameters("urn:mavedb:00000001-a-1",
		"UBE2I","ENSG00000103275",mafCutoff=mafCutoff,flip=flip,homozygous=FALSE
	)
	values <- ps$getParameters("urn:mavedb:00000001-a-1")
	expect_equal(values$mafCutoff,mafCutoff)
	expect_equal(values$flip,flip)
	expect_equal(ps$getStatus("urn:mavedb:00000001-a-1"),"pending")

	ps$close()
	teardown(dbfile)
	
})