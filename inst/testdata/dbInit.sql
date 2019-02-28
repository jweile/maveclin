/*
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
*/
CREATE TABLE IF NOT EXISTS scoresets (
	urn TEXT PRIMARY KEY, symbol TEXT, ensemblGeneID TEXT, mafCutoff REAL, 
	flip BOOLEAN, homozygous BOOLEAN, lastUpdated DATE, status TEXT
);
CREATE TABLE IF NOT EXISTS variants (
	accession TEXT PRIMARY KEY, scoreset TEXT, hgvs_pro TEXT, score REAL, sd REAL, se REAL,
	flippedScore REAL, clinsig TEXT, maf REAL, hom INTEGER, refset TEXT,
	llr REAL, llrCIleft REAL, llrCIright REAL,
	FOREIGN KEY(scoreset) REFERENCES scoresets(urn)
);
CREATE INDEX idx_var_ss ON variants (scoreset);
CREATE INDEX idx_var_hvgs ON variants (hgvs_pro);
