/* 
 * Copyright (C) 2018  Jochen Weile, Roth Lab
 *
 * This file is part of MaveClin.
 *
 * MaveClin is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MaveClin is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with MaveClin.  If not, see <https://www.gnu.org/licenses/>.
 */

$(document).ready(function(){

	var urn = $("#urn").text();

	//if there is no gold standard table to show, hide the box
	if ($("#tablebox").find("table").length == 0) {
		$("#imgbox").hide();
		$("#tablebox").hide();
	}
	$("#statusindicator").hide();

	$("#calibrateButton").click(function(){
		runCalibration();
	});

	$("#errordialog").dialog({
		autoOpen: false,
		buttons: {
			Close: function() {
				$(this).dialog("close");
			}
		}
	}).parent().addClass("ui-state-error");

	function showError(text) {
		$("#errormessage").text(text);
		$("#errordialog").dialog("open");
	}

	function statusMessage(msg) {
		if (msg == null) {
			$("#statusindicator").hide();
		} else {
			$("#statusmessage").text(msg);
			$("#statusindicator").show();
		}
	}

	function runCalibration() {

		var symbol = $("#genesymbol").val();
		var ensemblGeneID = $("#ensemblID").val();
		var flip = $("#flip").prop("checked");
		var mafCutoff = $("#mafCutoff").val();
		var homozygous = $("#homozygous").prop("checked");

		window.console&&console.log(
			"urn: "+urn+
			"\n symbol: "+symbol+
			"\n ensemblGeneID: "+ensemblGeneID+
			"\n mafCutoff: "+mafCutoff+
			"\n flip: "+flip+
			"\n homozygous: "+homozygous
		);

		statusMessage("pending");
		$("#calibrateButton").prop("disabled",true);

		$.post("runCalibration.R",{
			urn : urn,
			symbol : symbol,
			ensemblGeneID : ensemblGeneID,
			mafCutoff : mafCutoff,
			flip : flip,
			homozygous : homozygous
		}).done(function(data) {
			// data = JSON.parse(rawdata);
			switch(data.response) {
				case "busy":
					showError("Already processing existing request!")
					statusMessage(null);
				break;
				case "submitted":
					window.console&&console.log("Successfully submitted! Scheduling poll...");
					setTimeout(pollStatus,3000);
				break;
				default:
					showError("Unknown response:"+data.response);
					statusMessage(null);
			}
		}).fail(function(xhr,status,error) {
			showError(error);
			statusMessage(null);
			$("#calibrateButton").prop("disabled",false);
		});
	}

	function pollStatus() {
		window.console&&console.log("Polling status...");
		$.post("queryStatus.R",{
			urn : urn
		}).done(function(data) {
			// data = JSON.parse(rawdata);
			switch(data.status) {
				case "pending":
					//same as processing, so no 'break'
				case "processing":
					statusMessage(data.status);
					window.console&&console.log("Pending or processing. Scheduling another poll...");
					setTimeout(pollStatus,1000);
					break;
				case "calibrated":
					window.console&&console.log("Done. Loading table...");
					loadTable();
					statusMessage(null);
					$("#calibrateButton").prop("disabled",false);
					break;
				case "error":
					showError("Calibration failed!");
					//no 'break' here, so we also apply the 'new' condition
				case "new":
					showError(
						"Dataset status reported as new."+
						"\nPlease submit a bug report!"
					);
					statusMessage(null);
					$("#calibrateButton").prop("disabled",false);
					break;
				default:
					statusMessage(null);
					showError(
						"Unknown dataset status: "+data.status+
						"\nPlease submit a bug report!"
					);
			}
		}).fail(function(xhr,status,error) {
			showError(error);
			statusMessage(null);
			$("#calibrateButton").prop("disabled",false);
		});
	}

	function loadTable() {
		window.console&&console.log("Querying table...");
		$.post("getScoresetTable.R",{
			urn : urn
		}).done(function(data) {
			$("#tablebox").html(
				"<h3>Gold Standard Variants</h3>" + 
				data.table
			);
			$("#imgbox").find("img").attr("src",data.imgTarget);
			$("#tablebox").show();
			$("#imgbox").show();
		}).fail(function(xhr,status,error) {
			showError(error);
		});
	}

});