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

	//Log Likelihood Ratio
	var llr = parseFloat($("#llr").val());
	//Log Likelihood Ratio Confidence Interval Lower value
	var ciLeft = parseFloat($("#ciLeft").val());
	//Log Likelihood Ratio Confidence Interval Upper value
	var ciRight = parseFloat($("#ciRight").val());

	//Rebase LLR to log(2) and display nicely (with sig.digits etc)
	function formatLLR() {
		//dividing by log(2) to re-base logarithm
		var rebase = Math.log(2)
		$("#llrDisplay").text(
			(llr/rebase).toFixed(2) + 
			" CI: [ " + (ciLeft/rebase).toFixed(2) + 
			" ; " + (ciRight/rebase).toFixed(2) + " ]"
		);
	}
	formatLLR();

	//basic logit function to turn log odds into probability
	function logit(k) {
		return Math.exp(k)/(1 + Math.exp(k));
	}

	//calculates the posterior from the entered prior and LLR
	//then displays it nicely
	function calcPosterior() {
		var prior = parseFloat($("#prior").val());
		var loPrior  = Math.log(prior/(1-prior))
		var posterior = logit(prior + llr)
		var postLeft = logit(prior + ciLeft)
		var postRight = logit(prior + ciRight)
		$("#posteriorDisplay").text(
			posterior.toFixed(2) + 
			" CI: [ " + postLeft.toFixed(2) + 
			" ; " + postRight.toFixed(2) + " ]" +
			" (Prior: "+prior.toFixed(2) + ")"
		);
		$("#postImgContainer").html(
			"<img src='drawInterval.R"+
			"?lower="+postLeft+
			"&amp;mid="+posterior+
			"&amp;upper="+postRight+
			"&amp;type=p'/>"
		);
	}

	function checkPriorInput() {
		var priorVal = $("#prior").val();
		if (priorVal && !isNaN(parseFloat(priorVal))) {
			calcPosterior();
			$("#priorDialog").dialog("close");
		} else {
			$("#prior").animate({
				backgroundColor: "red"
			}).animate({
				backgroundColor: "white"
			});
		}
	}

	//define the prior/posterior dialog box
	$("#priorDialog").dialog({
		autoOpen: false,
		modal: true,
		buttons: {
			OK: function() {
				checkPriorInput();
			},
			Cancel: function() {
				$(this).dialog("close");
			}
		}
	});

	$('#priorDialog').on('keyup', function(e){
		if (e.keyCode == $.ui.keyCode.ENTER) {
			$(':button:contains("OK")').click();
		}
	});

	//clicking the posterior button opens a dialog
	$("#priorButton").click(function(){
		//open the dialog
		$("#priorDialog").dialog("open");
	});


});