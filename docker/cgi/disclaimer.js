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

$(document).ready(function() {

	var cookiename = "maveclinDisclaimerAgreed=";

	function rememberAgreed(daysToExpire) {
		var d = new Date();
		d.setTime(d.getTime() + (daysToExpire*24*60*60*1000));
		document.cookie = cookiename + 1 + ";expires=" + d.toUTCString() + ";path=/";
	}

	function hasAgreed() {
		var ca = decodeURIComponent(document.cookie).split(';');
		for(var i = 0; i < ca.length; i++) {
			var c = ca[i];
			while (c.charAt(0) == ' ') {
				c = c.substring(1);
			}
			if (c.indexOf(cookiename) == 0) {
				var value = c.substring(cookiename.length, c.length);
				return value == "1"
			}
		}
		return false;
	}

	$("#disclaimerDialog").html(`
		<p>This tool and the information provided within 
		are intended for information and research purposes only and are not certified 
		for clinical use. Healthcare professionals using this tool should exercise 
		their own clinical judgement as to the information provided. All users 
		employ this tool and the information within at their own risk. Consumers 
		using this tool with respect to a specific condition are should always seek 
		professional medical advice and never act on the information in this tool alone. </p>
		<p>While the information within this tool is based on experimental data, 
		we do not warrant the accuracy of the information, nor that it covers all 
		possible applications, precautions or adverse effects. We do not assume any 
		liability or responsibility for damage, injury or death to users, other 
		persons or property due to actions taken following the use of this tool.</p>
		<p>This website uses cookies, by using it, you declare your acceptance of this use.</p>`
	);

	$("#disclaimerDialog").dialog({
		autoOpen: !hasAgreed(),
		modal: true,
	    closeOnEscape: false,
	    width: 600,
	    open: function(event, ui) {
	        $(".ui-dialog-titlebar-close", ui.dialog | ui).hide();
	    },
		buttons: {
			"I agree": function() {
				rememberAgreed(7);
				$(this).dialog("close");
			}
		}
	});
	//clicking the posterior button opens a dialog
	$("#disclaimerLink").click(function(){
		$("#disclaimerDialog").dialog("open");
	});


});