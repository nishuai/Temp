

pephighlight = function(i, tooltip) {
    var pep = document.getElementById("pep" + i+".1");
    pep.setAttribute("stroke-width", "6px");
    pep = document.getElementById("pep" + i);
    rect = pep.getBoundingClientRect();
    // boxrange = "Element Position Is : " + rect.top + " - " + rect.right + " - " + rect.bottom + " - " + rect.left;
    drawToolTip(rect, "bottom", tooltip);
    // var peprect = document.getElementById("peprect" + i);
    // var peptt = document.getElementById("peptt" + i);
    // peprect.setAttribute("visibility", "visible");
    // peptt.setAttribute("visibility", "visible");
}

pepdim = function(i) {
    elementText = document.getElementById("tooltiptext");
    if (elementText != null) {
	elementText.setAttribute("visibility", "hidden");
    }
    elementRect = document.getElementById("tooltiprect");
    if (elementRect != null) {
	elementRect.setAttribute("visibility", "hidden");
    }
    var pep = document.getElementById("pep" + i + ".1");
    pep.setAttribute("stroke-width", "2.25px");
    // var peprect = document.getElementById("peprect" + i);
    // var peptt = document.getElementById("peptt" + i);
    // peprect.setAttribute("visibility", "hidden");
    // peptt.setAttribute("visibility", "hidden");
}

gfhighlight = function(i, tooltip) {
    // alert(tooltip.length);
    pepgroup = document.getElementById("gf." + i + ".1");
    pepgroup.setAttribute("stroke-width", "6px");
    pep = document.getElementById("gf." + i);
    rect = pep.getBoundingClientRect();
    // boxrange = "Element Position Is : " + rect.top + " - " + rect.right + " - " + rect.bottom + " - " + rect.left;
    drawToolTip(rect, "top", tooltip);
    // alert("Hewllo");
}

gfdim = function(i) {
    // alert(i);
    elementText = document.getElementById("tooltiptext");
    if (elementText != null) {
	elementText.setAttribute("visibility", "hidden");
    }
    elementRect = document.getElementById("tooltiprect");
    if (elementRect != null) {
	elementRect.setAttribute("visibility", "hidden");
    }
    pepgroup = document.getElementById("gf." + i + ".1");
    pepgroup.setAttribute("stroke-width", "2.25px");
}

drawToolTip = function(rect, where, tooltip) {
    dist = 6
    var svgNS = "http://www.w3.org/2000/svg";
    var xlinkNS = "http://www.w3.org/1999/xlink";
    element = document.getElementById("tooltiptext");
    if (element == null) {
	newText = document.createElementNS(svgNS,"text");
	newText.setAttributeNS(null,"x",rect.left+(rect.right-rect.left)/2);
	newText.setAttributeNS(null,"y",rect.top);
	newText.setAttributeNS(null,"id","tooltiptext");
	newText.setAttributeNS(null,"font-size","13px");
	newText.setAttributeNS(null,"font-family","Arial");
	newText.setAttributeNS(null,"text-anchor","middle");
	newText.setAttributeNS(null,"fill-opacity",1);
	newText.setAttributeNS(null,"alignment-baseline","central");
	newText.setAttributeNS(null,"fill","rgb(0,0,0)");
	var textNode = document.createTextNode(tooltip);
	newText.appendChild(textNode);
	document.getElementsByTagName("svg")[0].appendChild(newText);

	var bbox = newText.getBBox();
	w = bbox.width+4;
	h = bbox.height;

	newRect = document.createElementNS(svgNS,"rect");
	newRect.setAttributeNS(null,"x",rect.left+(rect.right-rect.left)/2 - w/2);
	if (where == "top") {
	    newRect.setAttributeNS(null,"y",rect.top-h-dist);
	} else {
	    newRect.setAttributeNS(null,"y",rect.top+dist);
	}
	newRect.setAttributeNS(null,"width",w);
	newRect.setAttributeNS(null,"height",h);
	newRect.setAttributeNS(null,"stroke-width",1);
	newRect.setAttributeNS(null,"stroke","rgb(0,0,0)");
	newRect.setAttributeNS(null,"id","tooltiprect");
	newRect.setAttributeNS(null,"fill","rgb(255,255,224)");
	document.getElementsByTagName("svg")[0].appendChild(newRect);
	document.getElementsByTagName("svg")[0].removeChild(newText);
	if (where == "top") {
	    newText.setAttributeNS(null,"y",rect.top-h/2-dist);
	} else {
	    newText.setAttributeNS(null,"y",rect.top+h/2+dist);
	}
	document.getElementsByTagName("svg")[0].appendChild(newText);
    } else {
	elementText.setAttributeNS(null,"x",rect.left+(rect.right-rect.left)/2);
	elementText.setAttributeNS(null,"y",rect.top);
	elementText.firstChild.nodeValue = tooltip

	var bbox = elementText.getBBox();
	w = bbox.width+6;
	h = bbox.height;

	elementRect = document.getElementById("tooltiprect");

	elementRect.setAttributeNS(null,"x",rect.left+(rect.right-rect.left)/2 - w/2);
	if (where == "top") {
	    elementRect.setAttributeNS(null,"y",rect.top-h-dist);
	} else {
	    elementRect.setAttributeNS(null,"y",rect.top+dist);
	}
	elementRect.setAttributeNS(null,"width",w);
	elementRect.setAttributeNS(null,"height",h);
	if (where == "top") {
	    elementText.setAttributeNS(null,"y",rect.top-h/2-dist);
	} else {
	    elementText.setAttributeNS(null,"y",rect.top+h/2+dist);
	}
	elementRect.setAttribute("visibility", "visible");
	elementText.setAttribute("visibility", "visible");
    }
}


// pepmessage = function(i,mytext) {
//     alert(mytext);
// }
