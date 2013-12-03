var numClusters = 20,
    clusterNames = {},
	w = 960/2,
	h = window.innerHeight - 120,
	margin = 50,
	startDate = 1, 
	endDate = 400,
	minExpression = 0,
	maxExpression = 1000,
	y = d3.scale.linear().domain([maxExpression, minExpression]).range([0 + margin, h - margin]),
	x = d3.scale.linear().domain([startDate, endDate]).range([0 + margin, w - margin]),
	timepoints = [1.0,4.0,21.0,116.0,185.0,186.0,255.0,289.0,290.0,292.0,294.0,297.0,301.0,307.0,311.0,322.0,329.0,369.0,380.0,400],
    maxAverageExpression = 0,
    normalize=false;

var xAxis = d3.svg.axis().scale(x).tickSize(0).tickSubdivide(true),
    yAxis = d3.svg.axis().scale(y).ticks(4).orient("left");
var geneCluster;
var currentDendrogram;
var nodes;
var color = d3.scale.category20();
var sheet = document.styleSheets[0];
for (var i = 0; i < numClusters; i++) {
    clusterNames['Cluster_' + i] = 'C' + i;
    d3.select("#filters")
        .append('a')
        .attr('id', 'C' + i)
        .html(i);
    sheet.addRule("#filters a.C" + i, "background: " + color(i) + "; color: #fff;", 1);
    sheet.addRule(".C" + i+".highlight", "stroke: " + color(i) + ";", 1);
    sheet.addRule(".C" + i+".highlight .average", "stroke: black;", 1);
}

var vis = d3.select("#vis")
    .attr("width", w)
    .attr("height", h)
    .append("svg:svg")
    .attr("width", w)
    .attr("height", h)
    .append("svg:g");
			
var line = d3.svg.line()
    .x(function(d,i) { return x(d.x); })
    .y(function(d) { return y(d.y); });
					
var genes_clusters = {}; // cluster1 : [path1, path2] ;
//currentclusters = cluster1 : [path1,path2], cluster2:...
// data = currentclusters.toarray
var allClusters = {};
var currentClusters = {};
var averages = {};
var geneDescriptions = {};
var maxExpressions = {};
getGeneDescriptions();
for (var k = 0; k < numClusters; k++) {
    getCluster(k);
}
var data = [];

// Add the x-axis.
vis.append("svg:g")
  .attr("class", "x axis")
  .attr("transform", "translate(0," + parseInt(h - margin)+ ")")
  .call(xAxis);

// Add the y-axis.
vis.append("svg:g")
  .attr("class", "y axis")
  .attr("transform", "translate(" + margin + ", 0)")
  .call(yAxis);

/* Dendrogram code */

var radius = 960 / 2 / 2;

var cluster = d3.layout.cluster()
    .size([360, radius - 250]);

var diagonal = d3.svg.diagonal.radial()
    .projection(function(d) { return [d.y, d.x / 180 * Math.PI]; });

var svg = d3.select("#cluster-info")
.append("svg")
    .attr("width", radius * 2)
    .attr("height", radius * 2)
  .append("g")
    .attr("transform", "translate(" + radius + "," + radius + ")");

d3.json("ProjectCode_v1/Output/5FPKMNormalizedAligned/geneCluster.json", function(error, root) {
  
  geneCluster = root;
  currentDendrogram = geneCluster;
  nodes = cluster.nodes(currentDendrogram);


});

d3.select(self.frameElement).style("height", radius * 2 + "px");


function onclick(d, i) {
    var currClass = d3.select(this).attr("class");
    if (d3.select(this).classed('selected')) {
        d3.select(this).attr("class", currClass.substring(0, currClass.length-9));
        
    } else {
        d3.select(this).classed('selected', true);
        
    }
}

function onmouseover(d, i) {
    var currClass = d3.select(this).attr("class");
    d3.select(this)
        .attr("class", currClass + " current");
    
    var geneName = $(this).attr("name").split(":")[0];
    var geneDescription = '';
    if (geneName in geneDescriptions) {
        geneDescription = geneDescriptions[geneName];
    }
    var blurb = '<h2>' + geneName + '</h2>';
    blurb += '<p>' + geneDescription + '</p>'    
    $("#default-blurb").hide();
    $("#blurb-content").html(blurb);
}

function onmouseout(d, i) {
    var currClass = d3.select(this).attr("class");
    var prevClass = currClass.substring(0, currClass.length-8);
    d3.select(this)
        .attr("class", prevClass);
    $("#default-blurb").show();
    $("#blurb-content").html('');
}

function showCluster(clusterCode) {
    var genes = d3.selectAll("path."+clusterCode);
    updateClusters(clusterCode);
    updateDendrogram();
    /*
    if (genes.classed('highlight')) {
        genes.classed('highlight', false);
        updateClusters(clusterCode);
    } else {
        genes.classed('highlight', true);
        updateClusters(clusterCode);
    }*/
    redraw();
    redrawDendrogram();
}

function getCluster(num) {
    d3.text('ProjectCode_v1/Output/0.05FPKMAgglomerativeClusters/OutputDirectoryHeirarchicalNoOfClusters_20/Cluster_' + num + '.csv', 'text/csv', function(text) {
    var clusterName = 'C' + num;
    allClusters[clusterName] = new Array();
    var genes = d3.csv.parseRows(text);
    var average = new Array(20);
    var allSignals = [];
    for (var a = 0; a < 20; a++) {
        average[a] = 0;
    }
    for (i=1; i < genes.length; i++) {
        var values = genes[i].slice(1, genes[i.length-1]);
        var currData = [];
        var currMaxExpression = 0;
        genes_clusters[genes[i][0]] = clusterName;
        var normCurrData = [];

        for (j=0; j < values.length; j++) {
            if (values[j] != '') {
                currData.push({ x: timepoints[j], y: parseFloat(values[j]) });  
                //average[j] += parseFloat(values[j])/genes.length;
                currMaxExpression = Math.max(currMaxExpression,parseFloat(values[j]));
                allSignals.push(toString(values));
            }
        }
        for (j=0; j < values.length; j++) {
            if (values[j] != '') {
                average[j] += parseFloat(values[j])/currMaxExpression/genes.length;
            }
        }
        maxExpressions[genes[i][0]] = currMaxExpression;
        allClusters[clusterName].push([genes[i][0], currData, clusterName])
    }
    var currMaxExpression = 0;
    //average = getAverageSignal(allSignals);
    for (l = 0; l< average.length; l++) {
        average[l] = {x: timepoints[l], y: average[l]};
        currMaxExpression = Math.max(average[l].y, currMaxExpression);
        //average[l] = {x: timepoints[l], y: average[l]};
        //currMaxExpression = Math.max(average[l].y, currMaxExpression);
    }
    maxAverageExpression = Math.max(maxAverageExpression, currMaxExpression);
    maxExpressions[clusterName + " average"] = currMaxExpression;
    averages[clusterName + " average"] = average;
    updateData();
    redraw();
}); 
}

function getGeneDescriptions() {
    d3.text('ProjectCode_v1/Output/5FPKMNormalizedAligned/geneDescriptions.txt', 'text/csv', function(text) {
        var genes = d3.csv.parseRows(text);
        for (i=1; i < genes.length; i++) {
            geneDescriptions[genes[i][0]] = genes[i][1];
        }
    });  
}

function updateClusters(clusterCode) {
    if (clusterCode in currentClusters) {
        delete currentClusters[clusterCode];
    } else {
        currentClusters[clusterCode] = allClusters[clusterCode];
    }

    updateData();
}

function updateData() {
    data = [];
    var maxValue = 0;
    for (c in currentClusters) {
        for (p in currentClusters[c]) {
            var pathinfo = currentClusters[c][p];
            var copiedpathdata = makeArrayCopy(pathinfo[1]);
            data.push({pathname:pathinfo[0], pathdata: copiedpathdata,clustername:pathinfo[2]});
            maxValue = Math.max(maxValue, maxExpressions[pathinfo[0]]);
        }
        var copiedpathdata = makeArrayCopy(averages[c + " average"]);
        data.push({pathname:c + " average", pathdata: copiedpathdata, clustername:c + " average"});
    }
    if (maxValue === 0) {
        maxExpression = maxAverageExpression;
        for (average in averages) {
            var copiedpathdata = makeArrayCopy(averages[average]);
            data.push({pathname:average, pathdata: copiedpathdata, clustername:average});
        }
    } else {
        maxExpression = maxValue;
    }

    if (normalize) {
        normalizeData();
    }
}

function redraw() {

    y.domain([maxExpression, minExpression]);
    vis.transition().duration(1000).select(".y.axis").call(yAxis);

    var paths = vis.selectAll("path")
        .data(data, function(d) { return d.pathname; });
    paths.enter().insert("svg:path", "g")
        .attr("name", function(d) { return d.pathname; })
        .attr("class", function(d) { return d.clustername.concat(" highlight"); })
        .attr("d", function(d) { return line(d.pathdata); })
        .on("mouseover", onmouseover)
        .on("mouseout", onmouseout);

    paths.transition().duration(1000)
        .attr("d", function(d) { return line(d.pathdata); })

    paths.exit().remove();
 
}

function changeNormalized() {
    normalize = !normalize;
    if (normalize) {
        $("#normalizebutton").text("Denormalize");
    } else {
        $("#normalizebutton").text("Normalize");
    }
    updateData();
    redraw();
}

function normalizeData() {
    for (i in data) {
        var pathMax = maxExpressions[data[i].pathname];
        if (data[i].pathname.split(" ").length !=2) {
            for (j in data[i].pathdata) {
                data[i].pathdata[j].y = data[i].pathdata[j].y/pathMax;
            }
        }
    }
    maxExpression = 1;
}
function makecopy(obj) {
    return jQuery.extend(true, {}, obj);
}

function makeArrayCopy(array) {
    var newArray = new Array(array.length);
    for (var i=0; i< newArray.length; i++) {
        newArray[i] = makecopy(array[i]);
    }
    return newArray;
}

function updateDendrogram() {
    var newObject = {name: "Gene Clusters"};
    var children = [];
    for (i in geneCluster.children) {
        var tempcluster = geneCluster.children[i];
        console.log(tempcluster);
        if (clusterNames[tempcluster.name] in currentClusters) {
            children.push(tempcluster);
        }
    }
    newObject["children"] = children;
    currentDendrogram = newObject;
    

}

function redrawDendrogram() {

    nodes = cluster.nodes(currentDendrogram);
    console.log(nodes);
    var link = svg.selectAll("path.link")
      .data(cluster.links(nodes));
    link.enter().append("path")
      .attr("class", "link")
      .attr("d", diagonal);

    link.transition().duration(1000)
        .attr("d", diagonal);

    link.exit().remove();

    var node = svg.selectAll("g.node")
        .data(nodes);

    var nodeEnter = node.enter().append("g")
        .attr("class", "node")
        .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; })

    nodeEnter.append("circle")
        .attr("r", 4.5);

    nodeEnter.append("text")
        .attr("dy", ".31em")
        .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
        .attr("transform", function(d) { return d.x < 180 ? "translate(8)" : "rotate(180)translate(-8)"; })
        .text(function(d) {
            if (typeof d.description != 'undefined') {
                return d.name + " ("+d.description+")";
            } else{
                return d.name;
            }
            
        });

    var nodeUpdate = node.transition().duration(1000).attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; })
    nodeUpdate.select("text").attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
        .attr("dy", ".31em")
        .attr("transform", function(d) { return d.x < 180 ? "translate(8)" : "rotate(180)translate(-8)"; })
        .text(function(d) {
            if (typeof d.description != 'undefined') {
                return d.name + " ("+d.description+")";
            } else{
                return d.name;
            }
            
        });

    node.exit().remove();
}

var thresholdPercent = 0.1

function toString(signal){
    floatSignal = toFloat(signal);
    thresholdValue = getThresholdValue(floatSignal);
    charArray = new Array(signal.length);
    charArray[0] = 'S';
    for (var i=1; i< signal.length; i++) {
        charArray[i] = getChar((floatSignal[i] - floatSignal[i-1]),thresholdValue);
    }
    return charArray;
}

function toFloat(signal){
    var floatArray = new Array(signal.length);
    for (var i=0; i< signal.length; i++) {
        floatArray[i] = parseFloat(signal[i]);
    }

    return floatArray;
}

function getChar(delta,threshold) {
    if (delta >= threshold) {
        return 'R'; // rise
    }       
    else if (delta <= -1 * threshold) {
        return 'D'; //drop
    }        
    return 'S'; //# steady
 }   

 function getThresholdValue(array) {
    maxValue = -Infinity;
    minValue = Infinity;
    for (var i = 0; i < array.length; i++ ) {
        maxValue = Math.max(maxValue, array[i]);
        minValue = Math.min(minValue, array[i]);
    }
    return thresholdPercent * (maxValue - minValue)
 }
var convert = {'R' : 0, 'S': 1, 'D' : 2};
var strength = {0:1, 1:0, 2:-1};
function getAverageSignal(signals) {
    var averagedSignal = new Array(signals[0].length);
    var PWM = [new Array(signals[0].length), new Array(signals[0].length), new Array(signals[0].length)];
    var minValue = Infinity;
    for (var i = 0; i< PWM[0].length; i++) {
        PWM[0][i] = 0;
        PWM[1][i] = 0;
        PWM[2][i] = 0;
    }
    for (var i = 0; i< signals.length; i++) {
        for (var j = 0; j < signals[0].length ; j++) {
            PWM[convert[signals[i][j]]][j] ++;
        }
    }
    for (var i = 0; i< PWM[0].length; i++) {
        var idx = 0;
        var maxFreq = -Infinity;
        for (var j = 0; j< PWM.length; j++) {
            if (maxFreq < PWM[j][i]) {
                idx = j;
                maxFreq = PWM[j][i];
            }
        }
        if (i === 0) {
            averagedSignal[i] = 0;
        } else {
            averagedSignal[i] = averagedSignal[i-1] + strength[idx];
        }
        minValue = Math.min(minValue, averagedSignal[i]);
        console.log(minValue);
    }
    for (var i = 0; i< averagedSignal.length; i++) {
        averagedSignal[i] = averagedSignal[i] - minValue;
    }
    console.log(averagedSignal);
    return averagedSignal;
}
